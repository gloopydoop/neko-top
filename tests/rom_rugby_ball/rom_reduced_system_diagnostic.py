#!/usr/bin/env python3
"""Replay an exported Neko-TOP ROM without starting Neko or DNS.

The exporter writes the exact Petrov BDF1 objects used by ``rom_model``:

    (H-dt*Li) a[n+1] = H a[n] + dt*(b + Le@a[n] + q(a[n]))

The optional ``--dns-coefficients`` input is a CSV with ``step,time,a_1,...``.
The optional ``--dns-fields`` input is a NumPy ``.npz`` containing either
``fields`` with shape ``(time, field, dof)`` or ``velocity`` with shape
``(time, 3, dof)`` and, optionally, ``scalar``.  Those fields are projected
onto the exported final trial basis using the exported trial mass matrix.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def read_matrix(path: Path) -> np.ndarray:
    value = np.loadtxt(path, delimiter=",")
    return np.atleast_2d(value)


def read_vector(path: Path) -> np.ndarray:
    return np.atleast_1d(np.loadtxt(path, delimiter=","))


def read_quadratic(path: Path, rank: int) -> np.ndarray:
    rows = np.genfromtxt(path, delimiter=",", names=True)
    rows = np.atleast_1d(rows)
    q = np.zeros((rank, rank, rank), dtype=float)
    for row in rows:
        q[int(row["i"]) - 1, int(row["j"]) - 1, int(row["k"]) - 1] = row["value"]
    return q


def read_coefficient_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    data = np.atleast_1d(data)
    names = data.dtype.names or ()
    modes = [name for name in names if name.startswith("a_")]
    if not modes:
        modes = [name for name in names if name.startswith("lambda_")]
    modes.sort(key=lambda name: int(name.split("_")[-1]))
    return np.asarray(data["time"], dtype=float), np.column_stack(
        [np.asarray(data[name], dtype=float) for name in modes]
    )


def read_export(directory: Path) -> tuple[dict, dict[str, np.ndarray]]:
    with (directory / "metadata.json").open() as stream:
        metadata = json.load(stream)
    files = metadata["files"]
    rank = int(metadata["retained_rank"])
    operators = {
        "H": read_matrix(directory / files["H"]),
        "trial_mass": read_matrix(directory / files["trial_mass"]),
        "b": read_vector(directory / files["b"]),
        "Le": read_matrix(directory / files["linear_explicit"]),
        "Li": read_matrix(directory / files["linear_implicit"]),
        "Q": read_quadratic(directory / files["Q"], rank),
        "a0": read_vector(directory / files["a0"]),
    }
    rows = int(metadata["trial_basis_rows"])
    trial = np.fromfile(directory / files["trial_basis"], dtype=np.float64)
    operators["trial_basis"] = trial.reshape((rows, rank), order="F")
    weights = np.fromfile(directory / files["weights"], dtype=np.float64)
    operators["weights"] = weights
    lift = np.fromfile(directory / files["lift_velocity"], dtype=np.float64)
    operators["lift_velocity"] = lift.reshape((3, weights.size))
    if "scalar_lift" in files:
        operators["scalar_lift"] = np.fromfile(
            directory / files["scalar_lift"], dtype=np.float64
        )
    return metadata, operators


def nonlinear_residual(q: np.ndarray, state: np.ndarray) -> np.ndarray:
    return np.einsum("j,k,ijk->i", state, state, q, optimize=True)


def residual_parts(operators: dict[str, np.ndarray], state: np.ndarray):
    constant = operators["b"]
    linear = operators["Le"] @ state
    nonlinear = nonlinear_residual(operators["Q"], state)
    implicit = operators["Li"] @ state
    return constant, linear, nonlinear, implicit


def replay(metadata: dict, operators: dict[str, np.ndarray]):
    dt = float(metadata["dt"])
    n_steps = int(metadata["n_steps"])
    h = operators["H"]
    li = operators["Li"]
    step_matrix = h - dt * li
    states = np.full((n_steps + 1, operators["a0"].size), np.nan)
    states[0] = operators["a0"]
    residual_norms = np.full((n_steps + 1, 5), np.nan)

    for step in range(n_steps + 1):
        if not np.all(np.isfinite(states[step])):
            break
        parts = residual_parts(operators, states[step])
        residual_norms[step, :4] = [np.linalg.norm(part) for part in parts]
        residual_norms[step, 4] = np.linalg.norm(parts[0] + parts[1] + parts[2])
        if step == n_steps:
            break
        rhs = h @ states[step] + dt * (parts[0] + parts[1] + parts[2])
        try:
            states[step + 1] = np.linalg.solve(step_matrix, rhs)
        except np.linalg.LinAlgError:
            break
    times = float(metadata["time_start"]) + dt * np.arange(n_steps + 1)
    return times, states, residual_norms, step_matrix


def eigenvalues(metadata: dict, operators: dict[str, np.ndarray]):
    a0 = operators["a0"]
    q = operators["Q"]
    q_jacobian = np.einsum(
        "ijk,k->ij", q + np.swapaxes(q, 1, 2), a0, optimize=True
    )
    j_explicit = operators["Le"] + q_jacobian
    continuous = np.linalg.solve(
        operators["H"], j_explicit + operators["Li"]
    )
    dt = float(metadata["dt"])
    discrete = np.linalg.solve(
        operators["H"] - dt * operators["Li"],
        operators["H"] + dt * j_explicit,
    )
    return continuous, discrete


def interpolate_states(source_time: np.ndarray, source_state: np.ndarray, time: np.ndarray):
    result = np.empty((time.size, source_state.shape[1]))
    for mode in range(source_state.shape[1]):
        result[:, mode] = np.interp(time, source_time, source_state[:, mode])
    return result


def project_dns_fields(path: Path, metadata: dict, operators: dict[str, np.ndarray]):
    data = np.load(path)
    if "fields" in data:
        fields = np.asarray(data["fields"], dtype=float)
    else:
        velocity = np.asarray(data["velocity"], dtype=float)
        if "scalar" in data:
            fields = np.concatenate((velocity, np.asarray(data["scalar"])[:, None, :]), axis=1)
        else:
            fields = velocity
    if fields.ndim != 3:
        raise ValueError("DNS fields must have shape (time, field, dof)")
    rank = int(metadata["retained_rank"])
    n = operators["weights"].size
    n_fields = fields.shape[1]
    if fields.shape[2] != n or fields.shape[0] < 1:
        raise ValueError("DNS field dimensions do not match exported basis")
    if operators["trial_basis"].shape[0] != n_fields * n:
        raise ValueError("DNS field count does not match exported coupled basis")
    lift = np.zeros_like(fields)
    lift[:, :3, :] = operators["lift_velocity"][None, :, :]
    if n_fields == 4:
        lift[:, 3, :] = operators["scalar_lift"][None, :]
    fluctuation = fields - lift
    metric = np.tile(operators["weights"], n_fields)
    if n_fields == 4:
        metric[3 * n :] *= float(metadata["scalar_metric_weight"])
    rhs = np.einsum(
        "nr,tn->tr", operators["trial_basis"], fluctuation.reshape(fields.shape[0], -1) * metric
    )
    coefficients = np.linalg.solve(operators["trial_mass"], rhs.T).T
    if "time" in data:
        time = np.asarray(data["time"], dtype=float)
    else:
        time = float(metadata["time_start"]) + float(metadata["dt"]) * np.arange(fields.shape[0])
    return time, coefficients


def write_trajectory(path: Path, time: np.ndarray, states: np.ndarray, residuals: np.ndarray):
    columns = [time, states]
    header = ["time"] + [f"a_{i + 1}" for i in range(states.shape[1])]
    columns.append(residuals)
    header += ["constant_norm", "linear_norm", "nonlinear_norm", "implicit_norm", "explicit_total_norm"]
    np.savetxt(path, np.column_stack(columns), delimiter=",", header=",".join(header), comments="")


def make_plot(path: Path, time: np.ndarray, states: np.ndarray, residuals: np.ndarray, eigen_data):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].plot(time, states)
    axes[0, 0].set_title("Reduced coefficients")
    axes[0, 0].set_xlabel("time")
    axes[0, 1].semilogy(time, np.maximum(residuals[:, :4], 1.0e-300))
    axes[0, 1].set_title("Residual components")
    axes[0, 1].legend(("constant", "linear", "nonlinear", "implicit"))
    axes[1, 0].semilogy(time, np.maximum(residuals[:, 4], 1.0e-300), label="explicit total")
    axes[1, 0].set_title("Explicit residual norm")
    axes[1, 0].legend()
    continuous, discrete = eigen_data
    axes[1, 1].scatter(continuous.real, continuous.imag, label="continuous")
    axes[1, 1].scatter(discrete.real, discrete.imag, marker="x", label="BDF1 map")
    axes[1, 1].axhline(0.0, color="0.7", linewidth=0.5)
    axes[1, 1].axvline(0.0, color="0.7", linewidth=0.5)
    axes[1, 1].set_title("Linearisation at a0")
    axes[1, 1].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_directory", type=Path)
    parser.add_argument("--actual-trajectory", type=Path, default=None)
    parser.add_argument("--dns-coefficients", type=Path, default=None)
    parser.add_argument("--dns-fields", type=Path, default=None)
    parser.add_argument("--output-directory", type=Path, default=None)
    args = parser.parse_args()

    metadata, operators = read_export(args.export_directory)
    output = args.output_directory or args.export_directory
    output.mkdir(parents=True, exist_ok=True)
    time, states, residuals, _ = replay(metadata, operators)
    continuous, discrete = eigenvalues(metadata, operators)
    write_trajectory(output / "replayed_trajectory.csv", time, states, residuals)
    np.savetxt(
        output / "linearised_eigenvalues.csv",
        np.column_stack((continuous.real, continuous.imag, discrete.real, discrete.imag)),
        delimiter=",",
        header="continuous_real,continuous_imag,bdf1_real,bdf1_imag",
        comments="",
    )

    summary = {
        "rank": int(metadata["retained_rank"]),
        "dt": float(metadata["dt"]),
        "n_steps": int(metadata["n_steps"]),
        "finite_replayed_steps": int(np.sum(np.all(np.isfinite(states), axis=1))),
        "max_replayed_coefficient": float(np.nanmax(np.abs(states))),
        "continuous_eigenvalue_max_real": float(np.max(continuous.real)),
        "bdf1_eigenvalue_max_abs": float(np.max(np.abs(discrete))),
    }

    actual_path = args.actual_trajectory
    if actual_path is None:
        candidate = args.export_directory / "reference_trajectory.csv"
        if candidate.exists():
            actual_path = candidate
    if actual_path is not None and actual_path.exists():
        actual_time, actual_states = read_coefficient_csv(actual_path)
        replay_at_actual = interpolate_states(time, states, actual_time)
        difference = replay_at_actual - actual_states
        summary["actual_trajectory_max_abs_error"] = float(np.nanmax(np.abs(difference)))
        summary["actual_trajectory_rms_error"] = float(
            np.sqrt(np.nanmean(difference * difference))
        )

    dns_path = args.dns_coefficients
    if dns_path is not None:
        dns_time, dns_states = read_coefficient_csv(dns_path)
    elif args.dns_fields is not None:
        dns_time, dns_states = project_dns_fields(args.dns_fields, metadata, operators)
    else:
        dns_time = dns_states = None
    if dns_time is not None:
        rom_at_dns = interpolate_states(time, states, dns_time)
        coefficient_error = rom_at_dns - dns_states
        dns_derivative = np.gradient(dns_states, dns_time, axis=0)
        predicted_derivative = np.vstack(
            [
                np.linalg.solve(
                    operators["H"],
                    sum(residual_parts(operators, state)),
                )
                for state in dns_states
            ]
        )
        summary["dns_projection_max_abs_coefficient_error"] = float(
            np.nanmax(np.abs(coefficient_error))
        )
        summary["dns_projection_rms_coefficient_error"] = float(
            np.sqrt(np.nanmean(coefficient_error * coefficient_error))
        )
        summary["dns_derivative_rms_error"] = float(
            np.sqrt(np.nanmean((predicted_derivative - dns_derivative) ** 2))
        )
        np.savetxt(
            output / "dns_projection_comparison.csv",
            np.column_stack((dns_time, dns_states, rom_at_dns, predicted_derivative, dns_derivative)),
            delimiter=",",
            header="time,dns_a,rom_a,predicted_derivative,dns_derivative",
            comments="",
        )

    summary["continuous_eigenvalues"] = [complex(value).__repr__() for value in continuous]
    summary["bdf1_eigenvalues"] = [complex(value).__repr__() for value in discrete]
    with (output / "diagnostic_summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2)
    plot_written = make_plot(output / "reduced_system_diagnostic.png", time, states, residuals, (continuous, discrete))

    print(json.dumps(summary, indent=2))
    if not plot_written:
        print("matplotlib is unavailable; CSV and eigenvalue outputs were still written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
