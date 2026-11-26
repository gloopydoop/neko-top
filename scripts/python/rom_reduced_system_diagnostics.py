#!/usr/bin/env python3
"""Integrate and diagnose a Neko exported reduced BDF1 system.

The tool reads only the small reduced CSV/JSON files.  It can reproduce the
exported BDF1 trajectory, refine its step, or compare it with SciPy's implicit
Radau/BDF solver.  It does not modify the ROM or claim that a reference
trajectory is an independent full-order validation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
from typing import Any

import numpy as np


def _load_matrix(path: Path, shape: tuple[int, int] | None = None) -> np.ndarray:
    array = np.loadtxt(path, delimiter=",")
    array = np.atleast_2d(np.asarray(array, dtype=np.float64))
    if shape is not None and array.shape != shape:
        raise ValueError(f"{path}: expected shape {shape}, found {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{path}: non-finite matrix data")
    return array


def load_system(directory: Path) -> dict[str, Any]:
    directory = directory.resolve()
    metadata_path = directory / "metadata.json"
    metadata_text = metadata_path.read_text(encoding="utf-8")
    metadata_warnings: list[str] = []
    try:
        metadata = json.loads(metadata_text)
    except json.JSONDecodeError:
        # Some Fortran ESw.dEe formatting of the largest finite value emits
        # `1.797...+308` instead of a JSON exponent `1.797...E+308`. Repair
        # only that unambiguous missing-exponent-marker form and report it.
        repaired, count = re.subn(
            r"(?<![eE])(?<=\d)([+-]\d{3})(?=\s*[,}\]])", r"E\1", metadata_text
        )
        if count == 0:
            raise
        metadata = json.loads(repaired)
        metadata_warnings.append(
            f"repaired {count} Fortran JSON number(s) missing the E before a three-digit exponent"
        )
    n = int(metadata["n_state"])
    if n < 1:
        raise ValueError("n_state must be positive")
    required = ("H", "b", "linear_explicit", "linear_implicit", "Q", "a0")
    files = metadata.get("files", {})
    missing = [key for key in required if key not in files or not (directory / files[key]).is_file()]
    if missing:
        raise FileNotFoundError(f"reduced system is missing required exports: {', '.join(missing)}")

    H = _load_matrix(directory / files["H"], (n, n))
    b = np.loadtxt(directory / files["b"], delimiter=",")
    b = np.asarray(b, dtype=np.float64).reshape(-1)
    Le = _load_matrix(directory / files["linear_explicit"], (n, n))
    Li = _load_matrix(directory / files["linear_implicit"], (n, n))
    a0 = np.loadtxt(directory / files["a0"], delimiter=",")
    a0 = np.asarray(a0, dtype=np.float64).reshape(-1)
    if b.size != n or a0.size != n:
        raise ValueError(f"vector size mismatch: n_state={n}, b={b.size}, a0={a0.size}")
    if np.max(np.abs(H - H.T)) > 1e-9 * max(1.0, np.max(np.abs(H))):
        raise ValueError("H is not symmetric to the reduced-system diagnostic tolerance")
    q_path = directory / files["Q"]
    Q = np.zeros((n, n, n), dtype=np.float64)
    with q_path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        expected = {"i", "j", "k", "value"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("Q.csv must have i,j,k,value columns")
        for row in reader:
            i, j, k = int(row["i"]) - 1, int(row["j"]) - 1, int(row["k"]) - 1
            if not (0 <= i < n and 0 <= j < n and 0 <= k < n):
                raise ValueError("Q.csv contains an index outside n_state")
            Q[i, j, k] += float(row["value"])
    if not all(np.all(np.isfinite(x)) for x in (H, b, Le, Li, Q, a0)):
        raise ValueError("reduced system contains non-finite data")
    return {
        "directory": directory,
        "metadata": metadata,
        "metadata_warnings": metadata_warnings,
        "H": H,
        "b": b,
        "Le": Le,
        "Li": Li,
        "Q": Q,
        "a0": a0,
    }


def _q_value(Q: np.ndarray, a: np.ndarray) -> np.ndarray:
    return np.einsum("ijk,j,k->i", Q, a, a, optimize=True)


def _q_jacobian(Q: np.ndarray, a: np.ndarray) -> np.ndarray:
    # d/da_k (a_j Q_ijk a_k), retaining the actual tensor without symmetry assumptions.
    return np.einsum("ilk,l->ik", Q, a, optimize=True) + np.einsum("ijl,j->il", Q, a, optimize=True)


def _symmetric_mass_summary(H: np.ndarray) -> dict[str, Any]:
    """Report spectrum and conditioning of the symmetric part of reduced H."""
    symmetric = 0.5 * (H + H.T)
    eigenvalues = np.linalg.eigvalsh(symmetric)
    spectral_scale = float(np.max(np.abs(eigenvalues)))
    pd_tolerance = float(
        np.finfo(np.float64).eps * max(1, H.shape[0]) * spectral_scale
    )
    smallest_abs = float(np.min(np.abs(eigenvalues)))
    condition = (
        float(spectral_scale / smallest_abs)
        if smallest_abs > pd_tolerance
        else float("inf")
    )
    return {
        "symmetric_eigenvalues": [float(value) for value in eigenvalues],
        "positive_definite": bool(eigenvalues[0] > pd_tolerance),
        "positive_definite_tolerance": pd_tolerance,
        "spectral_condition_number_2": condition,
    }


def integrate_bdf1(system: dict[str, Any], refinement: int = 1) -> tuple[np.ndarray, np.ndarray]:
    if refinement < 1:
        raise ValueError("refinement must be >= 1")
    meta = system["metadata"]
    base_dt = float(meta["dt"])
    n_steps = int(meta["n_steps"])
    dt = base_dt / refinement
    H, b, Le, Li, Q = (system[key] for key in ("H", "b", "Le", "Li", "Q"))
    step_matrix = H - dt * Li
    lu = None
    try:
        from scipy.linalg import lu_factor, lu_solve
        lu = lu_factor(step_matrix)
        solve = lambda rhs: lu_solve(lu, rhs)
    except ImportError:
        solve = lambda rhs: np.linalg.solve(step_matrix, rhs)
    states = np.empty((n_steps * refinement + 1, system["a0"].size), dtype=np.float64)
    states[0] = system["a0"]
    for k in range(1, states.shape[0]):
        a = states[k - 1]
        rhs = H @ a + dt * (b + Le @ a + _q_value(Q, a))
        states[k] = solve(rhs)
        if not np.all(np.isfinite(states[k])):
            raise FloatingPointError(f"non-finite reduced state at step {k}")
    times = float(meta["time_start"]) + np.arange(states.shape[0]) * dt
    return times, states


def integrate_scipy(system: dict[str, Any], method: str, rtol: float, atol: float):
    from scipy.integrate import solve_ivp

    meta = system["metadata"]
    n_steps = int(meta["n_steps"])
    dt = float(meta["dt"])
    t0 = float(meta["time_start"])
    t1 = float(meta["time_end"])
    H, b, Le, Li, Q = (system[key] for key in ("H", "b", "Le", "Li", "Q"))
    linear = Le + Li

    def rhs(_time, a):
        return np.linalg.solve(H, b + linear @ a + _q_value(Q, a))

    times = t0 + np.arange(n_steps + 1) * dt
    result = solve_ivp(rhs, (t0, t1), system["a0"], method=method,
                       t_eval=times, rtol=rtol, atol=atol)
    if not result.success:
        raise RuntimeError(f"SciPy {method} integration failed: {result.message}")
    return result.t, result.y.T


def _load_reference(system: dict[str, Any]) -> tuple[np.ndarray, np.ndarray] | None:
    file_name = system["metadata"].get("files", {}).get("reference_trajectory")
    if not file_name:
        return None
    path = system["directory"] / file_name
    if not path.is_file():
        return None
    table = np.genfromtxt(path, delimiter=",", names=True, dtype=np.float64, encoding="utf-8")
    names = table.dtype.names or ()
    coefficient_names = [name for name in names if name.startswith("a_")]
    if "time" not in names or not coefficient_names:
        raise ValueError("reference trajectory must contain time and a_i columns")
    return np.atleast_1d(table["time"]), np.column_stack([np.atleast_1d(table[name]) for name in coefficient_names])


def projected_fom_coordinate_availability(system: dict[str, Any]) -> dict[str, Any]:
    """Report whether exact projected-FOM coordinates accompany this system.

    The reduced ``reference_trajectory`` is produced by integrating the ROM;
    it is never treated as projected FOM data.  A future producer can declare
    actual coordinates with a ``projected_fom_coordinates`` metadata object,
    but it must bind them to the same candidate, basis, mesh/dof map, coupled
    metric, and time grid before the reader calls them comparable.
    """
    metadata = system["metadata"]
    files = metadata.get("files", {})
    outer_directory = system["directory"].parent
    basis_files = {
        key: files.get(key)
        for key in ("trial_basis", "lift_velocity", "scalar_lift", "weights")
    }
    available_basis_files = {
        key: str(system["directory"] / name)
        for key, name in basis_files.items()
        if name and (system["directory"] / name).is_file()
    }
    reference_name = files.get("reference_trajectory")
    reference_path = system["directory"] / reference_name if reference_name else None
    integrated_rom_reference = {
        "file": str(reference_path) if reference_path else None,
        "exists": bool(reference_path and reference_path.is_file()),
        "classification": "integrated_rom_coefficients",
    }
    descriptor = metadata.get("projected_fom_coordinates")
    required_metadata = (
        "candidate_design_id",
        "basis_coordinate_id",
        "mesh_dof_map_id",
    )
    missing = [
        f"reduced-system metadata field {name}"
        for name in required_metadata
        if not metadata.get(name)
    ]
    required_basis_files = ("trial_basis", "lift_velocity", "weights")
    missing.extend(
        f"basis coordinate export {key}"
        for key in required_basis_files
        if key not in available_basis_files
    )

    if not isinstance(descriptor, dict):
        missing.append("producer-declared projected_fom_coordinates metadata")
        missing.append("identified physical FOM snapshots or their projected coefficients")
        return {
            "status": "unavailable",
            "available": False,
            "reason": "no exact projected-FOM coordinate stream is declared",
            "missing": missing,
            "available_basis_files": available_basis_files,
            "integrated_rom_reference": integrated_rom_reference,
            "forward_field_stream_note": (
                "The outer archive's forward_fields stream is the ROM reconstruction; "
                "it is not identified as a physical FOM trajectory."
            ),
            "smallest_valid_producer_insertion_point": (
                "At each physical FOM snapshot output, while q_FOM, the candidate ROM "
                "basis V, lift q_L, and coupled mass metric are resident, compute "
                "H=V^T M V and solve H a_proj=V^T M(q_FOM-q_L)."
            ),
            "required_provenance": [
                "candidate design and output sample identifiers",
                "physical time for each snapshot",
                "basis coordinate identifier matching this reduced system",
                "mesh/dof-map identifier and velocity/scalar ordering",
                "coupled metric definition and scalar metric weight",
                "projected coefficients a_1..a_n with one row per FOM snapshot",
            ],
        }

    # A declaration is not enough by itself: require the coordinate identity
    # to match this reduced system and check the exported file/time grid.
    declared_file = descriptor.get("file")
    coordinate_path = outer_directory / declared_file if declared_file else None
    checks = {
        "candidate_design_id": descriptor.get("candidate_design_id")
        == metadata.get("candidate_design_id"),
        "basis_coordinate_id": descriptor.get("basis_coordinate_id")
        == metadata.get("basis_coordinate_id"),
        "mesh_dof_map_id": descriptor.get("mesh_dof_map_id")
        == metadata.get("mesh_dof_map_id"),
        "metric": descriptor.get("metric") == "coupled_H",
        "state_definition": descriptor.get("state_definition")
        == "solve(H, V^T M (q_FOM - q_L))",
        "coordinate_file": bool(coordinate_path and coordinate_path.is_file()),
    }
    coordinate_rows = 0
    projected_times: list[float] = []
    if checks["coordinate_file"]:
        with coordinate_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            names = set(reader.fieldnames or ())
            required_columns = {"sample_id", "time"} | {
                f"a_{i + 1}" for i in range(int(metadata["n_state"]))
            }
            checks["coordinate_columns"] = required_columns.issubset(names)
            for row in reader:
                coordinate_rows += 1
                try:
                    projected_times.append(float(row["time"]))
                    coefficients = [float(row[f"a_{i + 1}"]) for i in range(int(metadata["n_state"]))]
                except (KeyError, TypeError, ValueError):
                    checks["coordinate_columns"] = False
                    continue
                if not np.isfinite(projected_times[-1]) or not np.all(np.isfinite(coefficients)):
                    checks["coordinate_columns"] = False
            checks["coordinate_rows"] = coordinate_rows > 0
            checks["time_order"] = bool(
                projected_times
                and all(later > earlier for earlier, later in zip(projected_times, projected_times[1:]))
            )
            reference = _load_reference(system)
            if reference is None:
                checks["integrated_rom_reference"] = False
            else:
                reference_times, _ = reference
                checks["time_grid_match"] = all(
                    np.any(np.isclose(reference_times, value, rtol=1e-12, atol=1e-14))
                    for value in projected_times
                )
    else:
        checks["coordinate_columns"] = False
        checks["coordinate_rows"] = False
        checks["time_order"] = False
        checks["time_grid_match"] = False
    checks.setdefault("time_order", False)
    checks.setdefault("time_grid_match", False)
    checks.setdefault("integrated_rom_reference", integrated_rom_reference["exists"])
    failed_checks = [name for name, passed in checks.items() if not passed]
    return {
        "status": "coordinates_available" if not failed_checks else "unavailable",
        "available": not failed_checks,
        "comparison_performed": False,
        "reason": None if not failed_checks else "projected-FOM provenance contract is incomplete or mismatched",
        "failed_checks": failed_checks,
        "checks": checks,
        "file": str(coordinate_path) if coordinate_path else None,
        "coordinate_rows": coordinate_rows,
        "available_basis_files": available_basis_files,
        "integrated_rom_reference": integrated_rom_reference,
    }


def _trajectory_summary(system, times, states, reference):
    H, b, Le, Li, Q = (system[key] for key in ("H", "b", "Le", "Li", "Q"))
    linear = Le + Li
    energy = 0.5 * np.einsum("ti,ij,tj->t", states, H, states)
    energy_rate = np.empty(states.shape[0])
    nonlinear_production = np.empty(states.shape[0])
    for index, a in enumerate(states):
        q = _q_value(Q, a)
        energy_rate[index] = a @ (b + linear @ a + q)
        nonlinear_production[index] = a @ q
    report = {
        "steps": int(states.shape[0] - 1),
        "time_start": float(times[0]),
        "time_end": float(times[-1]),
        "max_abs_coefficient": float(np.max(np.abs(states))),
        "final_coefficient_norm": float(np.linalg.norm(states[-1])),
        "energy_initial": float(energy[0]),
        "energy_final": float(energy[-1]),
        "energy_min": float(np.min(energy)),
        "energy_max": float(np.max(energy)),
        "energy_rate_min": float(np.min(energy_rate)),
        "energy_rate_max": float(np.max(energy_rate)),
        "energy_rate_semantics": (
            "Continuous ODE derivative d(0.5*a^T H a)/dt evaluated at the states "
            "sampled by this BDF1 replay; it is not a discrete BDF1 energy increment."
        ),
        "nonlinear_energy_production_min": float(np.min(nonlinear_production)),
        "nonlinear_energy_production_max": float(np.max(nonlinear_production)),
    }
    if reference is not None:
        ref_times, ref_states = reference
        if ref_states.shape == states.shape and np.allclose(ref_times, times, rtol=1e-12, atol=1e-14):
            error = states - ref_states
            report["reference_trajectory_comparison"] = {
                "max_abs_error": float(np.max(np.abs(error))),
                "relative_l2_error": float(np.linalg.norm(error) / max(np.linalg.norm(ref_states), np.finfo(float).tiny)),
            }
        else:
            report["reference_trajectory_comparison"] = {
                "available": True,
                "comparable": False,
                "reason": "time grid or coefficient shape differs",
                "reference_steps": int(ref_states.shape[0] - 1),
            }
    else:
        report["reference_trajectory_comparison"] = None
    eigenvalues = np.linalg.eigvals(np.linalg.solve(H, linear))
    report["linearized_eigenvalues"] = [
        {"real": float(value.real), "imag": float(value.imag)}
        for value in sorted(eigenvalues, key=lambda z: (z.real, z.imag), reverse=True)
    ]
    report["max_linear_growth_rate"] = float(np.max(eigenvalues.real))
    report["reduced_mass_matrix"] = _symmetric_mass_summary(H)
    if states.shape[0] > 1:
        growth = []
        for a in states:
            jac = linear + _q_jacobian(Q, a)
            growth.append(np.max(np.linalg.eigvals(np.linalg.solve(H, jac)).real))
        sampled_abscissa = float(np.max(growth))
        report["max_sampled_frozen_jacobian_spectral_abscissa"] = sampled_abscissa
        # Preserve the legacy JSON key while stating its exact sampled meaning.
        report["trajectory_max_local_growth_rate"] = sampled_abscissa
        report["trajectory_max_local_growth_rate_semantics"] = (
            "Legacy alias for max_sampled_frozen_jacobian_spectral_abscissa; "
            "maximum real eigenvalue of the frozen reduced Jacobian sampled "
            "along this replay, not a general or transient growth bound."
        )
    report["transient_growth"] = {
        "measured": False,
        "reason": "No propagator-norm or non-normal transient-growth analysis is performed.",
    }
    return report


def diagnose(directory: Path, method: str = "bdf1", refinement: int = 1,
             rtol: float = 1e-8, atol: float = 1e-10) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    system = load_system(directory)
    if method == "bdf1":
        times, states = integrate_bdf1(system, refinement)
    else:
        if refinement != 1:
            raise ValueError("refinement applies to BDF1 only")
        times, states = integrate_scipy(system, method, rtol, atol)
    report = {
        "directory": str(system["directory"]),
        "method": method,
        "refinement": refinement,
        "basis_pair": [system["metadata"].get("basis_forward"), system["metadata"].get("basis_adjoint")],
        "candidate_rank": int(system["metadata"].get("trial_basis_rank", states.shape[1])),
        "metadata_warnings": system["metadata_warnings"],
        "projected_fom_coordinates": projected_fom_coordinate_availability(system),
        "report": _trajectory_summary(system, times, states, _load_reference(system)),
        "interpretation_note": (
            "Reduced ODE diagnostics only. Nonlinear energy production includes open-boundary and SEM effects; "
            "a positive value is not by itself proof of a divergence defect."
        ),
    }
    return report, times, states


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory containing exported reduced_system/metadata.json")
    parser.add_argument("--method", choices=("bdf1", "Radau", "BDF"), default="bdf1")
    parser.add_argument("--refinement", type=int, default=1, help="BDF1 steps per exported time step")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--output-dir", type=Path, help="optional directory for summary and coefficient trajectory")
    args = parser.parse_args()
    report, times, states = diagnose(args.directory, args.method, args.refinement, args.rtol, args.atol)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "rom_reduced_system_summary.json").write_text(rendered, encoding="utf-8")
        with (args.output_dir / "rom_reduced_system_trajectory.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["step", "time", *[f"a_{i + 1}" for i in range(states.shape[1])]])
            writer.writerows((i, time, *state) for i, (time, state) in enumerate(zip(times, states)))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
