#!/usr/bin/env python3

"""Configuration and numerical helpers for the separate DMD recovery path."""

from dataclasses import dataclass
import os
import re

import numpy as np


def rotate_existing_file(path: str) -> str | None:
    """Archive an existing DMD output using the POD naming convention."""
    if not os.path.exists(path):
        return None

    stem, ext = os.path.splitext(path)
    index = 1
    while True:
        candidate = f"{stem}_{index:03d}{ext}"
        if not os.path.exists(candidate):
            os.rename(path, candidate)
            return candidate
        index += 1



def as_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "t")
    return bool(value)


def strip_json_comments(text: str) -> str:
    out = []
    index = 0
    in_string = False
    escaped = False
    while index < len(text):
        char = text[index]
        if in_string:
            out.append(char)
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            index += 1
            continue
        if char == '"':
            in_string = True
            out.append(char)
            index += 1
            continue
        if char == "/" and index + 1 < len(text):
            if text[index + 1] == "/":
                index += 2
                while index < len(text) and text[index] != "\n":
                    index += 1
                continue
            if text[index + 1] == "*":
                index += 2
                while index + 1 < len(text):
                    if text[index] == "*" and text[index + 1] == "/":
                        index += 2
                        break
                    index += 1
                continue
        out.append(char)
        index += 1
    return "".join(out)


def load_case_config(case_path: str) -> dict:
    with open(case_path, "r", encoding="utf-8") as handle:
        text = strip_json_comments(handle.read())
    text = re.sub(r",\s*([}\]])", r"\1", text)
    import json
    return json.loads(text)


def count_enabled_scalars(case_cfg: dict) -> int:
    case_data = case_cfg.get("case", {})
    if not isinstance(case_data, dict):
        return 0
    scalar_cfg = case_data.get("scalar")
    if isinstance(scalar_cfg, dict):
        return int(as_bool(scalar_cfg.get("enabled"), True))
    scalar_cfgs = case_data.get("scalars")
    if not isinstance(scalar_cfgs, list):
        return 0
    return sum(
        1
        for entry in scalar_cfgs
        if not isinstance(entry, dict)
        or as_bool(entry.get("enabled"), True)
    )


@dataclass(frozen=True)
class DMDConfig:
    batch_size: int
    keep_modes: int
    rank: int
    i_stream: int
    dtype: type
    write_modes: bool
    mode_output_dtype: type
    mode_output_format: str
    mode_output_file_name: str
    include_scalar: bool
    debug: bool
    snapshot_dt: float
    amplitude_regularization: float
    conjugate_tolerance: float

    @property
    def n_fields(self) -> int:
        return 4 if self.include_scalar else 3

    @property
    def n_input_fields(self) -> int:
        return self.n_fields

    @property
    def closure_enabled(self) -> bool:
        return False

    @property
    def scalar_closure_enabled(self) -> bool:
        return False


def set_debug(enabled: bool) -> None:
    # Kept as a separate hook so the DMD driver does not import the standard
    # POD helper module, which also imports optional visualization backends.
    return None


def recv_field(ds, dtype: type) -> np.ndarray:
    from pysemtools.io.utils import get_fld_from_ndarray

    field = get_fld_from_ndarray(ds.recieve(), ds.lx, ds.ly, ds.lz, ds.nelv)
    return field.astype(dtype)


def recv_fields(ds, n_fields: int, dtype: type) -> list[np.ndarray]:
    return [recv_field(ds, dtype) for _ in range(n_fields)]


def make_pod(comm, bm: np.ndarray, cfg: DMDConfig):
    from pysemtools.rom.io_help import IoHelp
    from pysemtools.rom.pod import POD

    pod = POD(
        comm,
        number_of_modes_to_update=cfg.keep_modes,
        global_updates=True,
        auto_expand=False,
        bckend="numpy",
    )
    ioh = IoHelp(
        comm,
        number_of_fields=cfg.n_fields,
        batch_size=cfg.batch_size,
        field_size=bm.size,
        field_data_type=cfg.dtype,
        mass_matrix_data_type=cfg.dtype,
    )
    mass_fields = [np.copy(np.sqrt(bm)) for _ in range(cfg.n_fields)]
    ioh.copy_fieldlist_to_xi(mass_fields)
    ioh.bm1sqrt[:, :] = np.copy(ioh.xi[:, :])
    return pod, ioh


def write_modes_to_disk(comm, mesh, modes: np.ndarray, cfg: DMDConfig, output_index: int) -> int:
    if not cfg.write_modes:
        return output_index
    if cfg.mode_output_format not in ("fld", "nek5000"):
        raise ValueError(
            "DMD mode output currently supports only output_format=fld"
        )

    from pysemtools.io.wrappers import write_data

    file_name = os.path.normpath(cfg.mode_output_file_name)
    out_dir = os.path.dirname(file_name)
    stem = os.path.splitext(os.path.basename(file_name))[0]
    if comm.Get_rank() == 0 and out_dir:
        os.makedirs(out_dir, exist_ok=True)
    comm.Barrier()

    field_size = modes.shape[0] // cfg.n_fields
    names = []
    for mode_index in range(cfg.keep_modes):
        suffix = f"_mode{mode_index + 1:04d}"
        names.extend([f"u{suffix}", f"v{suffix}", f"w{suffix}"])
        if cfg.include_scalar:
            names.append(f"s{suffix}")
    data = {}
    for mode_index in range(cfg.keep_modes):
        for field_index in range(cfg.n_fields):
            start = field_index * field_size
            end = (field_index + 1) * field_size
            name = names[mode_index * cfg.n_fields + field_index]
            data[name] = modes[start:end, mode_index].astype(
                cfg.mode_output_dtype, copy=False
            )

    sample_base = os.path.join(out_dir, f"{stem}0")
    write_data(
        comm,
        fname=f"{sample_base}.f{output_index:05d}",
        data_dict=data,
        parallel_io=False,
        dtype=cfg.mode_output_dtype,
        msh=[mesh.x, mesh.y, mesh.z],
        write_mesh=output_index == 0,
    )
    comm.Barrier()
    if comm.Get_rank() == 0:
        with open(
            os.path.join(out_dir, f"{stem}0.nek5000"),
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write(f"filetemplate:         {stem}%01d.f%05d\n")
            handle.write("firsttimestep:     0\n")
            handle.write(f"numtimesteps:      {output_index + 1}\n")
    comm.Barrier()
    return output_index + 1


def load_dmd_config(case_path: str) -> tuple[dict, DMDConfig]:
    case = load_case_config(case_path)
    state_recovery = case.get("state_recovery", {})
    if not isinstance(state_recovery, dict):
        raise KeyError("case['state_recovery'] must be a dict")

    enabled_scalars = count_enabled_scalars(case)
    if enabled_scalars > 1:
        raise ValueError(
            "DMD state recovery currently supports exactly one enabled scalar."
        )
    include_scalar = enabled_scalars == 1

    dtype_name = str(state_recovery.get("dtype", "double")).strip().lower()
    if dtype_name == "single":
        dtype = np.float32
    elif dtype_name == "double":
        dtype = np.float64
    else:
        raise ValueError(f"Unsupported dtype '{dtype_name}'")

    output_precision = str(
        state_recovery.get("output_precision", "sp")
    ).strip().lower()
    if output_precision in ("sp", "single"):
        mode_output_dtype = np.float32
    elif output_precision in ("dp", "double"):
        mode_output_dtype = np.float64
    else:
        raise ValueError(f"Unsupported output_precision '{output_precision}'")

    n_modes = int(state_recovery["n_modes"])
    if n_modes < 1:
        raise ValueError("state_recovery.n_modes must be positive")
    rank = int(state_recovery.get("dmd_rank", n_modes))
    if rank < 1:
        raise ValueError("state_recovery.dmd_rank must be positive")
    if rank > n_modes:
        raise ValueError(
            "state_recovery.dmd_rank cannot exceed state_recovery.n_modes"
        )

    timestep = float(case["case"]["time"]["timestep"])
    cfg = DMDConfig(
        batch_size=int(state_recovery["batch_size"]),
        keep_modes=n_modes,
        rank=rank,
        i_stream=int(state_recovery["i_stream"]),
        dtype=dtype,
        write_modes=as_bool(state_recovery.get("write_modes")),
        mode_output_dtype=mode_output_dtype,
        mode_output_format=str(
            state_recovery.get("output_format", "fld")
        ).strip().lower(),
        mode_output_file_name=str(
            state_recovery.get("output_file_name", "DMD_modes")
        ).strip(),
        include_scalar=include_scalar,
        debug=as_bool(state_recovery.get("debug")),
        snapshot_dt=timestep * int(state_recovery["i_stream"]),
        amplitude_regularization=float(
            state_recovery.get("dmd_amplitude_regularization", 1.0e-12)
        ),
        conjugate_tolerance=float(
            state_recovery.get("dmd_conjugate_tolerance", 1.0e-8)
        ),
    )
    if cfg.batch_size < 1:
        raise ValueError("state_recovery.batch_size must be positive")
    if cfg.i_stream < 1:
        raise ValueError("state_recovery.i_stream must be positive")
    if cfg.amplitude_regularization < 0.0:
        raise ValueError("dmd_amplitude_regularization must be non-negative")
    if cfg.conjugate_tolerance <= 0.0:
        raise ValueError("dmd_conjugate_tolerance must be positive")
    return case, cfg


def add_snapshot(
    comm,
    pod,
    ioh,
    fields: list[np.ndarray],
) -> None:
    ioh.copy_fieldlist_to_xi(fields)
    ioh.load_buffer(scale_snapshot=True)
    if ioh.update_from_buffer:
        pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])
        ioh.buffer_index = 0
        ioh.update_from_buffer = False


def flush_buffer(comm, pod, ioh) -> None:
    if ioh.buffer_index > 0 and not ioh.update_from_buffer:
        pod.update(comm, buff=ioh.buff[:, :ioh.buffer_index])
        ioh.buffer_index = 0


def _allreduce(comm, value: np.ndarray) -> np.ndarray:
    result = np.empty_like(value)
    comm.Allreduce(value, result)
    return result


def _find_conjugate(eigenvalues, index, unused, tolerance):
    target = np.conjugate(eigenvalues[index])
    candidates = list(unused)
    if not candidates:
        return None
    distances = [abs(eigenvalues[j] - target) for j in candidates]
    candidate = candidates[int(np.argmin(distances))]
    scale = max(1.0, abs(target))
    if min(distances) <= tolerance * scale:
        return candidate
    return None


def compute_dmd(
    comm,
    pod_x,
    pod_y,
    bm: np.ndarray,
    initial_fields: list[np.ndarray],
    times: list[float],
    cfg: DMDConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute exact DMD modes from two mass-weighted streaming PODs.

    ``pod_x`` contains X=[x_0,...,x_(m-1)] and ``pod_y`` contains
    Y=[x_1,...,x_m].  The return values are local real modes, one spectral
    row per real mode, the selected complex eigenvalues, and complex DMD
    amplitudes.  The spectral rows contain

    ``[omega.real, omega.imag, b.real, b.imag, component]``.

    ``component == 0`` reconstructs ``real(b exp(omega t))`` and
    ``component == 1`` reconstructs ``imag(b exp(omega t))``.  This is the
    real-valued representation of a complex-conjugate DMD pair used by Neko.
    """
    if not times:
        raise RuntimeError("DMD requires at least one consecutive snapshot pair")
    if pod_x.u_1t is None or pod_y.u_1t is None:
        raise RuntimeError("DMD PODs did not receive any snapshots")

    # Convert the two weighted POD bases to physical fields. They are
    # mass-orthonormal after this operation.
    ux = np.asarray(pod_x.u_1t)
    uy = np.asarray(pod_y.u_1t)
    rx = min(cfg.rank, ux.shape[1], uy.shape[1])
    rx = min(rx, len(pod_x.d_1t), len(pod_y.d_1t))
    if rx < 1:
        raise RuntimeError("DMD rank collapsed to zero")

    ux = ux[:, :rx]
    uy = uy[:, :rx]
    sx = np.asarray(pod_x.d_1t[:rx], dtype=np.float64)
    sy = np.asarray(pod_y.d_1t[:rx], dtype=np.float64)
    vx = np.asarray(pod_x.vt_1t[:rx, :], dtype=np.float64)
    vy = np.asarray(pod_y.vt_1t[:rx, :], dtype=np.float64)
    if vx.shape != vy.shape:
        raise RuntimeError("The two streaming PODs have different snapshot counts")

    bm_flat = np.tile(np.asarray(bm).reshape(-1), cfg.n_fields)
    cross_local = (ux * bm_flat[:, None]).conj().T @ uy
    cross = _allreduce(comm, cross_local)

    safe_sx = np.maximum(sx, np.finfo(np.float64).eps)
    transfer = sy[:, None] * (vy @ vx.T) / safe_sx[None, :]
    reduced_operator = cross @ transfer
    eigenvalues, eigenvectors = np.linalg.eig(reduced_operator)
    exact_modes = uy @ transfer.astype(np.complex128) @ eigenvectors

    x0 = np.concatenate([np.asarray(field).reshape(-1) for field in initial_fields])
    weighted_modes = bm_flat[:, None] * exact_modes
    gram_local = exact_modes.conj().T @ weighted_modes
    rhs_local = exact_modes.conj().T @ (bm_flat * x0)
    gram = _allreduce(comm, gram_local)
    rhs = _allreduce(comm, rhs_local)
    scale = max(1.0, float(np.trace(gram).real / max(rx, 1)))
    gram += cfg.amplitude_regularization * scale * np.eye(rx)
    amplitudes = np.linalg.solve(gram, rhs)

    score = np.abs(amplitudes)
    if not np.any(score > 0.0):
        score = np.abs(eigenvalues)
    order = np.argsort(-score)
    unused = set(range(rx))
    real_modes = []
    spectral_rows = []
    selected_eigenvalues = []
    selected_amplitudes = []

    for index in order:
        index = int(index)
        if index not in unused:
            continue
        eigenvalue = eigenvalues[index]
        eig_scale = max(1.0, abs(eigenvalue))
        unused.remove(index)
        if abs(eigenvalue) <= np.finfo(np.float64).tiny:
            continue
        omega = np.log(complex(eigenvalue)) / cfg.snapshot_dt
        if abs(eigenvalue.imag) <= cfg.conjugate_tolerance * eig_scale:
            # A real eigenvector is only defined up to a complex phase by
            # numpy.linalg.eig. Remove that phase before storing a real mode,
            # and apply the inverse phase to its amplitude.
            mode = exact_modes[:, index]
            pivot = int(np.argmax(np.abs(mode)))
            phase = np.exp(-1j * np.angle(mode[pivot]))
            mode = mode * phase
            amplitude = amplitudes[index] / phase
            real_modes.append(mode.real)
            amplitude = complex(amplitude.real, 0.0)
            spectral_rows.append(
                [omega.real, omega.imag, amplitude.real, amplitude.imag, 0.0]
            )
            selected_eigenvalues.append(eigenvalue)
            selected_amplitudes.append(amplitude)
        else:
            conjugate = _find_conjugate(
                eigenvalues, index, unused, cfg.conjugate_tolerance
            )
            if conjugate is None:
                continue
            unused.remove(conjugate)
            if len(real_modes) + 2 > cfg.keep_modes:
                continue
            mode = exact_modes[:, index]
            amplitude = complex(amplitudes[index])
            real_modes.extend([2.0 * mode.real, -2.0 * mode.imag])
            spectral_rows.extend(
                [
                    [omega.real, omega.imag, amplitude.real, amplitude.imag, 0.0],
                    [omega.real, omega.imag, amplitude.real, amplitude.imag, 1.0],
                ]
            )
            selected_eigenvalues.append(eigenvalue)
            selected_amplitudes.append(amplitude)

        if len(real_modes) >= cfg.keep_modes:
            break

    if not real_modes:
        raise RuntimeError("DMD produced no real reconstructible modes")

    n_state = exact_modes.shape[0]
    modes = np.zeros((n_state, cfg.keep_modes), dtype=cfg.dtype)
    spectral_data = np.zeros((cfg.keep_modes, 5), dtype=np.float64)
    n_selected = min(len(real_modes), cfg.keep_modes)
    modes[:, :n_selected] = np.column_stack(real_modes[:n_selected]).real
    spectral_data[:n_selected, :] = np.asarray(
        spectral_rows[:n_selected], dtype=np.float64
    )
    return (
        modes,
        spectral_data,
        np.asarray(selected_eigenvalues, dtype=np.complex128),
        np.asarray(selected_amplitudes, dtype=np.complex128),
    )


def evaluate_spectral_coefficients(spectral_data, times) -> np.ndarray:
    """Evaluate the real DMD representation at the supplied physical times."""
    data = np.asarray(spectral_data, dtype=np.float64)
    times = np.asarray(times, dtype=np.float64)
    omega = data[:, 0] + 1j * data[:, 1]
    amplitude = data[:, 2] + 1j * data[:, 3]
    dynamics = np.exp(times[:, None] * omega[None, :]) * amplitude[None, :]
    coefficients = np.real(dynamics)
    imaginary = np.imag(dynamics)
    components = np.rint(data[:, 4]).astype(np.int64)
    coefficients[:, components == 1] = imaginary[:, components == 1]
    return coefficients


def write_spectral_data(
    comm, spectral_data, eigenvalues, amplitudes, times, rotate_existing
) -> None:
    if comm.Get_rank() != 0:
        return

    # Match the POD writer: when mode output is enabled, preserve each
    # completed decomposition before writing the next one. Keep the spectrum,
    # eigenvalue summary, and diagnostic dynamics from the same decomposition
    # together.
    if rotate_existing:
        for path in (
            "dmd_spectral.csv",
            "dmd_eigenvalues.csv",
            "dmd_modal_dynamics.csv",
        ):
            rotate_existing_file(path)

    header = "omega_real,omega_imag,amplitude_real,amplitude_imag,component"
    np.savetxt(
        "dmd_spectral.csv",
        np.asarray(spectral_data, dtype=np.float64),
        delimiter=",",
        header=header,
        comments="",
    )
    with open("dmd_eigenvalues.csv", "w", encoding="utf-8") as handle:
        handle.write("real,imag,amplitude_real,amplitude_imag\n")
        for eigenvalue, amplitude in zip(eigenvalues, amplitudes):
            handle.write(
                f"{eigenvalue.real:.18e},{eigenvalue.imag:.18e},"
                f"{amplitude.real:.18e},{amplitude.imag:.18e}\n"
            )
    np.savetxt(
        "dmd_modal_dynamics.csv",
        np.column_stack(
            [np.asarray(times, dtype=np.float64),
             evaluate_spectral_coefficients(spectral_data, times)]
        ),
        delimiter=",",
        header="t," + ",".join(
            f"a{index + 1}" for index in range(np.asarray(spectral_data).shape[0])
        ),
        comments="",
    )
