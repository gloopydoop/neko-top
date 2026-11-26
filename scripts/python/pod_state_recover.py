#!/usr/bin/env python3

import os
import sys
import time
import json
import subprocess
from dataclasses import replace

import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.utils import get_gradient
from pysemtools.io.adios2.stream import DataStreamer
from pysemtools.rom.io_help import IoHelp
from pysemtools.rom.pod import POD

from neko_communicator import CtrlClient
from neko_communicator import MODE_ADJOINT
from neko_communicator import MODE_FORWARD
from neko_communicator import MODE_STOP
from neko_communicator import PHASE_ADJ_DONE
from neko_communicator import PHASE_ADJ_RUNNING
from neko_communicator import PHASE_ADJ_SNAPSHOT
from neko_communicator import PHASE_FWD_DONE
from neko_communicator import PHASE_FWD_RUNNING
from neko_communicator import get_peer_root
from neko_communicator import make_local_comm
from pod_state_recover_tools import EnergyState
from pod_state_recover_tools import PODConfig
from pod_state_recover_tools import add_snapshot
from pod_state_recover_tools import available_modes
from pod_state_recover_tools import build_time_coefficients
from pod_state_recover_tools import flush_buffer
from pod_state_recover_tools import load_pod_config
from pod_state_recover_tools import make_pod
from pod_state_recover_tools import make_admissibility_projector
from pod_state_recover_tools import recv_field
from pod_state_recover_tools import recv_fields
from pod_state_recover_tools import report_energy_capture
from pod_state_recover_tools import set_debug
from pod_state_recover_tools import stream_mode_fields
from pod_state_recover_tools import write_modes_to_disk
from pod_state_recover_tools import write_singular_values
from pod_state_recover_tools import write_time_coefficients
from pod_state_recover_tools import WindowStatistics
from statistical_closure import fit_boussinesq_nut, fit_scalar_gammat
from recycled_pod_tools import (
    append_factor_history,
    common_basis_distributed,
    cross_mass_singular_values,
    factors_from_pod,
    recycle_distributed,
)
from pod_rank_mode_diagnostics import append_family_mode_diagnostics


def source_revision() -> str:
    """Return the repository revision used by this POD peer, if available."""

    try:
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
        revision = result.stdout.strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet"], cwd=root, check=False
        ).returncode != 0 or bool(subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=root, check=True, capture_output=True, text=True,
        ).stdout.strip())
        return revision + ("-dirty" if dirty else "")
    except (OSError, subprocess.SubprocessError):
        return os.environ.get("ROM_SOURCE_REVISION", "unknown")


ROM_OUTPUT_ROOT = "rom_output"


def outer_output_directory(outer: int) -> str:
    """Return the permanent POD output directory for one ROM outer index."""

    return os.path.join(ROM_OUTPUT_ROOT, f"outer_{outer:03d}")


def prepare_outer_output_directory(comm: MPI.Comm, outer: int) -> str:
    """Create one shared outer directory from the Python root rank."""

    path = outer_output_directory(outer)
    if comm.Get_rank() == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()
    return path


def append_pod_rank_provenance(
    comm: MPI.Comm,
    outer_dir: str,
    outer_iteration: int,
    family: str,
    requested_rank: int | None,
    actual_pod_rank: int | None,
    historical_candidate_rank: int | None,
    recycled_rank: int | None,
    rank_cap: int | None,
    status: str = "available",
) -> None:
    """Append per-family rank provenance without changing basis streaming."""

    if comm.Get_rank() != 0:
        return
    record = {
        "outer_iteration": int(outer_iteration),
        "family": family,
        "status": status,
        "requested_pod_rank": requested_rank,
        "actual_pod_rank": actual_pod_rank,
        "historical_candidate_rank": historical_candidate_rank,
        "recycled_rank": recycled_rank,
        "recycling_rank_cap": rank_cap,
    }
    path = os.path.join(outer_dir, "pod_rank_provenance.jsonl")
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
        handle.flush()


def update_latest_alias(comm: MPI.Comm, target: str, alias: str) -> None:
    """Maintain an explicit root-level symlink to the newest permanent file."""

    if comm.Get_rank() == 0:
        if os.path.lexists(alias):
            os.unlink(alias)
        relative_target = os.path.relpath(target, os.path.dirname(alias) or ".")
        os.symlink(relative_target, alias)
    comm.Barrier()


def write_basis_metadata(
    comm: MPI.Comm,
    outer: int,
    outer_dir: str,
    cfg: PODConfig,
    history_length: int,
    forward_history: list,
    adjoint_history: list,
    forward_candidate_modes: int,
    adjoint_candidate_modes: int,
    forward_rank: int,
    adjoint_rank: int,
    overlap_singular_values: np.ndarray | None,
    adjoint_enrichment: bool,
) -> None:
    """Write one explicit provenance manifest for a streamed ROM basis."""

    metadata = {
        "outer_iteration": outer,
        "basis_forward": cfg.basis_forward,
        "basis_adjoint": cfg.basis_adjoint,
        "formulation": cfg.formulation,
        "sensitivity_source": cfg.sensitivity_source,
        "basis_transform": cfg.basis_transform,
        "differentiate_basis": cfg.differentiate_basis,
        "basis_transform_tau": float(cfg.basis_transform_tau),
        "history_outer_iterations_requested": cfg.history_outer_iterations,
        "history_outer_iterations_included": history_length,
        "forward_history_iterations": [
            int(item.iteration) for item in forward_history
        ],
        "adjoint_history_iterations": [
            int(item.iteration) for item in adjoint_history
        ],
        "forward_candidate_modes": int(forward_candidate_modes),
        "adjoint_candidate_modes": int(adjoint_candidate_modes),
        "forward_retained_rank": int(forward_rank),
        "adjoint_retained_rank": int(adjoint_rank),
        "adjoint_enrichment": bool(adjoint_enrichment),
        "forward_pod_files": {
            "time_coefficients": "forward_pod_coefficients.csv",
            "singular_values": "forward_pod_singular_values.csv",
            "modes": "forward_pod_modes",
        },
        "forward_recycling_mode_directory": "forward_pod_modes",
        "forward_pod_contributors": [
            {
                "pod_iteration": int(item.iteration),
                "outer_iteration": int(item.iteration),
                "modes": os.path.join(
                    "..",
                    f"outer_{int(item.iteration):03d}",
                    "forward_pod_modes",
                ),
                "singular_values": os.path.join(
                    "..",
                    f"outer_{int(item.iteration):03d}",
                    "forward_pod_singular_values.csv",
                ),
            }
            for item in forward_history
        ],
        "adjoint_pod_files": (
            {
                "time_coefficients": "adjoint_pod_coefficients.csv",
                "singular_values": "adjoint_pod_singular_values.csv",
                "modes": "adjoint_pod_modes",
            }
            if adjoint_history
            else None
        ),
        "adjoint_pod_contributors": [
            {
                "pod_iteration": int(item.iteration),
                "outer_iteration": int(item.iteration),
                "modes": os.path.join(
                    "..",
                    f"outer_{int(item.iteration):03d}",
                    "adjoint_pod_modes",
                ),
                "singular_values": os.path.join(
                    "..",
                    f"outer_{int(item.iteration):03d}",
                    "adjoint_pod_singular_values.csv",
                ),
            }
            for item in adjoint_history
        ],
        "source_revision": source_revision() if comm.Get_rank() == 0 else "",
        "output_layout": "outer_NNN with standard fld counters and flat inner-indexed ROM CSV files",
        "paths_relative_to": "this basis_metadata.json",
    }
    if overlap_singular_values is not None:
        metadata["forward_adjoint_mass_overlap_singular_values"] = (
            np.asarray(overlap_singular_values, dtype=np.float64).tolist()
        )
    path = os.path.join(outer_dir, "basis_metadata.json")
    if comm.Get_rank() == 0:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)
    update_latest_alias(comm, path, "basis_metadata.json")


def stream_rom_basis(comm, ds, ctrl, common_basis, bm, cfg) -> None:
    """Stream the local common basis over the reverse ADIOS2 connection."""

    rank = comm.Get_rank()
    modes = np.asarray(common_basis.modes)
    field_size = np.asarray(bm).size
    expected_rows = cfg.n_input_fields * field_size
    if modes.ndim != 2 or modes.shape[0] != expected_rows:
        raise RuntimeError(
            "ROM basis rows do not match the streamed local input fields: "
            f"got {modes.shape[0] if modes.ndim == 2 else modes.ndim}, "
            f"expected {expected_rows}"
        )

    n_modes = int(modes.shape[1])
    n_fields = int(cfg.n_input_fields)
    comm.Barrier()
    if rank == 0:
        ctrl.send_basis_header(n_modes, n_fields)
    comm.Barrier()

    for mode_index in range(n_modes):
        for field_index in range(n_fields):
            first = field_index * field_size
            last = (field_index + 1) * field_size
            field = np.ascontiguousarray(
                modes[first:last, mode_index], dtype=cfg.dtype
            )
            ds.stream(field)

    comm.Barrier()


def stream_rom_bases(comm, ds, ctrl, trial_basis, test_basis, bm, cfg) -> None:
    """Stream independent forward-trial and adjoint-test basis blocks.

    Each block has its own header and payload.  The Fortran receiver consumes
    the first payload before waiting for the second header, so this preserves
    the existing control protocol without changing generic state recovery.
    """

    rank = comm.Get_rank()
    field_size = np.asarray(bm).size
    blocks = (trial_basis, test_basis)
    for block in blocks:
        modes = np.asarray(block.modes)
        expected_rows = cfg.n_input_fields * field_size
        if modes.ndim != 2 or modes.shape[0] != expected_rows:
            raise RuntimeError(
                "Petrov basis rows do not match the streamed local fields: "
                f"got {modes.shape[0] if modes.ndim == 2 else modes.ndim}, "
                f"expected {expected_rows}"
            )
        comm.Barrier()
        if rank == 0:
            ctrl.send_basis_header(int(modes.shape[1]), int(cfg.n_input_fields))
        comm.Barrier()
        for mode_index in range(modes.shape[1]):
            for field_index in range(cfg.n_input_fields):
                first = field_index * field_size
                last = (field_index + 1) * field_size
                ds.stream(np.ascontiguousarray(
                    modes[first:last, mode_index], dtype=cfg.dtype
                ))
        comm.Barrier()


def init_runtime(
    comm: MPI.Comm,
    cfg: PODConfig,
) -> tuple[
    DataStreamer,
    Mesh,
    Coef,
    np.ndarray,
    POD,
    IoHelp,
    list[np.ndarray],
]:
    ds = DataStreamer(comm)

    x = recv_field(ds, cfg.dtype)
    y = recv_field(ds, cfg.dtype)
    z = recv_field(ds, cfg.dtype)

    initial_fields = recv_fields(ds, cfg.n_input_fields, cfg.dtype)

    mesh = Mesh(comm, x=x, y=y, z=z, create_connectivity=False)
    coef = Coef(mesh, comm)
    pod, ioh = make_pod(comm, coef.B, cfg)
    return ds, mesh, coef, pod, ioh, initial_fields


def build_forward_only_basis(
    comm: MPI.Comm,
    ds: DataStreamer,
    ctrl: CtrlClient,
    mesh: Mesh,
    bm: np.ndarray,
    cfg: PODConfig,
    pod: POD,
    ioh: IoHelp,
    pod_iteration: int,
    forward_history: list,
    outer_dir: str,
) -> None:
    """Recycle and stream a forward-only common velocity/scalar basis.

    A forward-only experiment still uses the normal common-basis transport and
    Neko-side admissibility gate.  It simply supplies no adjoint columns to
    the recycling/orthogonalisation step.  This keeps zero adjoint enrichment
    distinct from silently omitting the basis handoff altogether.
    """

    append_factor_history(
        forward_history,
        factors_from_pod(
            pod, ioh, bm, cfg.n_input_fields, pod_iteration,
            cfg.scalar_metric_weight,
        ),
        cfg.history_outer_iterations,
    )
    local_field_weights = cfg.field_metric_weights[:cfg.n_input_fields]
    local_mass = np.concatenate([
        weight * np.asarray(bm, dtype=np.float64).reshape(-1)
        for weight in local_field_weights
    ])
    rank_cap = None
    recycled_forward = recycle_distributed(
        comm, forward_history, local_mass, rank_cap
    )
    append_pod_rank_provenance(
        comm,
        outer_dir,
        pod_iteration,
        "forward",
        cfg.keep_modes,
        available_modes(pod, cfg.keep_modes),
        sum(item.modes.shape[1] for item in forward_history),
        int(recycled_forward.modes.shape[1]),
        rank_cap,
    )
    append_pod_rank_provenance(
        comm,
        outer_dir,
        pod_iteration,
        "adjoint",
        None,
        None,
        0,
        0,
        None,
        status="not_requested",
    )
    common_basis, basis_diag = common_basis_distributed(
        comm,
        recycled_forward,
        None,
        local_mass,
        rank_cap,
        adjoint_projector=make_admissibility_projector(mesh, bm, cfg),
    )
    append_family_mode_diagnostics(
        comm,
        outer_dir,
        pod_iteration,
        "forward",
        recycled_forward.rank_diagnostics,
        basis_diag,
    )
    if comm.Get_rank() == 0:
        print(
            f"Forward basis: {cfg.basis_forward}\n"
            f"Adjoint basis: {cfg.basis_adjoint}\n"
            f"POD history requested: {cfg.history_outer_iterations}\n"
            f"POD history available: {len(forward_history)}\n"
            f"Forward candidate modes: {sum(item.modes.shape[1] for item in forward_history)}\n"
            "Adjoint candidate modes: 0\n"
            f"Forward streamed rank: {common_basis.modes.shape[1]}",
            flush=True,
        )
    stream_rom_basis(comm, ds, ctrl, common_basis, bm, cfg)
    history_path = os.path.join(outer_dir, "recycled_rom_history.json")
    if comm.Get_rank() == 0:
        with open(history_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "history_outer_iterations": cfg.history_outer_iterations,
                    "history_requested": cfg.history_outer_iterations,
                    "history_available": len(forward_history),
                    "forward_iterations": [item.iteration for item in forward_history],
                    "adjoint_iterations": [],
                    "forward_candidate_modes": int(sum(item.modes.shape[1] for item in forward_history)),
                    "adjoint_candidate_modes": 0,
                    "basis_rank": int(common_basis.modes.shape[1]),
                    "mass_gram_error": basis_diag["mass_gram_error"],
                    "adjoint_enrichment": False,
                    "basis_forward": cfg.basis_forward,
                    "basis_adjoint": cfg.basis_adjoint,
                    "basis_transform": cfg.basis_transform,
                    "differentiate_basis": cfg.differentiate_basis,
                    "basis_transform_tau": float(cfg.basis_transform_tau),
                    "formulation": cfg.formulation,
                    "sensitivity_source": cfg.sensitivity_source,
                    "streamed_forward_rank": int(common_basis.modes.shape[1]),
                    "source_revision": source_revision(),
                },
                handle,
                indent=2,
            )
    update_latest_alias(comm, history_path, "recycled_rom_history.json")
    write_basis_metadata(
        comm,
        pod_iteration,
        outer_dir,
        cfg,
        len(forward_history),
        forward_history,
        [],
        sum(item.modes.shape[1] for item in forward_history),
        0,
        int(common_basis.modes.shape[1]),
        0,
        None,
        False,
    )


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(
            f"Usage: {os.path.basename(sys.argv[0])} <case-file>"
        )

    world = MPI.COMM_WORLD
    comm = make_local_comm(world)
    rank = comm.Get_rank()
    peer_root = get_peer_root()

    case_path = os.path.abspath(sys.argv[1])
    if not os.path.isfile(case_path):
        raise FileNotFoundError(f"Case file not found: {case_path}")
    case_cfg, cfg = load_pod_config(case_path)
    set_debug(cfg.debug)

    ds, mesh, coef, pod, ioh, initial_fields = init_runtime(comm, cfg)
    bm = coef.B
    received_reference_design = False
    times = []
    energy = EnergyState()
    # The ordinary POD is always populated.  It is the independent
    # full-order-adjoint reconstruction basis and must not change when the
    # ROM-only transformed stream is enabled.
    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)
    window_stats = None
    if cfg.input_mode == "windowed_statistics":
        # The initial condition is not a statistical sample.  Window means are
        # inserted into the POD basis at FWD_DONE below.
        pod, ioh = make_pod(comm, bm, cfg)
        times = []
        energy = EnergyState()
        window_stats = WindowStatistics(
            cfg.window_duration,
            cfg.statistics_start_time,
            cfg.dtype,
            cfg.keep_partial_window,
        )
    mode_output_index = 0
    adjoint_mode_output_index = 0
    forward_pod_rank = None
    adjoint_pod_rank = None
    pod_iteration = 0
    singular_value_history_new = True
    adjoint_pod = None
    adjoint_ioh = None
    adjoint_times = []
    adjoint_energy = EnergyState()
    adjoint_cfg = replace(
        cfg,
        closure_enabled=False,
        scalar_closure="off",
        keep_modes=cfg.keep_modes,
        mode_output_file_name=cfg.adjoint_mode_output_file_name,
    )
    forward_history = []
    adjoint_history = []

    ctrl = None
    if rank == 0:
        ctrl = CtrlClient(world, peer_root, cfg.debug)

    my_count = getattr(ds, "py2f_field_my_count", None)
    if my_count is None:
        my_count = initial_fields[0].size
    zero_field = np.zeros(my_count, dtype=cfg.dtype)

    try:
        while True:
            if rank == 0:
                state = ctrl.read_state()
            else:
                state = (None, None, None, None)

            mode, phase, step, tcur = comm.bcast(state, root=0)

            if mode is None or mode == MODE_STOP:
                break

            if mode == MODE_FORWARD and phase == PHASE_FWD_RUNNING:
                if cfg.basis_transform != "none" and not received_reference_design:
                    # The first control event is the explicit mapped-design
                    # handoff sent by Fortran immediately before this batch.
                    # Keep it for provenance only: all forward POD snapshots
                    # and recycled factors remain ordinary physical fields.
                    _ = recv_field(ds, cfg.dtype)
                    received_reference_design = True
                raw_fields = recv_fields(ds, cfg.n_input_fields, cfg.dtype)
                if window_stats is None:
                    add_snapshot(comm, pod, ioh, bm, raw_fields, tcur, times, energy)
                else:
                    window_stats.add(raw_fields, tcur)
                continue

            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                # Reserve the permanent outer directory before writing any
                # window/closure diagnostics.  These diagnostics belong to
                # the POD basis produced by this forward lifecycle, not to a
                # moving root-level "latest" file.
                outer_iteration = pod_iteration
                outer_dir = prepare_outer_output_directory(
                    comm, outer_iteration
                )
                if window_stats is not None:
                    windows = window_stats.finish()
                    if not windows:
                        raise RuntimeError(
                            "windowed_statistics produced no valid windows; "
                            "check statistics_start_time, window_duration and run length"
                        )
                    if rank == 0:
                        with open(
                            os.path.join(
                                outer_dir, "window_statistics_manifest.json"
                            ),
                            "w",
                            encoding="utf-8",
                        ) as handle:
                            json.dump(
                                {
                                    "input_mode": cfg.input_mode,
                                    "window_duration": cfg.window_duration,
                                    "statistics_start_time": cfg.statistics_start_time,
                                    "design_id": case_cfg.get("design_id", None),
                                    "timestep": cfg.snapshot_dt / cfg.i_stream,
                                    "sample_stride": cfg.i_stream,
                                    "dtype": np.dtype(cfg.dtype).name,
                                    "mesh": {
                                        "lx": int(mesh.lx), "ly": int(mesh.ly),
                                        "lz": int(mesh.lz), "nelv": int(mesh.nelv),
                                    },
                                    "field_mapping": [
                                        "window_mean_u", "window_mean_v", "window_mean_w",
                                        "window_mean_phi" if cfg.include_scalar else None,
                                        "adjoint_nut" if cfg.closure_enabled else None,
                                        "adjoint_gammat" if cfg.scalar_closure_enabled else None,
                                    ],
                                    "windows": [{
                                        k: v for k, v in item.items()
                                        if k not in ("mean", "raw", "reynolds", "scalar_mean", "scalar_raw_flux", "scalar_flux")
                                    } for item in windows],
                                },
                                handle,
                                indent=2,
                            )
                        with open(
                            os.path.join(
                                outer_dir, "window_statistics_diagnostics.csv"
                            ),
                            "w",
                            encoding="utf-8",
                        ) as handle:
                            handle.write("window,center,samples,weight,mean_k_local_sum\n")
                            for idx, item in enumerate(windows):
                                reynolds = item["reynolds"]
                                local_k = 0.5 * sum(reynolds[j, j].sum() for j in range(3))
                                global_k = comm.allreduce(float(local_k), op=MPI.SUM)
                                handle.write(
                                    f"{idx},{item['center']:.16e},{item['samples']},"
                                    f"{item['weight']:.16e},{global_k:.16e}\n"
                                )
                    elif windows:
                        # Keep the collective aligned with rank zero's
                        # diagnostics write above.
                        for item in windows:
                            reynolds = item["reynolds"]
                            local_k = 0.5 * sum(reynolds[j, j].sum() for j in range(3))
                            comm.allreduce(float(local_k), op=MPI.SUM)
                    closure_rows = []
                    scalar_closure_rows = []
                    for item in windows:
                        pod_fields = item["mean"]
                        if cfg.closure_enabled:
                            gradients = get_gradient(mesh, coef, item["mean"][:3])
                            # pySEMTools returns (..., field, derivative).
                            gradients = np.moveaxis(gradients, (-2, -1), (0, 1))
                            nut, nut_diag = fit_boussinesq_nut(
                                item["reynolds"],
                                gradients + np.swapaxes(gradients, 0, 1),
                                regularisation=cfg.closure_regularisation,
                                maximum=cfg.closure_nut_max,
                                clip_negative=True,
                            )
                            pod_fields = list(pod_fields)
                            pod_fields.append(nut * cfg.closure_field_scale)
                            closure_rows.append((item["center"], nut_diag))
                        if cfg.scalar_closure_enabled:
                            scalar_gradient = get_gradient(
                                mesh, coef, [item["scalar_mean"]]
                            )
                            scalar_gradient = np.moveaxis(
                                scalar_gradient, (-2, -1), (0, 1)
                            )[0]
                            if cfg.scalar_closure == "flux_least_squares":
                                gammat, gammat_diag = fit_scalar_gammat(
                                    item["scalar_flux"],
                                    scalar_gradient,
                                    regularisation=cfg.gammat_regularisation,
                                    maximum=cfg.gammat_max,
                                    clip_negative=True,
                                )
                            else:
                                if not cfg.closure_enabled:
                                    raise RuntimeError(
                                        "scalar_closure=schmidt requires an eddy-viscosity closure"
                                    )
                                gammat = nut / cfg.turbulent_schmidt
                                gammat_diag = {
                                    "negative_fraction_raw": 0.0,
                                    "clipped_fraction": 0.0,
                                    "valid_fraction": 1.0,
                                    "mean_fit_relative_residual": 0.0,
                                    "max_gammat": float(np.max(gammat)),
                                }
                            pod_fields = list(pod_fields)
                            pod_fields.append(gammat * cfg.gammat_field_scale)
                            scalar_closure_rows.append((item["center"], gammat_diag))
                        add_snapshot(
                            comm, pod, ioh, bm, pod_fields, item["center"], times, energy,
                        )
                    if rank == 0:
                        print(
                            f"Windowed statistics: {len(windows)} windows, "
                            f"samples={[item['samples'] for item in windows]}",
                            flush=True,
                        )
                    if rank == 0 and closure_rows:
                        with open(
                            os.path.join(outer_dir, "closure_diagnostics.csv"),
                            "w",
                            encoding="utf-8",
                        ) as handle:
                            handle.write(
                                "center,valid_fraction,negative_fraction_raw,"
                                "clipped_fraction,mean_fit_relative_residual,max_nut\n"
                            )
                            for center, diag in closure_rows:
                                handle.write(
                                    f"{center:.16e},{diag['valid_fraction']:.16e},"
                                    f"{diag['negative_fraction_raw']:.16e},"
                                    f"{diag['clipped_fraction']:.16e},"
                                    f"{diag['mean_fit_relative_residual']:.16e},"
                                    f"{diag['max_nut']:.16e}\n"
                                )
                    if rank == 0 and scalar_closure_rows:
                        with open(
                            os.path.join(
                                outer_dir, "scalar_closure_diagnostics.csv"
                            ),
                            "w",
                            encoding="utf-8",
                        ) as handle:
                            handle.write(
                                "center,valid_fraction,negative_fraction_raw,"
                                "clipped_fraction,mean_fit_relative_residual,max_gammat\n"
                            )
                            for center, diag in scalar_closure_rows:
                                handle.write(
                                    f"{center:.16e},{diag['valid_fraction']:.16e},"
                                    f"{diag['negative_fraction_raw']:.16e},"
                                    f"{diag['clipped_fraction']:.16e},"
                                    f"{diag['mean_fit_relative_residual']:.16e},"
                                    f"{diag['max_gammat']:.16e}\n"
                                )
                flush_buffer(comm, pod, ioh)
                pod.scale_modes(comm, bm1sqrt=ioh.bm1sqrt, op="div")
                pod.rotate_local_modes_to_global(comm)

                coeffs, header, n_avail = build_time_coefficients(
                    comm, pod, cfg, times
                )
                forward_pod_rank = int(n_avail)
                forward_output_cfg = replace(
                    cfg,
                    mode_output_file_name=os.path.join(
                        outer_dir, "forward_pod_modes"
                    ),
                )
                write_time_coefficients(
                    comm,
                    coeffs,
                    header,
                    file_name=os.path.join(
                        outer_dir, "forward_pod_coefficients.csv"
                    ),
                    latest_alias="pod_time_coeffs.csv",
                )
                pod_iteration = outer_iteration
                write_singular_values(
                    comm,
                    pod,
                    cfg.keep_modes,
                    pod_iteration,
                    file_name=os.path.join(
                        outer_dir, "forward_pod_singular_values.csv"
                    ),
                    latest_alias="pod_singular_values.csv",
                )
                singular_value_history_new = False
                if rank == 0:
                    report_energy_capture(pod, energy, cfg.keep_modes)
                    ctrl.send_cmd(MODE_ADJOINT, PHASE_ADJ_RUNNING)

                stream_mode_fields(ds, ioh, pod, cfg, zero_field, n_avail)
                mode_output_index = write_modes_to_disk(
                    comm,
                    mesh,
                    ioh,
                    pod,
                    forward_output_cfg,
                    zero_field,
                    n_avail,
                    0,
                )
                if cfg.basis_adjoint != "forward":
                    adjoint_pod, adjoint_ioh = make_pod(comm, bm, adjoint_cfg)
                    adjoint_times = []
                    adjoint_energy = EnergyState()
                continue

            if mode == MODE_ADJOINT and phase == PHASE_ADJ_SNAPSHOT:
                fields = recv_fields(ds, cfg.n_input_fields, cfg.dtype)
                if cfg.basis_adjoint != "forward":
                    add_snapshot(
                        comm, adjoint_pod, adjoint_ioh, bm, fields, tcur,
                        adjoint_times, adjoint_energy
                    )
                continue

            if mode == MODE_ADJOINT and phase == PHASE_ADJ_DONE:
                if cfg.basis_adjoint == "forward":
                    build_forward_only_basis(
                        comm, ds, ctrl, mesh, bm, cfg, pod, ioh,
                        pod_iteration, forward_history,
                        outer_output_directory(pod_iteration),
                    )
                if cfg.basis_adjoint != "forward" and adjoint_pod is not None:
                    flush_buffer(comm, adjoint_pod, adjoint_ioh)
                    adjoint_pod.scale_modes(comm, bm1sqrt=adjoint_ioh.bm1sqrt, op="div")
                    adjoint_pod.rotate_local_modes_to_global(comm)
                    adjoint_coeffs, adjoint_header, adjoint_pod_rank = build_time_coefficients(
                        comm, adjoint_pod, adjoint_cfg, adjoint_times
                    )
                    adjoint_pod_rank = int(adjoint_pod_rank)
                    outer_dir = outer_output_directory(pod_iteration)
                    adjoint_output_cfg = replace(
                        adjoint_cfg,
                        mode_output_file_name=os.path.join(
                            outer_dir, "adjoint_pod_modes"
                        ),
                    )
                    write_time_coefficients(
                        comm,
                        adjoint_coeffs,
                        adjoint_header,
                        file_name=os.path.join(
                            outer_dir, "adjoint_pod_coefficients.csv"
                        ),
                        latest_alias="adjoint_pod_time_coeffs.csv",
                    )
                    write_singular_values(
                        comm,
                        adjoint_pod,
                        adjoint_cfg.keep_modes,
                        pod_iteration,
                        file_name=os.path.join(
                            outer_dir, "adjoint_pod_singular_values.csv"
                        ),
                        latest_alias="adjoint_pod_singular_values.csv",
                    )
                    if rank == 0:
                        print(
                            f"Adjoint POD: {len(adjoint_times)} snapshots, "
                            f"{adjoint_energy.count} accumulated samples",
                            flush=True,
                        )
                    adjoint_mode_output_index = write_modes_to_disk(
                        comm,
                        mesh,
                        adjoint_ioh,
                        adjoint_pod,
                        adjoint_output_cfg,
                        zero_field,
                        available_modes(adjoint_pod, adjoint_cfg.keep_modes),
                        0,
                    )
                    # Keep compact factor metadata for the next three
                    # outer iterations.  The local pySEMTools mode block is
                    # retained separately above for inspection and plotting.
                    # Retain per-iteration V Sigma factors, never a recycled
                    # aggregate.  The compact Gram-matrix SVD is distributed:
                    # only small rank-by-rank matrices are reduced over MPI.
                    append_factor_history(
                        forward_history,
                        factors_from_pod(
                            pod,
                            ioh,
                            bm,
                            cfg.n_input_fields,
                            pod_iteration,
                            cfg.scalar_metric_weight,
                        ),
                        cfg.history_outer_iterations,
                    )
                    append_factor_history(
                        adjoint_history,
                        factors_from_pod(
                            adjoint_pod, adjoint_ioh, bm,
                            adjoint_cfg.n_fields, pod_iteration,
                            cfg.scalar_metric_weight,
                        ),
                        cfg.history_outer_iterations,
                    )
                    history_length = min(len(forward_history), len(adjoint_history))
                    forward_candidate_modes = int(sum(
                        item.modes.shape[1] for item in forward_history
                    ))
                    adjoint_candidate_modes = int(sum(
                        item.modes.shape[1] for item in adjoint_history
                    ))
                    local_field_weights = cfg.field_metric_weights[:cfg.n_input_fields]
                    local_mass = np.concatenate([
                        weight * np.asarray(bm, dtype=np.float64).reshape(-1)
                        for weight in local_field_weights
                    ])
                    # Dual-space formulations always stream independent V_f and
                    # V_a blocks.  The Fortran handoff constructs V_E when the
                    # selected formulation requests it; Python never projects
                    # V_a into V_f or silently replaces it by a complement.
                    forward_cap = None
                    adjoint_cap = None
                    recycled_forward = recycle_distributed(
                        comm, forward_history, local_mass, forward_cap
                    )
                    recycled_adjoint = recycle_distributed(
                        comm, adjoint_history, local_mass, adjoint_cap
                    )
                    append_family_mode_diagnostics(
                        comm,
                        outer_dir,
                        pod_iteration,
                        "forward",
                        recycled_forward.rank_diagnostics,
                        None,
                    )
                    append_family_mode_diagnostics(
                        comm,
                        outer_dir,
                        pod_iteration,
                        "adjoint",
                        recycled_adjoint.rank_diagnostics,
                        None,
                    )
                    append_pod_rank_provenance(
                        comm,
                        outer_dir,
                        pod_iteration,
                        "forward",
                        cfg.keep_modes,
                        forward_pod_rank,
                        forward_candidate_modes,
                        int(recycled_forward.modes.shape[1]),
                        forward_cap,
                    )
                    append_pod_rank_provenance(
                        comm,
                        outer_dir,
                        pod_iteration,
                        "adjoint",
                        adjoint_cfg.keep_modes,
                        adjoint_pod_rank,
                        adjoint_candidate_modes,
                        int(recycled_adjoint.modes.shape[1]),
                        adjoint_cap,
                    )
                    if recycled_forward.modes.shape[1] < 1:
                        raise RuntimeError("ROM recycling produced no forward modes")
                    if recycled_adjoint.modes.shape[1] < 1:
                        raise RuntimeError("ROM recycling produced no adjoint modes")
                    if rank == 0:
                        print(
                            f"Forward basis: {cfg.basis_forward}\n"
                            f"Adjoint basis: {cfg.basis_adjoint}\n"
                            f"POD history requested: {cfg.history_outer_iterations}\n"
                            f"POD history available: {history_length}\n"
                            f"Forward candidate modes: {forward_candidate_modes}\n"
                            f"Adjoint candidate modes: {adjoint_candidate_modes}\n"
                            f"Forward streamed rank: {recycled_forward.modes.shape[1]}\n"
                            f"Adjoint streamed rank: {recycled_adjoint.modes.shape[1]}",
                            flush=True,
                        )
                    stream_rom_bases(
                        comm, ds, ctrl, recycled_forward, recycled_adjoint, bm, cfg
                    )
                    basis_rank = int(recycled_forward.modes.shape[1])
                    basis_diag = {"mass_gram_error": 0.0}
                    overlap_singular_values = cross_mass_singular_values(
                        comm,
                        recycled_forward.modes,
                        recycled_adjoint.modes,
                        local_mass,
                    )
                    history_path = os.path.join(
                        outer_dir, "recycled_rom_history.json"
                    )
                    if rank == 0:
                        with open(history_path, "w", encoding="utf-8") as handle:
                            json.dump(
                                {
                                    "history_outer_iterations": history_length,
                                    "history_requested": cfg.history_outer_iterations,
                                    "history_available": history_length,
                                    "forward_iterations": [item.iteration for item in forward_history],
                                    "adjoint_iterations": [item.iteration for item in adjoint_history],
                                    "forward_candidate_modes": forward_candidate_modes,
                                    "adjoint_candidate_modes": adjoint_candidate_modes,
                                    "basis_rank": basis_rank,
                                    "mass_gram_error": basis_diag["mass_gram_error"],
                                    "basis_forward": cfg.basis_forward,
                                    "basis_adjoint": cfg.basis_adjoint,
                                    "formulation": cfg.formulation,
                                    "sensitivity_source": cfg.sensitivity_source,
                                    "basis_transform": cfg.basis_transform,
                                    "differentiate_basis": cfg.differentiate_basis,
                                    "basis_transform_tau": float(cfg.basis_transform_tau),
                                    "streamed_forward_rank": int(recycled_forward.modes.shape[1]),
                                    "streamed_adjoint_rank": int(recycled_adjoint.modes.shape[1]),
                                    "streamed_enriched_candidate_rank": (
                                        int(recycled_forward.modes.shape[1] + recycled_adjoint.modes.shape[1])
                                        if cfg.formulation in ("forward_enriched_adjoint", "enriched_galerkin")
                                        else 0
                                    ),
                                    "forward_adjoint_mass_overlap_singular_values": overlap_singular_values.tolist(),
                                    "source_revision": source_revision(),
                                },
                                handle,
                                indent=2,
                            )
                    update_latest_alias(
                        comm, history_path, "recycled_rom_history.json"
                    )
                    write_basis_metadata(
                        comm,
                        pod_iteration,
                        outer_dir,
                        cfg,
                        history_length,
                        forward_history,
                        adjoint_history,
                        forward_candidate_modes,
                        adjoint_candidate_modes,
                        int(recycled_forward.modes.shape[1]),
                        int(recycled_adjoint.modes.shape[1]),
                        overlap_singular_values,
                        cfg.basis_adjoint != "forward",
                    )
                    adjoint_pod = None
                    adjoint_ioh = None
                # The next forward/adjoint lifecycle belongs to the next
                # physical design.  The current basis and all of its files
                # remain permanently under outer_NNN.
                pod_iteration += 1
                pod, ioh = make_pod(comm, bm, cfg)
                times = []
                energy = EnergyState()
                if window_stats is None:
                    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)
                else:
                    window_stats = WindowStatistics(
                        cfg.window_duration,
                        cfg.statistics_start_time,
                        cfg.dtype,
                        cfg.keep_partial_window,
                    )
                if cfg.basis_transform != "none":
                    received_reference_design = False
                continue

            time.sleep(0.001)

    finally:
        if rank == 0 and ctrl is not None:
            ctrl.close()
        comm.Barrier()
        ds.finalize()


if __name__ == "__main__":
    main()
