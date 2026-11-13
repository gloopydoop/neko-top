#!/usr/bin/env python3

import os
import sys
import time
import json

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
from neko_communicator import PHASE_FWD_DONE
from neko_communicator import PHASE_FWD_RUNNING
from neko_communicator import get_peer_root
from neko_communicator import make_local_comm
from urans_pod_state_recover_tools import EnergyState
from urans_pod_state_recover_tools import PODConfig
from urans_pod_state_recover_tools import add_snapshot
from urans_pod_state_recover_tools import build_time_coefficients
from urans_pod_state_recover_tools import flush_buffer
from urans_pod_state_recover_tools import load_pod_config
from urans_pod_state_recover_tools import make_pod
from urans_pod_state_recover_tools import recv_field
from urans_pod_state_recover_tools import recv_fields
from urans_pod_state_recover_tools import report_energy_capture
from urans_pod_state_recover_tools import set_debug
from urans_pod_state_recover_tools import stream_mode_fields
from urans_pod_state_recover_tools import write_modes_to_disk
from urans_pod_state_recover_tools import write_singular_values
from urans_pod_state_recover_tools import write_time_coefficients
from statistical_closure import fit_boussinesq_nut, fit_scalar_gammat


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
    times = []
    energy = EnergyState()
    # The initial condition is received only to establish the distributed
    # field layout. URANS-POD consumes one Fortran-computed statistics record
    # per completed window.
    times = []
    energy = EnergyState()
    window_metadata = []
    closure_rows = []
    scalar_closure_rows = []
    mode_output_index = 0
    pod_iteration = 0
    singular_value_history_new = True

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
                statistics = recv_fields(ds, cfg.n_statistics_fields, cfg.dtype)
                n_base = cfg.n_input_fields
                means = statistics[:n_base]
                reynolds_components = statistics[n_base:n_base + 6]
                reynolds = np.zeros(
                    (3, 3) + means[0].shape,
                    dtype=cfg.dtype,
                )
                reynolds[0, 0] = reynolds_components[0]
                reynolds[1, 1] = reynolds_components[1]
                reynolds[2, 2] = reynolds_components[2]
                reynolds[0, 1] = reynolds[1, 0] = reynolds_components[3]
                reynolds[0, 2] = reynolds[2, 0] = reynolds_components[4]
                reynolds[1, 2] = reynolds[2, 1] = reynolds_components[5]
                stats_offset = n_base + 6
                scalar_mean = means[3] if cfg.include_scalar else None
                scalar_flux = (
                    np.asarray(statistics[stats_offset:stats_offset + 3])
                    if cfg.include_scalar else None
                )

                local_k = 0.5 * sum(reynolds[j, j].sum() for j in range(3))
                global_k = comm.allreduce(float(local_k), op=MPI.SUM)
                pod_fields = list(means)
                nut = None
                if cfg.closure_enabled:
                    gradients = get_gradient(mesh, coef, means[:3])
                    gradients = np.moveaxis(gradients, (-2, -1), (0, 1))
                    nut, nut_diag = fit_boussinesq_nut(
                        reynolds,
                        gradients + np.swapaxes(gradients, 0, 1),
                        regularisation=cfg.closure_regularisation,
                        maximum=cfg.closure_nut_max,
                        clip_negative=True,
                    )
                    pod_fields.append(nut * cfg.closure_field_scale)
                    if rank == 0:
                        closure_rows.append((float(tcur), nut_diag))

                if cfg.scalar_closure_enabled:
                    scalar_gradient = get_gradient(mesh, coef, [scalar_mean])
                    scalar_gradient = np.moveaxis(
                        scalar_gradient, (-2, -1), (0, 1)
                    )[0]
                    if cfg.scalar_closure == "flux_least_squares":
                        gammat, gammat_diag = fit_scalar_gammat(
                            scalar_flux,
                            scalar_gradient,
                            regularisation=cfg.gammat_regularisation,
                            maximum=cfg.gammat_max,
                            clip_negative=True,
                        )
                    else:
                        if nut is None:
                            raise RuntimeError(
                                "scalar_closure=schmidt requires fit_eddy_viscosity"
                            )
                        gammat = nut / cfg.turbulent_schmidt
                        gammat_diag = {
                            "negative_fraction_raw": 0.0,
                            "clipped_fraction": 0.0,
                            "valid_fraction": 1.0,
                            "mean_fit_relative_residual": 0.0,
                            "max_gammat": float(np.max(gammat)),
                        }
                    pod_fields.append(gammat * cfg.gammat_field_scale)
                    if rank == 0:
                        scalar_closure_rows.append((float(tcur), gammat_diag))

                add_snapshot(
                    comm,
                    pod,
                    ioh,
                    bm,
                    pod_fields,
                    tcur,
                    times,
                    energy,
                    cfg.pod_field_weights,
                )
                if rank == 0:
                    window_metadata.append({
                        "window": len(window_metadata),
                        "center": float(tcur),
                        "end": float(tcur) + 0.5 * cfg.window_duration,
                        "turbulent_kinetic_energy_local_sum": float(global_k),
                    })
                continue

            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                if not times:
                    raise RuntimeError(
                        "URANS-POD produced no complete statistics windows."
                    )
                if rank == 0:
                    with open("window_statistics_manifest.json", "w", encoding="utf-8") as handle:
                        json.dump(
                            {
                                "type": "urans_pod",
                                "window_duration": cfg.window_duration,
                                "design_id": case_cfg.get("design_id", None),
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
                                "windows": window_metadata,
                            },
                            handle,
                            indent=2,
                        )
                    if closure_rows:
                        with open("closure_diagnostics.csv", "w", encoding="utf-8") as handle:
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
                    if scalar_closure_rows:
                        with open("scalar_closure_diagnostics.csv", "w", encoding="utf-8") as handle:
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
                write_time_coefficients(comm, coeffs, header, cfg.write_modes)
                pod_iteration += 1
                write_singular_values(
                    comm,
                    pod,
                    cfg.keep_modes,
                    pod_iteration,
                    singular_value_history_new,
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
                    cfg,
                    zero_field,
                    n_avail,
                    mode_output_index,
                )
                continue

            if mode == MODE_ADJOINT and phase == PHASE_ADJ_DONE:
                pod, ioh = make_pod(comm, bm, cfg)
                times = []
                energy = EnergyState()
                window_metadata = []
                closure_rows = []
                scalar_closure_rows = []
                continue

            time.sleep(0.001)

    finally:
        if rank == 0 and ctrl is not None:
            ctrl.close()
        comm.Barrier()
        ds.finalize()


if __name__ == "__main__":
    main()
