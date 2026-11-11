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
from pod_state_recover_tools import EnergyState
from pod_state_recover_tools import PODConfig
from pod_state_recover_tools import add_snapshot
from pod_state_recover_tools import build_time_coefficients
from pod_state_recover_tools import flush_buffer
from pod_state_recover_tools import load_pod_config
from pod_state_recover_tools import make_pod
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
                fields = recv_fields(ds, cfg.n_input_fields, cfg.dtype)
                if window_stats is None:
                    add_snapshot(comm, pod, ioh, bm, fields, tcur, times, energy)
                else:
                    window_stats.add(fields, tcur)
                continue

            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                if window_stats is not None:
                    windows = window_stats.finish()
                    if not windows:
                        raise RuntimeError(
                            "windowed_statistics produced no valid windows; "
                            "check statistics_start_time, window_duration and run length"
                        )
                    if rank == 0:
                        with open("window_statistics_manifest.json", "w", encoding="utf-8") as handle:
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
                        with open("window_statistics_diagnostics.csv", "w", encoding="utf-8") as handle:
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
                            comm, pod, ioh, bm, pod_fields, item["center"], times, energy
                        )
                    if rank == 0:
                        print(
                            f"Windowed statistics: {len(windows)} windows, "
                            f"samples={[item['samples'] for item in windows]}",
                            flush=True,
                        )
                    if rank == 0 and closure_rows:
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
                    if rank == 0 and scalar_closure_rows:
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
                write_time_coefficients(
                    comm, coeffs, header, cfg.write_modes
                )
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
                if window_stats is None:
                    add_snapshot(comm, pod, ioh, bm, initial_fields, 0.0, times, energy)
                else:
                    window_stats = WindowStatistics(
                        cfg.window_duration,
                        cfg.statistics_start_time,
                        cfg.dtype,
                        cfg.keep_partial_window,
                    )
                continue

            time.sleep(0.001)

    finally:
        if rank == 0 and ctrl is not None:
            ctrl.close()
        comm.Barrier()
        ds.finalize()


if __name__ == "__main__":
    main()
