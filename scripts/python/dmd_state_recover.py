#!/usr/bin/env python3

"""Streaming DMD state recovery for Neko-TOP."""

import os
import time

import numpy as np
from mpi4py import MPI

from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.adios2.stream import DataStreamer

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
from dmd_state_recover_tools import DMDConfig
from dmd_state_recover_tools import add_snapshot
from dmd_state_recover_tools import compute_dmd
from dmd_state_recover_tools import flush_buffer
from dmd_state_recover_tools import load_dmd_config
from dmd_state_recover_tools import make_pod
from dmd_state_recover_tools import recv_field
from dmd_state_recover_tools import recv_fields
from dmd_state_recover_tools import set_debug
from dmd_state_recover_tools import write_modes_to_disk
from dmd_state_recover_tools import write_spectral_data


def init_runtime(comm: MPI.Comm, cfg: DMDConfig):
    ds = DataStreamer(comm)
    x = recv_field(ds, cfg.dtype)
    y = recv_field(ds, cfg.dtype)
    z = recv_field(ds, cfg.dtype)
    initial_fields = recv_fields(ds, cfg.n_input_fields, cfg.dtype)
    mesh = Mesh(comm, x=x, y=y, z=z, create_connectivity=False)
    coef = Coef(mesh, comm)
    pod_x, ioh_x = make_pod(comm, coef.B, cfg)
    pod_y, ioh_y = make_pod(comm, coef.B, cfg)
    return ds, mesh, coef, pod_x, ioh_x, pod_y, ioh_y, initial_fields


def stream_modes(ds, modes: np.ndarray, cfg: DMDConfig) -> None:
    field_size = modes.shape[0] // cfg.n_fields
    for mode_index in range(cfg.keep_modes):
        flat_mode = modes[:, mode_index]
        for field_index in range(cfg.n_fields):
            start = field_index * field_size
            end = (field_index + 1) * field_size
            ds.stream(flat_mode[start:end].astype(cfg.dtype, copy=False))


def main() -> None:
    if len(os.sys.argv) < 2:
        raise SystemExit(
            f"Usage: {os.path.basename(os.sys.argv[0])} <case-file>"
        )

    world = MPI.COMM_WORLD
    comm = make_local_comm(world)
    rank = comm.Get_rank()
    peer_root = get_peer_root()
    case_path = os.path.abspath(os.sys.argv[1])
    if not os.path.isfile(case_path):
        raise FileNotFoundError(f"Case file not found: {case_path}")

    case_cfg, cfg = load_dmd_config(case_path)
    set_debug(cfg.debug)
    (
        ds,
        mesh,
        coef,
        pod_x,
        ioh_x,
        pod_y,
        ioh_y,
        initial_fields,
    ) = init_runtime(comm, cfg)

    ctrl = None
    if rank == 0:
        ctrl = CtrlClient(world, peer_root, cfg.debug)

    previous_fields = initial_fields
    previous_time = 0.0
    x_times = []
    my_count = getattr(ds, "py2f_field_my_count", None)
    if my_count is None:
        my_count = initial_fields[0].size
    zero_field = np.zeros(my_count, dtype=cfg.dtype)
    mode_output_index = 0

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
                add_snapshot(comm, pod_x, ioh_x, previous_fields)
                add_snapshot(comm, pod_y, ioh_y, fields)
                x_times.append(float(previous_time))
                previous_fields = fields
                previous_time = float(tcur)
                continue

            if mode == MODE_FORWARD and phase == PHASE_FWD_DONE:
                flush_buffer(comm, pod_x, ioh_x)
                flush_buffer(comm, pod_y, ioh_y)
                pod_x.scale_modes(comm, bm1sqrt=ioh_x.bm1sqrt, op="div")
                pod_y.scale_modes(comm, bm1sqrt=ioh_y.bm1sqrt, op="div")
                modes, spectral_data, eigenvalues, amplitudes = compute_dmd(
                    comm,
                    pod_x,
                    pod_y,
                    coef.B,
                    initial_fields,
                    x_times,
                    cfg,
                )
                write_spectral_data(
                    comm, spectral_data, eigenvalues, amplitudes, x_times
                )
                if rank == 0:
                    ctrl.send_cmd(MODE_ADJOINT, PHASE_ADJ_RUNNING)

                stream_modes(ds, modes, cfg)
                if cfg.write_modes:
                    mode_output_index = write_modes_to_disk(
                        comm,
                        mesh,
                        modes,
                        cfg,
                        mode_output_index,
                    )
                continue

            if mode == MODE_ADJOINT and phase == PHASE_ADJ_DONE:
                pod_x, ioh_x = make_pod(comm, coef.B, cfg)
                pod_y, ioh_y = make_pod(comm, coef.B, cfg)
                previous_fields = initial_fields
                previous_time = 0.0
                x_times = []
                continue

            time.sleep(0.001)
    finally:
        if rank == 0 and ctrl is not None:
            ctrl.close()
        comm.Barrier()
        ds.finalize()


if __name__ == "__main__":
    main()
