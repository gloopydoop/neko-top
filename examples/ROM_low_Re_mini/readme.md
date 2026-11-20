# ROM_low_Re_mini

This is the minimum coupled velocity/passive-scalar ROM experiment. It uses
the real Neko-TOP fluid, scalar transport, viscous dissipation, Brinkman
dissipation, scalar-mixing objective, design mapping and `rom_mma` lifecycle.
The POD basis is forward-only (`stream_adjoint=false`); no adjoint columns are
used for enrichment. The Neko-side admissibility diagnostics remain active.

From the repository root:

```bash
source prepare.env
./run.sh ROM_low_Re_mini -c
```

The normal launcher prepares `mini_box.nmsh`, builds the user driver and
starts one Neko rank plus one Python ADIOS2 rank by default. The Python
environment is taken from `prepare.env` when it is present.

The main experimental controls are in `ROM_low_Re_mini.case`:

* `state_recovery.rank_forward` and `rom.rank_forward`: forward POD rank;
* `rom.inner_iterations`: reduced MMA iterations per refresh;
* `optimization.solver.mma.move_limit`: MMA move limit;
* `state_recovery.history_outer_iterations`: recycling history length;
* `rom.scalar_metric_weight`: scalar contribution to the coupled state metric.

The run writes the standard optimizer CSV plus ROM forward/adjoint trajectory,
design and sensitivity fields in the archived result directory. This is an
engineering experiment unless the numerical validation records explicitly
show that the admissibility, objective, gradient and full-order gates pass.
