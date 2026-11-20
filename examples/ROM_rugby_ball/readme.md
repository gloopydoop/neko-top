# ROM Rugby Ball

## Current status

There are two deliberately separate paths:

* The normal launcher is the real Neko/Python integration path. It keeps the
  admissibility and numerical-correctness diagnostics active. By default the
  case uses an explicit engineering warning fallback when the available basis
  fails admissibility; set `rom.strict_admissibility` to `true` to make that
  gate fail closed.
* The engineering wrapper below runs the reduced lifecycle with the supplied
  forward POD field so that the trajectory, real objectives, MMA update,
  full-order refresh and output contract can be inspected. It intentionally
  does not call the production basis-admissibility handoff and is **not
  scientifically validated**.

For the available forward POD field, the independent continuity comparison is
recorded in
`results/ROM_rugby_ball/forward_continuity_diagnostics.csv`. The original
snapshots have assembled weak-continuity residuals from `8.85395e-3` down to
`6.56370e-3`; after lift homogenisation they are `7.86043e-3`--`7.80413e-3`.
The restricted POD modes are smaller but nonzero (`5.22491e-5`--`7.65418e-4`).
Thus the basis does not introduce the whole discrepancy, but it also does not
pass the strict assembled constraint. This is distinct from the SEM strong
divergence diagnostic and remains a scientific blocker.

ROM_rugby_ball is the Neko-TOP recycled adjoint-enriched POD--Galerkin
example. It follows the same full-order/streaming execution model as
POD_rugby_ball: Neko runs the rugby-ball forward and adjoint solves while an
MPI Python peer receives the fields through ADIOS2.

The ROM case deliberately selects native GLL quadrature (`dealias=false`) and
the explicitly configured Neko-compatible BDF1 reduced time map
(`rom.time_scheme="bdf1"`, matching `case.numerics.time_order=1`). The
physical velocity lift is `[1,0,0]`; snapshots are
stored in the homogeneous lifted space and reconstructed with the lift. The
ROM rejects scalar cases, variable timesteps and dealiased objectives.
Basis diagnostics always report Neko's discrete divergence, weak continuity,
BC and pressure residuals; strict mode rejects a failed gate, while default
engineering mode continues only with an explicit warning and marks the run
uncertified.

The peer leaves physical snapshots unchanged during basis construction. The
common basis is checked against Neko's native essential-boundary and discrete
divergence constraints after transfer; no boundary-node zeroing surrogate is
used.

When the native projector is available, the peer performs:

1. streamed forward and objective-forced adjoint snapshot collection;
2. separate forward and adjoint PODs;
3. FIFO recycling over the configured positive history length;
4. the configured essential-velocity admissibility projection and one common
   mass-orthonormal basis;
5. ROM basis and validation diagnostics.

Neko then runs `rom_mma`: `rom_problem_t` constructs the reduced operators,
runs the forward and backward reduced trajectories, and computes the real
objectives and mapped sensitivities. `rom_optimizer_t` owns the MMA updates,
full-order validation diagnostics and logging. Each outer iteration is a full
forward/adjoint refresh, and `rom.inner_iterations` fixed-basis reduced MMA
updates are made between refreshes. Fixed Galerkin operators are assembled
once per basis; mapped Brinkman operators and exact objective/constraint
correlations are assembled once per candidate. The proposed candidate is
evaluated again by the full-order problem and is always retained. Validation
still reports whether the ROM prediction is scientifically trustworthy, but it
does not roll the engineering trajectory back.

Run it from the repository root:

    source prepare.env
    ./run.sh -c -r ROM_rugby_ball

The `-c` option cleans launcher logs and `-r` forces execution when
`results/ROM_rugby_ball` already exists. On a fresh checkout,
`./run.sh ROM_rugby_ball -c` is sufficient.

The main case controls are in `ROM_rugby_ball.case`:

* `state_recovery.rank_forward` — forward POD rank;
* `rom.inner_iterations` — reduced MMA iterations per refresh;
* `optimization.solver.mma.move_limit` — MMA's per-inner-update move limit;
* `state_recovery.history_outer_iterations` — positive recycling history
  length.

Do not treat the warning-fallback run as scientifically validated. The strict
independent check can be repeated without the lifecycle:

    source prepare.env
    mpirun -n 1 ./build/tests/rom_rugby_ball/ROM_rugby_ball_projection_diagnostic \
      results/ROM_rugby_ball/ROM_rugby_ball.case \
      results/ROM_rugby_ball/POD_modes/POD_modes0.fld

The forward-snapshot versus POD weak-continuity comparison can be repeated
without rebuilding the basis or rerunning DNS/adjoint:

    source prepare.env
    mpirun -n 1 ./build/tests/rom_rugby_ball/ROM_rugby_ball_forward_continuity_diagnostic \
      results/ROM_rugby_ball/ROM_rugby_ball.case \
      logs/POD_rugby_ball/forward_fields0.fld \
      results/ROM_rugby_ball/POD_modes/POD_modes0.fld \
      17 results/ROM_rugby_ball/forward_continuity_diagnostics.csv

## Engineering integration test (not validation)

Run the reproducible forward-only engineering lifecycle through the normal
repository launcher with:

    source prepare.env
    ROM_ENGINEERING=1 ./run.sh ROM_rugby_ball -r -c

The `-r` is useful when a previous result directory exists; on a fresh
checkout, `./run.sh ROM_rugby_ball -c` is sufficient. The direct ROM-specific
wrapper is also available:

    source prepare.env
    ./tests/rom_rugby_ball/run_engineering.sh

The wrapper uses the real case, forward POD field, Neko objective routines,
PDE-filter/RAMP design mapping, reduced forward/adjoint, configured inner MMA
iterations, full-order refresh and validation diagnostics. It runs one short
outer step by default (`HORIZON=0.02`, `OUTER_STEPS=1`) and writes to
`results/ROM_rugby_ball_engineering_forward_only_verified/`. The output is kept
outside `results/ROM_rugby_ball/` because the normal repository launcher
cleans that case directory before archiving.

The successful reference run proposed a genuine MMA move of `8.0e-2`. Even if
the full-order validation reports failed objective or gradient agreement, the
candidate remains active and is used by the next outer iteration. This is
evidence of an executable engineering path, not evidence of an admissible or
scientifically validated ROM.

Useful overrides are shell variables, for example:

    HORIZON=0.02 OUTER_STEPS=2 MPI_RANKS=1 \
      OUT_DIR=results/ROM_rugby_ball/engineering_forward_only_trial \
      ./tests/rom_rugby_ball/run_engineering.sh

Change `state_recovery.rank_forward`, `rom.inner_iterations`,
`optimization.solver.mma.move_limit`, and
`state_recovery.history_outer_iterations` in the case file before running.
The supplied diagnostic field contains one mode, so the engineering hook caps
its supplied-basis rank at one; a larger configured rank requires a
corresponding multi-mode field and the scientific path still requires
admissibility.

The engineering output contract is:

* `engineering_summary.csv` — proposed and retained design changes and elapsed
  time;
* `rom_validation.csv` — predicted/full-order objective and constraint values,
  four validation errors, validation status and candidate-retained status;
* `engineering_rom_status.txt` — explicit non-certification and admissibility
  status;
* `ROM_forward_fields_iter_10.f*` — reduced forward trajectory samples;
* `ROM_design_iter_10.f*` — initial/proposed design fields;
* `ROM_sensitivity_iter_10.f*` — reduced sensitivity fields.

The full-order refresh is evaluated internally for validation and sensitivity;
the CSV is the authoritative comparison. The normal scientific command
remains:

    source prepare.env
    ./run.sh -c -r ROM_rugby_ball

It must not be replaced by the engineering wrapper when making a scientific
claim. The normal path reports the active admissibility gate in
`rom_projection_diagnostics.csv`; with `strict_admissibility=true` it stops at
that gate until the forward/adjoint basis handoff is repaired.

The diagnostic build targets are registered under `tests/rom_rugby_ball`:

    source prepare.env
    cmake --build build --target ROM_rugby_ball_diagnostics -j6

The executable paths are under `build/tests/rom_rugby_ball/`. The repository
launcher reports the status of every case in `logs/`; in this development
checkout its overall status can remain nonzero because of the pre-existing
unrelated `low_Re` error. The `ROM_rugby_ball` engineering case itself is
complete when its output directory and `engineering_summary.csv` are present;
the `low_Re` log is intentionally not modified.

The standard runner builds examples/ROM_rugby_ball/neko, generates box.nmsh
with prepare.sh, and launches the Neko and Python peers with the same MPMD
helper used by POD_rugby_ball. Results are archived under
results/ROM_rugby_ball/.

The monitor files are:

* `results/ROM_rugby_ball/optimization_data.csv`
* `results/ROM_rugby_ball/rom_validation.csv`
* outer full-order outputs: `forward_fields0.f*`, `design0.f*`, and
  `sensitivity0.f*`; the sample counter is the outer iteration
* inner ROM outputs: `ROM_forward_fields_iter_<outer>_0.f*`,
  `ROM_design_iter_<outer>_0.f*`, and `ROM_sensitivity_iter_<outer>_0.f*`
* `results/ROM_rugby_ball/recycled_rom_history.json`
* `results/ROM_rugby_ball/adjoint_pod_time_coeffs*.csv`
* `results/ROM_rugby_ball/adjoint_pod_singular_values*.csv`
* `results/ROM_rugby_ball/adjoint_POD_modes/adjoint_POD_modes0.f*`
