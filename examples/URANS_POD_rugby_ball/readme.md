# URANS/POD rugby ball comparison

This is an opt-in statistical-POD comparison case using the separate
`type: "urans_pod"` recovery implementation. The standard `type: "pod"`
implementation is unchanged. Fortran accumulates the velocity/scalar means,
Reynolds-stress moments, and scalar-flux moments over each fixed window and
sends exactly one statistics record to Python at the end of that window.
Python fits the closure fields from that record and submits one POD snapshot
containing the window means plus the fitted closure fields.

Run it through the repository root runner:

```bash
./run.sh URANS_POD_rugby_ball -c
```

The legacy instantaneous comparison remains available as
`examples/POD_rugby_ball`. Window metadata is written to
`window_statistics_manifest.json`, and closure-fit diagnostics to
`closure_diagnostics.csv` in the run directory. The recovered eddy-viscosity
mode is registered as `adjoint_nut` and is consumed by the full-stress adjoint
operator. Closure fields use an epsilon-sized inner-product weight, so they
remain reconstructible without materially contributing to POD energy; the
velocity/scalar fields use unit weight. The fluid and optimization setup is
otherwise aligned with the standard POD rugby-ball case.

For a scalar-enabled case, add `scalar_closure: "flux_least_squares"` under
`state_recovery`, together with `gammat_regularisation` and optional
`gammat_max`; use `scalar_closure: "schmidt"` to select the explicitly
labelled `Gamma_T = nu_T/Sc_T` fallback and set `turbulent_schmidt`.
