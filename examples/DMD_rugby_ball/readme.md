# DMD rugby ball

This is the streaming dynamic-mode-decomposition example. It uses the
separate `type: "dmd"` state-recovery implementation and forms consecutive
snapshot matrices from the initial state and the streamed forward states.

Run it through the repository root runner:

```bash
./run.sh DMD_rugby_ball -c
```

The Python driver maintains two streaming PODs, computes exact DMD modes from
their mass-weighted reduced data, and sends real-valued modes back to Neko.
Complex-conjugate DMD pairs are represented by paired real and imaginary modes.
The `dmd_rank` value must not exceed `n_modes`; `n_modes` is the number of real
mode fields sent back to Neko.

The DMD spectrum is written to `dmd_spectral.csv`; the eigenvalue summary is
written to `dmd_eigenvalues.csv`, and the evaluated modal dynamics are written
to `dmd_modal_dynamics.csv`. Before each new decomposition, existing files are
archived with the POD-style suffixes `_001`, `_002`, and so on. For example,
`dmd_spectral_001.csv` is the first archived spectrum.
