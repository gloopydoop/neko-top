# Reduced-system replay diagnostic

`ROM_low_Re_mini` can export the exact reduced system used by its BDF1
integrator. The case enables this with:

```json
"export_reduced_system": true,
"reduced_system_export_directory": "rom_reduced_system"
```

After the first ROM evaluation at the unchanged reference design, the export
directory contains `metadata.json`, `H.csv`, `b.csv`,
`linear_explicit.csv`, `linear_implicit.csv`, `Q.csv`, `a0.csv`, the split
lift/Brinkman files, and the final trial/test bases as NumPy-readable raw
`float64` binaries. If that ROM pass reaches the end of its trajectory,
`reference_trajectory.csv` is also written.

Replay it without starting Neko or DNS:

```bash
python tests/rom_rugby_ball/rom_reduced_system_diagnostic.py \
    /path/to/rom_reduced_system
```

The script writes `replayed_trajectory.csv`,
`linearised_eigenvalues.csv`, `diagnostic_summary.json`, and (when
Matplotlib is installed) `reduced_system_diagnostic.png`.

Optional DNS comparisons are supported without a new DNS run:

```bash
python tests/rom_rugby_ball/rom_reduced_system_diagnostic.py \
    /path/to/rom_reduced_system \
    --dns-coefficients dns_projected_coefficients.csv
```

The coefficient CSV must contain `step,time,a_1,...,a_r`. Alternatively,
`--dns-fields dns_snapshots.npz` accepts `fields` with shape
`(time,field,dof)`, or `velocity` with shape `(time,3,dof)` plus an optional
`scalar` array. Those fields are lift-subtracted and projected with the
exported final trial basis and trial mass matrix.

The replay uses the exact sign convention and indexing exported by Neko:

```text
(H - dt*linear_implicit) a[n+1]
    = H a[n] + dt*(b + linear_explicit*a[n] + q(a[n]))
q_i(a) = a[j] * Q[i,j,k] * a[k]
```
