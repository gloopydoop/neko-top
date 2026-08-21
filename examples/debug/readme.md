Debug batch example for the split POD mixer branch.

Run it with:

```bash
./run.sh --submit LUMI-G debug
```

The batch job executes a sequence of MPI, Python, and Neko smoke tests and
writes per-step logs plus `summary.txt` under `debug_artifacts/`.

After the framework cleanup completes, the archived files are in:

```text
results/debug/
```
