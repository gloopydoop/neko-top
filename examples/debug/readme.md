Debug batch example for the split POD mixer branch.

Run it with:

```bash
./run.sh --submit LUMI-G debug
```

This variant is trimmed to the remaining high-value checks:

- Python-only MPMD sanity.
- Current `./neko` binary fingerprint and launch behavior.
- Optional control binaries run under the same batch harness.
- LUMI-style 8-rank GPU launch variants inspired by the Neko LUMI discussion.
- Shared MPMD tests on short smoke cases for `no_pod` and, when that passes,
  the real POD driver.

The runner writes per-step logs plus `summary.txt` under `debug_artifacts/`.
It generates short smoke cases at runtime so healthy runs can finish instead of
always timing out.

Optional environment variables:

```bash
export DEBUG_WORKING_NEKO=/path/to/known-good/neko
export DEBUG_NODEVICE_NEKO=/path/to/rebuild/without/device_mpi/neko
export DEBUG_BAD_NEKO=/path/to/known-bad/neko
export DEBUG_SMOKE_STEPS=1
export DEBUG_CPUS_PER_TASK=7
export DEBUG_LUMI_GPU_MAP=4,5,2,3,6,7,0,1
export DEBUG_LUMI_CPU_BIND=cores
```

If `DEBUG_WORKING_NEKO` or `DEBUG_NODEVICE_NEKO` are set and accessible from
the batch node, the debug run executes the same focused launch checks with
those binaries too.

After the framework cleanup completes, the archived files are in:

```text
results/debug/
```
