# Low-Re POD, Two Nodes {#low-re-pod-2nodes}

This example is the two-node POD state-recovery low-Re mixer case used on
LUMI-G. It keeps the original 3:1:1 mesh aspect ratio but uses a 60x20x20
mesh, giving 24,000 elements versus 12,288 in `low_Re`.

For manual testing outside Slurm, launch it from the repository root with
`./run.sh low_Re_2nodes`, overriding `RUN_NODES`,
`NEKO_RANKS_PER_NODE`, and `PY_RANKS_PER_NODE` if needed.

For a LUMI submission, run `./run.sh --submit LUMI-G low_Re_2nodes`. The
example-specific job script under `scripts/jobscripts/LUMI-G/low_Re_2nodes`
sets up the two-node full-node layout used by the coupled run:

- 8 Neko ranks per node, one GPU-backed rank per MI250x GCD
- 48 Python ranks per node, using the remaining CPU-only tasks

That job path generates `select_gpu` and `mpmd.conf` automatically before
calling `srun --multi-prog`, so the cluster-specific placement stays with the
example rather than in the generic MPMD helpers.
