# POD Mixer {#pod-mixer}

This example is the POD state-recovery mixer case derived from the
`low_Re` setup used on LUMI. It is the cluster-oriented POD example in this
split stack.

For manual testing outside Slurm, launch it from the repository root with
`./run.sh POD_mixer`, overriding `NEKO_RANKS` and `PY_RANKS` if needed.

For a LUMI submission, run `./run.sh --submit LUMI-G POD_mixer`. The
example-specific job script under `scripts/jobscripts/LUMI-G/POD_mixer` sets up
the full-node layout used by the coupled run:

- 8 Neko ranks per node, one GPU-backed rank per MI250x GCD
- 48 Python ranks per node, using the remaining CPU-only tasks

That job path generates `select_gpu` and `mpmd.conf` automatically before
calling `srun --multi-prog`, so the cluster-specific placement stays with the
example rather than in the generic MPMD helpers.
