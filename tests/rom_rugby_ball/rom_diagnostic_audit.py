#!/usr/bin/env python3
"""Inspect a ROM diagnostic audit without rerunning Neko.

The Fortran collector writes gathered little-endian float64 vectors in
``vectors/``.  The vector order is the concatenation of local nodal vectors in
``NEKO_COMM`` rank order; ``vector_partition.csv`` records the offsets.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def vector(root: Path, name: str) -> np.ndarray:
    path = root / "vectors" / name
    return np.fromfile(path, dtype="<f8")


def finite_float(value: str) -> float:
    return float(value.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit_directory", type=Path)
    parser.add_argument("--plot", action="store_true", help="write audit.png")
    args = parser.parse_args()
    root = args.audit_directory.resolve()

    config = json.loads((root / "configuration.json").read_text())
    iteration = rows(root / "iteration_diagnostics.csv")
    gradients = rows(root / "gradient_comparison.csv")
    finite_difference = rows(root / "directional_fd.csv")
    objectives = rows(root / "objective_comparison.csv")

    selected_errors = []
    for row in iteration:
        selected = vector(root, row["gradient_after"])
        projected = vector(root, row["gradient_projected_after"])
        selected_errors.append(float(np.linalg.norm(selected - projected)) /
                              max(1.0, float(np.linalg.norm(projected))))

    summary = {
        "configuration": config,
        "iteration_rows": len(iteration),
        "gradient_rows": len(gradients),
        "finite_difference_rows": len(finite_difference),
        "objective_rows": len(objectives),
        "max_selected_projected_vector_error": max(selected_errors, default=None),
        "rom_objective_changes": [
            finite_float(row["rom_objective_change"]) for row in iteration
        ],
        "mma_directional_products": [
            finite_float(row["g_mma_dot_update"]) for row in iteration
        ],
        "gradient_comparison": gradients,
        "directional_fd": finite_difference,
        "objective_comparison": objectives,
    }
    (root / "diagnostic_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"audit: {root}")
    print(f"source: {config.get('source_revision')} ({config.get('source_worktree')})")
    print(f"formulation/source: {config.get('formulation')} / "
          f"{config.get('sensitivity_source')}")
    print(f"evaluated MMA transitions: {len(iteration)}")
    print("ROM objective changes:", summary["rom_objective_changes"])
    print("g_MMA dot update:", summary["mma_directional_products"])
    print("max selected-vs-projected gradient error:",
          summary["max_selected_projected_vector_error"])
    for row in objectives:
        print(row["point"], "ROM", row["rom_total"], "FOM", row["full_order_total"],
              "ROM-FOM", row["rom_minus_full"])

    if not args.plot:
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; summary JSON was still written")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    inner = np.array([int(row["inner_iteration"]) for row in iteration])
    before = np.array([finite_float(row["rom_objective_before"]) for row in iteration])
    after = np.array([finite_float(row["rom_objective_after"]) for row in iteration])
    axes[0, 0].plot(inner, before, "o-", label="before MMA")
    axes[0, 0].plot(inner, after, "s-", label="after MMA")
    axes[0, 0].set(xlabel="inner iteration", ylabel="ROM objective")
    axes[0, 0].legend()

    mma_dot = np.array([finite_float(row["g_mma_dot_update"]) for row in iteration])
    exact_dot = np.array([finite_float(row["g_exact_dot_update"]) for row in iteration])
    projected_dot = np.array([finite_float(row["g_projected_dot_update"]) for row in iteration])
    actual = after - before
    axes[0, 1].plot(inner, actual, "o-", label="actual ROM change")
    axes[0, 1].plot(inner, mma_dot, "s-", label=r"$g_{MMA}^T\Delta x$")
    axes[0, 1].plot(inner, exact_dot, "^-", label=r"$g_G^T\Delta x$")
    axes[0, 1].plot(inner, projected_dot, "v-", label=r"$g_P^T\Delta x$")
    axes[0, 1].set(xlabel="inner iteration", ylabel="directional value")
    axes[0, 1].legend(fontsize=8)

    points = [row["point"] for row in gradients]
    exact_error = [finite_float(row["relative_exact_full"]) for row in gradients]
    projected_error = [finite_float(row["relative_projected_full"]) for row in gradients]
    axes[1, 0].bar(np.arange(len(points)) - 0.18, exact_error, 0.36,
                   label="exact Galerkin vs FOM")
    axes[1, 0].bar(np.arange(len(points)) + 0.18, projected_error, 0.36,
                   label="projected adjoint vs FOM")
    axes[1, 0].set_xticks(np.arange(len(points)), points)
    axes[1, 0].set_ylabel("relative gradient error")
    axes[1, 0].legend(fontsize=8)

    fd = [row for row in finite_difference if row["feasible"].upper() == "T"]
    for name in sorted({row["direction"] for row in fd}):
        subset = [row for row in fd if row["direction"] == name]
        eps = np.array([finite_float(row["epsilon"]) for row in subset])
        err_g = np.abs([finite_float(row["fd_minus_exact"]) for row in subset])
        err_p = np.abs([finite_float(row["fd_minus_projected"]) for row in subset])
        axes[1, 1].loglog(eps, err_g, "o-", label=f"{name}: exact")
        axes[1, 1].loglog(eps, err_p, "s--", label=f"{name}: projected")
    axes[1, 1].set(xlabel="epsilon", ylabel="absolute FD derivative error")
    axes[1, 1].legend(fontsize=7)
    fig.savefig(root / "audit.png", dpi=160)


if __name__ == "__main__":
    main()
