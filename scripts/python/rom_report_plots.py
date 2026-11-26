#!/usr/bin/env python3
"""Plot recorded ROM/FOM archive diagnostics without inventing comparisons.

Plots are emitted only for scalar or coordinate data actually present in the
archive. The JSON manifest records the data scope, caveats, and unavailable
diagnostics alongside the figures.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return [dict(row) for row in csv.DictReader(stream, skipinitialspace=True)]


def _float(row: dict[str, str], name: str) -> float | None:
    raw = row.get(name)
    if raw is None or not raw.strip():
        return None
    try:
        number = float(raw)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _save(fig, path: Path) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def _plot_objectives(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    points = []
    for run in run_dirs:
        diagnostic = run / "rom_gradient_diagnostics.csv"
        if not diagnostic.is_file():
            continue
        rows = _read_rows(diagnostic)
        for index, row in enumerate(rows):
            row_id = row.get("design_id", "").strip()
            ref_id = row.get("reference_design_id", "").strip()
            candidate_id = row.get("candidate_design_id", "").strip()
            if ref_id and candidate_id:
                identity = "same design" if ref_id == candidate_id else "different designs"
            elif row_id:
                identity = f"design {row_id}"
            else:
                identity = "design identity unverified"
            label = run.name + f" / row {index + 1} ({identity})"
            values = []
            for field in ("rom_reference_objective", "rom_candidate_objective",
                          "full_order_reference_objective", "full_order_candidate_objective"):
                value = _float(row, field)
                values.append(value)
            if any(value is not None for value in values):
                points.append((label, values))
    if not points:
        manifest["unavailable"]["objective_comparison_plot"] = "No scalar ROM/FOM reference/candidate objective rows found."
        return

    labels = [point[0] for point in points]
    x = np.arange(len(points))
    fig, ax = plt.subplots(figsize=(max(8, len(points) * 1.1), 5))
    series = [
        ("ROM reference", 0, "o", "tab:blue"),
        ("ROM candidate", 1, "s", "tab:cyan"),
        ("FOM reference", 2, "o", "tab:orange"),
        ("FOM candidate", 3, "s", "tab:red"),
    ]
    for title, slot, marker, color in series:
        vals = [values[slot] for _, values in points]
        mask = [value is not None for value in vals]
        if any(mask):
            ax.scatter(x[np.asarray(mask)], np.asarray([v for v in vals if v is not None]),
                       label=title, marker=marker, color=color)
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.set_ylabel("objective")
    ax.set_title("Recorded ROM/FOM objective values (each row kept in its source archive)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    path = output / "objective_comparisons.png"
    manifest["plots"].append({
        "type": "objective_comparisons",
        "path": _save(fig, path),
        "data_scope": "Recorded reference/candidate scalar rows; no cross-run join was performed.",
        "identity": "row labels use exported design IDs when available; otherwise identity is unverified",
    })


def _plot_components(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    available = False
    for run in run_dirs:
        paths = sorted(run.glob("rom_output/outer_*/objective_diagnostics_inner_*.csv"))
        rows = []
        for path in paths:
            rows.extend(_read_rows(path))
        fields = ["objective", "scalar_mixing", "viscous_dissipation", "brinkman_dissipation"]
        fields = [field for field in fields if any(_float(row, field) is not None for row in rows)]
        if not rows or not fields:
            continue
        available = True
        iters = [_float(row, "inner_iteration") for row in rows]
        x = np.asarray([i if i is not None else idx for idx, i in enumerate(iters)])
        fig, ax = plt.subplots(figsize=(8, 5))
        for field in fields:
            ys = [_float(row, field) for row in rows]
            mask = np.asarray([y is not None for y in ys])
            if np.any(mask):
                ax.plot(x[mask], np.asarray([y for y in ys if y is not None]), marker="o", label=field)
        ax.set_xlabel("inner iteration")
        ax.set_ylabel("recorded objective value")
        ax.set_title(f"Objective components: {run.name} (recorded rows)")
        ax.grid(True, alpha=0.25)
        ax.legend()
        path = output / f"{run.name}_objective_components.png"
        manifest["plots"].append({
            "type": "objective_components",
            "path": _save(fig, path),
            "data_scope": "Per-run inner objective diagnostics; rows are not joined to FOM designs.",
        })
    if not available:
        manifest["unavailable"]["objective_component_plot"] = "No inner objective component CSV rows found."


def _plot_performance(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    available_timing = False
    available_collectives = False
    for run in run_dirs:
        groups: list[tuple[str, list[Path], list[dict[str, str]]]] = []
        detailed = run / "rom_performance.csv"
        if detailed.is_file():
            groups.append(("rom_performance.csv", [detailed], _read_rows(detailed)))
        inner_paths = sorted(run.glob("rom_output/outer_*/objective_diagnostics_inner_*.csv"))
        if inner_paths:
            inner_rows = []
            for path in inner_paths:
                inner_rows.extend(_read_rows(path))
            groups.append(("objective_diagnostics_inner_*.csv", inner_paths, inner_rows))
        opt = run / "optimization_data.csv"
        if opt.is_file():
            groups.append(("optimization_data.csv", [opt], _read_rows(opt)))

        for file_index, (group_name, source_paths, rows) in enumerate(groups):
            if not rows:
                continue
            fields = list(rows[0])
            timing = [f for f in fields if f.endswith("_seconds")]
            if timing:
                available_timing = True
                x_values = [_float(row, "inner_iteration") for row in rows]
                x = np.asarray([v if v is not None else i for i, v in enumerate(x_values)])
                fig, ax = plt.subplots(figsize=(10, 5))
                for field in timing:
                    ys = [_float(row, field) for row in rows]
                    valid = np.asarray([v is not None for v in ys])
                    if np.any(valid):
                        ax.plot(x[valid], np.asarray([v for v in ys if v is not None]), marker="o", label=field)
                ax.set_xlabel("inner iteration or recorded row")
                ax.set_ylabel("seconds")
                ax.set_title(f"Recorded timing: {run.name} / {group_name}")
                ax.grid(True, alpha=0.25)
                ax.legend(fontsize="small", ncol=2)
                name = f"{run.name}_timing_{file_index:02d}.png"
                manifest["plots"].append({
                    "type": "timing",
                    "path": _save(fig, output / name),
                    "data_scope": [str(path) for path in source_paths],
                    "note": "Timing scopes are kept separate by CSV source.",
                })

            count_fields = [f for f in fields if f.endswith("_calls") or f.endswith("_elements") or f.endswith("_bytes")]
            if count_fields:
                available_collectives = True
                groups = {
                    "calls": [f for f in count_fields if f.endswith("_calls")],
                    "elements": [f for f in count_fields if f.endswith("_elements")],
                    "bytes": [f for f in count_fields if f.endswith("_bytes")],
                }
                active_groups = [(name, fs) for name, fs in groups.items() if fs]
                fig, axes = plt.subplots(len(active_groups), 1, figsize=(10, 3.0 * len(active_groups)), squeeze=False)
                x_values = [_float(row, "inner_iteration") for row in rows]
                x = np.asarray([v if v is not None else i for i, v in enumerate(x_values)])
                for axis, (unit, unit_fields) in zip(axes[:, 0], active_groups):
                    for field in unit_fields:
                        ys = [_float(row, field) for row in rows]
                        valid = np.asarray([v is not None for v in ys])
                        if np.any(valid):
                            axis.plot(x[valid], np.asarray([v for v in ys if v is not None]), marker="o", label=field)
                    axis.set_ylabel(unit)
                    axis.grid(True, alpha=0.25)
                    axis.legend(fontsize="small")
                axes[-1, 0].set_xlabel("inner iteration or recorded row")
                fig.suptitle(f"Recorded collective telemetry: {run.name} / {group_name}")
                name = f"{run.name}_collectives_{file_index:02d}.png"
                manifest["plots"].append({
                    "type": "collective_telemetry",
                    "path": _save(fig, output / name),
                    "data_scope": [str(path) for path in source_paths],
                })
    if not available_timing:
        manifest["unavailable"]["timing_plot"] = "No measured *_seconds columns found."
    if not available_collectives:
        manifest["unavailable"]["collective_plot"] = "No recorded *_calls, *_elements or *_bytes columns found."


def _plot_ranks(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    labels, values = [], []
    for run in run_dirs:
        for path in sorted(run.glob("rom_output/outer_*/basis_metadata.json")):
            try:
                metadata = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            label = run.name + "/" + path.parent.parent.name
            row = [metadata.get(key) for key in (
                "forward_candidate_modes", "forward_retained_rank",
                "adjoint_candidate_modes", "adjoint_retained_rank")]
            if any(isinstance(value, (int, float)) for value in row):
                labels.append(label)
                values.append(row)
    if not values:
        manifest["unavailable"]["rank_transition_plot"] = "No recorded candidate/retained rank pairs found."
        manifest["unavailable"]["complete_rank_transition_plot"] = "Intermediate stage ranks are not available in current metadata."
        return
    x = np.arange(len(labels))
    width = 0.19
    names = ["F candidate", "F retained", "A candidate", "A retained"]
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.2), 5))
    for col, name in enumerate(names):
        ys = [row[col] for row in values]
        mask = np.asarray([isinstance(value, (int, float)) for value in ys])
        if np.any(mask):
            offsets = x[mask] + (col - 1.5) * width
            ax.bar(offsets, np.asarray([value for value in ys if isinstance(value, (int, float))]), width, label=name)
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylabel("rank")
    ax.set_title("Recorded candidate and retained basis ranks (intermediate stages unavailable)")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    manifest["plots"].append({
        "type": "rank_transition_subset",
        "path": _save(fig, output / "recorded_rank_pairs.png"),
        "data_scope": "Metadata candidate/retained ranks only; not a complete Python-to-Neko transition.",
    })
    manifest["unavailable"]["complete_rank_transition_plot"] = (
        "Requested POD, actual POD, historical candidate, Python recycled, Neko incoming, boundary/continuity, and final orthonormal ranks are not all exported."
    )


def _plot_divergence(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    produced = 0
    for run in run_dirs:
        for path in sorted(run.glob("rom_output/outer_*/basis_diagnostics/*.csv")):
            rows = _read_rows(path)
            if not rows or not any("strong_divergence_before" in row for row in rows):
                continue
            columns = ["strong_divergence_before", "strong_divergence_after",
                       "weak_continuity_before", "weak_continuity_after",
                       "bc_residual_before", "bc_residual_after"]
            fields = [name for name in columns if any(_float(row, name) is not None for row in rows)]
            if not fields:
                continue
            mode = [_float(row, "mode") for row in rows]
            x = np.asarray([value if value is not None else i + 1 for i, value in enumerate(mode)])
            fig, ax = plt.subplots(figsize=(9, 5))
            for field in fields:
                ys = [_float(row, field) for row in rows]
                mask = np.asarray([value is not None for value in ys])
                if np.any(mask):
                    ax.plot(x[mask], np.asarray([value for value in ys if value is not None]), marker="o", label=field)
            ax.set_xlabel("mode row (family not labelled in source CSV)")
            ax.set_ylabel("recorded residual norm")
            ax.set_title(f"Per-mode divergence/continuity/BC: {run.name} / {path.parent.parent.name}")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize="small", ncol=2)
            safe = f"{run.name}_divergence_{produced:02d}.png"
            manifest["plots"].append({
                "type": "divergence_and_boundary_residuals",
                "path": _save(fig, output / safe),
                "data_scope": str(path),
                "note": "Per-file modes are plotted only; the existing filename does not identify forward versus adjoint family.",
            })
            produced += 1
    if produced == 0:
        manifest["unavailable"]["divergence_plot"] = "No per-mode divergence/continuity/BC CSV values found."


def _plot_coefficients(run_dirs: list[Path], output: Path, manifest: dict[str, Any]) -> None:
    produced = 0
    for run in run_dirs:
        for path in sorted(run.glob("rom_output/outer_*/forward_coefficients_inner_*.csv")) + sorted(
            run.glob("rom_output/outer_*/adjoint_coefficients_inner_*.csv")
        ):
            rows = _read_rows(path)
            if not rows or "time" not in rows[0]:
                continue
            mode_columns = [name for name in rows[0] if name.startswith("a_") or name.startswith("lambda_")]
            if not mode_columns:
                continue
            times = [_float(row, "time") for row in rows]
            if any(value is None for value in times):
                continue
            fig, ax = plt.subplots(figsize=(10, 5))
            plotted = 0
            for column in mode_columns:
                values = [_float(row, column) for row in rows]
                if any(value is None for value in values):
                    continue
                ax.plot(times, values, label=column)
                plotted += 1
            if not plotted:
                plt.close(fig)
                continue
            ax.set_xlabel("time")
            ax.set_ylabel("integrated ROM coefficient")
            ax.set_title(f"{run.name}: {path.name} (single candidate coordinate system)")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize="small", ncol=2)
            out = output / f"{run.name}_coefficients_{produced:03d}.png"
            manifest["plots"].append({
                "type": "integrated_rom_coefficient_trajectory",
                "path": _save(fig, out),
                "data_scope": str(path),
                "note": "Each candidate trajectory is plotted separately; coefficients are not compared across changing bases.",
            })
            produced += 1
    if produced == 0:
        manifest["unavailable"]["coefficient_trajectory_plot"] = "No compatible integrated ROM coefficient trajectories found."
    else:
        manifest["unavailable"]["projected_fom_coefficient_plot"] = (
            "No projected-FOM coefficient trajectories found; integrated ROM coefficients are not substitutes."
        )


def create_report(run_paths: list[Path], output_dir: Path) -> dict[str, Any]:
    runs = [path.resolve() for path in run_paths]
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "run_directories": [str(path) for path in runs],
        "plots": [],
        "unavailable": {
            "gradient_vector_alignment_plot": "No per-design gradient vectors or alignment matrix were found by the archive audit.",
            "projected_fom_trajectory_plot": "No projected-FOM coefficient coordinates were found; integrated ROM coefficients are not used as a substitute.",
        },
        "comparison_policy": "No data are joined across runs without matching exported design IDs. Candidate rows without IDs are plotted as reported and labelled unverified.",
    }
    _plot_objectives(runs, output_dir, manifest)
    _plot_components(runs, output_dir, manifest)
    _plot_performance(runs, output_dir, manifest)
    _plot_ranks(runs, output_dir, manifest)
    _plot_divergence(runs, output_dir, manifest)
    _plot_coefficients(runs, output_dir, manifest)
    (output_dir / "rom_report_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="+", type=Path, help="one or more ROM result directories")
    parser.add_argument("--output-dir", required=True, type=Path, help="directory for PNG plots and manifest")
    args = parser.parse_args()
    manifest = create_report(args.run, args.output_dir)
    print(f"Wrote {len(manifest['plots'])} plots and {args.output_dir / 'rom_report_manifest.json'}")


if __name__ == "__main__":
    main()
