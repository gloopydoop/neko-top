#!/usr/bin/env python3
"""Extract supported objective and performance comparisons from ROM archives.

This reader summarizes recorded scalar diagnostics. It does not reconstruct
spatial gradients, assume two coefficient trajectories share coordinates, or
infer projected-FOM values from integrated-ROM coefficients.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        return [dict(row) for row in reader]


def _number(row: dict[str, str], name: str) -> float | None:
    value = row.get(name)
    if value is None or not value.strip():
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _truth(value: str | None) -> bool | None:
    if value is None or not value.strip():
        return None
    value = value.strip().casefold()
    if value in {"true", "t", ".true.", "1", "yes"}:
        return True
    if value in {"false", "f", ".false.", "0", "no"}:
        return False
    return None


def _objective_comparisons(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"status": "unavailable", "reason": "rom_gradient_diagnostics.csv not found", "rows": []}
    output = []
    for row in _rows(path):
        rom_ref = _number(row, "rom_reference_objective")
        rom_cand = _number(row, "rom_candidate_objective")
        fom_ref = _number(row, "full_order_reference_objective")
        fom_cand = _number(row, "full_order_candidate_objective")
        fom_ref_ok = _truth(row.get("full_order_reference_available")) is True
        if fom_ref is not None and fom_ref < 0.0:
            fom_ref = None
        if fom_cand is not None and fom_cand < 0.0:
            fom_cand = None
        if not fom_ref_ok:
            fom_ref = None
        output.append({
            "outer_iteration": _number(row, "outer_iteration"),
            "sensitivity_source": row.get("sensitivity_source", "").strip() or None,
            "rom_reference_objective": rom_ref,
            "rom_candidate_objective": rom_cand,
            "rom_candidate_minus_reference": (
                rom_cand - rom_ref if rom_cand is not None and rom_ref is not None else None
            ),
            "full_order_reference_objective": fom_ref,
            "full_order_candidate_objective": fom_cand,
            "full_order_candidate_minus_reference": (
                fom_cand - fom_ref if fom_cand is not None and fom_ref is not None else None
            ),
            "reported_predicted_rom_improvement": _number(row, "predicted_rom_improvement"),
            "reported_actual_full_order_improvement": _number(row, "actual_full_order_improvement"),
            "same_design_identity_verified": False,
            "design_identity_missing": ["reference_design_id", "candidate_design_id"],
        })
    return {"status": "available", "path": str(path), "rows": output}


def _inner_objectives(run_dir: Path) -> dict[str, Any]:
    paths = sorted(run_dir.glob("rom_output/outer_*/objective_diagnostics_inner_*.csv"))
    if not paths:
        return {"status": "unavailable", "reason": "inner objective diagnostics not found", "rows": []}
    output = []
    for path in paths:
        for row in _rows(path):
            output.append({
                "path": str(path),
                "outer_iteration": _number(row, "outer_iteration"),
                "inner_iteration": _number(row, "inner_iteration"),
                "objective": _number(row, "objective"),
                "scalar_mixing": _number(row, "scalar_mixing"),
                "viscous_dissipation": _number(row, "viscous_dissipation"),
                "brinkman_dissipation": _number(row, "brinkman_dissipation"),
                "volume_constraint": _number(row, "volume_constraint"),
                "update_status": row.get("update_status", "").strip() or None,
                "design_file": row.get("design_file", "").strip() or None,
                "forward_coefficients_file": row.get("forward_coefficients", "").strip() or None,
                "adjoint_coefficients_file": row.get("adjoint_coefficients", "").strip() or None,
            })
    return {"status": "available", "rows": output}


def _timing_summary(path: Path) -> dict[str, Any]:
    rows = _rows(path)
    fields = list(rows[0]) if rows else []
    values: dict[str, list[float]] = {}
    for row in rows:
        for name in fields:
            if name in {"outer_iteration", "inner_iteration", "phase", "source", "status"}:
                continue
            value = _number(row, name)
            if value is not None:
                values.setdefault(name, []).append(value)
    metrics = {}
    for name, numbers in values.items():
        if name.endswith("_seconds") or name.endswith("_calls") or name.endswith("_elements") or name.endswith("_bytes"):
            metrics[name] = {
                "samples": len(numbers),
                "min": min(numbers),
                "mean": sum(numbers) / len(numbers),
                "max": max(numbers),
            }
    return {"path": str(path), "rows": len(rows), "metrics": metrics}


def _performance(run_dir: Path) -> dict[str, Any]:
    paths: list[Path] = []
    preferred = run_dir / "rom_performance.csv"
    if preferred.is_file():
        paths.append(preferred)
    # Preserve each telemetry scope as its own record instead of pooling
    # unrelated optimizer and inner-candidate timers.
    for path in sorted(run_dir.glob("rom_output/outer_*/objective_diagnostics_inner_*.csv")):
        paths.append(path)
    opt = run_dir / "optimization_data.csv"
    if opt.is_file():
        paths.append(opt)
    if not paths:
        return {
            "status": "unavailable",
            "reason": "no rom_performance.csv, inner timing rows, or optimization_data.csv",
            "measurements": [],
            "detailed_assembly_and_collective_counts": "unavailable",
        }
    measurements = [_timing_summary(path) for path in paths]
    detailed = any(
        "host_collective_calls" in summary["metrics"] for summary in measurements
    )
    return {
        "status": "available",
        "measurements": measurements,
        "detailed_assembly_and_collective_counts": "available" if detailed else "unavailable",
    }


def summarize_run(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    metadata_path = run_dir / "rom_output" / "outer_000" / "basis_metadata.json"
    metadata = None
    if metadata_path.is_file():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            metadata = {"status": "invalid_json", "path": str(metadata_path)}

    gradient_files = sorted(run_dir.glob("*gradient*.csv")) + sorted(
        run_dir.glob("rom_output/outer_*/**/*gradient*.csv")
    )
    projected_files = sorted(run_dir.glob("rom_output/outer_*/**/*projected*fom*"))
    return {
        "run_directory": str(run_dir),
        "basis_metadata": metadata,
        "objective_comparisons": _objective_comparisons(run_dir / "rom_gradient_diagnostics.csv"),
        "inner_rom_objective_history": _inner_objectives(run_dir),
        "performance": _performance(run_dir),
        "gradient_vector_alignment": {
            "status": "unavailable",
            "reason": "archive contains no exported g_FOM, g_frozen, g_exact, g_proj vectors or dot-product alignment matrix",
            "gradient_diagnostic_files_found": [str(path) for path in gradient_files],
            "scalar_gradient_error_fields_are_not_vectors": True,
        },
        "projected_fom_coordinates": {
            "status": "available" if projected_files else "unavailable",
            "files_found": [str(path) for path in projected_files],
            "reason": None if projected_files else (
                "no projected-FOM coefficient data found; integrated forward/adjoint coefficient files "
                "do not identify projections of synchronized FOM snapshots"
            ),
        },
        "fom_rom_same_design_comparison": {
            "status": "partial" if (run_dir / "rom_gradient_diagnostics.csv").is_file() else "unavailable",
            "note": "Recorded reference/candidate objective scalars are summarized when present; design IDs, synchronized trajectory identity, and projected-FOM basis coordinates are absent.",
        },
    }


def summarize(paths: list[Path]) -> dict[str, Any]:
    return {"run_count": len(paths), "runs": [summarize_run(path) for path in paths]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="+", type=Path, help="one or more result archive/run directories")
    parser.add_argument("-o", "--output", type=Path, help="optional JSON output path")
    args = parser.parse_args()
    rendered = json.dumps(summarize(args.run), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
