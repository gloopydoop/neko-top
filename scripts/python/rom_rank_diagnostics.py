#!/usr/bin/env python3
"""Summarise ROM basis-rank evidence without opening large field files.

The current run metadata records some rank transitions, while several Python
POD and handoff stages are not exported yet. This reader reports recorded
values verbatim and lists absent transitions explicitly; it never infers that
two different rank fields refer to the same stage.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


FAMILY_FIELDS = {
    "forward": {
        "candidate_rank": "forward_candidate_modes",
        "retained_rank": "forward_retained_rank",
        "history_iterations": "forward_history_iterations",
        "requested_snapshot_pod_rank": None,
        "actual_snapshot_pod_rank": None,
        "historical_candidate_rank": None,
        "python_recycled_rank": None,
        "neko_incoming_rank": None,
        "rank_after_boundary_restriction": None,
        "rank_after_continuity_treatment": None,
        "final_orthonormal_rank": None,
    },
    "adjoint": {
        "candidate_rank": "adjoint_candidate_modes",
        "retained_rank": "adjoint_retained_rank",
        "history_iterations": "adjoint_history_iterations",
        "requested_snapshot_pod_rank": None,
        "actual_snapshot_pod_rank": None,
        "historical_candidate_rank": None,
        "python_recycled_rank": None,
        "neko_incoming_rank": None,
        "rank_after_boundary_restriction": None,
        "rank_after_continuity_treatment": None,
        "final_orthonormal_rank": None,
    },
}

NEKO_RANK_STAGES = (
    "requested_block_rank",
    "received_rank",
    "post_boundary_restriction_rank",
    "post_continuity_treatment_rank",
    "post_linear_dependence_rank",
    "final_retained_rank",
)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _csv_summary(path: Path) -> dict[str, Any]:
    """Read scalar columns only; large basis and field arrays are never read."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        columns = list(reader.fieldnames or [])
        rows = list(reader)

    numeric: dict[str, list[float]] = {}
    true_counts: dict[str, int] = {}
    for row in rows:
        for name, raw in row.items():
            if raw is None or not raw.strip():
                continue
            value = raw.strip()
            if value.casefold() in {"true", "t", ".true."}:
                true_counts[name] = true_counts.get(name, 0) + 1
                continue
            if value.casefold() in {"false", "f", ".false."}:
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            if math.isfinite(number):
                numeric.setdefault(name, []).append(number)

    metrics: dict[str, dict[str, float]] = {}
    for name, values in numeric.items():
        metrics[name] = {
            "min": min(values),
            "max": max(values),
        }
    return {
        "path": str(path),
        "columns": columns,
        "rows": len(rows),
        "numeric_ranges": metrics,
        "true_counts": true_counts,
        "per_mode_records": (
            [
                {"mode_index": row.get("mode", row.get("mode_index")), "values": row}
                for row in rows
            ]
            if any(name.strip().casefold() in {"mode", "mode_index"} for name in columns)
            else None
        ),
        "rank_evidence_columns": [
            name for name in columns
            if any(token in name.casefold() for token in ("rank", "mode_index", "singular"))
        ],
    }


def _family_summary(metadata: dict[str, Any], family: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    missing: list[str] = []
    for output_name, source_name in FAMILY_FIELDS[family].items():
        if source_name is None:
            result[output_name] = None
            missing.append(output_name)
        elif source_name in metadata:
            result[output_name] = metadata[source_name]
        else:
            result[output_name] = None
            missing.append(output_name)
    if family == "adjoint":
        result["enrichment_enabled"] = metadata.get("adjoint_enrichment")
    result["missing_transitions"] = missing
    return result


def _pod_provenance_records(outer_dir: Path, outer_iteration: Any) -> dict[str, Any]:
    path = outer_dir / "pod_rank_provenance.jsonl"
    records: dict[str, list[dict[str, Any]]] = {"forward": [], "adjoint": []}
    warnings: list[str] = []
    ignored_outer = 0
    if not path.is_file():
        return {
            "path": str(path),
            "status": "not_emitted",
            "records": records,
            "warnings": warnings,
            "ignored_other_outer_records": ignored_outer,
        }

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                warnings.append(f"line {line_number}: invalid JSON ({error.msg})")
                continue
            if not isinstance(record, dict):
                warnings.append(f"line {line_number}: record is not a JSON object")
                continue
            family = str(record.get("family", "")).casefold()
            if family not in records:
                warnings.append(f"line {line_number}: unknown family {family!r}")
                continue
            try:
                record_outer = int(record["outer_iteration"])
                expected_outer = int(outer_iteration)
            except (KeyError, TypeError, ValueError):
                warnings.append(f"line {line_number}: missing or invalid outer_iteration")
                continue
            if record_outer != expected_outer:
                ignored_outer += 1
                continue
            records[family].append(record)

    status = "available" if any(records.values()) else "no_matching_outer_family_records"
    return {
        "path": str(path),
        "status": status,
        "records": records,
        "warnings": warnings,
        "ignored_other_outer_records": ignored_outer,
    }


def _python_mode_diagnostics(outer_dir: Path, outer_iteration: Any, family: str) -> dict[str, Any]:
    """Read append-only Python per-mode evidence without inferring missing data."""
    path = outer_dir / f"pod_basis_mode_diagnostics_{family}.jsonl"
    records: list[dict[str, Any]] = []
    warnings: list[str] = []
    try:
        expected_outer = int(outer_iteration)
    except (TypeError, ValueError):
        return {
            "path": str(path),
            "status": "outer_iteration_unknown",
            "records": [],
            "warnings": ["cannot select records without a valid outer_iteration"],
        }
    if not path.is_file():
        return {"path": str(path), "status": "not_emitted", "records": [], "warnings": warnings}
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                warnings.append(f"line {line_number}: invalid JSON ({error.msg})")
                continue
            if not isinstance(record, dict):
                warnings.append(f"line {line_number}: record is not a JSON object")
                continue
            if str(record.get("family", "")).casefold() != family:
                continue
            try:
                record_outer = int(record["outer_iteration"])
            except (KeyError, TypeError, ValueError):
                warnings.append(f"line {line_number}: missing or invalid outer_iteration")
                continue
            if record_outer == expected_outer:
                records.append(record)
    status = "not_emitted_for_outer_family" if not records else "available" if len(records) == 1 else "multiple_records"
    return {"path": str(path), "status": status, "records": records, "warnings": warnings}


def _normalize_pod_record(record: dict[str, Any]) -> dict[str, Any]:
    """Keep configured rank cap distinct from the observed POD rank."""
    return {
        "status": record.get("status"),
        "configured_requested_pod_max_rank": record.get("requested_pod_rank"),
        "requested_rank_semantics": "configured POD keep_modes maximum; not an observed rank",
        "actual_snapshot_pod_rank": record.get("actual_pod_rank"),
        "actual_rank_semantics": (
            "producer-recorded available POD modes after the keep_modes cap"
        ),
        "historical_candidate_rank": record.get("historical_candidate_rank"),
        "python_recycled_rank": record.get("recycled_rank"),
        "recycling_rank_cap": record.get("recycling_rank_cap"),
        "raw_record": record,
    }


def _neko_rank_stages(outer_dir: Path, outer_iteration: Any, family: str) -> dict[str, Any]:
    path = outer_dir / "basis_diagnostics" / f"rom_rank_provenance_{family}.csv"
    matches: dict[str, list[dict[str, Any]]] = {stage: [] for stage in NEKO_RANK_STAGES}
    warnings: list[str] = []
    try:
        expected_outer = int(outer_iteration)
    except (TypeError, ValueError):
        return {
            "path": str(path),
            "file_status": "outer_iteration_unknown",
            "stages": {
                stage: {"status": "outer_iteration_unknown", "records": []}
                for stage in NEKO_RANK_STAGES
            },
            "warnings": ["cannot join without a valid outer_iteration"],
        }
    if not path.is_file():
        return {
            "path": str(path),
            "file_status": "not_emitted",
            "stages": {
                stage: {"status": "not_emitted", "records": []}
                for stage in NEKO_RANK_STAGES
            },
            "warnings": warnings,
        }

    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        required = {"outer_iteration", "family", "stage", "rank", "status", "rejected_count", "rejection_reason"}
        if not required.issubset(set(reader.fieldnames or ())):
            warnings.append("rank CSV is missing one or more required provenance columns")
        else:
            for row_number, row in enumerate(reader, 2):
                try:
                    row_outer = int(row["outer_iteration"])
                except (TypeError, ValueError):
                    warnings.append(f"row {row_number}: invalid outer_iteration")
                    continue
                if row_outer != expected_outer:
                    continue
                if str(row.get("family", "")).strip().casefold() != family:
                    continue
                stage = str(row.get("stage", "")).strip()
                if stage not in matches:
                    warnings.append(f"row {row_number}: unknown stage {stage!r}")
                    continue
                try:
                    reported_rank = int(row["rank"])
                except (TypeError, ValueError):
                    reported_rank = None
                try:
                    rejected_count = int(row["rejected_count"])
                except (TypeError, ValueError):
                    rejected_count = None
                status = str(row.get("status", "")).strip() or "unknown"
                matches[stage].append({
                    "rank": reported_rank if reported_rank is not None and reported_rank >= 0 else None,
                    "reported_rank": reported_rank,
                    "status": status,
                    "rejected_count": rejected_count if rejected_count is not None and rejected_count >= 0 else None,
                    "reported_rejected_count": rejected_count,
                    "rejection_reason": row.get("rejection_reason"),
                })

    any_matches = any(matches.values())
    stages: dict[str, Any] = {}
    for stage, records in matches.items():
        if not records:
            stages[stage] = {"status": "not_emitted_for_outer_family", "records": []}
        elif len(records) == 1:
            stages[stage] = {"status": records[0]["status"], "records": records}
        else:
            # Repeated candidate evaluations share outer/family keys in this
            # CSV. Preserve them rather than guessing which row pairs with POD.
            stages[stage] = {"status": "multiple_records", "records": records}
    return {
        "path": str(path),
        "file_status": "available" if any_matches else "no_matching_outer_family_rows",
        "stages": stages,
        "warnings": warnings,
    }


def _rank_provenance_join(outer_dir: Path, outer_iteration: Any) -> dict[str, Any]:
    pod = _pod_provenance_records(outer_dir, outer_iteration)
    families: dict[str, Any] = {}
    for family in ("forward", "adjoint"):
        pod_records = pod["records"][family]
        families[family] = {
            "python_pod": {
                "status": (
                    "not_emitted" if pod["status"] == "not_emitted"
                    else "missing_for_outer_family" if not pod_records
                    else str(pod_records[0].get("status", "available")) if len(pod_records) == 1
                    else "multiple_records"
                ),
                "records": [_normalize_pod_record(record) for record in pod_records],
            },
            "python_mode_diagnostics": _python_mode_diagnostics(
                outer_dir, outer_iteration, family
            ),
            "neko": _neko_rank_stages(outer_dir, outer_iteration, family),
            "candidate_level_alignment": "unavailable",
            "alignment_note": (
                "Join key is outer_iteration+family only. The current Neko rank CSV "
                "does not carry a candidate/evaluation identifier, and repeated rows "
                "are retained without pairing them to an individual Python record."
            ),
        }
    return {
        "join_key": ["outer_iteration", "family"],
        "outer_iteration": outer_iteration,
        "python_source": {
            "path": pod["path"],
            "status": pod["status"],
            "warnings": pod["warnings"],
            "ignored_other_outer_records": pod["ignored_other_outer_records"],
        },
        "families": families,
    }


def summarize_outer(outer_dir: Path) -> dict[str, Any]:
    """Summarise one ``outer_NNN`` result directory."""
    outer_dir = outer_dir.resolve()
    metadata_path = outer_dir / "basis_metadata.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"basis_metadata.json not found in {outer_dir}")
    metadata = _read_json(metadata_path)
    outer_iteration = metadata.get("outer_iteration")

    diagnostic_files = sorted((outer_dir / "basis_diagnostics").glob("*.csv"))
    forward_files: list[dict[str, Any]] = []
    adjoint_files: list[dict[str, Any]] = []
    shared_files: list[dict[str, Any]] = []
    for path in diagnostic_files:
        summary = _csv_summary(path)
        name = path.name.casefold()
        if "forward" in name:
            forward_files.append(summary)
        elif "adjoint" in name:
            adjoint_files.append(summary)
        else:
            # Existing files such as rom_direct_pod_diagnostics.csv do not
            # identify which independently constructed basis they describe.
            shared_files.append(summary)

    missing_shared = [] if diagnostic_files else ["basis_diagnostics/*.csv"]
    return {
        "outer_directory": str(outer_dir),
        "outer_iteration": metadata.get("outer_iteration"),
        "rank_provenance": _rank_provenance_join(outer_dir, outer_iteration),
        "basis_pair": {
            "forward": metadata.get("basis_forward"),
            "adjoint": metadata.get("basis_adjoint"),
        },
        "basis_ranks": {
            "forward": _family_summary(metadata, "forward"),
            "adjoint": _family_summary(metadata, "adjoint"),
            "enriched_rank": metadata.get("enriched_rank"),
            "forward_adjoint_overlap": metadata.get("forward_adjoint_overlap"),
            "missing_enriched_fields": [
                field for field in ("enriched_rank", "forward_adjoint_overlap")
                if field not in metadata
            ],
        },
        "diagnostics": {
            "forward_files": forward_files,
            "adjoint_files": adjoint_files,
            "shared_or_unlabelled_files": shared_files,
            "missing": missing_shared,
            "note": (
                "Unlabelled CSVs are not assigned to forward or adjoint. "
                "Rank stages absent from metadata remain null."
            ),
        },
        "metadata_path": str(metadata_path),
    }


def discover_outer_dirs(path: Path) -> list[Path]:
    path = path.resolve()
    if path.is_file():
        if path.name != "basis_metadata.json":
            raise ValueError("file input must be a basis_metadata.json")
        return [path.parent]
    if not path.is_dir():
        raise FileNotFoundError(path)
    return sorted({p.parent for p in path.rglob("basis_metadata.json")})


def summarize(path: Path) -> dict[str, Any]:
    directories = discover_outer_dirs(path)
    return {
        "input": str(path.resolve()),
        "outer_count": len(directories),
        "outer_summaries": [summarize_outer(directory) for directory in directories],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="result archive, outer directory, or basis_metadata.json")
    parser.add_argument("-o", "--output", type=Path, help="write JSON report to this path (default: stdout)")
    args = parser.parse_args()
    report = summarize(args.path)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
