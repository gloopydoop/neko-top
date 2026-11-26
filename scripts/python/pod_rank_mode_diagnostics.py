"""Append-only per-mode evidence for Python POD recycling stages.

The diagnostics are derived from the small Gram matrices and mode vectors
already computed by ``recycle_distributed`` and ``common_basis_distributed``.
They do not infer physical causes for a rejected direction.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


FAMILIES = ("forward", "adjoint")


def _candidate_status(diagnostics: dict[str, Any]) -> tuple[dict[int, dict], dict[int, dict]]:
    accepted = {
        int(item["candidate_index_zero_based"]): item
        for item in diagnostics.get("accepted_modes", [])
        if "candidate_index_zero_based" in item
    }
    rejected = {
        int(item["candidate_index_zero_based"]): item
        for item in diagnostics.get("rejected_modes", [])
        if "candidate_index_zero_based" in item
    }
    return accepted, rejected


def build_family_mode_record(
    outer_iteration: int,
    family: str,
    recycling_diagnostics: dict[str, Any] | None,
    common_basis_diagnostics: dict[str, Any] | None,
    *,
    status: str = "available",
) -> dict[str, Any]:
    """Create one JSON-safe family record from computed Python diagnostics.

    Norm fields explicitly name the Python mass metric and projector stage.
    They must not be interpreted as Neko's boundary-restriction or pressure
    continuity measurements.
    """

    family = str(family).casefold()
    if family not in FAMILIES:
        raise ValueError(f"family must be one of {FAMILIES}, got {family!r}")
    recycling = recycling_diagnostics or {}
    basis = common_basis_diagnostics or {}
    sources = list(basis.get("input_modes", []))
    pre = list(basis.get("pre_projection_norms", []))
    post = list(basis.get("post_projection_norms", []))
    residual = list(basis.get("pre_normalisation_orthogonal_residual_norms", []))
    if any(len(values) not in (0, len(sources)) for values in (pre, post, residual)):
        raise ValueError("common-basis per-mode arrays do not align with input_modes")

    accepted, rejected = _candidate_status(basis)
    retained_eigen = list(recycling.get("retained_gram_eigen_indices_descending", []))
    mode_records: list[dict[str, Any]] = []
    for index, source in enumerate(sources):
        if str(source.get("source_family", "")).casefold() != family:
            continue
        accepted_item = accepted.get(index)
        rejected_item = rejected.get(index)
        source_mode = source.get("source_mode_index_zero_based")
        recycled_eigen_index = None
        if isinstance(source_mode, int) and 0 <= source_mode < len(retained_eigen):
            recycled_eigen_index = int(retained_eigen[source_mode])
        record = {
            **source,
            "candidate_index_zero_based": index,
            "recycled_source_gram_eigen_index_descending": recycled_eigen_index,
            "pre_projection_mass_norm": pre[index] if pre else None,
            "post_python_admissibility_projection_mass_norm": post[index] if post else None,
            "pre_normalisation_orthogonal_residual_mass_norm": residual[index] if residual else None,
            "status": "accepted" if accepted_item is not None else "rejected" if rejected_item else "not_emitted",
            "retained_basis_index_zero_based": (
                next((j for j, item in enumerate(basis.get("accepted_modes", []))
                      if int(item.get("candidate_index_zero_based", -1)) == index), None)
            ),
            "rejection_reason": rejected_item.get("reason") if rejected_item else None,
        }
        mode_records.append(record)

    family_recycling_modes = [
        item for item in recycling.get("input_modes", [])
        if isinstance(item, dict)
    ]
    return {
        "outer_iteration": int(outer_iteration),
        "family": family,
        "status": status,
        "recycling_status": "available" if recycling_diagnostics is not None else "not_emitted",
        "common_basis_status": "available" if common_basis_diagnostics is not None else "not_emitted",
        "producer_stage": "python_historical_recycling_and_common_basis",
        "mode_norm_metric": "global coupled mass metric supplied to Python recycling",
        "projection_semantics": (
            "post_python_admissibility_projection_mass_norm is after the optional "
            "Python projector; it is not Neko boundary restriction or continuity treatment"
        ),
        "recycling": {
            "candidate_rank": recycling.get("candidate_rank"),
            "recycled_rank": recycling.get("recycled_rank"),
            "rank_cap": recycling.get("rank_cap"),
            "tolerance": recycling.get("tolerance"),
            "gram_eigenvalue_threshold": recycling.get("gram_eigenvalue_threshold"),
            "gram_eigenvalues_descending": recycling.get("gram_eigenvalues_descending", []),
            "gram_singular_values_descending": recycling.get("singular_values_descending", []),
            "retained_gram_eigen_indices_descending": retained_eigen,
            "input_modes": family_recycling_modes,
            "rejected_modes": recycling.get("rejected_modes", []),
        },
        "common_basis": {
            "candidate_rank": basis.get("candidate_rank"),
            "retained_rank": basis.get("retained_rank"),
            "rank_cap": basis.get("rank_cap"),
            "admissibility_projected": basis.get("admissibility_projected"),
            "mass_gram_error": basis.get("mass_gram_error"),
            "gram_singular_values": basis.get("gram_singular_values", []),
            "modes": mode_records,
        },
        "cause_attribution": "not_inferred_from_rank_or_norm_diagnostics",
    }


def append_family_mode_diagnostics(
    comm,
    outer_dir: str | os.PathLike[str],
    outer_iteration: int,
    family: str,
    recycling_diagnostics: dict[str, Any] | None,
    common_basis_diagnostics: dict[str, Any] | None,
    *,
    status: str = "available",
) -> Path | None:
    """Append one family record on the Python root rank; preserve old files."""

    family = str(family).casefold()
    record = build_family_mode_record(
        outer_iteration,
        family,
        recycling_diagnostics,
        common_basis_diagnostics,
        status=status,
    )
    if comm.Get_rank() != 0:
        return None
    path = Path(outer_dir) / f"pod_basis_mode_diagnostics_{family}.jsonl"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
        handle.flush()
    return path
