#!/usr/bin/env python3
"""Project physical FOM snapshots into a provenance-matched coupled ROM basis.

Each input is a NumPy ``.npz`` artifact. Snapshot archives contain
``snapshots``, ``times``, and ``snapshot_ids``; the other archives contain
``basis``, ``lift``, or ``metric`` respectively. Every archive must also
contain scalar string metadata ``design_id``, ``basis_coordinate_id``,
``mesh_dof_map_id``, and ``field_order`` (a JSON list of field names).
The metric archive additionally declares ``metric_definition='coupled_H'``.
The metric may be a positive diagonal vector or a symmetric positive-definite
dense matrix. Array row/column order is preserved exactly as supplied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


MAX_GRAM_CONDITION = 1.0e12
PROVENANCE_KEYS = ("design_id", "basis_coordinate_id", "mesh_dof_map_id")


@dataclass(frozen=True)
class ProjectionResult:
    gram: np.ndarray
    coefficients: np.ndarray
    times: np.ndarray
    snapshot_ids: tuple[str, ...]
    provenance: dict[str, Any]
    gram_condition: float


def _load_npz(path: str | Path) -> dict[str, np.ndarray]:
    path = Path(path)
    try:
        with np.load(path, allow_pickle=False) as archive:
            return {name: np.array(archive[name], copy=True) for name in archive.files}
    except (OSError, ValueError) as exc:
        raise ValueError(f"cannot read NumPy archive {path}: {exc}") from exc


def _required_array(archive: dict[str, np.ndarray], key: str, path: str | Path) -> np.ndarray:
    if key not in archive:
        raise ValueError(f"{path}: missing required array {key!r}")
    return archive[key]


def _scalar_text(archive: dict[str, np.ndarray], key: str, path: str | Path) -> str:
    array = _required_array(archive, key, path)
    if array.shape != () or array.dtype.kind not in "US":
        raise ValueError(f"{path}: {key!r} must be a scalar string")
    value = str(array.item()).strip()
    if not value:
        raise ValueError(f"{path}: {key!r} must not be empty")
    return value


def _metadata(archive: dict[str, np.ndarray], path: str | Path) -> dict[str, Any]:
    result = {key: _scalar_text(archive, key, path) for key in PROVENANCE_KEYS}
    field_order_text = _scalar_text(archive, "field_order", path)
    try:
        field_order = json.loads(field_order_text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path}: field_order must be a JSON list") from exc
    if (
        not isinstance(field_order, list)
        or not field_order
        or any(not isinstance(name, str) or not name.strip() for name in field_order)
        or len(set(field_order)) != len(field_order)
    ):
        raise ValueError(f"{path}: field_order must contain unique nonempty field names")
    result["field_order"] = field_order
    return result


def project_snapshots(
    snapshots_path: str | Path,
    basis_path: str | Path,
    lift_path: str | Path,
    metric_path: str | Path,
    *,
    max_condition: float = MAX_GRAM_CONDITION,
) -> ProjectionResult:
    """Compute coupled-metric least-squares coefficients for every snapshot.

    The basis and lift columns/rows are consumed in their stored order. No
    orthonormality assumption is made: the Gram system ``H a = V^T M(q-lift)``
    is solved explicitly after SPD and conditioning checks.
    """
    if not math.isfinite(max_condition) or max_condition <= 1.0:
        raise ValueError("max_condition must be finite and greater than one")

    paths = {
        "snapshots": Path(snapshots_path),
        "basis": Path(basis_path),
        "lift": Path(lift_path),
        "metric": Path(metric_path),
    }
    artifacts = {name: _load_npz(path) for name, path in paths.items()}
    metadata = {name: _metadata(artifacts[name], paths[name]) for name in artifacts}

    reference = metadata["snapshots"]
    for artifact_name in ("basis", "lift", "metric"):
        candidate = metadata[artifact_name]
        for key in (*PROVENANCE_KEYS, "field_order"):
            if candidate[key] != reference[key]:
                raise ValueError(
                    f"{artifact_name} {key} mismatch: {candidate[key]!r} != {reference[key]!r}"
                )

    snapshots = _required_array(artifacts["snapshots"], "snapshots", paths["snapshots"])
    basis = _required_array(artifacts["basis"], "basis", paths["basis"])
    lift = _required_array(artifacts["lift"], "lift", paths["lift"])
    metric = _required_array(artifacts["metric"], "metric", paths["metric"])
    times = _required_array(artifacts["snapshots"], "times", paths["snapshots"])
    snapshot_ids_array = _required_array(
        artifacts["snapshots"], "snapshot_ids", paths["snapshots"]
    )

    if snapshots.ndim != 2 or basis.ndim != 2:
        raise ValueError("snapshots and basis must be rank-2 arrays")
    sample_count, dof_count = snapshots.shape
    if sample_count < 1 or dof_count < 1 or basis.shape[0] != dof_count or basis.shape[1] < 1:
        raise ValueError("snapshot/basis dimensions are empty or incompatible")
    if lift.shape != (dof_count,):
        raise ValueError(f"lift must have shape ({dof_count},), got {lift.shape}")
    if dof_count % len(reference["field_order"]) != 0:
        raise ValueError("dof count is not divisible by declared field_order length")
    if metric.ndim == 1:
        if metric.shape != (dof_count,):
            raise ValueError(f"diagonal metric must have shape ({dof_count},), got {metric.shape}")
        if np.any(metric <= 0.0):
            raise ValueError("diagonal coupled metric must be strictly positive")
        weighted_basis = metric[:, None] * basis
        metric_definition = "coupled_H_diagonal"
    elif metric.ndim == 2:
        if metric.shape != (dof_count, dof_count):
            raise ValueError(
                f"dense metric must have shape ({dof_count},{dof_count}), got {metric.shape}"
            )
        if not np.allclose(metric, metric.T, rtol=1.0e-12, atol=1.0e-14):
            raise ValueError("dense coupled metric must be symmetric")
        try:
            np.linalg.cholesky(metric)
        except np.linalg.LinAlgError as exc:
            raise ValueError("dense coupled metric must be positive definite") from exc
        weighted_basis = metric @ basis
        metric_definition = "coupled_H_dense"
    else:
        raise ValueError("metric must be a diagonal vector or dense matrix")
    if not np.isfinite(snapshots).all() or not np.isfinite(basis).all() or not np.isfinite(lift).all():
        raise ValueError("snapshots, basis, and lift must be finite")
    if not np.isfinite(metric).all() or not np.isfinite(weighted_basis).all():
        raise ValueError("metric and metric-weighted basis must be finite")
    declared_metric = _scalar_text(artifacts["metric"], "metric_definition", paths["metric"])
    if declared_metric != "coupled_H":
        raise ValueError("metric_definition must be exactly 'coupled_H'")

    times = np.asarray(times, dtype=np.float64)
    if times.shape != (sample_count,) or not np.isfinite(times).all():
        raise ValueError("times must be a finite vector with one value per snapshot")
    if sample_count > 1 and np.any(np.diff(times) <= 0.0):
        raise ValueError("snapshot times must be strictly increasing")
    if snapshot_ids_array.shape != (sample_count,) or snapshot_ids_array.dtype.kind not in "US":
        raise ValueError("snapshot_ids must be a string vector with one ID per snapshot")
    snapshot_ids = tuple(str(item).strip() for item in snapshot_ids_array.tolist())
    if any(not item for item in snapshot_ids) or len(set(snapshot_ids)) != sample_count:
        raise ValueError("snapshot_ids must be unique nonempty strings")

    gram = basis.T @ weighted_basis
    gram = 0.5 * (gram + gram.T)
    if not np.isfinite(gram).all():
        raise ValueError("Gram matrix contains nonfinite values")
    try:
        np.linalg.cholesky(gram)
    except np.linalg.LinAlgError as exc:
        raise ValueError("basis Gram matrix must be positive definite") from exc
    condition = float(np.linalg.cond(gram))
    if not math.isfinite(condition) or condition > max_condition:
        raise ValueError(
            f"basis Gram condition number {condition:.6e} exceeds limit {max_condition:.6e}"
        )

    differences = snapshots - lift[None, :]
    rhs = differences @ weighted_basis
    try:
        coefficients = np.linalg.solve(gram, rhs.T).T
    except np.linalg.LinAlgError as exc:
        raise ValueError("failed to solve coupled projection Gram system") from exc
    if not np.isfinite(rhs).all() or not np.isfinite(coefficients).all():
        raise ValueError("projection produced nonfinite RHS or coefficients")

    provenance = dict(reference)
    provenance.update(
        {
            "metric_definition": declared_metric,
            "metric_storage": metric_definition,
            "basis_rank": int(basis.shape[1]),
            "dof_count": int(dof_count),
            "snapshot_count": int(sample_count),
            "basis_column_order_preserved": True,
        }
    )
    return ProjectionResult(
        gram=gram,
        coefficients=coefficients,
        times=times,
        snapshot_ids=snapshot_ids,
        provenance=provenance,
        gram_condition=condition,
    )


def write_projection(result: ProjectionResult, output_prefix: str | Path) -> tuple[Path, Path, Path]:
    """Write coefficients, Gram matrix, and provenance/conditioning report."""
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    coefficients_path = prefix.with_name(prefix.name + "_coefficients.csv")
    gram_path = prefix.with_name(prefix.name + "_gram.csv")
    report_path = prefix.with_name(prefix.name + "_report.json")
    mode_count = result.coefficients.shape[1]
    with coefficients_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["snapshot_id", "time", "design_id", "basis_coordinate_id", "mesh_dof_map_id"]
            + [f"a_{index + 1}" for index in range(mode_count)]
        )
        for row, (snapshot_id, time_value) in enumerate(zip(result.snapshot_ids, result.times)):
            writer.writerow(
                [
                    snapshot_id,
                    f"{time_value:.17g}",
                    result.provenance["design_id"],
                    result.provenance["basis_coordinate_id"],
                    result.provenance["mesh_dof_map_id"],
                ]
                + [f"{value:.17g}" for value in result.coefficients[row]]
            )
    np.savetxt(gram_path, result.gram, delimiter=",", fmt="%.17e")
    report = {
        **result.provenance,
        "gram_condition_number": result.gram_condition,
        "max_gram_condition_number": MAX_GRAM_CONDITION,
        "snapshot_ids": list(result.snapshot_ids),
        "times": result.times.tolist(),
        "coefficient_coordinate": "actual ROM basis column order",
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return coefficients_path, gram_path, report_path


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshots", type=Path, help="NPZ with physical snapshot matrix and times")
    parser.add_argument("basis", type=Path, help="NPZ with actual ROM basis matrix")
    parser.add_argument("lift", type=Path, help="NPZ with matching affine lift vector")
    parser.add_argument("metric", type=Path, help="NPZ with coupled metric vector/matrix")
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args(argv)
    try:
        result = project_snapshots(args.snapshots, args.basis, args.lift, args.metric)
        outputs = write_projection(result, args.output_prefix)
    except ValueError as exc:
        parser.error(str(exc))
    print(
        f"projected {result.provenance['snapshot_count']} snapshots onto rank "
        f"{result.provenance['basis_rank']} basis; Gram condition="
        f"{result.gram_condition:.6e}; outputs={', '.join(map(str, outputs))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
