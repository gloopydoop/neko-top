#!/usr/bin/env python3
"""Validate and compare indexed CPU ROM assembly diagnostic CSVs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


RELATIVE_TOLERANCE = 5.0e-13
IMPLEMENTATIONS = ("packed", "legacy")
CSV_COLUMNS = ("operator", "implementation", "i", "j", "k", "value")
OUTPUT_SHAPES = {
    "mass": (3, 3),
    "diffusion": (3, 3),
    "brinkman": (3, 3),
    "brinkman_lift": (3,),
    "convection": (3, 3, 3),
    "scalar_linear": (3, 3),
    "scalar_constant": (3,),
    "scalar_quadratic": (3, 3, 3),
    "scalar_implicit": (3, 3),
    "petrov_cross_mass": (3, 3),
    "petrov_brinkman": (3, 3),
    "petrov_brinkman_lift": (3,),
}

Key = tuple[str, str, int, int, int]


def _expected_keys() -> set[Key]:
    expected: set[Key] = set()
    for operator, shape in OUTPUT_SHAPES.items():
        if len(shape) == 1:
            indices = tuple((i, 0, 0) for i in range(1, shape[0] + 1))
        elif len(shape) == 2:
            indices = tuple(
                (i, j, 0)
                for i in range(1, shape[0] + 1)
                for j in range(1, shape[1] + 1)
            )
        else:
            indices = tuple(
                (i, j, k)
                for i in range(1, shape[0] + 1)
                for j in range(1, shape[1] + 1)
                for k in range(1, shape[2] + 1)
            )
        for implementation in IMPLEMENTATIONS:
            expected.update(
                (operator, implementation, i, j, k) for i, j, k in indices
            )
    return expected


EXPECTED_KEYS = _expected_keys()


def read_output_csv(path: str | Path) -> dict[Key, float]:
    """Read a diagnostic dump, requiring every expected key exactly once."""
    path = Path(path)
    values: dict[Key, float] = {}
    try:
        stream = path.open(newline="", encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"cannot open {path}: {exc}") from exc

    with stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != CSV_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {CSV_COLUMNS}, got {reader.fieldnames}"
            )
        try:
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(row.get(column) is None for column in CSV_COLUMNS):
                    raise ValueError(f"{path}:{line_number}: malformed CSV row")
                operator = row["operator"]
                implementation = row["implementation"]
                if operator not in OUTPUT_SHAPES:
                    raise ValueError(f"{path}:{line_number}: unexpected operator {operator!r}")
                if implementation not in IMPLEMENTATIONS:
                    raise ValueError(
                        f"{path}:{line_number}: unexpected implementation {implementation!r}"
                    )
                try:
                    i, j, k = (int(row[name]) for name in ("i", "j", "k"))
                    value = float(row["value"])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{path}:{line_number}: invalid index/value") from exc
                key: Key = (operator, implementation, i, j, k)
                if key in values:
                    raise ValueError(f"{path}:{line_number}: duplicate key {key}")
                if key not in EXPECTED_KEYS:
                    raise ValueError(f"{path}:{line_number}: unexpected key {key}")
                if not math.isfinite(value):
                    raise ValueError(f"{path}:{line_number}: nonfinite value for {key}")
                values[key] = value
        except csv.Error as exc:
            raise ValueError(f"{path}: invalid CSV: {exc}") from exc

    missing = EXPECTED_KEYS - values.keys()
    if missing:
        sample = sorted(missing)[:5]
        raise ValueError(f"{path}: missing {len(missing)} expected keys; examples: {sample}")
    return values


def compare_outputs(
    rank1_path: str | Path,
    rank2_path: str | Path,
    output_path: str | Path,
) -> list[dict[str, str | int]]:
    """Validate both complete dumps and write scale-aware rank comparisons."""
    rank1 = read_output_csv(rank1_path)
    rank2 = read_output_csv(rank2_path)
    if rank1.keys() != rank2.keys():  # Defensive; each reader enforces same schema.
        raise ValueError("rank output key sets differ")

    grouped: dict[tuple[str, str], list[Key]] = defaultdict(list)
    for key in sorted(EXPECTED_KEYS):
        grouped[key[:2]].append(key)

    results: list[dict[str, str | int]] = []
    for (operator, implementation), keys in sorted(grouped.items()):
        values1 = [rank1[key] for key in keys]
        values2 = [rank2[key] for key in keys]
        differences = [a - b for a, b in zip(values1, values2)]
        norm1 = math.sqrt(sum(value * value for value in values1))
        norm2 = math.sqrt(sum(value * value for value in values2))
        difference_norm = math.sqrt(sum(value * value for value in differences))
        scale = max(norm1, norm2, math.ulp(0.0))
        relative = difference_norm / scale
        if not all(math.isfinite(value) for value in (norm1, norm2, difference_norm, relative)):
            raise ValueError(f"nonfinite comparison metric for {operator}/{implementation}")
        results.append(
            {
                "operator": operator,
                "implementation": implementation,
                "entries": len(keys),
                "max_abs_rank1_vs_rank2": f"{max(map(abs, differences)):.16e}",
                "relative_frobenius_rank1_vs_rank2": f"{relative:.16e}",
                "rank1_frobenius": f"{norm1:.16e}",
                "rank2_frobenius": f"{norm2:.16e}",
                "relative_tolerance": f"{RELATIVE_TOLERANCE:.1e}",
                "status": "PASS" if relative <= RELATIVE_TOLERANCE else "FAIL",
            }
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = (
        "operator",
        "implementation",
        "entries",
        "max_abs_rank1_vs_rank2",
        "relative_frobenius_rank1_vs_rank2",
        "rank1_frobenius",
        "rank2_frobenius",
        "relative_tolerance",
        "status",
    )
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(results)
    return results


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rank1_csv", type=Path)
    parser.add_argument("rank2_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args(argv)
    try:
        results = compare_outputs(args.rank1_csv, args.rank2_csv, args.output_csv)
    except ValueError as exc:
        parser.error(str(exc))
    max_relative = max(float(row["relative_frobenius_rank1_vs_rank2"]) for row in results)
    failures = sum(row["status"] != "PASS" for row in results)
    print(
        f"validated {len(EXPECTED_KEYS)} indexed values per file across "
        f"{len(results)} groups; max relative Frobenius={max_relative:.16e}; "
        f"gate={RELATIVE_TOLERANCE:.1e}; failures={failures}; output={args.output_csv}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
