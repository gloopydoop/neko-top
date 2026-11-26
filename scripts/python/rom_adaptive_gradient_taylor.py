#!/usr/bin/env python3
"""Summarize an optimizer directional-FD audit as a first-order Taylor check.

This tool only postprocesses ``mma_directional_fd.csv``.  It does not alter
the optimizer, select a gradient, or decide scientific acceptance.  The
reported remainders are meaningful only when the audit's exact-algebraic
directional dot is enabled and the rows come from the adaptive mode.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def summarize(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no audit rows in {path}")
    required = {
        "mode", "h", "total_fd", "exact_algebraic_audit_enabled",
        "exact_algebraic_unverified_dot_step", "total_at_0", "total_at_h",
        "total_at_2h",
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"missing audit columns: {', '.join(sorted(missing))}")
    if any(row["exact_algebraic_audit_enabled"].strip().upper() != "T" for row in rows):
        raise ValueError("exact-algebraic directional audit was not enabled for every row")

    adaptive = [row for row in rows if row["mode"].strip() == "adaptive"]
    if not adaptive:
        raise ValueError("CSV contains no adaptive-basis rows")
    results: list[dict[str, float | str]] = []
    for row in adaptive:
        h, fd, grad = (float(row[key]) for key in
                       ("h", "total_fd", "exact_algebraic_unverified_dot_step"))
        f0, f1, f2 = (float(row[key]) for key in
                      ("total_at_0", "total_at_h", "total_at_2h"))
        if not all(math.isfinite(v) for v in (h, fd, grad, f0, f1, f2)) or h <= 0:
            raise ValueError("non-finite audit value or nonpositive h")
        rem_h = abs(f1 - f0 - h * grad)
        rem_2h = abs(f2 - f0 - 2.0 * h * grad)
        results.append({
            "outer": row.get("outer", ""), "h": h, "fd": fd, "gradient_dot": grad,
            "abs_fd_error": abs(fd - grad), "remainder_h": rem_h,
            "remainder_2h": rem_2h,
            "remainder_ratio": rem_2h / rem_h if rem_h else math.nan,
        })
    return sorted(results, key=lambda item: float(item["h"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="optimizer mma_directional_fd.csv")
    args = parser.parse_args()
    print("outer,h,adaptive_fd,exact_algebraic_dot,abs_fd_error,R(h),R(2h),R(2h)/R(h)")
    for row in summarize(args.csv):
        print("{outer},{h:.8g},{fd:.12g},{gradient_dot:.12g},{abs_fd_error:.6g},"
              "{remainder_h:.6g},{remainder_2h:.6g},{remainder_ratio:.8g}".format(**row))
    print("Scope: directional adaptive-basis audit only; no acceptance gate is applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
