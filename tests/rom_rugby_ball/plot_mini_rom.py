#!/usr/bin/env python3
"""Plot the inspectable histories from ROM_low_Re_mini.

This is intentionally a small test-side utility: it reads the CSV files
written by the ROM experiment and never changes the optimisation result.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    import matplotlib.pyplot as plt

    optimization = _read(args.result_dir / "optimization_data.csv")
    validation_path = args.result_dir / "rom_validation.csv"
    validation = _read(validation_path) if validation_path.exists() else []
    coefficient_files = sorted(args.result_dir.glob("ROM_forward_coefficients_iter_*.csv"))
    coefficients = _read(coefficient_files[-1]) if coefficient_files else []
    x = list(range(len(optimization)))
    labels = [row["phase"] for row in optimization]

    figure, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=False)
    axes[0].plot(x, [float(row["total_objective"]) for row in optimization], "o-", label="total")
    axes[0].plot(x, [float(row["scalar_mixing"]) for row in optimization], "o-", label="scalar mixing")
    axes[0].plot(x, [float(row["viscous_dissipation"]) for row in optimization], "o-", label="viscous")
    axes[0].plot(x, [float(row["brinkman_dissipation"]) for row in optimization], "o-", label="Brinkman")
    axes[0].set_ylabel("objective")
    axes[0].legend()
    axes[1].plot(x, [float(row["volume_constraint"]) for row in optimization], "o-", label="volume")
    axes[1].plot(x, [float(row["max_design_change"]) for row in optimization], "o-", label="max design change")
    axes[1].set_ylabel("constraint/change")
    axes[1].legend()
    if validation:
        vx = [int(row["outer_iteration"]) for row in validation]
        axes[2].plot(vx, [float(row["predicted_objective"]) for row in validation], "o-", label="ROM")
        axes[2].plot(vx, [float(row["actual_objective"]) for row in validation], "o-", label="full order")
        axes[2].set_ylabel("refresh objective")
        axes[2].legend()
    if coefficients:
        coefficient_names = [name for name in coefficients[0] if name.startswith("a_")]
        time = [float(row["time"]) for row in coefficients]
        for name in coefficient_names:
            axes[3].plot(time, [float(row[name]) for row in coefficients], "o-", label=name)
        axes[3].set_ylabel("ROM coefficient")
        axes[3].set_xlabel("time")
        axes[3].legend()
        axes[3].set_title(f"{coefficient_files[-1].name} (common velocity/scalar basis)")
    else:
        axes[3].text(0.02, 0.5, "No ROM coefficient CSV found", transform=axes[3].transAxes)
        axes[3].set_axis_off()
    axes[2].set_xlabel("outer refresh")
    axes[0].set_title("ROM_low_Re_mini: " + ", ".join(labels))
    figure.tight_layout()
    output = args.output or args.result_dir / "rom_mini_history.png"
    figure.savefig(output, dpi=150)
    print(output)


if __name__ == "__main__":
    main()
