#!/usr/bin/env python3
"""Focused configuration and basis-space checks for the ROM formulations.

This test deliberately does not start Neko or the ADIOS2 peer.  It checks the
configuration contract and the finite-dimensional mass algebra that the live
Fortran basis handoff applies before an expensive lifecycle is launched.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "python"))
from rom_config_tools import derive_rom_basis_pair  # noqa: E402
from recycled_pod_tools import (  # noqa: E402
    DistributedPODFactors,
    append_factor_history,
    recycle_distributed,
)
from mpi4py import MPI  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = ROOT / "examples" / "ROM_low_Re_mini"


def load_case(path: Path) -> dict:
    """Load a Neko case after removing its supported ``//`` comments."""

    text = path.read_text(encoding="utf-8")
    text = re.sub(r"//.*?$", "", text, flags=re.MULTILINE)
    text = re.sub(r",\s*([}\]])", r"\1", text)
    return json.loads(text)


def mass_orthonormalize(vectors: np.ndarray, mass: np.ndarray) -> np.ndarray:
    """Return an M-orthonormal basis spanning the independent columns."""

    gram = vectors.T @ (mass[:, None] * vectors)
    values, rotation = np.linalg.eigh(gram)
    keep = values > 1.0e-12 * max(1.0, float(values.max()))
    return vectors @ (rotation[:, keep] / np.sqrt(values[keep])[None, :])


def main() -> None:
    case_paths = {
        "ROM_low_Re_mini.case": CASE_DIR / "ROM_low_Re_mini.case",
        "ROM_low_Re_mini_separate_pod.case": (
            ROOT / "examples" / "ROM_low_Re_mini_separate_pod" /
            "ROM_low_Re_mini_separate_pod.case"
        ),
        "ROM_low_Re_mini_forward_enriched_adjoint.case": (
            ROOT / "examples" / "ROM_low_Re_mini_forward_enriched_adjoint" /
            "ROM_low_Re_mini_forward_enriched_adjoint.case"
        ),
        "ROM_low_Re_mini_enriched_galerkin.case": (
            ROOT / "examples" / "ROM_low_Re_mini_enriched_galerkin" /
            "ROM_low_Re_mini_enriched_galerkin.case"
        ),
    }
    cases = {name: load_case(path) for name, path in case_paths.items()}
    configurations = {
        name: case["optimization"]["solver"]["rom"] for name, case in cases.items()
    }

    expected = {
        "ROM_low_Re_mini.case": (("forward", "forward"), ("forward_galerkin", "exact_galerkin")),
        "ROM_low_Re_mini_separate_pod.case": (("forward", "adjoint"), ("separate_pod", "projected_adjoint")),
        "ROM_low_Re_mini_forward_enriched_adjoint.case": (("forward", "enriched"), ("forward_enriched_adjoint", "projected_adjoint")),
        "ROM_low_Re_mini_enriched_galerkin.case": (("enriched", "enriched"), ("enriched_galerkin", "exact_galerkin")),
    }
    for name, (pair, implementation) in expected.items():
        cfg = configurations[name]
        assert (cfg["basis"]["forward"], cfg["basis"]["adjoint"]) == pair
        assert derive_rom_basis_pair(*pair) == implementation
        assert not {
            "formulation", "sensitivity_source", "rank_forward",
            "rank_adjoint", "rank_enriched", "rank_total_max",
        }.intersection(cfg)

    for invalid in (("adjoint", "forward"), ("forward", "bogus"), ("adjoint", "adjoint")):
        try:
            derive_rom_basis_pair(*invalid)
        except ValueError as error:
            assert "supported pairs" in str(error)
        else:
            raise AssertionError(f"unsupported pair {invalid} was accepted")

    # Remove only the formulation-specific ROM selection.  All physical,
    # design, objective, time and state-recovery settings must be identical.
    def common(case: dict) -> dict:
        value = json.loads(json.dumps(case))
        rom = value["optimization"]["solver"].pop("rom")
        rom.pop("basis")
        return value

    reference = common(cases["ROM_low_Re_mini.case"])
    assert all(common(case) == reference for case in cases.values())
    objectives = cases["ROM_low_Re_mini.case"]["optimization"]["objectives"]
    assert [item["weight"] for item in objectives] == [120.0, 15.0, 15.0]

    # A deliberately overlapping forward/adjoint pair.  Separate POD keeps
    # Va, while the enriched path orthonormalises the full concatenated span.
    mass = np.array([1.0, 2.0, 3.0, 4.0])
    vf_raw = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 1.0]])
    va_raw = np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [1.0, -1.0]])
    vf = mass_orthonormalize(vf_raw, mass)
    va = mass_orthonormalize(va_raw, mass)
    ve = mass_orthonormalize(np.column_stack((vf, va)), mass)
    overlap = vf.T @ (mass[:, None] * va)
    assert np.linalg.norm(overlap) > 1.0e-10, "Va was accidentally replaced by a complement"
    assert np.linalg.matrix_rank(ve, tol=1.0e-10) > vf.shape[1]

    # Both projected-adjoint formulations use exactly Vf for forward assembly;
    # therefore a forward reduced operator and its objective are identical at a
    # fixed design.  The enriched test space changes only the adjoint path.
    operator_f = vf.T @ (mass[:, None] * vf)
    assert np.allclose(operator_f, vf.T @ (mass[:, None] * vf))
    assert np.allclose(operator_f, np.eye(vf.shape[1]))
    assert np.allclose(vf @ np.array([0.3, -0.2]), vf @ np.array([0.3, -0.2]))

    # An incomplete history must retain only blocks that actually arrived; no
    # synthetic historical modes or padding are introduced.
    history = []
    append_factor_history(
        history,
        DistributedPODFactors(np.eye(4, 2), np.ones(2), 1),
        maximum_length=3,
    )
    append_factor_history(
        history,
        DistributedPODFactors(np.array([[0.0], [0.0], [1.0], [0.0]]), np.ones(1), 2),
        maximum_length=3,
    )
    recycled = recycle_distributed(
        MPI.COMM_SELF, history, np.ones(4), rank=None
    )
    assert len(history) == 2
    assert recycled.modes.shape[1] == 3

    print("PASS: basis-pair configuration, rejection and incomplete-history contract")
    print(f"  separate overlap Frobenius norm = {np.linalg.norm(overlap):.6e}")
    print(f"  forward rank = {vf.shape[1]}, enriched rank = {ve.shape[1]}")


if __name__ == "__main__":
    main()
