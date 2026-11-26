"""Small, dependency-free ROM basis-pair configuration helpers."""

from __future__ import annotations

import math


SUPPORTED_BASIS_PAIRS = (
    ("forward", "forward"),
    ("forward", "adjoint"),
    ("forward", "enriched"),
    ("enriched", "enriched"),
)
SUPPORTED_BASIS_TRANSFORMS = ("none", "algebraic", "projected")


def derive_rom_basis_pair(forward: str, adjoint: str) -> tuple[str, str]:
    """Derive the existing ROM implementation from public basis names."""

    pair = (str(forward).strip().lower(), str(adjoint).strip().lower())
    implementation = {
        ("forward", "forward"): ("forward_galerkin", "exact_galerkin"),
        ("forward", "adjoint"): ("separate_pod", "projected_adjoint"),
        ("forward", "enriched"): (
            "forward_enriched_adjoint", "projected_adjoint"
        ),
        ("enriched", "enriched"): ("enriched_galerkin", "exact_galerkin"),
    }
    if pair not in implementation:
        supported = ", ".join(
            f"({left},{right})" for left, right in SUPPORTED_BASIS_PAIRS
        )
        raise ValueError(
            f"unsupported ROM basis pair ({pair[0]},{pair[1]}); "
            f"supported pairs are {supported}"
        )
    return implementation[pair]


def derive_rom_transform(
    rom_cfg: dict, forward: str, adjoint: str
) -> tuple[str, float, bool, str]:
    """Validate the two ROM mode/sensitivity choices and derive the source.

    Forward/adjoint space selection is handled independently by
    :func:`derive_rom_basis_pair`. This helper keeps transformed-POD source and
    legacy sensitivity variants out of the public configuration.
    """

    removed = {
        "basis_transform_model", "algebraic_sensitivity", "basis_pod_source",
        "formulation", "sensitivity_source", "rank_forward", "rank_adjoint",
        "rank_enriched", "rank_total_max",
    }.intersection(rom_cfg)
    if removed:
        raise ValueError(
            "removed ROM options: " + ", ".join(sorted(removed))
            + "; use basis.forward, basis.adjoint, state_recovery.n_modes, "
            "basis_transform and differentiate_basis"
        )
    transform = str(rom_cfg.get("basis_transform", "none")).strip().lower()
    if transform not in SUPPORTED_BASIS_TRANSFORMS:
        raise ValueError(
            "basis_transform must be none, algebraic or projected; "
            "vector_potential is not implemented yet"
        )
    tau = float(rom_cfg.get("basis_transform_tau", 1.0))
    if not math.isfinite(tau) or tau < 0.0:
        raise ValueError("basis_transform_tau must be finite and nonnegative")
    differentiate = rom_cfg.get("differentiate_basis", False)
    if isinstance(differentiate, str):
        normalized = differentiate.strip().lower()
        if normalized in ("true", "yes", "1", "on"):
            differentiate = True
        elif normalized in ("false", "no", "0", "off", ""):
            differentiate = False
        else:
            raise ValueError("differentiate_basis must be a boolean")
    elif not isinstance(differentiate, bool):
        raise ValueError("differentiate_basis must be a boolean")

    pair = (str(forward).strip().lower(), str(adjoint).strip().lower())
    _, sensitivity = derive_rom_basis_pair(*pair)
    if transform == "projected" and pair != ("forward", "forward"):
        raise ValueError(
            "projected basis transformation currently requires forward:forward"
        )
    if transform == "algebraic":
        if differentiate:
            sensitivity = "exact_algebraic"
        elif pair in (("forward", "forward"), ("enriched", "enriched")):
            sensitivity = "lagged_basis"
    elif transform == "projected" and differentiate:
        # Always differentiate the forward Galerkin objective with its
        # consistent adjoint. Independent adjoint spaces remain diagnostic.
        sensitivity = "exact_galerkin"
    return transform, tau, differentiate, sensitivity
