"""Window-statistics closure algebra used by the statistical POD pipeline.

The routines in this module are deliberately independent of pySEMTools.  They
operate on local field arrays, which makes the sign, tensor convention and
regularisation testable without starting an MPI/Neko run.
"""

from __future__ import annotations

import numpy as np


def reynolds_stress(raw_moments: np.ndarray, mean_velocity: np.ndarray) -> np.ndarray:
    """Return centered Reynolds stresses from raw second moments.

    Components are ordered ``uu,vv,ww,uv,uw,vw`` and the returned tensor has
    shape ``(3, 3, ...)``.
    """
    raw = np.asarray(raw_moments)
    mean = np.asarray(mean_velocity)
    if raw.shape[0] != 6 or mean.shape[0] != 3:
        raise ValueError("expected six raw moments and three mean velocities")
    out = np.zeros((3, 3) + raw.shape[1:], dtype=np.result_type(raw, mean))
    out[0, 0] = raw[0] - mean[0] * mean[0]
    out[1, 1] = raw[1] - mean[1] * mean[1]
    out[2, 2] = raw[2] - mean[2] * mean[2]
    out[0, 1] = out[1, 0] = raw[3] - mean[0] * mean[1]
    out[0, 2] = out[2, 0] = raw[4] - mean[0] * mean[2]
    out[1, 2] = out[2, 1] = raw[5] - mean[1] * mean[2]
    return out


def fit_boussinesq_nut(
    reynolds: np.ndarray,
    gradient_sum: np.ndarray,
    regularisation: float = 1.0e-12,
    maximum: float | None = None,
    clip_negative: bool = True,
) -> tuple[np.ndarray, dict[str, float]]:
    """Fit kinematic eddy viscosity using ``-R_d = nut * A``.

    ``gradient_sum`` is the full symmetric gradient A, with shape ``(3,3,...)``.
    The returned diagnostics are fractions over the final array.
    """
    r = np.asarray(reynolds)
    a = np.asarray(gradient_sum)
    if r.shape[0:2] != (3, 3) or a.shape != r.shape:
        raise ValueError("reynolds and gradient_sum must both have shape (3,3,...)")
    k = 0.5 * (r[0, 0] + r[1, 1] + r[2, 2])
    rd = r.copy()
    rd[0, 0] -= 2.0 * k / 3.0
    rd[1, 1] -= 2.0 * k / 3.0
    rd[2, 2] -= 2.0 * k / 3.0
    aa = np.sum(a * a, axis=(0, 1))
    numerator = np.sum(rd * a, axis=(0, 1))
    nut = -numerator / (aa + regularisation * regularisation)
    valid = aa > regularisation * regularisation
    raw_nut = nut.copy()
    negative = raw_nut < 0.0
    clipped = np.zeros_like(negative, dtype=bool)
    if clip_negative:
        nut = np.maximum(nut, 0.0)
        clipped |= raw_nut < 0.0
    if maximum is not None:
        if maximum < 0.0:
            raise ValueError("maximum eddy viscosity must be non-negative")
        clipped |= nut > maximum
        nut = np.minimum(nut, maximum)
    denom = np.maximum(np.sum(rd * rd, axis=(0, 1)), regularisation**2)
    residual = np.sqrt(np.sum((rd + nut * a) ** 2, axis=(0, 1)) / denom)
    diagnostics = {
        "valid_fraction": float(np.mean(valid)),
        "negative_fraction_raw": float(np.mean(negative)),
        "clipped_fraction": float(np.mean(clipped)),
        "mean_fit_relative_residual": float(np.mean(residual)),
        "max_nut": float(np.max(nut)) if nut.size else 0.0,
    }
    return nut, diagnostics


def fit_scalar_gammat(
    flux: np.ndarray,
    scalar_gradient: np.ndarray,
    regularisation: float = 1.0e-12,
    maximum: float | None = None,
    clip_negative: bool = True,
) -> tuple[np.ndarray, dict[str, float]]:
    """Fit turbulent scalar diffusivity from ``q = -Gamma_T grad(Phi)``."""
    q = np.asarray(flux)
    grad = np.asarray(scalar_gradient)
    if q.shape[0] != 3 or grad.shape != q.shape:
        raise ValueError("flux and scalar_gradient must have shape (3,...)")
    gg = np.sum(grad * grad, axis=0)
    gamma = -np.sum(q * grad, axis=0) / (gg + regularisation**2)
    raw_gamma = gamma.copy()
    negative = raw_gamma < 0.0
    clipped = np.zeros_like(negative, dtype=bool)
    if clip_negative:
        gamma = np.maximum(gamma, 0.0)
        clipped |= raw_gamma < 0.0
    if maximum is not None:
        clipped |= gamma > maximum
        gamma = np.minimum(gamma, maximum)
    return gamma, {
        "negative_fraction_raw": float(np.mean(negative)),
        "valid_fraction": float(np.mean(gg > regularisation**2)),
        "clipped_fraction": float(np.mean(clipped)),
        "mean_fit_relative_residual": float(
            np.mean(np.abs(q + gamma[None, ...] * grad)) /
            max(float(np.mean(np.abs(q))), regularisation)
        ),
        "max_gammat": float(np.max(gamma)) if gamma.size else 0.0,
    }
