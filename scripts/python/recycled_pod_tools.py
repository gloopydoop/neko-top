"""Distributed compact-factor recycling for the pySEMTools POD driver."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from mpi4py import MPI


@dataclass
class DistributedPODFactors:
    modes: np.ndarray
    singular_values: np.ndarray
    iteration: int
    # Small, JSON-safe evidence explaining numerical rank decisions made
    # while constructing this factor block.
    rank_diagnostics: dict = field(default_factory=dict)


def append_factor_history(
    history: list[DistributedPODFactors],
    factors: DistributedPODFactors,
    maximum_length: int = 3,
) -> None:
    """Append one outer-iteration factor block and retain a bounded window.

    The ROM contract is a three-outer-iteration window.  Keeping this small
    state transition here makes it independently testable and prevents the
    forward and adjoint histories from acquiring subtly different trimming
    rules in the streaming driver.
    """

    if maximum_length < 1:
        raise ValueError("maximum_length must be positive")
    history.append(factors)
    if len(history) > maximum_length:
        del history[:-maximum_length]


def factors_from_pod(
    pod,
    ioh,
    bm: np.ndarray,
    n_fields: int,
    iteration: int,
    scalar_metric_weight: float = 1.0,
):
    local_modes = np.asarray(getattr(pod, "u_1t", np.empty((0, 0))))
    singular = np.asarray(getattr(pod, "d_1t", np.empty(0)), dtype=np.float64)
    # A legacy POD case may include closure fields after the physical state.
    # Recycled ROM factors contain only coupled velocity/scalar fields.
    field_size = np.asarray(bm).size
    if local_modes.ndim == 2 and field_size > 0 and local_modes.shape[0] % field_size == 0:
        total_fields = local_modes.shape[0] // field_size
        if n_fields < total_fields:
            local_modes = np.concatenate(
                [local_modes[i * field_size:(i + 1) * field_size, :]
                 for i in range(n_fields)], axis=0
            )
    n_modes = min(local_modes.shape[1] if local_modes.ndim == 2 else 0, singular.size)
    if n_modes == 0:
        return DistributedPODFactors(
            np.empty((0, 0), dtype=np.float64), np.empty(0), iteration
        )
    # Store the local metric-weighted vectors expected by the global Gram
    # factorisation.  The scalar metric is applied by the caller's `mass`
    # vector exactly once; applying it here would double-weight the factors.
    if n_fields >= 4 and scalar_metric_weight <= 0.0:
        raise ValueError("scalar_metric_weight must be positive")
    local_modes = local_modes[:, :n_modes]
    return DistributedPODFactors(
        np.asarray(local_modes, dtype=np.float64),
        singular[:n_modes].copy(),
        iteration,
    )


def _global_gram(comm: MPI.Comm, left: np.ndarray, right: np.ndarray, mass: np.ndarray):
    mass = np.asarray(mass, dtype=np.float64).reshape(-1)
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.ndim != 2 or right.ndim != 2:
        raise ValueError("POD factors must be two-dimensional local matrices")
    if left.shape[0] != mass.size or right.shape[0] != mass.size:
        raise ValueError(
            "POD factor/mass local sizes disagree: "
            f"left={left.shape}, right={right.shape}, mass={mass.shape}"
        )
    local = left.T @ (mass[:, None] * right)
    return comm.allreduce(local, op=MPI.SUM)


def cross_mass_singular_values(
    comm: MPI.Comm,
    left: np.ndarray,
    right: np.ndarray,
    mass: np.ndarray,
) -> np.ndarray:
    """Return singular values of the global cross mass product ``left^T M right``."""

    gram = _global_gram(comm, left, right, mass)
    return np.linalg.svd(gram, compute_uv=False) if gram.size else np.empty(0)


def recycle_distributed(
    comm: MPI.Comm,
    factors: list[DistributedPODFactors],
    mass: np.ndarray,
    rank: int | None = None,
    tolerance: float = 1.0e-12,
) -> DistributedPODFactors:
    """Recompress ``V Sigma`` factors using only global small Gram matrices."""

    factors = [item for item in factors if item.modes.shape[1]]
    if not factors:
        return DistributedPODFactors(
            np.empty((mass.size, 0)), np.empty(0), 0,
            rank_diagnostics={
                "candidate_rank": 0,
                "recycled_rank": 0,
                "rank_cap": rank,
                "tolerance": float(tolerance),
                "gram_eigenvalues": [],
                "input_modes": [],
                "rejected_modes": [],
            },
        )
    compact = np.concatenate(
        [item.modes * item.singular_values[None, :] for item in factors], axis=1
    )
    gram = _global_gram(comm, compact, compact, mass)
    gram = 0.5 * (gram + gram.T)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    threshold = tolerance * max(1.0, values[0] if values.size else 1.0)
    keep = values > threshold
    if rank is not None:
        keep_indices = np.flatnonzero(keep)[: int(rank)]
    else:
        keep_indices = np.flatnonzero(keep)
    column_gram = _global_gram(comm, compact, compact, mass)
    input_modes = []
    column = 0
    for factor in factors:
        for mode_index in range(factor.modes.shape[1]):
            input_modes.append({
                "source_iteration": int(factor.iteration),
                "source_mode_index_zero_based": int(mode_index),
                "source_singular_value": float(factor.singular_values[mode_index]),
                "incoming_mass_norm": float(np.sqrt(max(0.0, column_gram[column, column]))),
            })
            column += 1
    rejected_modes = []
    kept_set = set(int(idx) for idx in keep_indices)
    numerically_admissible = set(int(idx) for idx in np.flatnonzero(keep))
    for eig_index in range(values.size):
        if eig_index in kept_set:
            continue
        if eig_index not in numerically_admissible:
            reason = "gram_eigenvalue_below_tolerance"
        elif rank is not None and len(keep_indices) < len(numerically_admissible):
            reason = "rank_cap"
        else:
            reason = "not_selected"
        rejected_modes.append({
            "gram_eigen_index_descending": int(eig_index),
            "gram_eigenvalue": float(values[eig_index]),
            "reason": reason,
        })
    rank_diagnostics = {
        "candidate_rank": int(compact.shape[1]),
        "recycled_rank": int(keep_indices.size),
        "rank_cap": None if rank is None else int(rank),
        "tolerance": float(tolerance),
        "gram_eigenvalue_threshold": float(threshold),
        "gram_eigenvalues_descending": [float(value) for value in values],
        "singular_values_descending": [
            float(np.sqrt(max(0.0, value))) for value in values
        ],
        "input_modes": input_modes,
        "retained_gram_eigen_indices_descending": [int(idx) for idx in keep_indices],
        "rejected_modes": rejected_modes,
    }
    if keep_indices.size == 0:
        return DistributedPODFactors(
            np.empty((mass.size, 0)), np.empty(0), factors[-1].iteration,
            rank_diagnostics=rank_diagnostics,
        )
    values = values[keep_indices]
    vectors = vectors[:, keep_indices]
    modes = compact @ (vectors / np.sqrt(values)[None, :])
    return DistributedPODFactors(
        modes, np.sqrt(values), factors[-1].iteration,
        rank_diagnostics=rank_diagnostics,
    )


def common_basis_distributed(
    comm: MPI.Comm,
    forward: DistributedPODFactors,
    adjoint: DistributedPODFactors | None,
    mass: np.ndarray,
    total_rank: int | None = None,
    tolerance: float = 1.0e-10,
    adjoint_projector: Callable[[np.ndarray], np.ndarray] | None = None,
) -> tuple[DistributedPODFactors, dict]:
    """Combine factors with global M-orthogonalisation.

    ``adjoint_projector`` is the case's explicit physical admissibility
    projector.  When supplied it is applied to both forward and adjoint
    columns before the global Gram-Schmidt pass.  Applying it to the forward
    columns also handles a streamed initial snapshot that predates Neko's
    essential-boundary update.  Neko still validates the resulting velocity
    basis against its divergence/pressure admissibility requirements.
    """

    def project_columns(columns: np.ndarray) -> np.ndarray:
        if adjoint_projector is None or columns.shape[1] == 0:
            return columns
        projected = np.column_stack([
            np.asarray(adjoint_projector(columns[:, idx])).reshape(-1)
            for idx in range(columns.shape[1])
        ])
        if projected.shape != columns.shape:
            raise ValueError(
                "the physical admissibility projector changed the POD factor shape: "
                f"before={columns.shape}, after={projected.shape}"
            )
        return projected

    q = []
    if forward.modes.shape[1]:
        q.append(project_columns(forward.modes))
    if adjoint is not None and adjoint.modes.shape[1]:
        if adjoint_projector is None:
            raise ValueError(
                "an explicit physical admissibility projector is required "
                "before adjoint enrichment"
            )
        residual = project_columns(adjoint.modes)
        if q:
            fwd = np.concatenate(q, axis=1)
            residual -= fwd @ _global_gram(comm, fwd, residual, mass)
        q.append(residual)
    if not q:
        return DistributedPODFactors(np.empty((mass.size, 0)), np.empty(0), 0), {
            "mass_gram_error": 0.0,
            "admissibility_projected": float(adjoint_projector is not None),
            "candidate_rank": 0,
            "retained_rank": 0,
            "input_modes": [],
            "rejected_modes": [],
            "gram_singular_values": [],
        }
    vectors = np.concatenate(q, axis=1)
    provenance = []
    for family, factors in (("forward", forward), ("adjoint", adjoint)):
        if factors is None:
            continue
        for mode_index in range(factors.modes.shape[1]):
            provenance.append({
                "source_family": family,
                "source_iteration": int(factors.iteration),
                "source_mode_index_zero_based": int(mode_index),
            })
    accepted = []
    accepted_provenance = []
    rejected_modes = []
    pre_projection_norms = []
    post_projection_norms = []
    orthogonal_residual_norms = []
    for index in range(vectors.shape[1]):
        vector = vectors[:, index].copy()
        pre_projection_norm_sq = float(comm.allreduce(
            np.dot(vector, mass * vector), op=MPI.SUM
        ))
        pre_projection_norms.append(float(np.sqrt(max(0.0, pre_projection_norm_sq))))
        if adjoint_projector is not None:
            vector = np.asarray(adjoint_projector(vector)).reshape(-1)
        post_projection_norm_sq = float(comm.allreduce(
            np.dot(vector, mass * vector), op=MPI.SUM
        ))
        post_projection_norms.append(float(np.sqrt(max(0.0, post_projection_norm_sq))))
        if accepted:
            basis = np.column_stack(accepted)
            vector -= basis @ _global_gram(comm, basis, vector[:, None], mass).ravel()
            vector -= basis @ _global_gram(comm, basis, vector[:, None], mass).ravel()
        norm_sq = float(comm.allreduce(np.dot(vector, mass * vector), op=MPI.SUM))
        orthogonal_residual_norms.append(float(np.sqrt(max(0.0, norm_sq))))
        if norm_sq <= tolerance * tolerance:
            rejected_modes.append({
                **provenance[index],
                "candidate_index_zero_based": int(index),
                "reason": (
                    "admissibility_projection_zeroed_mode"
                    if post_projection_norm_sq <= tolerance * tolerance
                    else "dependent_after_mass_orthogonalisation"
                ),
                "pre_projection_norm": pre_projection_norms[-1],
                "post_projection_norm": post_projection_norms[-1],
                "orthogonal_residual_norm": orthogonal_residual_norms[-1],
            })
            continue
        accepted.append(vector / np.sqrt(norm_sq))
        accepted_provenance.append({
            **provenance[index],
            "candidate_index_zero_based": int(index),
        })
        if total_rank is not None and len(accepted) >= int(total_rank):
            for rejected_index in range(index + 1, vectors.shape[1]):
                rejected_modes.append({
                    **provenance[rejected_index],
                    "candidate_index_zero_based": int(rejected_index),
                    "reason": "rank_cap",
                })
            break
    basis = np.column_stack(accepted) if accepted else np.empty((mass.size, 0))
    gram = _global_gram(comm, basis, basis, mass) if accepted else np.empty((0, 0))
    error = float(np.max(np.abs(gram - np.eye(len(accepted))))) if accepted else 0.0
    basis_values = _global_gram(comm, basis, basis, mass) if accepted else np.empty((0, 0))
    diagnostics = {
        "mass_gram_error": error,
        "admissibility_projected": float(adjoint_projector is not None),
        "candidate_rank": int(vectors.shape[1]),
        "retained_rank": int(len(accepted)),
        "rank_cap": None if total_rank is None else int(total_rank),
        "input_modes": provenance,
        "accepted_modes": accepted_provenance,
        "rejected_modes": rejected_modes,
        "pre_projection_norms": pre_projection_norms,
        "post_projection_norms": post_projection_norms,
        "pre_normalisation_orthogonal_residual_norms": orthogonal_residual_norms,
        "gram_singular_values": [
            float(value) for value in np.linalg.svd(basis_values, compute_uv=False)
        ],
    }
    return DistributedPODFactors(
        basis, np.ones(len(accepted)), forward.iteration,
        rank_diagnostics=diagnostics,
    ), diagnostics
