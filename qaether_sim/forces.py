from __future__ import annotations

import numpy as np

from .geometry import minimum_image


def phase_stiffness(theta: np.ndarray, i: int, j: int, k_core: float, lambda_phase: float) -> float:
    if theta is None or lambda_phase == 0.0:
        return k_core
    return float(k_core * (1.0 + lambda_phase * np.cos(theta[i] - theta[j])))


def harmonic_forces(state, k_core: float = 1.0, lambda_phase: float = 0.0):
    n = state.n
    forces = np.zeros_like(state.positions)
    if n < 2:
        return forces, 0.0, 0.0

    i_idx, j_idx = np.triu_indices(n, k=1)
    displacements = minimum_image(
        state.positions[j_idx] - state.positions[i_idx],
        state.box,
    )
    distances = np.linalg.norm(displacements, axis=1)
    active = (distances > 1.0e-12) & (distances < state.ell_q)
    if not np.any(active):
        return forces, 0.0, 0.0

    i_active = i_idx[active]
    j_active = j_idx[active]
    rij = displacements[active]
    dist = distances[active]
    overlap = state.ell_q - dist
    if state.theta is None or lambda_phase == 0.0:
        stiffness = np.full_like(overlap, k_core)
    else:
        stiffness = k_core * (
            1.0 + lambda_phase * np.cos(state.theta[i_active] - state.theta[j_active])
        )

    force_vectors = -(stiffness * overlap / dist)[:, None] * rij
    np.add.at(forces, i_active, force_vectors)
    np.add.at(forces, j_active, -force_vectors)

    energy = 0.5 * np.sum(stiffness * overlap * overlap)
    pressure_like = np.sum(stiffness * overlap * dist)
    return forces, float(energy), float(pressure_like / max(state.volume, 1.0e-12))
