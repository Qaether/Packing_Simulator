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
    energy = 0.0
    pressure_like = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            rij = minimum_image(state.positions[j] - state.positions[i], state.box)
            dist = float(np.linalg.norm(rij))
            if dist <= 1.0e-12 or dist >= state.ell_q:
                continue
            overlap = state.ell_q - dist
            kij = phase_stiffness(state.theta, i, j, k_core, lambda_phase)
            fmag = kij * overlap / dist
            fij = -fmag * rij
            forces[i] += fij
            forces[j] -= fij
            energy += 0.5 * kij * overlap * overlap
            pressure_like += float(np.linalg.norm(fij) * dist)
    return forces, float(energy), float(pressure_like / max(state.volume, 1.0e-12))
