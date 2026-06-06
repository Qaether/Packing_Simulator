from __future__ import annotations

import numpy as np


def assign_phases(state, seed: int, omega_q: float = 1.0):
    rng = np.random.default_rng(seed)
    state.theta = rng.uniform(0.0, 2.0 * np.pi, size=state.n)
    state.metadata.update({"phase_seed": seed, "omega_q": omega_q})
    return state
