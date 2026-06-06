from __future__ import annotations

import numpy as np


def random_micro_jitter(state, amplitude: float, seed: int):
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=state.positions.shape)
    norms = np.linalg.norm(directions, axis=1)
    directions = directions / np.maximum(norms[:, None], 1.0e-12)
    state.positions += amplitude * directions
    state.wrap()
    state.metadata.update({"perturbation": "random_micro_jitter", "amplitude": amplitude, "seed": seed})
    return state
