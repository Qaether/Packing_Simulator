from __future__ import annotations

import numpy as np

from .geometry import packing_box_length
from .state import QaetherState


def random_gas(n: int, phi: float, seed: int, radius: float = 0.5, ell_q: float = 1.0) -> QaetherState:
    rng = np.random.default_rng(seed)
    length = packing_box_length(n, radius, phi)
    positions = rng.random((n, 3)) * length
    velocities = np.zeros((n, 3), dtype=float)
    return QaetherState(
        positions=positions,
        velocities=velocities,
        box=np.array([length, length, length], dtype=float),
        radius=radius,
        ell_q=ell_q,
        metadata={"kind": "random_gas", "seed": seed, "phi": phi},
    )


def fcc_lattice(n: int, phi: float, radius: float = 0.5, ell_q: float = 1.0) -> QaetherState:
    basis = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
        dtype=float,
    )
    cells = int(np.ceil((n / 4.0) ** (1.0 / 3.0)))
    pts = []
    for i in range(cells):
        for j in range(cells):
            for k in range(cells):
                for b in basis:
                    pts.append(np.array([i, j, k], dtype=float) + b)
                    if len(pts) == n:
                        break
                if len(pts) == n:
                    break
            if len(pts) == n:
                break
        if len(pts) == n:
            break
    pts = np.asarray(pts)
    length = packing_box_length(n, radius, phi)
    pts = pts / cells * length
    return QaetherState(
        positions=pts,
        velocities=np.zeros((n, 3), dtype=float),
        box=np.array([length, length, length], dtype=float),
        radius=radius,
        ell_q=ell_q,
        metadata={"kind": "fcc", "phi": phi},
    )


def hcp_lattice(n: int, phi: float, radius: float = 0.5, ell_q: float = 1.0) -> QaetherState:
    # Compact periodic proxy: AB layers with triangular in-plane offsets, rescaled to target phi.
    rows = int(np.ceil(n ** (1.0 / 3.0))) + 1
    pts = []
    dz = np.sqrt(2.0 / 3.0)
    for k in range(rows):
        layer_offset = np.array([0.5, np.sqrt(3.0) / 6.0, 0.0]) if k % 2 else np.zeros(3)
        for j in range(rows):
            for i in range(rows):
                p = np.array([i + 0.5 * (j % 2), j * np.sqrt(3.0) / 2.0, k * dz]) + layer_offset
                pts.append(p)
                if len(pts) == n:
                    break
            if len(pts) == n:
                break
        if len(pts) == n:
            break
    pts = np.asarray(pts)
    mins = pts.min(axis=0)
    spans = pts.max(axis=0) - mins + 1.0
    length = packing_box_length(n, radius, phi)
    pts = (pts - mins) / spans.max() * length
    return QaetherState(
        positions=pts,
        velocities=np.zeros((n, 3), dtype=float),
        box=np.array([length, length, length], dtype=float),
        radius=radius,
        ell_q=ell_q,
        metadata={"kind": "hcp", "phi": phi},
    )
