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


def _find_grid_dimensions(n_target: int, cell_size: int = 4):
    if n_target <= 0 or n_target % cell_size != 0:
        raise ValueError(f"lattice size N must be a positive multiple of {cell_size}")

    cell_count = n_target // cell_size
    candidates = []
    for nx in range(1, cell_count + 1):
        if cell_count % nx:
            continue
        remaining = cell_count // nx
        for ny in range(1, remaining + 1):
            if remaining % ny:
                continue
            nz = remaining // ny
            dims = tuple(sorted((nx, ny, nz), reverse=True))
            candidates.append((dims[0] / dims[-1], dims))
    return min(candidates)[1]


def fcc_lattice(n: int, phi: float, radius: float = 0.5, ell_q: float = 1.0) -> QaetherState:
    """
    Construct FCC lattice using an orthogonal 4-particle unit cell basis, scaled to target phi.
    If nearest-neighbor distance a_nn is used:
    lattice constant a = sqrt(2) * a_nn
    box = [nx * a, ny * a, nz * a]
    """
    # Orthogonal basis coordinates normalized to [0, 1) unit box
    basis = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
        dtype=float,
    )
    nx, ny, nz = _find_grid_dimensions(n)

    sphere_volume = 4.0 * np.pi * radius**3 / 3.0
    a = (4.0 * sphere_volume / phi) ** (1.0 / 3.0)

    pts = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for b in basis:
                    pt = (np.array([i, j, k], dtype=float) + b) * a
                    pts.append(pt)
    pts = np.asarray(pts)
    if len(pts) != n:
        raise RuntimeError("FCC replication did not generate exactly N vertices")
    box = np.array([nx * a, ny * a, nz * a], dtype=float)

    return QaetherState(
        positions=pts,
        velocities=np.zeros((n, 3), dtype=float),
        box=box,
        radius=radius,
        ell_q=ell_q,
        metadata={"kind": "fcc", "phi": phi},
    )


def hcp_lattice(n: int, phi: float, radius: float = 0.5, ell_q: float = 1.0) -> QaetherState:
    """
    Construct HCP lattice using an orthogonal 4-particle unit cell basis, scaled to target phi.
    Periodicity scale: [a, a * sqrt(3), a * sqrt(8/3)]
    """
    # Relative basis positions normalized to [1.0, 1.0, 1.0] relative periodic dimensions
    basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 1.0/3.0, 0.5],
        [0.5, 5.0/6.0, 0.5]
    ], dtype=float)

    nx, ny, nz = _find_grid_dimensions(n)

    # Unit dimensions proportions: [1, sqrt(3), sqrt(8/3)]
    aspect_x = 1.0
    aspect_y = np.sqrt(3.0)
    aspect_z = np.sqrt(8.0 / 3.0)

    sphere_volume = 4.0 * np.pi * radius**3 / 3.0
    cell_aspect_volume = aspect_x * aspect_y * aspect_z
    a = (4.0 * sphere_volume / (phi * cell_aspect_volume)) ** (1.0 / 3.0)

    scale_vector = np.array([a * aspect_x, a * aspect_y, a * aspect_z], dtype=float)

    pts = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for b in basis:
                    pt = (np.array([i, j, k], dtype=float) + b) * scale_vector
                    pts.append(pt)
    pts = np.asarray(pts)
    if len(pts) != n:
        raise RuntimeError("HCP replication did not generate exactly N vertices")
    box = scale_vector * np.array([nx, ny, nz], dtype=float)

    return QaetherState(
        positions=pts,
        velocities=np.zeros((n, 3), dtype=float),
        box=box,
        radius=radius,
        ell_q=ell_q,
        metadata={"kind": "hcp", "phi": phi},
    )
