from __future__ import annotations

import numpy as np


def minimum_image(displacement: np.ndarray, box: np.ndarray) -> np.ndarray:
    return displacement - box * np.round(displacement / box)


def pair_displacement(positions: np.ndarray, i: int, j: int, box: np.ndarray) -> np.ndarray:
    return minimum_image(positions[j] - positions[i], box)


def distance(positions: np.ndarray, i: int, j: int, box: np.ndarray) -> float:
    return float(np.linalg.norm(pair_displacement(positions, i, j, box)))


def packing_box_length(n: int, radius: float, phi: float) -> float:
    sphere_volume = 4.0 * np.pi * radius**3 / 3.0
    return float((n * sphere_volume / phi) ** (1.0 / 3.0))


def tetra_volume(points: np.ndarray) -> float:
    a, b, c, d = points
    return float(abs(np.dot(b - a, np.cross(c - a, d - a))) / 6.0)


def centered_points(indices, positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Unwrap candidate vertices into a local minimum-image patch."""
    indices = list(indices)
    anchor = positions[indices[0]]
    pts = [anchor]
    for idx in indices[1:]:
        pts.append(anchor + minimum_image(positions[idx] - anchor, box))
    return np.asarray(pts)


def overlap_metrics(state) -> dict:
    numerical_tolerance = 1.0e-12
    overlaps = []
    overlapping_pairs = 0
    pair_count = state.n * (state.n - 1) // 2
    for i in range(state.n - 1):
        for j in range(i + 1, state.n):
            d = distance(state.positions, i, j, state.box)
            overlap = state.ell_q - d
            if overlap > numerical_tolerance:
                overlapping_pairs += 1
                overlaps.append(overlap)
    if overlaps:
        arr = np.asarray(overlaps, dtype=float)
        max_overlap = float(arr.max())
        mean_overlap = float(arr.mean())
        energy = float(0.5 * np.sum(arr * arr))
    else:
        max_overlap = 0.0
        mean_overlap = 0.0
        energy = 0.0
    return {
        "max_overlap": max_overlap,
        "mean_overlap": mean_overlap,
        "overlap_energy_per_vertex": energy / max(state.n, 1),
        "fraction_overlapping_pairs": overlapping_pairs / max(pair_count, 1),
    }
