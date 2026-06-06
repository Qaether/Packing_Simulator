from __future__ import annotations

from itertools import combinations

from .geometry import centered_points, distance, tetra_volume


def detect_t_motifs(graph, state, triangles=None, epsilon_cycle: float = 0.10, epsilon_volume: float = 1.0e-4):
    motifs = []
    triangle_set = set(tuple(sorted(t)) for t in triangles) if triangles is not None else None
    if triangle_set is None:
        triangle_set = {
            tuple(sorted((a, b, c)))
            for a, b in graph.edges
            for c in set(graph.neighbors(a)) & set(graph.neighbors(b))
        }

    candidates = set()
    for triangle in triangle_set:
        common = set(graph.neighbors(triangle[0]))
        common &= set(graph.neighbors(triangle[1]))
        common &= set(graph.neighbors(triangle[2]))
        for fourth in common:
            candidates.add(tuple(sorted((*triangle, fourth))))

    for nodes in sorted(candidates):
        if any(abs(distance(state.positions, a, b, state.box) - state.ell_q) > epsilon_cycle for a, b in combinations(nodes, 2)):
            continue
        pts = centered_points(nodes, state.positions, state.box)
        if tetra_volume(pts) <= epsilon_volume:
            continue
        if all(tuple(sorted(face)) in triangle_set for face in combinations(nodes, 3)):
            motifs.append(tuple(sorted(nodes)))
    return motifs
