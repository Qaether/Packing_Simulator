from __future__ import annotations

from itertools import combinations

from .geometry import centered_points, distance, tetra_volume


def detect_t_motifs(graph, state, triangles=None, epsilon_cycle: float = 0.10, epsilon_volume: float = 1.0e-4):
    motifs = []
    triangle_set = set(tuple(sorted(t)) for t in (triangles or []))
    for nodes in combinations(graph.nodes, 4):
        if graph.subgraph(nodes).number_of_edges() != 6:
            continue
        if any(abs(distance(state.positions, a, b, state.box) - state.ell_q) > epsilon_cycle for a, b in combinations(nodes, 2)):
            continue
        pts = centered_points(nodes, state.positions, state.box)
        if tetra_volume(pts) <= epsilon_volume:
            continue
        faces_ok = True
        if triangle_set:
            for face in combinations(nodes, 3):
                if tuple(sorted(face)) not in triangle_set:
                    faces_ok = False
                    break
        if faces_ok:
            motifs.append(tuple(sorted(nodes)))
    return motifs
