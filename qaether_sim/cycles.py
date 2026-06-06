from __future__ import annotations

from itertools import combinations

import numpy as np

from .geometry import centered_points, distance


def primitive_triangles(graph, state, epsilon_cycle: float = 0.10):
    out = []
    for nodes in combinations(graph.nodes, 3):
        if graph.subgraph(nodes).number_of_edges() != 3:
            continue
        pts = centered_points(nodes, state.positions, state.box)
        area = np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0])) / 2.0
        if area <= 1.0e-8:
            continue
        if all(abs(distance(state.positions, a, b, state.box) - state.ell_q) <= epsilon_cycle for a, b in combinations(nodes, 2)):
            out.append(tuple(sorted(nodes)))
    return out


def primitive_squares(graph, state, epsilon_cycle: float = 0.10, epsilon_planar: float = 0.12):
    out = set()
    for nodes in combinations(graph.nodes, 4):
        sub = graph.subgraph(nodes)
        if sub.number_of_edges() != 4:
            continue
        if sorted(dict(sub.degree()).values()) != [2, 2, 2, 2]:
            continue
        pts = centered_points(nodes, state.positions, state.box)
        _, s, _ = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
        if len(s) >= 3 and s[2] > epsilon_planar:
            continue
        ok = True
        for a, b in sub.edges:
            if abs(distance(state.positions, a, b, state.box) - state.ell_q) > epsilon_cycle:
                ok = False
                break
        if ok:
            out.add(tuple(sorted(nodes)))
    return sorted(out)
