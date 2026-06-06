from __future__ import annotations

from itertools import combinations

import numpy as np

from .geometry import centered_points, distance


def primitive_triangles(graph, state, epsilon_cycle: float = 0.10):
    out = []
    for a, b in graph.edges:
        if b < a:
            a, b = b, a
        common = set(graph.neighbors(a)) & set(graph.neighbors(b))
        for c in common:
            if c <= b:
                continue
            nodes = (a, b, c)
            pts = centered_points(nodes, state.positions, state.box)
            area = np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0])) / 2.0
            if area <= 1.0e-8:
                continue
            if all(
                abs(distance(state.positions, u, v, state.box) - state.ell_q) <= epsilon_cycle
                for u, v in combinations(nodes, 2)
            ):
                out.append(nodes)
    return sorted(set(out))


def primitive_squares(graph, state, epsilon_cycle: float = 0.10, epsilon_planar: float = 0.12):
    out = set()
    for a, b in combinations(graph.nodes, 2):
        if graph.has_edge(a, b):
            continue
        common = set(graph.neighbors(a)) & set(graph.neighbors(b))
        for c, d in combinations(common, 2):
            if graph.has_edge(c, d):
                continue
            nodes = tuple(sorted((a, b, c, d)))
            if nodes in out:
                continue
            pts = centered_points(nodes, state.positions, state.box)
            _, singular_values, _ = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
            if len(singular_values) >= 3 and singular_values[2] > epsilon_planar:
                continue
            sub = graph.subgraph(nodes)
            if all(
                abs(distance(state.positions, u, v, state.box) - state.ell_q) <= epsilon_cycle
                for u, v in sub.edges
            ):
                out.add(nodes)
    return sorted(out)
