from __future__ import annotations

from itertools import combinations
from collections import Counter

import numpy as np

from .geometry import centered_points


def _disjoint(pairs):
    seen = set()
    for a, b in pairs:
        if a in seen or b in seen:
            return False
        seen.add(a)
        seen.add(b)
    return True


def _square_cycle_edges(square_nodes, subgraph):
    sub = subgraph.subgraph(square_nodes)
    if sub.number_of_edges() != 4:
        return None
    if sorted(dict(sub.degree()).values()) != [2, 2, 2, 2]:
        return None
    return {tuple(sorted(edge)) for edge in sub.edges}


def detect_o_motifs(
    graph,
    state,
    triangles=None,
    squares=None,
    epsilon_center: float = 0.15,
    epsilon_perp: float = 0.25,
):
    motifs = set()
    triangle_set = set(tuple(sorted(t)) for t in (triangles or []))
    square_set = set(tuple(sorted(s)) for s in (squares or []))
    candidate_nodes = [node for node, degree in graph.degree if degree >= 4]
    checked = 0
    for nodes in combinations(candidate_nodes, 6):
        checked += 1
        # Keep the reference detector bounded; large production runs should swap in
        # a dedicated octahedral candidate generator before this exact validator.
        if checked > 50000:
            break
        sub = graph.subgraph(nodes)
        if sub.number_of_edges() != 12:
            continue
        pair_triple = [tuple(sorted(p)) for p in combinations(nodes, 2) if not graph.has_edge(*p)]
        if len(pair_triple) != 3 or not _disjoint(pair_triple):
            continue
        node_pos = {node: centered_points(nodes, state.positions, state.box)[idx] for idx, node in enumerate(nodes)}
        centers = []
        axes = []
        for a, b in pair_triple:
            pa, pb = node_pos[a], node_pos[b]
            centers.append((pa + pb) / 2.0)
            axis = pb - pa
            norm = np.linalg.norm(axis)
            if norm <= 1.0e-12:
                break
            axes.append(axis / norm)
        else:
            centers = np.asarray(centers)
            if np.max(np.linalg.norm(centers - centers.mean(axis=0), axis=1)) > epsilon_center:
                continue
            if any(abs(float(np.dot(a, b))) > epsilon_perp for a, b in combinations(axes, 2)):
                continue
            if triangle_set:
                motif_triangles = [tuple(sorted(t)) for t in combinations(nodes, 3) if tuple(sorted(t)) in triangle_set]
                tri_count = len(motif_triangles)
                if tri_count < 8:
                    continue
                edge_triangle_count = Counter()
                for tri in motif_triangles:
                    for edge in combinations(tri, 2):
                        edge_triangle_count[tuple(sorted(edge))] += 1
                if any(edge_triangle_count[edge] != 2 for edge in sub.edges):
                    continue
            if square_set:
                motif_squares = [tuple(sorted(s)) for s in combinations(nodes, 4) if tuple(sorted(s)) in square_set]
                if len(motif_squares) < 3:
                    continue
                square_edge_counts = Counter()
                used_square_count = 0
                for square in motif_squares:
                    cycle_edges = _square_cycle_edges(square, sub)
                    if cycle_edges is None:
                        continue
                    used_square_count += 1
                    square_edge_counts.update(cycle_edges)
                if used_square_count < 3:
                    continue
                if set(square_edge_counts) != {tuple(sorted(edge)) for edge in sub.edges}:
                    continue
                if any(count != 1 for count in square_edge_counts.values()):
                    continue
            motifs.add(nodes)
    return sorted(motifs)
