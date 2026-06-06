from __future__ import annotations

import pandas as pd

from .contact_graph import build_contact_graph, graph_summary
from .cycles import primitive_squares, primitive_triangles
from .motifs_O import detect_o_motifs
from .motifs_T import detect_t_motifs
from .geometry import overlap_metrics


def analyze_graph(phi, state, cfg, graph):
    triangles = primitive_triangles(graph, state, cfg.epsilon_cycle)
    squares = primitive_squares(graph, state, cfg.epsilon_cycle, cfg.epsilon_planar)
    t_motifs = detect_t_motifs(graph, state, triangles, cfg.epsilon_cycle, cfg.epsilon_volume)
    o_motifs = detect_o_motifs(graph, state, triangles, squares, cfg.epsilon_center, cfg.epsilon_perp)
    summary = graph_summary(graph)
    summary.update(
        {
            "phi": phi,
            "phi_achieved": state.phi,
            "triangles": len(triangles),
            "squares": len(squares),
            "t_motifs": len(t_motifs),
            "o_motifs": len(o_motifs),
            "t_density": len(t_motifs) / max(state.n, 1),
            "o_density": len(o_motifs) / max(state.n, 1),
            **overlap_metrics(state),
        }
    )
    return {
        "graph": graph,
        "triangles": triangles,
        "squares": squares,
        "t_motifs": t_motifs,
        "o_motifs": o_motifs,
        "summary": summary,
    }


def analyze_snapshot(phi, state, cfg):
    graph = build_contact_graph(state, cfg.epsilon_contact)
    return analyze_graph(phi, state, cfg, graph)


def build_atlas(snapshots, cfg):
    analyses = {}
    rows = []
    for phi, state in snapshots.items():
        item = analyze_snapshot(phi, state, cfg)
        analyses[phi] = item
        rows.append(item["summary"])
    return analyses, pd.DataFrame(rows).sort_values("phi").reset_index(drop=True)
