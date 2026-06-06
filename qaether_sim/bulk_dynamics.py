from __future__ import annotations

import pandas as pd

from .contact_graph import build_contact_graph
from .forces import harmonic_forces
from .graph_atlas import analyze_snapshot


def pressure_off_dynamics(state, cfg, steps=None, lambda_phase=None):
    steps = cfg.dynamics_steps if steps is None else steps
    lambda_phase = cfg.lambda_phase if lambda_phase is None else lambda_phase
    initial_graph = build_contact_graph(state, cfg.epsilon_contact)
    initial_edges = set(tuple(sorted(e)) for e in initial_graph.edges)
    initial_analysis = analyze_snapshot(0.0, state.copy(), cfg)
    initial_t = set(initial_analysis["t_motifs"])
    initial_o = set(initial_analysis["o_motifs"])
    rows = []
    for step in range(steps + 1):
        if step % cfg.snapshot_stride == 0 or step == steps:
            analysis = analyze_snapshot(float(state.phi), state.copy(), cfg)
            edges = set(tuple(sorted(e)) for e in analysis["graph"].edges)
            t_motifs = set(analysis["t_motifs"])
            o_motifs = set(analysis["o_motifs"])
            rows.append(
                {
                    "step": step,
                    "time": step * cfg.dt,
                    "phi": state.phi,
                    "energy": harmonic_forces(state, cfg.k_core, lambda_phase)[1],
                    "edges": len(edges),
                    "triangles": len(analysis["triangles"]),
                    "squares": len(analysis["squares"]),
                    "t_motifs": len(t_motifs),
                    "o_motifs": len(o_motifs),
                    "S_E": len(initial_edges & edges) / max(len(initial_edges), 1),
                    "S_T": len(initial_t & t_motifs) / max(len(initial_t), 1),
                    "S_O": len(initial_o & o_motifs) / max(len(initial_o), 1),
                }
            )
        forces, _, _ = harmonic_forces(state, cfg.k_core, lambda_phase)
        state.velocities = forces / max(cfg.damping, 1.0e-12)
        state.positions += cfg.dt * state.velocities
        state.wrap()
    return pd.DataFrame(rows)
