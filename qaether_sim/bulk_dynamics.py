from __future__ import annotations

import h5py
import numpy as np
import pandas as pd

from .contact_graph import build_hysteretic_contact_graph
from .forces import harmonic_forces
from .geometry import minimum_image, overlap_metrics
from .graph_atlas import analyze_graph


def _survival(initial: set, current: set) -> float:
    if not initial:
        return float("nan")
    return len(initial & current) / len(initial)


def _classify_dynamics(state, metrics: dict) -> str:
    is_stressed = (
        metrics["max_overlap"] > 1.0e-3
        or metrics["overlap_energy_per_vertex"] > 1.0e-6
    )
    is_perturbed = state.metadata.get("perturbation") is not None
    if is_stressed:
        return "A. Residual-stress relaxation"
    if is_perturbed:
        return "B. Perturbation-response"
    return "none"


def pressure_off_dynamics(state, cfg, steps=None, lambda_phase=None, traj_path: str = None):
    steps = cfg.dynamics_steps if steps is None else steps
    lambda_phase = cfg.lambda_phase if lambda_phase is None else lambda_phase

    metrics = overlap_metrics(state)
    dynamics_type = _classify_dynamics(state, metrics)
    state.metadata["dynamics_type"] = dynamics_type

    current_graph = build_hysteretic_contact_graph(
        state,
        previous_graph=None,
        epsilon_on=cfg.epsilon_contact_on,
        epsilon_off=cfg.epsilon_contact_off,
    )
    initial_analysis = analyze_graph(float(state.phi), state, cfg, current_graph)
    initial_edges = {tuple(sorted(edge)) for edge in current_graph.edges}
    initial_triangles = set(initial_analysis["triangles"])
    initial_squares = set(initial_analysis["squares"])
    initial_t = set(initial_analysis["t_motifs"])
    initial_o = set(initial_analysis["o_motifs"])

    previous_snapshot = {
        "edges": initial_edges,
        "triangles": initial_triangles,
        "squares": initial_squares,
        "t": initial_t,
        "o": initial_o,
    }
    cumulative_displacement = np.zeros_like(state.positions)

    if traj_path:
        with h5py.File(traj_path, "w") as h5:
            h5.attrs["dynamics_type"] = dynamics_type
            h5.attrs["radius"] = state.radius
            h5.attrs["ell_q"] = state.ell_q

    rows = []
    for step in range(steps + 1):
        if step > 0:
            old_positions = state.positions.copy()
            forces, _, _ = harmonic_forces(state, cfg.k_core, lambda_phase)
            state.velocities = forces / max(cfg.damping, 1.0e-12)
            state.positions += cfg.dt * state.velocities
            state.wrap()
            cumulative_displacement += minimum_image(
                state.positions - old_positions, state.box
            )
            current_graph = build_hysteretic_contact_graph(
                state,
                previous_graph=current_graph,
                epsilon_on=cfg.epsilon_contact_on,
                epsilon_off=cfg.epsilon_contact_off,
            )

        if step % cfg.snapshot_stride != 0 and step != steps:
            continue

        analysis = analyze_graph(float(state.phi), state, cfg, current_graph)
        edges = {tuple(sorted(edge)) for edge in current_graph.edges}
        triangles = set(analysis["triangles"])
        squares = set(analysis["squares"])
        t_motifs = set(analysis["t_motifs"])
        o_motifs = set(analysis["o_motifs"])

        rows.append(
            {
                "step": step,
                "time": step * cfg.dt,
                "phi": state.phi,
                "energy": harmonic_forces(state, cfg.k_core, lambda_phase)[1],
                "msd": float(
                    np.mean(np.sum(cumulative_displacement**2, axis=1))
                ),
                "edges": len(edges),
                "triangles": len(triangles),
                "squares": len(squares),
                "t_motifs": len(t_motifs),
                "o_motifs": len(o_motifs),
                "initial_edges": len(initial_edges),
                "initial_triangles": len(initial_triangles),
                "initial_squares": len(initial_squares),
                "initial_t_motifs": len(initial_t),
                "initial_o_motifs": len(initial_o),
                "S_E": _survival(initial_edges, edges),
                "S_T": _survival(initial_t, t_motifs),
                "S_O": _survival(initial_o, o_motifs),
                "edge_edit_distance": len(initial_edges ^ edges),
                "edge_birth": len(edges - previous_snapshot["edges"]),
                "edge_death": len(previous_snapshot["edges"] - edges),
                "tri_birth": len(triangles - previous_snapshot["triangles"]),
                "tri_death": len(previous_snapshot["triangles"] - triangles),
                "square_birth": len(squares - previous_snapshot["squares"]),
                "square_death": len(previous_snapshot["squares"] - squares),
                "t_birth": len(t_motifs - previous_snapshot["t"]),
                "t_death": len(previous_snapshot["t"] - t_motifs),
                "o_birth": len(o_motifs - previous_snapshot["o"]),
                "o_death": len(previous_snapshot["o"] - o_motifs),
                "dynamics_type": dynamics_type,
                "contact_rule": "hysteretic",
            }
        )

        previous_snapshot = {
            "edges": edges,
            "triangles": triangles,
            "squares": squares,
            "t": t_motifs,
            "o": o_motifs,
        }

        if traj_path:
            with h5py.File(traj_path, "a") as h5:
                group = h5.create_group(f"step_{step}")
                group.create_dataset("positions", data=state.positions)
                group.create_dataset("velocities", data=state.velocities)
                group.create_dataset("box", data=state.box)
                group.create_dataset(
                    "cumulative_displacement", data=cumulative_displacement
                )
                if state.theta is not None:
                    group.create_dataset("theta", data=state.theta)

    return pd.DataFrame(rows)
