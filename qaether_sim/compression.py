from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

from .config import ExperimentConfig
from .forces import harmonic_forces
from .geometry import minimum_image, overlap_metrics, packing_box_length
from .contact_graph import build_contact_graph
from .jamming import compute_z_rattler_removed


def rescale_to_phi(state, phi: float):
    new_length = packing_box_length(state.n, state.radius, phi)
    old_box = state.box.copy()
    state.positions = state.positions / old_box * new_length
    state.box = np.array([new_length, new_length, new_length], dtype=float)
    state.wrap()


def relax(state, cfg: ExperimentConfig, steps: int = None, lambda_phase: float = None):
    steps = cfg.relax_steps if steps is None else steps
    lambda_phase = cfg.lambda_phase if lambda_phase is None else lambda_phase
    relax_dt = cfg.dt if cfg.relax_dt is None else cfg.relax_dt
    rows = []
    converged_steps = 0

    # Evaluate initial forces
    forces, energy, pressure = harmonic_forces(state, cfg.k_core, lambda_phase)
    state.velocities = forces / max(cfg.damping, 1.0e-12)
    prev_positions = state.positions.copy()

    for step in range(steps):
        # Update coordinates using existing velocities
        state.positions += relax_dt * state.velocities
        state.wrap()
        disp = minimum_image(state.positions - prev_positions, state.box)
        msd = float(np.mean(np.sum(disp * disp, axis=1)))
        prev_positions = state.positions.copy()

        # Re-evaluate forces/energy at the new positions
        forces, energy, pressure = harmonic_forces(state, cfg.k_core, lambda_phase)
        state.velocities = forces / max(cfg.damping, 1.0e-12)
        force_residual = float(np.sqrt(np.mean(np.sum(forces * forces, axis=1))))

        rows.append({
            "step": step,
            "energy": energy,
            "pressure": pressure,
            "msd_step": msd,
            "force_balance_residual": force_residual,
        })

        if (
            step + 1 >= min(cfg.relax_min_steps, steps)
            and msd <= cfg.relax_msd_tolerance
            and force_residual <= cfg.relax_force_tolerance
        ):
            converged_steps += 1
            if converged_steps >= cfg.relax_convergence_window:
                break
        else:
            converged_steps = 0
    return pd.DataFrame(rows)


def compression_sweep(state, cfg: ExperimentConfig, phi_targets: Iterable[float] = None):
    phi_targets = list(cfg.phi_targets if phi_targets is None else phi_targets)
    snapshots: Dict[float, object] = {}
    rows: List[dict] = []
    for phi in phi_targets:
        rescale_to_phi(state, phi)
        trace = relax(state, cfg)
        final = trace.iloc[-1].to_dict()

        # Calculate force balance residual on final relaxed state
        force_res = float(final["force_balance_residual"])

        # Calculate Z_rattler_removed
        graph = build_contact_graph(state, cfg.epsilon_contact)
        z_rattler = compute_z_rattler_removed(graph)

        overlap = overlap_metrics(state)
        hard_like = overlap["max_overlap"] <= 1.0e-3 and overlap["overlap_energy_per_vertex"] <= 1.0e-6
        final.update(
            {
                "phi_target": phi,
                "phi_achieved": state.phi,
                "phi": state.phi,
                "n": state.n,
                "relaxation_steps_used": len(trace),
                "state_label": "hard_sphere_like" if hard_like else "soft_overcompressed_or_stressed",
                "force_balance_residual": force_res,
                "z_rattler_removed": z_rattler,
                **overlap,
            }
        )
        rows.append(final)
        snap = state.copy()
        snap.metadata.update(
            {
                "phi_target": phi,
                "phi_achieved": state.phi,
                "stage": "compression",
                "state_label": final["state_label"],
                "force_balance_residual": force_res,
                "z_rattler_removed": z_rattler,
                **overlap,
            }
        )
        snapshots[phi] = snap
    return snapshots, pd.DataFrame(rows)
