from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

from .config import ExperimentConfig
from .forces import harmonic_forces
from .geometry import overlap_metrics, packing_box_length


def rescale_to_phi(state, phi: float):
    new_length = packing_box_length(state.n, state.radius, phi)
    old_box = state.box.copy()
    state.positions = state.positions / old_box * new_length
    state.box = np.array([new_length, new_length, new_length], dtype=float)
    state.wrap()


def relax(state, cfg: ExperimentConfig, steps: int = None, lambda_phase: float = None):
    steps = cfg.relax_steps if steps is None else steps
    lambda_phase = cfg.lambda_phase if lambda_phase is None else lambda_phase
    rows = []
    prev_positions = state.positions.copy()
    for step in range(steps):
        forces, energy, pressure = harmonic_forces(state, cfg.k_core, lambda_phase)
        state.velocities = forces / max(cfg.damping, 1.0e-12)
        state.positions += cfg.dt * state.velocities
        state.wrap()
        disp = state.positions - prev_positions
        msd = float(np.mean(np.sum(disp * disp, axis=1)))
        rows.append({"step": step, "energy": energy, "pressure": pressure, "msd_step": msd})
        prev_positions = state.positions.copy()
    return pd.DataFrame(rows)


def compression_sweep(state, cfg: ExperimentConfig, phi_targets: Iterable[float] = None):
    phi_targets = list(cfg.phi_targets if phi_targets is None else phi_targets)
    snapshots: Dict[float, object] = {}
    rows: List[dict] = []
    for phi in phi_targets:
        rescale_to_phi(state, phi)
        trace = relax(state, cfg)
        final = trace.iloc[-1].to_dict()
        overlap = overlap_metrics(state)
        hard_like = overlap["max_overlap"] <= 1.0e-3 and overlap["overlap_energy_per_vertex"] <= 1.0e-6
        final.update(
            {
                "phi_target": phi,
                "phi_achieved": state.phi,
                "phi": state.phi,
                "n": state.n,
                "state_label": "hard_sphere_like" if hard_like else "soft_overcompressed_or_stressed",
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
                **overlap,
            }
        )
        snapshots[phi] = snap
    return snapshots, pd.DataFrame(rows)
