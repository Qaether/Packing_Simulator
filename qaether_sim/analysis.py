from __future__ import annotations

import json
import os

import pandas as pd

from .bulk_dynamics import pressure_off_dynamics
from .compression import compression_sweep
from .config import ExperimentConfig
from .graph_atlas import analyze_snapshot, build_atlas
from .initial_conditions import random_gas
from .jamming import classify_compression_trace
from .lattice_benchmarks import build_lattice_benchmarks, close_packing_phi
from .perturbation import random_micro_jitter
from .phase import assign_phases


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def run_smoke_pipeline(out_dir: str, n: int = 32, seed: int = 0) -> dict:
    ensure_dir(out_dir)
    cfg = ExperimentConfig(
        phi_targets=[0.20, 0.35, 0.50],
        relax_steps=20,
        dynamics_steps=30,
        snapshot_stride=10,
        dt=0.01,
    )
    state = random_gas(n=n, phi=0.10, seed=seed, radius=cfg.radius, ell_q=cfg.ell_q)
    snapshots, compression_df = compression_sweep(state, cfg)
    compression_df.to_csv(os.path.join(out_dir, "energy_curve.csv"), index=False)
    atlas, atlas_df = build_atlas(snapshots, cfg)
    atlas_df.to_csv(os.path.join(out_dir, "compression_graph_motif_atlas.csv"), index=False)
    pd.DataFrame(classify_compression_trace(compression_df)).to_csv(
        os.path.join(out_dir, "state_selection_report.csv"), index=False
    )

    selected_phi = max(snapshots)
    selected = snapshots[selected_phi].copy()
    selected.save_h5(os.path.join(out_dir, "selected_state_random.h5"))
    dyn_df = pressure_off_dynamics(selected.copy(), cfg)
    dyn_df.to_csv(os.path.join(out_dir, "pressure_off_summary.csv"), index=False)

    phase_state = assign_phases(selected.copy(), seed=seed + 100, omega_q=cfg.omega_q)
    phase_df = pressure_off_dynamics(phase_state, cfg, lambda_phase=0.3)
    phase_df.to_csv(os.path.join(out_dir, "phase_sweep_summary.csv"), index=False)

    perturbed = random_micro_jitter(selected.copy(), amplitude=1.0e-3, seed=seed + 200)
    pert_df = pressure_off_dynamics(perturbed, cfg)
    pert_df.to_csv(os.path.join(out_dir, "perturbed_pressure_off_summary.csv"), index=False)

    lattice_rows = []
    lattice_reference_phi = close_packing_phi() * 0.99
    for name, lattice_state in build_lattice_benchmarks(n, lattice_reference_phi, cfg.radius, cfg.ell_q).items():
        lattice_state.metadata.update(
            {
                "benchmark_type": "ideal_close_packing_reference",
                "comparison_phi": selected_phi,
                "phi_target": lattice_reference_phi,
                "phi_achieved": lattice_state.phi,
            }
        )
        item = analyze_snapshot(lattice_state.phi, lattice_state, cfg)
        row = dict(item["summary"])
        row["kind"] = name
        row["benchmark_type"] = "ideal_close_packing_reference"
        row["comparison_phi"] = selected_phi
        lattice_rows.append(row)
        lattice_state.save_h5(os.path.join(out_dir, f"{name}_state_phi_{lattice_state.phi:.3f}.h5"))
    pd.DataFrame(lattice_rows).to_csv(os.path.join(out_dir, "hcp_fcc_motif_summary.csv"), index=False)

    summary = {
        "n": n,
        "seed": seed,
        "phi_targets": cfg.phi_targets,
        "atlas_rows": int(len(atlas_df)),
        "pressure_off_rows": int(len(dyn_df)),
        "phase_rows": int(len(phase_df)),
        "perturbed_rows": int(len(pert_df)),
    }
    with open(os.path.join(out_dir, "smoke_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    return summary
