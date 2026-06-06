from __future__ import annotations

import json
import os
import shutil
from typing import List

import pandas as pd

from .bulk_dynamics import pressure_off_dynamics
from .compression import compression_sweep
from .config import ExperimentConfig
from .graph_atlas import analyze_snapshot, build_atlas
from .initial_conditions import random_gas
from .jamming import classify_compression_trace
from .lattice_benchmarks import build_lattice_benchmarks, build_toto_labeled_benchmarks, close_packing_phi
from .perturbation import random_micro_jitter
from .phase import assign_phases


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def _write_state_selection(seed_dir: str, snapshots: dict, compression_df: pd.DataFrame) -> list:
    trace_labels = classify_compression_trace(compression_df)
    labels_df = pd.DataFrame(trace_labels)
    labels_df.to_csv(os.path.join(seed_dir, "state_selection_report.csv"), index=False)

    candidates = {
        "jammed": labels_df[labels_df["label"] == "jammed_candidate"],
        "frustrated": labels_df[labels_df["label"] == "frustrated_candidate"],
    }
    candidates["jammed"].to_csv(os.path.join(seed_dir, "jammed_state_candidates.csv"), index=False)
    candidates["frustrated"].to_csv(os.path.join(seed_dir, "frustrated_state_candidates.csv"), index=False)

    for label in ("jammed", "frustrated", "random"):
        path = os.path.join(seed_dir, f"selected_state_{label}.h5")
        if os.path.exists(path):
            os.remove(path)

    selected_manifest = []
    for label in ("jammed", "frustrated"):
        if candidates[label].empty:
            continue
        row = candidates[label].iloc[-1]
        phi_target = float(row["phi_target"])
        selected_state = snapshots[phi_target].copy()
        filename = f"selected_state_{label}.h5"
        selected_state.save_h5(os.path.join(seed_dir, filename))
        selected_manifest.append({"label": label, "phi": row["phi"], "path": filename})

    flowing_candidates = labels_df[labels_df["label"] == "flowing"]
    if not flowing_candidates.empty:
        random_row = flowing_candidates.iloc[len(flowing_candidates) // 2]
        random_phi_target = float(random_row["phi_target"])
        random_phi = random_row["phi"]
    else:
        random_phi_target = max(snapshots)
        random_phi = snapshots[random_phi_target].phi
    random_state = snapshots[random_phi_target].copy()
    random_state.save_h5(os.path.join(seed_dir, "selected_state_random.h5"))
    selected_manifest.append({"label": "random", "phi": random_phi, "path": "selected_state_random.h5"})

    pd.DataFrame(selected_manifest).to_csv(os.path.join(seed_dir, "selected_states_manifest.csv"), index=False)
    return selected_manifest


def refresh_state_selection(out_dir: str) -> None:
    from .state import QaetherState

    energy_path = os.path.join(out_dir, "energy_curve.csv")
    compression_df = pd.read_csv(energy_path)
    for seed, seed_df in compression_df.groupby("seed", sort=True):
        seed_dir = os.path.join(out_dir, f"seed_{int(seed)}")
        seed_df = seed_df.sort_values("phi_target").reset_index(drop=True)
        snapshots = {
            float(row.phi_target): QaetherState.load_h5(
                os.path.join(seed_dir, f"snapshot_phi_{row.phi_target:.3f}.h5")
            )
            for row in seed_df.itertuples()
        }
        _write_state_selection(seed_dir, snapshots, seed_df)


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


def run_production_pipeline(
    out_dir: str,
    n: int = 64,
    seeds: List[int] = None,
    phi_targets: List[float] = None,
    smoke: bool = False,
    max_stage: int = 10,
) -> dict:
    if max_stage not in (7, 10):
        raise ValueError("max_stage must be 7 (geometry atlas) or 10 (dynamics pilots)")

    ensure_dir(out_dir)
    seeds = [0] if seeds is None else list(seeds)

    # Scale config options depending on smoke flag
    if smoke:
        cfg = ExperimentConfig(
            phi_targets=[0.20, 0.35, 0.50],
            relax_steps=20,
            dynamics_steps=30,
            snapshot_stride=10,
            dt=0.01,
        )
    else:
        cfg = ExperimentConfig(
            phi_targets=phi_targets if phi_targets is not None else [
                0.20, 0.30, 0.40, 0.50, 0.55, 0.58, 0.60,
                0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.735, 0.740
            ],
            relax_dt=0.2,
            relax_steps=3000,
            dynamics_steps=120,
            snapshot_stride=10,
            dt=0.01,
        )

    # Stage 0 outputs
    schema_dir = os.path.join(os.path.dirname(__file__), "schema")
    for schema_name in ("config_schema.json", "metadata_schema.json"):
        shutil.copyfile(
            os.path.join(schema_dir, schema_name),
            os.path.join(out_dir, schema_name),
        )

    with open(os.path.join(out_dir, "protocol_definitions.md"), "w", encoding="utf-8") as f:
        f.write(
            "# Protocol Definitions\n\n"
            "- Exclusion radius R_Q = 0.5\n"
            "- Unit length ell_Q = 1.0\n"
            "- Contact graph uses minimum-image distances under PBC\n"
            "- Geometry atlas runs with phase, perturbation, and dynamics disabled\n"
            "- Hard-sphere-like overlap thresholds: max_overlap <= 1e-3 and "
            "overlap_energy_per_vertex <= 1e-6\n"
        )

    if max_stage >= 9:
        with open(os.path.join(out_dir, "phase_model_spec.md"), "w", encoding="utf-8") as f:
            f.write("# Phase Model Spec\n\n- k_ij = k_core * (1 + lambda * cos(theta_i - theta_j))\n")

    all_compression_rows = []

    for seed in seeds:
        seed_dir = os.path.join(out_dir, f"seed_{seed}")
        ensure_dir(seed_dir)

        # Stage 1: Initial Random Gas
        state = random_gas(n=n, phi=0.10, seed=seed, radius=cfg.radius, ell_q=cfg.ell_q)
        state.save_h5(os.path.join(seed_dir, f"initial_state_seed_{seed}.h5"))

        # Stage 2: Compression Sweep
        snapshots, compression_df = compression_sweep(state, cfg)
        compression_df["seed"] = seed
        all_compression_rows.append(compression_df)

        # Save per-phi snapshots
        for phi, snap in snapshots.items():
            snap.save_h5(os.path.join(seed_dir, f"snapshot_phi_{phi:.3f}.h5"))

        # Stage 3-6: Atlas & Graph Analysis
        analyses, atlas_df = build_atlas(snapshots, cfg)
        atlas_df["seed"] = seed
        atlas_df.to_csv(os.path.join(seed_dir, "compression_graph_motif_atlas.csv"), index=False)

        # Save detailed graph / cycle / motif files per phi
        for phi, item in analyses.items():
            graph = item["graph"]
            triangles = item["triangles"]
            squares = item["squares"]
            t_motifs = item["t_motifs"]
            o_motifs = item["o_motifs"]

            pd.DataFrame(list(graph.edges), columns=["source", "target"]).to_csv(
                os.path.join(seed_dir, f"graph_edges_phi_{phi:.3f}.csv"), index=False
            )
            pd.DataFrame(triangles, columns=["v0", "v1", "v2"]).to_csv(
                os.path.join(seed_dir, f"primitive_triangles_phi_{phi:.3f}.csv"), index=False
            )
            pd.DataFrame(squares, columns=["v0", "v1", "v2", "v3"]).to_csv(
                os.path.join(seed_dir, f"primitive_squares_phi_{phi:.3f}.csv"), index=False
            )
            pd.DataFrame(t_motifs, columns=["v0", "v1", "v2", "v3"]).to_csv(
                os.path.join(seed_dir, f"T_motifs_phi_{phi:.3f}.csv"), index=False
            )
            pd.DataFrame(o_motifs, columns=["v0", "v1", "v2", "v3", "v4", "v5"]).to_csv(
                os.path.join(seed_dir, f"O_motifs_phi_{phi:.3f}.csv"), index=False
            )

        # State Selection (Stage 7)
        selected_manifest = _write_state_selection(seed_dir, snapshots, compression_df)
        highest_phi = max(snapshots.keys())
        highest_state = snapshots[highest_phi].copy()

        if max_stage >= 8:
            # Stage 8: Pressure-off Dynamics for all representative states
            from .state import QaetherState
            for selected_item in selected_manifest:
                label = selected_item["label"]
                state_path = os.path.join(seed_dir, selected_item["path"])
                sel_state = QaetherState.load_h5(state_path)

                dyn_df = pressure_off_dynamics(
                    sel_state.copy(), cfg,
                    traj_path=os.path.join(seed_dir, f"pressure_off_trajectory_{label}.h5")
                )
                dyn_df.to_csv(os.path.join(seed_dir, f"pressure_off_summary_{label}.csv"), index=False)

                # Compatibility fallback
                if label == "random":
                    dyn_df.to_csv(os.path.join(seed_dir, "pressure_off_summary.csv"), index=False)

        if max_stage >= 9:
            # Stage 9: Phase Coupling Sweep
            phase_rows = []
            for lp in [0.0, 0.1, 0.3, 0.5]:
                phase_state = assign_phases(highest_state.copy(), seed=seed + 100, omega_q=cfg.omega_q)
                p_df = pressure_off_dynamics(phase_state, cfg, lambda_phase=lp)
                p_df["lambda_phase"] = lp
                phase_rows.append(p_df)
            pd.concat(phase_rows).to_csv(os.path.join(seed_dir, "phase_sweep_by_lambda.csv"), index=False)

        if max_stage >= 10:
            # Stage 10: Perturbation Sweep
            pert_rows = []
            pert_meta = []
            for amp in [1.0e-4, 1.0e-3, 1.0e-2]:
                pert_state = random_micro_jitter(highest_state.copy(), amplitude=amp, seed=seed + 200)
                pert_state.save_h5(os.path.join(seed_dir, f"perturbed_initial_state_amp_{amp:.1e}.h5"))

                p_df = pressure_off_dynamics(pert_state, cfg)
                p_df["amplitude"] = amp
                pert_rows.append(p_df)

                pert_meta.append({"amplitude": amp, "seed": seed + 200, "state_label": "random_micro_jitter"})

            pd.concat(pert_rows).to_csv(os.path.join(seed_dir, "perturbed_pressure_off_summary.csv"), index=False)
            pd.DataFrame(pert_meta).to_csv(os.path.join(seed_dir, "perturbation_metadata.csv"), index=False)

    # FCC / HCP lattice benchmarks at out_dir level
    lattice_rows = []
    ideal_lattice_phi = close_packing_phi() * 0.99

    # Ideal Reference Lattice
    for name, l_state in build_lattice_benchmarks(n, ideal_lattice_phi, cfg.radius, cfg.ell_q).items():
        l_state.metadata.update({"benchmark_type": "ideal_close_packing_reference"})
        item = analyze_snapshot(l_state.phi, l_state, cfg)
        row = dict(item["summary"])
        row["kind"] = name
        row["benchmark_type"] = "ideal_close_packing_reference"
        row["hard_sphere_valid"] = (
            row["max_overlap"] <= 1.0e-3
            and row["overlap_energy_per_vertex"] <= 1.0e-6
        )
        row["toto_valid"] = False
        lattice_rows.append(row)
        l_state.save_h5(os.path.join(out_dir, f"{name}_state_phi_{l_state.phi:.3f}.h5"))

    # TOTO-labeled Constructed Lattice Reference
    for name, l_state in build_toto_labeled_benchmarks(n, cfg.radius, cfg.ell_q).items():
        item = analyze_snapshot(l_state.phi, l_state, cfg)
        row = dict(item["summary"])
        hard_sphere_valid = bool(l_state.metadata["hard_sphere_valid"])
        toto_valid = hard_sphere_valid and bool(item["t_motifs"]) and bool(item["o_motifs"])
        l_state.metadata["toto_valid"] = toto_valid
        row["kind"] = name
        row["benchmark_type"] = "toto_labeled_constructed"
        row["hard_sphere_valid"] = hard_sphere_valid
        row["toto_valid"] = toto_valid
        lattice_rows.append(row)
        l_state.save_h5(os.path.join(out_dir, f"{name}_constructed_state.h5"))

    pd.DataFrame(lattice_rows).to_csv(os.path.join(out_dir, "hcp_fcc_motif_summary.csv"), index=False)

    # Write graph_motif_atlas_report.md
    with open(os.path.join(out_dir, "graph_motif_atlas_report.md"), "w", encoding="utf-8") as f:
        f.write(
            "# Compression Graph Motif Atlas Report\n\n"
            f"- Completed run for Stages 0-{max_stage}.\n"
            f"- N = {n}\n"
            f"- Seeds = {seeds}\n"
            f"- Phi targets = {cfg.phi_targets}\n"
        )

    # Master summary files
    if all_compression_rows:
        pd.concat(all_compression_rows).to_csv(os.path.join(out_dir, "energy_curve.csv"), index=False)

    return {"status": "success", "seeds": seeds, "n": n, "max_stage": max_stage}
