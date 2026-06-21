import os
import tempfile
import unittest

import numpy as np

from qaether_sim.analysis import run_production_pipeline, run_smoke_pipeline, run_stage8_pressure_off
from qaether_sim.bulk_dynamics import pressure_off_dynamics
from qaether_sim.compression import relax
from qaether_sim.config import ExperimentConfig
from qaether_sim.contact_graph import build_contact_graph, build_hysteretic_contact_graph, graph_summary
from qaether_sim.cycles import primitive_squares, primitive_triangles
from qaether_sim.geometry import overlap_metrics
from qaether_sim.initial_conditions import fcc_lattice, hcp_lattice, random_gas
from qaether_sim.jamming import classify_compression_trace
from qaether_sim.lattice_benchmarks import close_packing_phi
from qaether_sim.motifs_O import detect_o_motifs
from qaether_sim.motifs_T import detect_t_motifs
from qaether_sim.state import QaetherState


class QaetherSmokeTests(unittest.TestCase):
    def test_hdf5_roundtrip(self):
        state = random_gas(8, 0.1, seed=1)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "state.h5")
            state.save_h5(path)
            loaded = QaetherState.load_h5(path)
        self.assertEqual(loaded.n, state.n)
        self.assertTrue(np.allclose(loaded.positions, state.positions))
        self.assertTrue(np.allclose(loaded.box, state.box))

    def test_cycle_and_t_motif_detection_on_fcc(self):
        cfg = ExperimentConfig(epsilon_contact=0.08, epsilon_cycle=0.20)
        state = fcc_lattice(32, 0.60)
        graph = build_contact_graph(state, cfg.epsilon_contact)
        triangles = primitive_triangles(graph, state, cfg.epsilon_cycle)
        squares = primitive_squares(graph, state, cfg.epsilon_cycle, cfg.epsilon_planar)
        t_motifs = detect_t_motifs(graph, state, triangles, cfg.epsilon_cycle, cfg.epsilon_volume)
        self.assertGreater(graph.number_of_edges(), 0)
        self.assertIsInstance(triangles, list)
        self.assertIsInstance(squares, list)
        self.assertIsInstance(t_motifs, list)

    def test_ideal_tetrahedron_and_octahedron_motifs(self):
        cfg = ExperimentConfig(epsilon_contact=0.06, epsilon_cycle=0.08, epsilon_center=0.08, epsilon_perp=0.08)
        tetra_points = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, np.sqrt(3.0) / 2.0, 0.0],
                [0.5, np.sqrt(3.0) / 6.0, np.sqrt(2.0 / 3.0)],
            ]
        ) + 2.0
        tetra = QaetherState(tetra_points, np.zeros_like(tetra_points), np.array([6.0, 6.0, 6.0]))
        graph = build_contact_graph(tetra, cfg.epsilon_contact)
        triangles = primitive_triangles(graph, tetra, cfg.epsilon_cycle)
        self.assertEqual(len(detect_t_motifs(graph, tetra, triangles, cfg.epsilon_cycle, cfg.epsilon_volume)), 1)

        a = 1.0 / np.sqrt(2.0)
        octa_points = np.array(
            [
                [a, 0.0, 0.0],
                [-a, 0.0, 0.0],
                [0.0, a, 0.0],
                [0.0, -a, 0.0],
                [0.0, 0.0, a],
                [0.0, 0.0, -a],
            ]
        ) + 2.0
        octa = QaetherState(octa_points, np.zeros_like(octa_points), np.array([6.0, 6.0, 6.0]))
        graph = build_contact_graph(octa, cfg.epsilon_contact)
        triangles = primitive_triangles(graph, octa, cfg.epsilon_cycle)
        squares = primitive_squares(graph, octa, cfg.epsilon_cycle, cfg.epsilon_planar)
        self.assertEqual(len(detect_o_motifs(graph, octa, triangles, squares, cfg.epsilon_center, cfg.epsilon_perp)), 1)

    def test_hysteretic_contact_rule(self):
        state = QaetherState(
            positions=np.array([[0.0, 0.0, 0.0], [1.015, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            box=np.array([4.0, 4.0, 4.0]),
        )
        graph = build_hysteretic_contact_graph(state, epsilon_on=0.01, epsilon_off=0.02)
        self.assertEqual(graph.number_of_edges(), 0)

        previous = build_contact_graph(state, epsilon_contact=0.02)
        graph = build_hysteretic_contact_graph(state, previous, epsilon_on=0.01, epsilon_off=0.02)
        self.assertEqual(graph.number_of_edges(), 1)

    def test_graph_summary_reports_open_and_cycle_structure(self):
        import networkx as nx

        graph = nx.Graph()
        graph.add_nodes_from(range(8))
        graph.add_edges_from([(0, 1), (1, 2), (3, 4), (4, 5), (5, 3), (5, 6)])
        summary = graph_summary(graph)
        self.assertEqual(summary["isolated_nodes"], 1)
        self.assertEqual(summary["degree1_nodes"], 3)
        self.assertEqual(summary["degree2_nodes"], 3)
        self.assertEqual(summary["bridge_edges"], 3)
        self.assertEqual(summary["cycle_nodes"], 3)
        self.assertEqual(summary["noncycle_nodes"], 5)
        self.assertEqual(summary["chain_components"], 1)
        self.assertEqual(summary["tree_components"], 1)

    def test_force_free_state_is_flowing(self):
        import pandas as pd

        trace = pd.DataFrame(
            [
                {
                    "phi": 0.10,
                    "phi_target": 0.10,
                    "pressure": 0.0,
                    "energy": 0.0,
                    "msd_step": 0.0,
                    "force_balance_residual": 0.0,
                    "z_rattler_removed": 0.0,
                    "max_overlap": 0.0,
                    "overlap_energy_per_vertex": 0.0,
                }
            ]
        )
        labels = classify_compression_trace(trace)
        self.assertEqual(labels[0]["label"], "flowing")

        trace.loc[0, "z_rattler_removed"] = 6.0
        labels = classify_compression_trace(trace)
        self.assertEqual(labels[0]["label"], "flowing")

        trace.loc[0, "pressure"] = 2.0e-6
        trace.loc[0, "energy"] = 2.0e-8
        labels = classify_compression_trace(trace)
        self.assertEqual(labels[0]["label"], "flowing")

    def test_periodic_lattices_are_complete_and_non_overlapping(self):
        phi = close_packing_phi()
        for n in (24, 32, 64, 128):
            for builder in (fcc_lattice, hcp_lattice):
                state = builder(n, phi)
                self.assertEqual(state.n, n)
                self.assertAlmostEqual(state.phi, phi)
                self.assertLessEqual(overlap_metrics(state)["max_overlap"], 1.0e-12)

        with self.assertRaises(ValueError):
            fcc_lattice(30, phi)

    def test_relaxation_uses_minimum_image_msd_and_post_update_energy(self):
        state = QaetherState(
            positions=np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            box=np.array([5.0, 5.0, 5.0]),
        )
        trace = relax(state, ExperimentConfig(relax_steps=1, dt=0.1))
        self.assertAlmostEqual(trace.iloc[0]["msd_step"], 0.0025)
        self.assertAlmostEqual(trace.iloc[0]["energy"], 0.08)

    def test_relaxation_converges_early_with_separate_timestep(self):
        state = random_gas(24, 0.10, seed=3)
        cfg = ExperimentConfig(
            relax_dt=0.2,
            relax_steps=3000,
            relax_min_steps=20,
            relax_convergence_window=5,
        )
        trace = relax(state, cfg)
        self.assertLess(len(trace), cfg.relax_steps)
        self.assertLessEqual(overlap_metrics(state)["max_overlap"], 1.0e-3)

    def test_pressure_off_uses_one_hysteretic_graph_and_nan_survival(self):
        cfg = ExperimentConfig(
            dynamics_steps=1,
            snapshot_stride=1,
            epsilon_contact=0.06,
            epsilon_contact_on=0.01,
            epsilon_contact_off=0.02,
        )
        state = QaetherState(
            positions=np.array([[0.0, 0.0, 0.0], [1.04, 0.0, 0.0]]),
            velocities=np.zeros((2, 3)),
            box=np.array([5.0, 5.0, 5.0]),
        )
        result = pressure_off_dynamics(state, cfg)
        self.assertTrue((result["edges"] == 0).all())
        self.assertTrue((result["edge_edit_distance"] == 0).all())
        self.assertTrue(result["S_E"].isna().all())
        self.assertTrue(result["S_T"].isna().all())
        self.assertTrue(result["S_O"].isna().all())
        self.assertTrue((result["dynamics_type"] == "none").all())

        phase_state = state.copy()
        phase_state.theta = np.array([0.0, 1.0])
        phase_result = pressure_off_dynamics(phase_state, cfg, lambda_phase=0.3)
        self.assertTrue((phase_result["dynamics_type"] == "none").all())

    def test_full_smoke_pipeline_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = run_smoke_pipeline(tmp, n=24, seed=2)
            expected = [
                "energy_curve.csv",
                "compression_graph_motif_atlas.csv",
                "pressure_off_summary.csv",
                "phase_sweep_summary.csv",
                "perturbed_pressure_off_summary.csv",
                "hcp_fcc_motif_summary.csv",
                "smoke_summary.json",
            ]
            for name in expected:
                self.assertTrue(os.path.exists(os.path.join(tmp, name)), name)
        self.assertEqual(summary["atlas_rows"], 3)
        self.assertGreater(summary["pressure_off_rows"], 0)

    def test_stage7_pipeline_stops_before_dynamics(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = run_production_pipeline(tmp, n=24, seeds=[0], smoke=True, max_stage=7)
            seed_dir = os.path.join(tmp, "seed_0")
            self.assertTrue(os.path.exists(os.path.join(tmp, "config_schema.json")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "metadata_schema.json")))
            self.assertTrue(os.path.exists(os.path.join(seed_dir, "compression_graph_motif_atlas.csv")))
            self.assertTrue(os.path.exists(os.path.join(seed_dir, "selected_states_manifest.csv")))
            self.assertFalse(os.path.exists(os.path.join(seed_dir, "pressure_off_summary.csv")))
            self.assertFalse(os.path.exists(os.path.join(seed_dir, "phase_sweep_by_lambda.csv")))
            self.assertFalse(os.path.exists(os.path.join(seed_dir, "perturbed_pressure_off_summary.csv")))
        self.assertEqual(summary["max_stage"], 7)

    def test_stage8_reuses_phase_a_states(self):
        with tempfile.TemporaryDirectory() as tmp:
            phase_a = os.path.join(tmp, "phase_a")
            stage8 = os.path.join(tmp, "stage8")
            run_production_pipeline(phase_a, n=24, seeds=[0], smoke=True, max_stage=7)
            summary = run_stage8_pressure_off(phase_a, stage8, steps=1, snapshot_stride=1)
            final = os.path.join(stage8, "pressure_off_final_summary.csv")
            manifest = os.path.join(stage8, "stage8_run_manifest.csv")
            self.assertTrue(os.path.exists(final))
            self.assertTrue(os.path.exists(manifest))
            self.assertEqual(summary["runs"], 6)
            self.assertEqual(len(__import__("pandas").read_csv(final)), 6)


if __name__ == "__main__":
    unittest.main()
