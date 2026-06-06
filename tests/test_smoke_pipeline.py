import os
import tempfile
import unittest

import numpy as np

from qaether_sim.analysis import run_smoke_pipeline
from qaether_sim.config import ExperimentConfig
from qaether_sim.contact_graph import build_contact_graph, build_hysteretic_contact_graph
from qaether_sim.cycles import primitive_squares, primitive_triangles
from qaether_sim.initial_conditions import fcc_lattice, random_gas
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


if __name__ == "__main__":
    unittest.main()
