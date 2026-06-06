from __future__ import annotations

import numpy as np

from .initial_conditions import fcc_lattice, hcp_lattice


def close_packing_phi() -> float:
    return float(np.pi / (3.0 * np.sqrt(2.0)))


def build_lattice_benchmarks(n: int, phi: float = None, radius: float = 0.5, ell_q: float = 1.0):
    reference_phi = close_packing_phi() * 0.99 if phi is None else phi
    return {
        "hcp": hcp_lattice(n, reference_phi, radius, ell_q),
        "fcc": fcc_lattice(n, reference_phi, radius, ell_q),
    }
