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


def build_toto_labeled_benchmarks(n: int, radius: float = 0.5, ell_q: float = 1.0):
    """
    Construct HCP/FCC structures at the exact close packing scale, ensuring
    nearest-neighbor distances are exactly ell_q, guaranteeing clean TOTO incidence conditions.
    """
    from .geometry import overlap_metrics

    ideal_phi = close_packing_phi()
    hcp = hcp_lattice(n, ideal_phi, radius, ell_q)
    fcc = fcc_lattice(n, ideal_phi, radius, ell_q)

    hcp_metrics = overlap_metrics(hcp)
    fcc_metrics = overlap_metrics(fcc)

    # Strictly validate that overlap conforms to our hard-sphere-like tolerances
    hcp_valid = (
        hcp_metrics["max_overlap"] <= 1.0e-3
        and hcp_metrics["overlap_energy_per_vertex"] <= 1.0e-6
    )
    fcc_valid = (
        fcc_metrics["max_overlap"] <= 1.0e-3
        and fcc_metrics["overlap_energy_per_vertex"] <= 1.0e-6
    )

    hcp.metadata.update({
        "benchmark_type": "toto_labeled_constructed",
        "phi_target": ideal_phi,
        "phi_achieved": hcp.phi,
        "hard_sphere_valid": bool(hcp_valid),
        "toto_valid": None,
        **hcp_metrics
    })
    fcc.metadata.update({
        "benchmark_type": "toto_labeled_constructed",
        "phi_target": ideal_phi,
        "phi_achieved": fcc.phi,
        "hard_sphere_valid": bool(fcc_valid),
        "toto_valid": None,
        **fcc_metrics
    })
    return {
        "hcp_toto": hcp,
        "fcc_toto": fcc
    }
