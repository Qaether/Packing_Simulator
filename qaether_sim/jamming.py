from __future__ import annotations

import networkx as nx
import numpy as np


def compute_z_rattler_removed(graph: nx.Graph) -> float:
    """Recursively remove nodes with degree <= 3 (rattlers) and compute Z on the core."""
    core = graph.copy()
    while True:
        rattlers = [node for node, deg in core.degree() if deg <= 3]
        if not rattlers:
            break
        core.remove_nodes_from(rattlers)

    n_core = core.number_of_nodes()
    if n_core == 0:
        return 0.0
    return float(2.0 * core.number_of_edges() / n_core)


def classify_compression_trace(summary_df, pressure_quantile: float = 0.8):
    """
    Classify each state in the compression trace as flowing, jammed_candidate, or frustrated_candidate.
    Uses Z_rattler_removed, force_balance_residual, and dP/dphi if available.
    """
    rows = []
    if summary_df.empty:
        return rows

    p_cut = summary_df["pressure"].quantile(pressure_quantile)

    # Calculate dP/dphi (bulk modulus proxy) in trace
    phi_vals = summary_df["phi"].values
    p_vals = summary_df["pressure"].values
    dp_dphi = np.zeros_like(phi_vals)
    for i in range(1, len(phi_vals)):
        dphi = phi_vals[i] - phi_vals[i-1]
        if dphi > 1e-8:
            dp_dphi[i] = (p_vals[i] - p_vals[i-1]) / dphi

    msd_cut = max(float(summary_df["msd_step"].median()), 1.0e-12)
    for position, (_, row) in enumerate(summary_df.iterrows()):
        label = "flowing"

        # Check coordinates and forces
        z_rattler = row.get("z_rattler_removed", 0.0)
        force_res = row.get("force_balance_residual", 0.0)
        dp_dphi_val = dp_dphi[position]

        # Geometric near-contacts and numerical-scale overlaps can raise Z or
        # pressure without carrying meaningful load. Use the same overlap gate
        # as the hard-sphere/soft-overcompressed state classification.
        if "max_overlap" in row and "overlap_energy_per_vertex" in row:
            has_active_contacts = (
                row["max_overlap"] > 1.0e-3
                or row["overlap_energy_per_vertex"] > 1.0e-6
            )
        else:
            has_active_contacts = row["pressure"] > 1.0e-6 or row["energy"] > 1.0e-8

        if has_active_contacts:
            low_mobility = row["msd_step"] <= msd_cut
            constrained = z_rattler >= 4.0
            load_bearing = row["pressure"] >= max(float(p_cut), 1.0e-6) or dp_dphi_val > 10.0
            force_balanced = force_res < 1.0e-4

            if constrained and load_bearing and force_balanced and low_mobility:
                label = "jammed_candidate"
            elif low_mobility and (constrained or load_bearing):
                label = "frustrated_candidate"

        rows.append({
            "phi": row["phi"],
            "phi_target": row.get("phi_target", row["phi"]),
            "label": label,
            "z_rattler_removed": z_rattler,
            "force_balance_residual": force_res,
            "dp_dphi": dp_dphi_val
        })
    return rows
