#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os

import pandas as pd


METRICS = [
    "isolated_nodes",
    "degree1_nodes",
    "degree2_nodes",
    "low_degree_0_2_nodes",
    "bridge_edges",
    "cycle_nodes",
    "noncycle_nodes",
    "t_or_o_motif_nodes",
    "non_t_o_cycle_nodes",
]


def main():
    parser = argparse.ArgumentParser(description="Compare topology fractions across Phase A system sizes.")
    parser.add_argument("result_dirs", nargs="+", help="Phase A result directories")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-report", required=True)
    args = parser.parse_args()

    frames = []
    for result_dir in args.result_dirs:
        energy = pd.read_csv(os.path.join(result_dir, "energy_curve.csv"))
        n = int(energy["n"].iloc[0])
        topology = pd.read_csv(os.path.join(result_dir, "structure_topology_by_phi.csv"))
        for metric in METRICS:
            topology[f"{metric}_fraction"] = topology[metric] / n
        grouped = topology.groupby("phi")[[f"{metric}_fraction" for metric in METRICS]].mean().reset_index()
        grouped.insert(0, "N", n)
        frames.append(grouped)

    comparison = pd.concat(frames, ignore_index=True).sort_values(["N", "phi"])
    comparison.to_csv(args.out_csv, index=False)

    selected_phi = [0.20, 0.30, 0.40, 0.50, 0.55, 0.58, 0.64, 0.70, 0.74]
    table = comparison[comparison["phi"].round(3).isin(selected_phi)].copy()
    columns = [
        "N",
        "phi",
        "isolated_nodes_fraction",
        "low_degree_0_2_nodes_fraction",
        "bridge_edges_fraction",
        "cycle_nodes_fraction",
        "t_or_o_motif_nodes_fraction",
        "non_t_o_cycle_nodes_fraction",
    ]
    with open(args.out_report, "w", encoding="utf-8") as report:
        report.write("# Structure Topology Scaling: N=64 and N=128\n\n")
        report.write(
            "Fractions are seed means. `non_t_o_cycle_nodes` is a measurable proxy for "
            "cycle-rich structure not assigned to a validated T/O motif; it is not a direct "
            "incomplete-cage detector.\n\n"
        )
        report.write("| " + " | ".join(columns) + " |\n")
        report.write("| " + " | ".join(["---"] * len(columns)) + " |\n")
        for values in table[columns].itertuples(index=False, name=None):
            formatted = [
                str(int(value)) if column == "N" else f"{float(value):.4f}"
                for column, value in zip(columns, values)
            ]
            report.write("| " + " | ".join(formatted) + " |\n")


if __name__ == "__main__":
    main()
