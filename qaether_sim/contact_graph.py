from __future__ import annotations

import networkx as nx

from .geometry import distance


def build_contact_graph(state, epsilon_contact: float = 0.06) -> nx.Graph:
    cutoff = state.ell_q * (1.0 + epsilon_contact)
    graph = nx.Graph()
    graph.add_nodes_from(range(state.n))
    for i in range(state.n - 1):
        for j in range(i + 1, state.n):
            d = distance(state.positions, i, j, state.box)
            if d <= cutoff:
                graph.add_edge(i, j, distance=d)
    graph.graph["phi"] = state.phi
    return graph


def build_hysteretic_contact_graph(
    state,
    previous_graph: nx.Graph | None = None,
    epsilon_on: float = 0.01,
    epsilon_off: float = 0.02,
) -> nx.Graph:
    on_cutoff = state.ell_q * (1.0 + epsilon_on)
    off_cutoff = state.ell_q * (1.0 + epsilon_off)
    graph = nx.Graph()
    graph.add_nodes_from(range(state.n))
    previous_edges = set()
    if previous_graph is not None:
        previous_edges = {tuple(sorted(edge)) for edge in previous_graph.edges}
    for i in range(state.n - 1):
        for j in range(i + 1, state.n):
            d = distance(state.positions, i, j, state.box)
            edge = (i, j)
            if edge in previous_edges:
                keep = d < off_cutoff
            else:
                keep = d <= on_cutoff
            if keep:
                graph.add_edge(i, j, distance=d)
    graph.graph["phi"] = state.phi
    graph.graph["contact_rule"] = "hysteretic"
    graph.graph["epsilon_on"] = epsilon_on
    graph.graph["epsilon_off"] = epsilon_off
    return graph


def graph_summary(graph: nx.Graph) -> dict:
    degrees = [d for _, d in graph.degree()]
    components = [len(c) for c in nx.connected_components(graph)]
    return {
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "mean_degree": sum(degrees) / max(len(degrees), 1),
        "max_degree": max(degrees) if degrees else 0,
        "components": len(components),
        "largest_component": max(components) if components else 0,
        "clustering": nx.average_clustering(graph) if graph.number_of_nodes() else 0.0,
    }
