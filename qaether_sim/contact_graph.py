from __future__ import annotations

import networkx as nx

from .geometry import distance


def cycle_node_set(graph: nx.Graph) -> set:
    nodes = set()
    for component in nx.biconnected_components(graph):
        if len(component) >= 3:
            nodes.update(component)
    return nodes


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
    component_nodes = list(nx.connected_components(graph))
    bridge_edges = list(nx.bridges(graph))
    cycle_nodes = cycle_node_set(graph)

    chain_components = 0
    tree_components = 0
    for nodes in component_nodes:
        subgraph = graph.subgraph(nodes)
        if len(nodes) >= 2 and subgraph.number_of_edges() == len(nodes) - 1:
            tree_components += 1
            if max(dict(subgraph.degree()).values(), default=0) <= 2:
                chain_components += 1

    return {
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "mean_degree": sum(degrees) / max(len(degrees), 1),
        "max_degree": max(degrees) if degrees else 0,
        "components": len(components),
        "largest_component": max(components) if components else 0,
        "clustering": nx.average_clustering(graph) if graph.number_of_nodes() else 0.0,
        "isolated_nodes": sum(degree == 0 for degree in degrees),
        "degree1_nodes": sum(degree == 1 for degree in degrees),
        "degree2_nodes": sum(degree == 2 for degree in degrees),
        "low_degree_0_2_nodes": sum(degree <= 2 for degree in degrees),
        "bridge_edges": len(bridge_edges),
        "cycle_nodes": len(cycle_nodes),
        "noncycle_nodes": graph.number_of_nodes() - len(cycle_nodes),
        "chain_components": chain_components,
        "tree_components": tree_components,
    }
