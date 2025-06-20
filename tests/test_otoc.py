"""
    Tests for the average out-of-time-order correlator.
"""
import pytest
import numpy as np
import networkx as nx
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.application.otoc import average_otoc
from paulie.common.pauli_string_factory import (
    get_pauli_string as p,
    get_identity
)

generators_list = [
    ["I"], ["X"], ["Y"], ["Z"],
    ["XX", "YY", "ZZ"],
    ["ZI", "IZ", "XX"],
    ["XI", "IX", "XX", "YY"],
    ["XY", "YX", "YZ", "ZY"],
    ["XI", "IX", "YI", "IY", "ZZ"],
    ["XI", "ZZ", "YI", "IY", "XY", "YX"],
    ["XI", "IX", "YI", "IY", "ZI", "IZ", "XX"],
]

def naive_otoc(generators: PauliStringCollection,
                 v: PauliString, w: PauliString) -> float:
    """
    Computes the Haar averaged out-of-time-order correlator between two Pauli strings
    using NetworkX.

    We can compute this as
    1 - 2 * |{W, P} = 0 : P in connected component of V| / |connected component of V|
    where we take the commutator graph. (arXiV:2502.16404)

    Args:
        generators: Generating set of the Pauli string DLA.
        v: Pauli string V
        w: Pauli string W
    """
    # Generate commutator graph
    vertices, edges = generators.get_commutator_graph()
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    # Get connected component of V
    v_connected_component = nx.node_connected_component(graph, str(v))
    # Count the number of elements t in the connected component of V
    # that anticommute with W
    anticommute_count = sum(not w | p(t) for t in v_connected_component)
    return 1 - 2 * anticommute_count / len(v_connected_component)

@pytest.mark.parametrize("generators", generators_list)
def test_average_otoc_matches_naive(generators: list[str]) -> None:
    """
    Test that average_otoc(g, v, w) == naive_otoc(g, v, w).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    for v in all_paulis:
        for w in all_paulis:
            assert average_otoc(g, v, w) == pytest.approx(naive_otoc(g, v, w))

@pytest.mark.parametrize("generators", generators_list)
def test_average_otoc_is_symmetric(generators: list[str]) -> None:
    """
    Test that average_otoc(g, v, w) == average_otoc(g, w, v).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    for v in all_paulis:
        for w in all_paulis:
            assert average_otoc(g, v, w) == pytest.approx(average_otoc(g, w, v))

@pytest.mark.parametrize("generators", generators_list)
def test_average_eq_initial_otoc_for_commutants(generators: list[str]) -> None:
    """
    Test that the average_otoc(g, v, w) == tr[W @ V @ W @ V] / d if V is a
    commutant of the DLA (since it would commute with the evolution unitary by
    the Baker-Campbell-Haussdorf formula).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    commutants = g.get_commutants()
    d = 2 ** len(generators[0])
    for v in commutants:
        vmat = v.get_matrix()
        for w in all_paulis:
            wmat = w.get_matrix()
            analytical_value = np.trace(wmat @ vmat @ wmat @ vmat) / d
            assert average_otoc(g, v, w) == pytest.approx(analytical_value)
