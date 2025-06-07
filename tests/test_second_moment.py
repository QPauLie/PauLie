"""
    Test second-order twirl routine.
"""
import pytest
import numpy as np
import networkx as nx
from paulie.common.pauli_string_factory import get_pauli_string as p, get_identity
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.application.second_moment import second_moment

def create_swap_matrix(n: int) -> np.ndarray:
    """
    Construct SWAP matrix on n qubits.

    The SWAP matrix is given by sum_{i} sum_{j} |i, j> <j, i|.
    """
    swap = np.zeros((4 ** n, 4 ** n))
    for i in range(4 ** n):
        # Take the binary repr of i, split it in half and swap two halves
        t = np.binary_repr(i, width=2*n)
        j = int(t[n:] + t[:n], base=2)
        swap[i, j] = 1
    return swap

def average_otoc(generators: PauliStringCollection, v: PauliString, w: PauliString) -> float:
    """
    Compute the OTOC using the commutator graph.
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
    anticommute_count = sum([not w | p(t) for t in v_connected_component])
    return 1 - 2 * anticommute_count / len(v_connected_component)

def otoc_from_twirl(generators: PauliStringCollection, v: PauliString, w: PauliString) -> float:
    """
    Compute the OTOC using Corollary 3 in Appendix A of arXiV:2502.16404.
    """
    vkronv = p([(1, str(v) + str(v))])
    wkronw = p([(1, str(w) + str(w))])
    wkronw_twirled = second_moment(generators, wkronw)
    VVop = vkronv.get_matrix()
    WWtop = wkronw_twirled.get_matrix()
    swap = create_swap_matrix(len(v))
    return np.trace(VVop @ WWtop @ swap) / 2 ** len(v)

def incorrect_otoc_from_twirl(generators: PauliStringCollection, v: PauliString, w: PauliString) -> float:
    """
    Compute the OTOC using an incorrect SWAP matrix.
    """
    vkronv = p([(1, str(v) + str(v))])
    wkronw = p([(1, str(w) + str(w))])
    wkronw_twirled = second_moment(generators, wkronw)
    VVop = vkronv.get_matrix()
    WWtop = wkronw_twirled.get_matrix()
    swap = np.eye(4 ** len(v))
    return np.trace(VVop @ WWtop @ swap) / 2 ** len(v)

generator_list = [
    ["I"], ["X"], ["Y"], ["Z"],
    ["ZI", "IZ", "XX"],
    ["XI", "IX", "XY"],
    ["XX", "YY", "ZZ"],
]

@pytest.mark.parametrize("generators", generator_list)
def test_otoc_correct(generators: list[str]):
    """
    Test that the OTOC computed analytically and from the twirl match.
    """
    gens = p(generators)
    i = get_identity(len(generators[0]))
    analytical_values = []
    computed_values = []
    for v in i.gen_all_pauli_strings():
        for w in i.gen_all_pauli_strings():
            analytical_values.append(average_otoc(gens, v, w))
            computed_values.append(otoc_from_twirl(gens, v, w))
    assert all(a_val == pytest.approx(c_val) for a_val, c_val in zip(analytical_values, computed_values))

@pytest.mark.parametrize("generators", generator_list)
def test_otoc_incorrect(generators: list[str]):
    """
    Test that the OTOC computed analytically and from the twirl using
    incorrect SWAP matrix do not match.
    """
    gens = p(generators)
    i = get_identity(len(generators[0]))
    analytical_values = []
    computed_values = []
    for v in i.gen_all_pauli_strings():
        for w in i.gen_all_pauli_strings():
            analytical_values.append(average_otoc(gens, v, w))
            computed_values.append(incorrect_otoc_from_twirl(gens, v, w))
    assert any(a_val != pytest.approx(c_val) for a_val, c_val in zip(analytical_values, computed_values))
