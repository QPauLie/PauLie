"""
Module for computing the average Out-of-Time-Order Correlator (OTOC).
"""

import networkx as nx
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection


def average_otoc(
    op_v: PauliString, op_w: PauliString, system_generators: PauliStringCollection
) -> float:
    """
    Computes the average Out-of-Time-Order Correlator (OTOC) between two
    Pauli strings, op_v and op_w, for a system defined by its generators.

    The calculation is based on Corollary 3 of arXiv:2502.16404, which relates
    the average OTOC to a counting problem on the commutator graph.

    Args:
        op_v (PauliString): The first Pauli string operator.
        op_w (PauliString): The second Pauli string operator.
        system_generators (PauliStringCollection): The generators of the system.

    Returns:
        The average OTOC value, a float between -1 and 1.
    """
    num_qubits = system_generators.get_size()
    if len(op_v) != num_qubits or len(op_w) != num_qubits:
        raise ValueError(
            "Operators V, W, and system generators must act on the same number of qubits."
        )

    # Step 1: Find the connected components of both V and W.
    nodes, edges = system_generators.get_commutator_graph()
    comm_graph = nx.Graph(edges)
    comm_graph.add_nodes_from(nodes)

    component_of_v = None
    for component in nx.connected_components(comm_graph):
        if str(op_v) in component:
            component_of_v = component
            break

    # Step 2: Apply the correct formula based on Corollary 3
    # First, check if W is in the same component as V
    if str(op_w) in component_of_v:
        # Case 1: C(V) = C(W). Average OTOC is 1.0.
        return 1.0

    # Case 2: C(V) != C(W). Use the counting formula.
    component_size = len(component_of_v)

    anticommuting_count = 0
    for t_str in component_of_v:
        op_t = PauliString(pauli_str=t_str)
        if not op_w.commutes_with(op_t):
            anticommuting_count += 1

    if component_size == 0:
        return 1.0

    return 1.0 - 2.0 * (anticommuting_count / component_size)
