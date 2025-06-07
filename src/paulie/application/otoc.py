"""
Module for computing the average Out-of-Time-Order Correlator (OTOC).
"""
import networkx as nx

from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection


def average_otoc(
    static_op: PauliString,
    initial_evolved_op: PauliString,
    system_generators: PauliStringCollection,
) -> float:
    """
    Computes the average Out-of-Time-Order Correlator (OTOC) between a
    static operator and a time-evolved operator for a system ensemble.

    The calculation is based on Corollary 3 of arXiv:2502.16404. The OTOC
    is defined as F(W, V_t) where W is the static operator and V_t is the
    time-evolved one. This function computes E[F(static_op, (initial_evolved_op)_t)].

    The formula is:
    1 - 2 * |{T in C(V) : {W, T} = 0}| / |C(V)|
    where V is the initial_evolved_op, W is the static_op, and {W, T}
    is the anticommutator. If both operators are in the same component,
    the average OTOC is 1.0.

    Args:
        static_op (PauliString): The operator that is not time-evolved.
        initial_evolved_op (PauliString): The initial state of the operator
                                           that is time-evolved.
        system_generators (PauliStringCollection): The generators of the system.

    Returns:
        The average OTOC value.
    """
    num_qubits = system_generators.get_size()
    if len(static_op) != num_qubits or len(initial_evolved_op) != num_qubits:
        raise ValueError(
            "All operators and system generators must act on the same "
            "number of qubits."
        )

    # 3.1 Apply the logic from Corollary 3
    # Case 1: If the static op is in the same component, OTOC is 1.
    if static_op == initial_evolved_op:
        return 1.0

    # 1. Build the commutator graph
    all_paulis = list(PauliString(n=num_qubits).gen_all_pauli_strings())
    nodes = [str(p) for p in all_paulis]
    edges = set()

    for p in all_paulis:
        for g in system_generators:
            q = p.adjoint_map(g)
            if q is not None:
                edges.add(tuple(sorted((str(p), str(q)))))

    comm_graph = nx.Graph(list(edges))
    comm_graph.add_nodes_from(nodes)

    # 3.2 Determine the component of the initial_evolved_op
    # 2. Find the component of the operator that is being time-evolved.
    try:
        component_of_evolved_op = nx.node_connected_component(
            comm_graph, str(initial_evolved_op)
        )
    except nx.NetworkXError:
        return 1.0  # Should be unreachable if op is valid

    if str(static_op) in component_of_evolved_op:
        return 1.0

    # 3.3 Count the number of operators in the component that anticommute with the static_op
    # Case 2: Different components. Use the counting formula.
    component_size = len(component_of_evolved_op)
    if component_size == 0:
        return 1.0  # Should be unreachable

    # Count operators T in the component that anticommute with the static_op
    anticommuting_count = 0
    for t_str in component_of_evolved_op:
        op_t = PauliString(pauli_str=t_str)
        if not static_op.commutes_with(op_t):
            anticommuting_count += 1

    return 1.0 - 2.0 * (anticommuting_count / component_size)
