"""
Function for constructing a canonical graph
"""
import networkx as nx
from paulie import get_pauli_string
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.classifier.classification import Classification
from paulie.classifier.connected_canonicalizer import ConnectedCanonicalizer
from paulie.classifier.tracked_connected_canonicalizer import TrackedConnectedCanonicalizer


def canonical_graph(gens: PauliStringCollection, detect_independent_subset: bool=False
    )->Classification:
    """
    Function for constructing a canonical graph
    Args:
        gens (PauliStringCollection): Algebra Generator Collection
        detect_independent_subset (bool): Detect subset independent generators
    Returns:
        (Classification): Algebra classifier 
    """
    verts, edges, _ = gens.get_graph()
    g = nx.Graph()
    g.add_nodes_from(verts)
    g.add_edges_from(edges)
    ccs = nx.connected_components(g)
    classification = Classification()
    for cc in ccs:
        vertex_stack = [get_pauli_string(s) for s in nx.dfs_preorder_nodes(g.subgraph(cc))]
        vertex_stack.reverse()
        if detect_independent_subset:
            conn_canon = TrackedConnectedCanonicalizer()
        else:
            conn_canon = ConnectedCanonicalizer()
        classification.add(conn_canon.build_canonical_graph(vertex_stack.copy()))
    return classification
