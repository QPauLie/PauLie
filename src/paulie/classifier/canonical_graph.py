import networkx as nx
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.classification import Classification, Morph
from paulie.classifier.connected_canonicalizer import ConnectedCanonicalizer
from paulie.classifier.tracked_connected_canonicalizer import TrackedConnectedCanonicalizer


def canonical_graph(gens: PauliStringCollection, detect_independent_subset: bool=False
    )->Classification:
    verts, edges, _ = gens.get_graph()
    g = nx.Graph()
    g.add_nodes_from(verts)
    g.add_edges_from(edges)
    ccs = nx.connected_components(g)
    if detect_independent_subset:
        conn_canon = TrackedConnectedCanonicalizer()
    else:
        conn_canon = ConnectedCanonicalizer()
    classification = Classification()
    independent_subset = []
    for cc in ccs:
        vertex_stack = [get_pauli_string(s) for s in nx.dfs_preorder_nodes(g.subgraph(cc))]
        vertex_stack.reverse()
        conn_canon.build_canonical_graph(vertex_stack.copy())
        classification.add(conn_canon.get_morph())
    return classification
