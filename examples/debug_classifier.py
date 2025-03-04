
from paulie.classifier.classify import classify
from paulie.common.ext_k_local import get_k_local_algebra_generators
from paulie.graphs.graph_view import get_graph_view
from paulie.helpers.drawing import plot_graph



def debug_classification(n, name, debug = True):
    print(f"Debugging transform for algebra {name} size {n}")
    print("--------------------------------------------------")
    generators = get_k_local_algebra_generators(n, name)
    classification = classify(generators, debug=debug)
    print("--------------------------------------------------")
    print(f"algebra = {classification.get_algebra()}")
    vertices, edges, edge_labels =  get_graph_view(generators)
    #plot_graph(vertices, edges, edge_labels)

    cannonics = classification.get_vertices()
    vertices, edges, edge_labels =  get_graph_view(cannonics)
    plot_graph(vertices, edges, edge_labels)


if __name__ == '__main__':
    #debug_transform(600, "a1", True)
    debug_classification(10, "a22", True)






