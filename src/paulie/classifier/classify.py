from paulie.classifier.classification import Classification
from paulie.classifier.morph_factory import MorphFactory
from paulie.graphs.subgraphs import get_subgraphs



def classify(vertices, debug=False, record=None):
    """Split the original graph into connected subgraphs and transform them to a cononical type."""
    subgraphs = get_subgraphs(vertices)
    classification = Classification()

    for subgraph in subgraphs:
        morph_factory = MorphFactory(debug, record)
        classification.add(morph_factory.build(subgraph).get_morph())
  
    return classification

