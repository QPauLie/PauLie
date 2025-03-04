from paulie.helpers.classify_generators import classify_generators
from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph


# Animation building transformation anti-commutation graph
# generators - list of generators
# size - Generator extensions to size size
def animation_anti_commutation_graph(generators, size=0, storage=None, interval=1000, initGraph=False):
    record = RecordGraph()
    classify_generators(generators, size=size, record=record, initGraph=initGraph)
    animation_graph(record, storage=storage, interval=interval)



