from paulie.helpers.transform import app_transform_to_canonics
from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph


# Animation building transformation anti-commutation graph
# generators - list of generators
# size - Generator extensions to size size
def animation_anti_commutation_graph(generators, size=0, storage=None, interval=1000, initGraph=False):
    record = RecordGraph()
    app_transform_to_canonics(generators, size=size, record=record, initGraph=initGraph)
    animation_graph(record, storage=storage, interval=interval)



