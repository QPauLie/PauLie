from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection


# Animation building transformation anti-commutation graph
# generators - list of generators
def animation_anti_commutation_graph(generators: PauliStringCollection, storage=None, interval:int =1000, initGraph:bool=False):
    record = RecordGraph()
    generators.set_record(record)
    generators.get_class()
    animation_graph(record, interval=interval, storage=storage)
