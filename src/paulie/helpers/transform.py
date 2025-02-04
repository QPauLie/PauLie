from paulie.classifier.transform import transform_to_canonics
from paulie.helpers.recording import recording_graph
from paulie.common.pauli import get_array_pauli_arrays
from paulie.common.ext_k_local import get_k_local_generators


#Partial function for constructing a canonical graph
# generators - list of generators
# size - Generator extensions to size size
# debug - debbuging
# record - recording 
def app_transform_to_canonics(generators, size = 0, debug=False, record=None, initGraph=False):
    if size == 0:
        bitGenerators = get_array_pauli_arrays(generators)
    else:
        bitGenerators = get_k_local_generators(size, generators)
    if record is not None and initGraph:
        recording_graph(record, bitGenerators)
    return transform_to_canonics(bitGenerators, debug, record)
