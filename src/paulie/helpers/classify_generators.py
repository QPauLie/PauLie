from paulie.classifier.classify import classify
from paulie.helpers.recording import recording_graph
from paulie.common.pauli import get_array_pauli_arrays
from paulie.common.ext_k_local import get_k_local_generators


#Partial function for constructing a canonical graph
# generators - list of generators
# size - Generator extensions to size size
# debug - debbuging
# record - recording 
def classify_generators(generators, size = 0, debug=False, record=None, initGraph=False):
    if size == 0:
        bit_generators = get_array_pauli_arrays(generators)
    else:
        bit_generators = get_k_local_generators(size, generators)

    if record is not None and initGraph:
        recording_graph(record, bit_generators)
    return classify(bit_generators, debug, record)
