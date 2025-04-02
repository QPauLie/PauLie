from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.algebras import get_lie_algebra
from time import perf_counter
from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph



def recording_classification(n, name):
    print(f"Debugging classification for algebra {name} size {n}")
    print("--------------------------------------------------")

    generators = p(get_lie_algebra(name), n = n)
    generators.set_debug(False)

    record = RecordGraph()
    generators.set_record(record)

    print("--------------------------------------------------")
    print(f"algebra = {generators.get_class().get_algebra()}")
    print(f"record = {record}")
    animation_graph(record, storage={"filename": "data/recording.gif", "writer": "pillow"})

if __name__ == '__main__':
    start_time = perf_counter()
    recording_classification(9, "a22")
    end_time = perf_counter()
    print(f"time {end_time - start_time: 0.4f} sec.")
