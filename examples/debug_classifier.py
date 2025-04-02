from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.algebras import get_lie_algebra
from time import perf_counter
from paulie.helpers.recording import RecordGraph



def debug_classification(n, name):
    print(f"Debugging classification for algebra {name} size {n}")
    print("--------------------------------------------------")

    generators = p(get_lie_algebra(name), n = n)
    generators.set_debug(False)

    print("--------------------------------------------------")
    print(f"algebra = {generators.get_class().get_algebra()}")

if __name__ == '__main__':
    start_time = perf_counter()
    debug_classification(10, "a22")
    end_time = perf_counter()
    print(f"time {end_time - start_time: 0.4f} sec.")
