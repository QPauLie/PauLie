from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.all_pauli_strings import get_random_list
from time import perf_counter
from itertools import combinations
from pauliarray import PauliArray


def perfomence_mult(g1):
    for a, b in combinations(g1, 2):
        a @ b

def perfomence_check_comm(g1):
    for a, b in combinations(g1, 2):
        a ^ b

def perfomence_commutator_pauliarray(g1):
    for a, b in combinations(g1, 2):
        a.commute_with(b)

if __name__ == "__main__":
#    set_factory(PauliStringType.BITARRAY)
#    exit()
    g = get_random_list(100, 1000)
    g1 = p(g)
    start_time = perf_counter()
    perfomence_mult(g1)
    end_time = perf_counter()
    print(f"@ time {end_time - start_time: 0.4f} sec.")

    start_time = perf_counter()
    perfomence_check_comm(g1)
    end_time = perf_counter()
    print(f"| time {end_time - start_time: 0.4f} sec.")

    paulis = PauliArray.from_labels(g)

    start_time = perf_counter()
    perfomence_commutator_pauliarray(paulis)
    end_time = perf_counter()
    print(f"pauliarray | time {end_time - start_time: 0.4f} sec.")
