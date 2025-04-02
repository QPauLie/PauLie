from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.algebras import get_lie_algebra



def print_pauli_string(pauli_strings):
    print(f"count = {len(pauli_strings)}")
    print(f"{pauli_strings}")

def print_k_local(generators, n: int):
    print(f"k local pauli strings for size = {n}")
    pauli_strings = p(generators, n = n)
    print_pauli_string(pauli_strings)

def print_k_local_by_algebra(name: str, n: int):
    print(f"k local pauli strings for algebra = {name} size = {n} ")
    print_k_local(get_lie_algebra(name), n)

if __name__ == "__main__":

    v = p("XY",n=6) + "ZI"
    print(f"{v}")


    print_k_local(["XXZ", "YX"], 4)
    print_k_local_by_algebra("a5", 4)
    print_k_local(["IX", "ZXI"], 5)
    print_k_local_by_algebra("b3", 4)

    v = p("IXY", n = 7) + "XZI"
    print(f"{v}")

    g = p(["XX", "XY"], 4)
    print(f"{g}")

    print(f"{g + p('Z_2')}")

    print(f"X | Y {p('X')|p('Y')}")
    print(f"I | X {p('I')|p('X')}")

    print(f"X ^ Y {p('X')^p('Y')}")
    print(f"I ^ X {p('I')^p('X')}")

    print(f"X @ Y {p('X')@p('Y')}")
    print(f"I @ X {p('I')@p('X')}")
