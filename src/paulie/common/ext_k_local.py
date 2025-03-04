from paulie.common.pauli import (
    get_I,
    get_Y,
    get_pauli_array,
    get_pauli_string,
    inc_pauli_array,
    is_sub_in_array,
)
from paulie.common.algebras import get_lie_algebra
from paulie.common.nested import get_nested_nodes_in_algebra
from typing import Generator


def gen_k_local(n: int, p: list[int], converter=None, used=None) -> Generator[list[int], None, None]:
    """Generates k-local Pauli operators."""
    if n < len(p) // 2:
        raise ValueError(f"Size must be greater than {len(p)//2}")

    used = used or []
    np = n - len(p) // 2

    for k in range(np + 1):
        left = get_I(k) + p + get_I(np - k)
        if left in used:
            continue

        used.append(left)
        yield left if converter is None else converter(left)


def gen_k_local_ext(n: int, p: list[int], converter=None, used=None) -> Generator[list[int], None, None]:
    """Generates extended k-local Pauli operators."""
    if n < 2:
        raise ValueError("Size must be greater than 1")

    used = used or []
    np = n - len(p) // 2

    for k in range(np + 1):
        left, left_one = get_I(k), get_Y(k)
        right, right_one = get_I(np - k), get_Y(np - k)

        while True:
            gen = get_I(0) + left + p + right

            if right == right_one:
                if left == left_one:
                    break
                left = inc_pauli_array(left)
            else:
                right = inc_pauli_array(right)

            if any(is_sub_in_array(u, gen) for u in used):
                continue

            used.append(gen)
            yield gen if converter is None else converter(gen)
 

def gen_k_local_by_string(n: int, pauli_string: str, used=None) -> Generator[list[int], None, None]:
    """Wrapper for gen_k_local using a Pauli string."""
    yield from gen_k_local(n, get_pauli_array(pauli_string), used=used)


def gen_k_local_generators(n: int, generators: list[str]) -> Generator[list[int], None, None]:
    """Generates k-local operators for a set of generators."""
    used = []
    for g in generators:
        yield from gen_k_local_by_string(n, g, used=used)


def gen_k_local_algebra_generators(n: int, name: str) -> Generator[list[int], None, None]:
    """Generates k-local algebra generators."""
    used = []
    for g in get_lie_algebra()[name]:
        yield from gen_k_local_by_string(n, g, used=used)


def gen_k_local_string(n: int, pauli_string: str, used=None) -> Generator[str, None, None]:
    """Wrapper for gen_k_local converting results to strings."""
    yield from gen_k_local(n, get_pauli_array(pauli_string), get_pauli_string, used=used)


def gen_k_local_ext_by_string(n: int, pauli_string: str) -> Generator[list[int], None, None]:
    """Wrapper for gen_k_local_ext using a Pauli string."""
    yield from gen_k_local_ext(n, get_pauli_array(pauli_string))


def gen_k_local_ext_string(n: int, pauli_string: str) -> Generator[str, None, None]:
    """Wrapper for gen_k_local_ext converting results to strings."""
    yield from gen_k_local_ext(n, get_pauli_array(pauli_string), get_pauli_string)


def get_k_local_generators(n: int, generators: list[str]) -> list[list[int]]:
    """Returns k-local operators as a list."""
    return list(gen_k_local_generators(n, generators))

               

def gen_k_local_string_generators(n: int, generators: list[str]) -> Generator[str, None, None]:
    """Generates k-local operators as strings."""
    used = []
    for g in generators:
        yield from gen_k_local_string(n, g, used=used)


def get_k_local_string_generators(n: int, generators: list[str]) -> list[str]:
    """Returns k-local operators as strings in a list."""
    return list(gen_k_local_string_generators(n, generators))


def get_k_local_algebra_generators(n: int, name: str) -> list[list[int]]:
    """Returns k-local algebra generators as a list."""
    return list(gen_k_local_algebra_generators(n, name))


def gen_k_local_string_algebra_generators(n: int, name: str) -> Generator[str, None, None]:
    """Generates k-local algebra generators as strings."""
    used = []
    for g in get_lie_algebra()[name]:
        yield from gen_k_local_string(n, g, used=used)


def get_k_local_string_algebra_generators(n: int, name: str) -> list[str]:
    """Returns k-local algebra generators as strings in a list."""
    return list(gen_k_local_string_algebra_generators(n, name))


def gen_k_local_nested_nodes_in_algebra(n: int, name: str) -> Generator[list[int], None, None]:
    """Generates k-local nested nodes in an algebra."""
    used = []
    for node in get_nested_nodes_in_algebra(name):
        yield from gen_k_local_ext(n, node, used=used)
        used.append(node)


def get_k_local_nested_nodes_in_algebra(n: int, name: str) -> list[list[int]]:
    """Returns k-local nested nodes in an algebra as a list."""
    return list(gen_k_local_nested_nodes_in_algebra(n, name))