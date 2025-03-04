from itertools import combinations
from paulie.common.algebras import get_lie_algebra
from paulie.common.pauli import (
    gen_all_nodes,
    get_array_pauli_strings,
    get_pauli_array,
    is_commutative,
    multiply_pauli_arrays,
)
from bitarray import bitarray


def get_nested(p: bitarray) -> list[list[bitarray]]:
    """Finds pairs of Pauli operators that multiply to `p` and do not commute."""
    n = len(p) // 2
    return [
        [a.copy(), b.copy()]
        for a, b in combinations([g.copy() for g in gen_all_nodes(n)], 2)
        if a != p and b != p and not is_commutative(a, b) and multiply_pauli_arrays(a, b) == p
    ]


def get_nested_strings(pauli_string: str) -> list[list[str]]:
    """Returns nested pairs as Pauli strings."""
    return [
        get_array_pauli_strings(pair) 
        for pair in get_nested(get_pauli_array(pauli_string))
    ]


def get_nested_nodes_in_algebra(name: str) -> list[bitarray]:
    """Returns unique nodes in a nested algebra."""
    return list({
        node for g in get_lie_algebra()[name]
        for pair in get_nested(get_pauli_array(g)) 
        for node in pair
    })
