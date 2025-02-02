

from paulie.common.algebras import get_algebra_generators
from paulie.common.pauli import (
    gen_all_nodes,
    get_array_pauli_strings,
    get_pauli_array,
    is_commutate,
    multi_pauli_arrays,
)


def get_nested(p):
    nested = []
    for a in gen_all_nodes(len(p)//2):
        if a != p:
            for b in gen_all_nodes(len(p)//2):
                if b != p and a < b:
                    if is_commutate(a, b) is False:
                        if multi_pauli_arrays(a, b) == p:
                            nested.append([a.copy(), b.copy()])
    return nested


def get_nested_by_string(pauli_string):
    return get_nested(get_pauli_array(pauli_string))


def get_nested_strings(pauli_string):
    nested = get_nested_by_string(pauli_string)
    return  list(map(get_array_pauli_strings, nested))


def get_nested_algebra(name):
    generators = get_algebra_generators(name)
    nested = []
    for g in generators:
        nested += get_nested_by_string(g)
    return nested


def get_nested_string_algebra(name):
    generators = get_algebra_generators(name)
    nested = []
    for g in generators:
        nested += get_nested_strings(g)
    return nested


def get_nodes_in_nested(nested):
    nodes = []
    for pair in nested:
        for node in pair:
            if node not in nodes:
                nodes.append(node)
    return nodes


def get_nested_nodes_in_algebra(name):
    return get_nodes_in_nested(get_nested_algebra(name))