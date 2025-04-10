"""
A string representation of all Pauli strings of a given size. Needed to search for commutators
A bitarray is used for this as a simple way to sequentially increase Pauli strings
"""

from random import randint, choice


def get_random(n: int):
    """ Get random Pauli String lenght n """
    return''.join([choice("IXYZ") for _ in range(n)])

def get_random_k_local(k:int, n:int):
    """ Get random k local Pauli String lenght n """
    if k > n:
        raise ValueError("Invalid args: first arg grater than second")
    pauli_string = get_random(k)
    if k < n:
        position = randint(0, n-k)
        pauli_string = "".join("I" for _ in range(position)) + pauli_string + "".join(["I" for _ in range(position+k, n)])
    return pauli_string

def get_random_list(n:int, size: int):
    """ Get random list of Pauli String lenght n """
    return [get_random(n) for _ in range(size)]
