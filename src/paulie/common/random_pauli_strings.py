"""
    Module to generate random Pauli strings.
"""
from random import randint, choice


def get_random(n: int) -> str:
    """
    Get random Pauli string of length n.

    Args:
        n (int): Length of Pauli string.
    Returns: str
      Random Pauli string.
    """
    return ''.join([choice("IXYZ") for _ in range(n)])

def get_random_k_local(k: int, n: int) -> str:
    """
    Get random k-local Pauli string of length n.

    Args:
        k (int): Length of locality.
        n (int): Length of Pauli string.
    Returns: str
        Random k-local Pauli string of length n.

    Raises:
         ValueError: If k > n.
    """
    if k > n:
        raise ValueError("Invalid args: first arg greater than second")
    pauli_string = get_random(k)
    if k < n:
        position = randint(0, n - k)
        pauli_string = "I" * position + pauli_string + "I" * (n - position - k)
    return pauli_string

def get_random_list(n: int, size: int) -> list[str]:
    """
    Get list of random Pauli strings of length n.

    Args:
         n (int): Length of each Pauli string.
         size (int): Size of list.
    Returns: list[str]
        Random list of Pauli strings of length n.
    """
    return [get_random(n) for _ in range(size)]
