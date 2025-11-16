"""
Pauli's String Factory. Responsible for creating instances
of Pauli strings of various implementations
"""
from typing import Generator, Union
from bitarray import bitarray
from paulie.common.pauli_string_linear import PauliStringLinear
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection

def get_identity(n: int) -> PauliString:
    """
    Get an identity of a given length
    Args: n - length of Pauli string
    Returns: identity
    """
    return PauliString(n=n)

def get_single(n: int, i: int, label: str) -> PauliString:
    """
    Get a Pauli string with a single value at position
    Args: n - length of Pauli string
          i - position in Pauli string
          label - PauliString in single position
    Returns: PauliString with set label
    """
    p = get_identity(n)
    p[i] = label
    return p

def get_last(n:int) -> PauliString:
    """
    Get an all Y of a given length
    Args: n - length of Pauli string
    Returns: all Y
    """
    return PauliString(bits = bitarray([1] * (2 * n)))

def get_pauli_string(o, n:int = None) -> PauliString|PauliStringCollection:
    """
    Get Pauli strings in their current representation.
    Args: 
         o - a Pauli string or a collection of Pauli strings.
         n - length of Pauli strings
    Returns: if o is a Pauli string, then its current instantiation value n is created
             otherwise PauliStringCollection is created - a collection of Pauli strings.
             Given n, the collection is expanded as k-local.
             Where k is the maximum length of a Pauli string in a given collection
    """
    if isinstance(o, str):
        return PauliString(pauli_str=o, n=n)
    if isinstance(o, PauliString):
        return o

    if isinstance(o, list):
        if len(o) > 0 and  isinstance(o[0], tuple):
            return PauliStringLinear(o)
    generators = PauliStringCollection([PauliString(pauli_str=p) if isinstance(p, str)
                 else PauliString(pauli_str=str(p)) for p in o])
    if n is not None:
        return PauliStringCollection(list(gen_k_local_generators(n, generators.get())))
    return generators


class Used:
    """
    Helper class for monitoring previously created Pauli strings
    """
    def __init__(self) -> None:
        """
        init class
        """
        self.clear()

    def clear(self) -> None:
        """
        Clear set
        Args: empty
        Returns: None
        """
        self.used = set()

    def append(self, p: PauliString) -> None:
        """Append to set
        Args:
            p: Pauli string to append in the set
        Returns:
            None
        """
        self.used.add(p)

    def is_used(self, p: PauliString) -> bool:
        """
        Checking a Pauli string in a set
        Args:
            p: Pauli string for checking in the used
        Returns:
            True if Pauli string in a set
        """
        return p in self.used


def gen_k_local(n: int, p: PauliString, used:Used=None) -> Generator[list[PauliString], None, None]:
    """
    Generates k-local Pauli strings.
    Args:
      n: a length of Paulistring
      p: local Pauli string 
      used: a repository of already generated strings
    Returns:
      k-local string generator
    Raises:
          ValueError: 
               n < len(p)
    """
    if n < len(p):
        raise ValueError(f"Size must be greater than {len(p)}")

    used = used or Used()
    n_p = n - len(p)
    for k in range(n_p + 1):
        left = get_identity(k) + p + get_identity(n_p - k)
        if used.is_used(left):
            continue

        used.append(left)
        yield left


def gen_k_local_generators(n: int,
                           generators: 'Union[list[str], list[PauliString], PauliStringCollection]',
                           used: Used = None) -> Generator[list[PauliString], None, None]:
    """
    Generates k-local operators for a set of generators.
    Args:
      n: a length of Paulistring
      generators: a collection of Pauli strings
      used: a repository of already generated strings
    Returns:
      k-local string generator
    """
    used = used or Used()
    longest = max(generators, key=len)
    for g in generators:
        if isinstance(g, str):
            g = get_pauli_string(g, n = len(longest))
        yield from gen_k_local(n, g, used=used)
