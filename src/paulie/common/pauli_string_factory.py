"""
Module for creating instances of Pauli strings of various implementations.
"""

from collections.abc import Generator
from typing import overload
from bitarray import bitarray
from paulie.common.pauli_string_linear import PauliStringLinear
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection


def get_identity(n: int) -> PauliString:
    """
    Get an identity of a given length.

    Args:
        n (int): Length of Pauli string.
    Returns:
        PauliString: Identity of the given length.
    """
    return PauliString(n=n)


def get_single(n: int, i: int, label: str) -> PauliString:
    """
    Get a Pauli string with a single value at position.

    Args:
        n (int): Length of Pauli string.
        i (int): Position in Pauli string.
        label (str): Pauli at position.
    Returns:
        PauliString: PauliString with set label.
    """
    p = get_identity(n)
    p[i] = label
    return p


def get_last(n: int) -> PauliString:
    """
    Get a Pauli string with `Y` at all positions of a given length.

    Args:
        n (int): Length of Pauli string.
    Returns:
        PauliString: Pauli string of given length with `Y` at every position.
    """
    return PauliString(bits=bitarray([1] * (2 * n)))


@overload
def get_pauli_string(o: str, n: int | None = None) -> PauliString: ...

@overload
def get_pauli_string(o: PauliString, n: int | None = None) -> PauliString: ...

@overload
def get_pauli_string(
    o: list[tuple[float, str | PauliString]], n: int | None = None
) -> PauliStringLinear: ...

@overload
def get_pauli_string(
    o: list[str] | list[PauliString] | PauliStringCollection, n: int | None = None
) -> PauliStringCollection: ...

def get_pauli_string(
    o: (
        str
        | PauliString
        | list[tuple[float, str | PauliString]]
        | list[str]
        | list[PauliString]
        | PauliStringCollection
    ),
    n: int | None = None,
) -> PauliString | PauliStringLinear | PauliStringCollection:
    """
    Get Pauli strings in their current representation.

    Args:
         o (PauliString|PauliStringCollection): A Pauli string or a collection of Pauli strings.
         n (int, optional): Length of Pauli strings
    Returns:
        PauliString|PauliStringLinear|PauliStringCollection:
        If `o` is a Pauli string, then it is tensored with identities at the end until its
        length is `n`. Otherwise a collection of Pauli strings is created. Given n, the
        collection is expanded as a k-local set where k is the maximum length of a Pauli
        string in the given collection.
    """
    if isinstance(o, str):
        return PauliString(pauli_str=o, n=n)
    if isinstance(o, PauliString):
        return o

    if isinstance(o, list):
        if len(o) > 0 and isinstance(o[0], tuple):
            return PauliStringLinear(o)
    generators = PauliStringCollection(
        [
            (
                PauliString(pauli_str=p)
                if isinstance(p, str)
                else PauliString(pauli_str=str(p))
            )
            for p in o
        ]
    )
    if n is not None:
        return PauliStringCollection(list(gen_k_local_generators(n, generators.get())))
    return generators


def gen_k_local(
    n: int, p: PauliString, used: set[PauliString] = None
) -> Generator[PauliString, None, None]:
    """
    Generates k-local Pauli strings.

    Examples:
        >>> from paulie.common.pauli_string_factory import gen_k_local
        >>> from paulie import get_pauli_string as p
        >>> print([s for s in gen_k_local(4, p('XY'))])
        [PauliString(XYII), PauliString(IXYI), PauliString(IIXY)]

    Args:
        n (int): Length of Pauli string.
        p (PauliString): Local Pauli string.
        used (set[PauliString], optional): Set of already generated strings.
    Yields:
        k-local strings.
    Raises:
        ValueError: If the desired length is less than the length of the local Pauli string.
    """
    if n < len(p):
        raise ValueError(f"Size must be greater than {len(p)}")

    used = used or set()
    n_p = n - len(p)
    for k in range(n_p + 1):
        left = get_identity(k) + p + get_identity(n_p - k)
        if left in used:
            continue
        used.add(left)
        yield left


def get_all_k_local(
    n: int, generators: list[str] | list[PauliString] | PauliStringCollection
) -> PauliStringCollection:
    """
    Get all k-local Pauli strings for a set of generators.

    Examples:
        >>> from paulie.common.pauli_string_factory import get_all_k_local
        >>> from paulie import get_pauli_string as p
        >>> print(get_all_k_local(4, p(['XY', 'Z'])))
        PauliStringCollection([PauliString(XYII), PauliString(IXYI), PauliString(IIXY),
        PauliString(ZIII), PauliString(IZII), PauliString(IIZI), PauliString(IIIZ)])

    Args:
        n (int): Length of Pauli string.
        generators (list[str] | list[PauliString] | PauliStringCollection): Collection of Pauli
            strings.
    Returns:
        PauliStringCollection: All k-local Pauli strings.
    """
    return PauliStringCollection(list(gen_k_local_generators(n, generators)))


def gen_k_local_generators(
    n: int,
    generators: list[str] | list[PauliString] | PauliStringCollection,
    used: set[PauliString] = None,
) -> Generator[PauliString, None, None]:
    """
    Generates k-local operators for a set of generators.

    Examples:
        >>> from paulie.common.pauli_string_factory import gen_k_local_generators
        >>> from paulie import get_pauli_string as p
        >>> print([s for s in gen_k_local_generators(4, p(['XY', 'Z']))])
        [PauliString(XYII), PauliString(IXYI), PauliString(IIXY), PauliString(ZIII),
        PauliString(IZII), PauliString(IIZI)]

    Args:
        n (int): Length of Pauli string.
        generators (list[str] | list[PauliString] | PauliStringCollection): Collection of Pauli
            strings.
        used (set[PauliString]): Set of already generated strings.
    Yields:
        k-local strings.
    """
    used = used or set()
    generators_list = []
    if isinstance(generators, PauliStringCollection):
        generators_list = generators.get()
    else:
        generators_list = list(generators)

    longest = max(generators_list, key=len)
    for g in generators_list:
        if isinstance(g, str):
            g = get_pauli_string(g, n=len(longest))
        yield from gen_k_local(n, g, used=used)


def gen_all_pauli_strings(n: int) -> Generator[PauliString, None, None]:
    """
    Generates all Pauli strings of a given length.

    Args:
        n (int): Length of Pauli string.
    Yields:
        Generator[PauliString, None, None]: Generator of all Pauli strings.
    """
    pauli_string = get_identity(n)
    last = get_last(n)

    while pauli_string != last:
        yield pauli_string.copy()
        pauli_string.inc()
    yield pauli_string.copy()
