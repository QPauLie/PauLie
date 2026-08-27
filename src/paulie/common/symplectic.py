"""
    Symplectic (GF(2)) representation of Pauli strings.

    A Pauli string on ``n`` qubits is a vector ``v = [x | z]`` over GF(2) of length
    ``2n``, ignoring phases. Two Pauli strings commute iff their symplectic form

    .. math::

        \\langle v, w \\rangle = x_v \\cdot z_w + z_v \\cdot x_w \\bmod 2

    vanishes, and the product of two Pauli strings is ``v XOR w`` up to phase. In this
    picture questions about commutation become linear algebra over GF(2), which is what
    makes them polynomial rather than exponential in ``n``.

    ``PauliString.bits`` is interleaved per qubit as ``[x_0, z_0, x_1, z_1, ...]``, so
    ``bits[::2]`` is the X part and ``bits[1::2]`` the Z part.
"""

from collections.abc import Iterable

import numpy as np
from pauliebits import pauliebits

from paulie.common.pauli_string_bitarray import PauliString


def to_symplectic(generators: Iterable[PauliString]) -> np.ndarray:
    """
    Pack Pauli strings into their symplectic representation.

    Args:
        generators (Iterable[PauliString]): The Pauli strings. They must all act on the
            same number of qubits.
    Returns:
        numpy.ndarray: An ``(m, 2n)`` array of 0/1, laid out as ``[X | Z]``.
    """
    gens = list(generators)
    if not gens:
        return np.zeros((0, 0), dtype=np.uint8)
    n = len(gens[0])
    v = np.empty((len(gens), 2 * n), dtype=np.uint8)
    for i, p in enumerate(gens):
        bits = p.bits
        v[i, :n] = np.frombuffer(bits[::2].unpack(), dtype=np.uint8) & 1
        v[i, n:] = np.frombuffer(bits[1::2].unpack(), dtype=np.uint8) & 1
    return v


def from_symplectic(v: np.ndarray) -> list[PauliString]:
    """
    Unpack a symplectic representation back into Pauli strings.

    Inverse of :func:`to_symplectic`; phases are not represented and so are not restored.

    Args:
        v (numpy.ndarray): An ``(m, 2n)`` array of 0/1, laid out as ``[X | Z]``.
    Returns:
        list[PauliString]: The corresponding Pauli strings.
    """
    n = v.shape[1] // 2
    strings = []
    for row in v:
        x = pauliebits()
        x.pack(np.ascontiguousarray(row[:n], dtype=np.uint8).tobytes())
        z = pauliebits()
        z.pack(np.ascontiguousarray(row[n:], dtype=np.uint8).tobytes())
        interleaved = pauliebits(2 * n)
        interleaved[::2] = x
        interleaved[1::2] = z
        p = PauliString(n=n)
        p.bits = interleaved
        strings.append(p)
    return strings


def symplectic_gram(v: np.ndarray) -> np.ndarray:
    """
    Anticommutation graph of a symplectic representation.

    Args:
        v (numpy.ndarray): An ``(m, 2n)`` array as returned by :func:`to_symplectic`.
    Returns:
        numpy.ndarray: An ``(m, m)`` adjacency matrix with ``a[i, j] = 1`` iff ``P_i``
        and ``P_j`` anticommute. The diagonal is zero.
    """
    n = v.shape[1] // 2
    x = v[:, :n].astype(np.int32)
    z = v[:, n:].astype(np.int32)
    a = (x @ z.T + z @ x.T) & 1
    np.fill_diagonal(a, 0)
    return a.astype(np.int32)


def row_reduce(a: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """
    Reduced row echelon form over GF(2).

    Args:
        a (numpy.ndarray): A 0/1 matrix.
    Returns:
        tuple[numpy.ndarray, list[int]]: The reduced matrix and its pivot columns. The
        rank is the number of pivots.
    """
    m = a.copy().astype(np.uint8) % 2
    pivots: list[int] = []
    row = 0
    for col in range(m.shape[1]):
        if row >= m.shape[0]:
            break
        candidates = np.nonzero(m[row:, col])[0]
        if candidates.size == 0:
            continue
        src = row + int(candidates[0])
        if src != row:
            m[[row, src]] = m[[src, row]]
        others = np.nonzero(m[:, col])[0]
        others = others[others != row]
        m[others] ^= m[row]
        pivots.append(col)
        row += 1
    return m, pivots


def null_space(a: np.ndarray, width: int | None = None) -> np.ndarray:
    """
    Basis of the null space of a matrix over GF(2).

    Args:
        a (numpy.ndarray): A 0/1 matrix with ``width`` columns.
        width (int, optional): Number of columns, needed only when ``a`` has no rows.
    Returns:
        numpy.ndarray: A ``(width - rank, width)`` array whose rows are a basis of
        ``{v : a v = 0}``.
    """
    if width is None:
        width = a.shape[1]
    if a.shape[0] == 0:
        return np.eye(width, dtype=np.uint8)

    reduced, pivots = row_reduce(a)
    free = [c for c in range(width) if c not in set(pivots)]
    basis = np.zeros((len(free), width), dtype=np.uint8)
    for k, col in enumerate(free):
        basis[k, col] = 1
        for r, p in enumerate(pivots):
            basis[k, p] = reduced[r, col]
    return basis


def span(basis: np.ndarray) -> np.ndarray:
    """
    Every vector in the span of a GF(2) basis.

    The output has ``2 ** len(basis)`` rows and is therefore only usable for small
    bases; prefer working with the basis itself where possible.

    Args:
        basis (numpy.ndarray): A ``(k, width)`` array of basis vectors.
    Returns:
        numpy.ndarray: A ``(2 ** k, width)`` array of all linear combinations, starting
        with the zero vector.
    """
    k, width = basis.shape
    out = np.zeros((1 << k, width), dtype=np.uint8)
    for i in range(k):
        half = 1 << i
        out[half:2 * half] = out[:half] ^ basis[i]
    return out
