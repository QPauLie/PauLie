"""
    Fast Pauli basis matrix decomposition algorithm (arXiV:2403.11644).
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import get_pauli_string
from collections import deque

def matrix_decomposition(matrix: np.ndarray) -> list[tuple[complex, PauliString]]:
    # Shape checks
    if matrix.ndim != 2:
        raise ValueError("matrix must be a 2D ndarray")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"expected square matrix but matrix dimensions \
                         are ({matrix.shape[0]}, {matrix.shape[1]})")
    if matrix.shape[0] == 1:
        raise ValueError(f"input must be a matrix, not a scalar")
    if int(matrix.shape[0]).bit_count() != 1:
        raise ValueError(f"expected square matrix with power of two \
                         dimensions but matrix dimensions are \
                         ({matrix.shape[0]}, {matrix.shape[1]})")
    # Numpy dtype checks?
    # General decomposition (future additions: diagonal/tridiagonal/band-diagonal cases?)
    return _decomp_general(matrix)

def _compute_coeff(k0: int, k: np.ndarray, m: np.ndarray, nY: int, matrix: np.ndarray) -> complex:
    n = matrix.shape[0]
    coeff = 0
    for j in range(n):
        coeff += m[j] * matrix[k0 + k[j]][j] / n
    return (-1j) ** (nY % 4) * coeff

def _compute_k0(paulistr: str) -> int:
    k0, p = 0, 1
    for c in reversed(paulistr):
        if c == 'X' or c == 'Y':
            k0 += p
        p <<= 1
    return k0

def _update_general(k: np.ndarray, m: np.ndarray, l: int, node: str) -> None:
    # The nature of the updates forces us to walk the tree in IXYZ order
    if node == 'I':
        k[2 ** l : 2 ** (l + 1)] = k[0 : 2 ** l] + 2 ** l
        m[2 ** l : 2 ** (l + 1)] = m[0 : 2 ** l]
    elif node == 'X':
        k[2 ** l : 2 ** (l + 1)] -= 2 ** (l + 1)
    elif node == 'Y':
        m[2 ** l : 2 ** (l + 1)] *= -1
    elif node == 'Z':
        k[2 ** l : 2 ** (l + 1)] += 2 ** (l + 1)

def _decomp_general(matrix: np.ndarray) -> list[tuple[complex, PauliString]]:
    n = matrix.shape[0]
    k = np.zeros(n, dtype=np.int64)
    m = np.zeros(n, dtype=np.int64)
    nY = 0
    coeffs = []
    # Initialize m[0] (k[0] is computed at leaves)
    m[0] = 1
    # Walk the tree
    dfs_stack = deque([('I', 0), ('X', 0), ('Y', 0), ('Z', 0)])
    current_string = deque()
    while dfs_stack:
        node, depth = dfs_stack.popleft()
        while len(current_string) > depth:
            rnode = current_string.popleft()
            if rnode == 'Y':
                nY -= 1
        if node == 'Y':
            nY += 1
        current_string.appendleft(node)
        _update_general(k, m, depth, node)
        if 2 ** len(current_string) == n:
            paulistr = ''.join(current_string)
            k0 = _compute_k0(paulistr)
            coeff = _compute_coeff(k0, k, m, nY, matrix)
            if not np.isclose(np.abs(coeff), 0):
                coeffs.append((coeff, get_pauli_string(paulistr)))
        else:
            dfs_stack.appendleft(('I', 1 + depth))
            dfs_stack.appendleft(('X', 1 + depth))
            dfs_stack.appendleft(('Y', 1 + depth))
            dfs_stack.appendleft(('Z', 1 + depth))
    return coeffs
