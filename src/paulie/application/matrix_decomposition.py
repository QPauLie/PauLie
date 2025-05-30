"""
    Fast Pauli basis matrix decomposition algorithm (arXiV:2403.11644).
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString, CODEC
from paulie.common.pauli_string_factory import get_pauli_string, get_identity
from collections import deque, Counter

def matrix_decomposition(matrix: np.ndarray, tol: float=1e-8) -> dict[PauliString, complex]:
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
    n = matrix.shape[0]
    log2n = int(n).bit_length() - 1
    j = np.arange(n)
    k = np.zeros(n, dtype=np.int64)
    m = np.zeros(n, dtype=np.int8)
    k0 = 0
    m[0] = 1
    coeffs = Counter()
    dfs_stack = deque([('Z', 0), ('Y', 0), ('X', 0), ('I', 0)])
    phase = [1, -1j, -1, 1j]
    pstr = get_identity(log2n)
    while dfs_stack:
        node, depth = dfs_stack.pop()
        if node == 'I':
            k[2 ** depth : 2 ** (depth + 1)] = k[0 : 2 ** depth] + 2 ** depth
            m[2 ** depth : 2 ** (depth + 1)] = m[0 : 2 ** depth]
        elif node == 'X':
            k0 += 2 ** depth
            k[2 ** depth : 2 ** (depth + 1)] -= 2 ** (depth + 1)
        elif node == 'Y':
            m[2 ** depth : 2 ** (depth + 1)] *= -1
        elif node == 'Z':
            k0 -= 2 ** depth
            k[2 ** depth : 2 ** (depth + 1)] += 2 ** (depth + 1)
        # Set bits manually for performance
        nodebits = CODEC[node]
        pstr.bits[2 * (log2n - depth - 1)] = nodebits[0]
        pstr.bits[2 * (log2n - depth - 1) + 1] = nodebits[1]
        pstr.bits_even[log2n - depth - 1] = nodebits[0]
        pstr.bits_odd[log2n - depth - 1] = nodebits[1]
        if depth + 1 == log2n:
            coeff = np.dot(m, matrix[k0 + k, j]) / n
            if np.abs(coeff) > tol:
                coeffs[pstr.copy()] += phase[(pstr.bits_even & pstr.bits_odd).count() % 4] * coeff
        else:
            dfs_stack.extend([('Z', depth + 1), ('Y', depth + 1), ('X', depth + 1), ('I', depth + 1)])
    return dict(coeffs)

def matrix_decomposition_diagonal(diag: np.ndarray, tol: float=1e-8) -> dict[PauliString, complex]:
    # Shape checks
    if diag.ndim != 1:
        raise ValueError("matrix must be a 1D ndarray")
    if diag.shape[0] == 1:
        raise ValueError(f"input cannot be scalar")
    if int(diag.shape[0]).bit_count() != 1:
        raise ValueError(f"expected 1D ndarray with power of two \
                         length but length is {diag.shape[0]}")
    n = diag.shape[0]
    log2n = int(n).bit_length() - 1
    m = np.zeros(n, dtype=np.int8)
    m[0] = 1
    coeffs = Counter()
    dfs_stack = deque([('Z', 0), ('I', 0)])
    pstr = get_identity(log2n)
    while dfs_stack:
        node, depth = dfs_stack.pop()
        if node == 'I':
            m[2 ** depth : 2 ** (depth + 1)] = m[0 : 2 ** depth]
        elif node == 'Z':
            m[2 ** depth : 2 ** (depth + 1)] *= -1
        # Set bits manually for performance
        nodebits = CODEC[node]
        pstr.bits[2 * (log2n - depth - 1)] = nodebits[0]
        pstr.bits[2 * (log2n - depth - 1) + 1] = nodebits[1]
        pstr.bits_even[log2n - depth - 1] = nodebits[0]
        pstr.bits_odd[log2n - depth - 1] = nodebits[1]
        if depth + 1 == log2n:
            coeff = np.dot(m, diag) / n
            if np.abs(coeff) > tol:
                coeffs[pstr.copy()] += coeff
        else:
            dfs_stack.extend([('Z', depth + 1), ('I', depth + 1)])
    return dict(coeffs)
