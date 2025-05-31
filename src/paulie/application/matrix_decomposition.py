"""
    Fast Pauli basis matrix decomposition algorithm (arXiV:2403.11644).
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString, CODEC
from paulie.common.pauli_string_factory import get_identity
from collections import deque, Counter

def _mat_to_vec(A: np.ndarray) -> list[complex]:
    if A.shape[0] == 2:
        return [A[0, 0], A[1, 1], A[0, 1], A[1, 0]]
    n = A.shape[0] // 2
    return _mat_to_vec(A[:n,:n]) + _mat_to_vec(A[n:,n:]) + _mat_to_vec(A[:n,n:]) + _mat_to_vec(A[n:,:n])

def _fast_pauli_transform(a: list) -> None:
    h = 1
    while h < len(a):
        for i in range(0, len(a), h << 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                z = a[j + 2 * h]
                w = a[j + 3 * h]
                a[j] = (x + y) / 2
                a[j + h] = (x - y) / 2
                a[j + 2 * h] = (z + w) / 2
                a[j + 3 * h] = 1j * (z - w) / 2
        h <<= 2

def matrix_decomposition(matrix: np.ndarray, tol: float=1e-8) -> dict[PauliString, np.complex128]:
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
    B = _mat_to_vec(matrix)
    _fast_pauli_transform(B)
    pstr = get_identity(log2n)
    res = dict()
    for i in range(4 ** log2n):
        if np.abs(B[i]) > tol:
            res[pstr.copy()] = B[i]
        pstr.inc()
    return res

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
