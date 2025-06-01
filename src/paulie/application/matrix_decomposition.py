"""
    Fast Pauli basis matrix decomposition algorithm (arXiV:2403.11644).
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import get_identity

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

def matrix_decomposition(matrix: np.ndarray, tol: float=1e-8) -> dict[PauliString, complex]:
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
    pind = 0
    if np.abs(B[0]) > tol:
        res[pstr.copy()] = B[0]
    for i in range(1, 4 ** log2n):
        # Iterate in Gray code order and manually set bits for performance
        ind = i ^ (i >> 1)
        bpos = (ind ^ pind).bit_length() - 1
        pstr.bits[2 * log2n - 1 - bpos] ^= 1
        if bpos & 1:
            pstr.bits_even[log2n - 1 - (bpos >> 1)] ^= 1
        else:
            pstr.bits_odd[log2n - 1 - (bpos >> 1)] ^= 1
        if np.abs(B[ind]) > tol:
            res[pstr.copy()] = B[ind]
        pind = ind
    return res

def _fast_walsh_hadamard_transform(a: np.ndarray) -> None:
    h = 1
    while h < len(a):
        for i in range(0, len(a), h << 1):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = (x + y) / 2
                a[j + h] = (x - y) / 2
        h <<= 1

def matrix_decomposition_diagonal(diag: np.ndarray, tol: float=1e-8) -> dict[PauliString, complex]:
    if diag.ndim != 1:
        raise ValueError("matrix must be a 1D ndarray")
    if diag.shape[0] == 1:
        raise ValueError(f"input cannot be scalar")
    if int(diag.shape[0]).bit_count() != 1:
        raise ValueError(f"expected 1D ndarray with power of two \
                         length but length is {diag.shape[0]}")
    n = diag.shape[0]
    log2n = int(n).bit_length() - 1
    B = diag.copy()
    _fast_walsh_hadamard_transform(B)
    pstr = get_identity(log2n)
    res = dict()
    pind = 0
    if np.abs(B[0]) > tol:
        res[pstr.copy()] = B[0]
    for i in range(1, 2 ** log2n):
        # Iterate in Gray code order and manually set bits for performance
        ind = i ^ (i >> 1)
        bpos = (ind ^ pind).bit_length() - 1
        pstr.bits_even[log2n - 1 - bpos] ^= 1
        pstr.bits[2 * (log2n - bpos) - 1] ^= 1
        if np.abs(B[ind]) > tol:
            res[pstr.copy()] = B[ind]
        pind = ind
    return res
