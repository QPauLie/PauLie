"""
    Fast Pauli basis matrix decomposition algorithm.
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import get_identity
from bitarray import bitarray

def _pauli_ord(row: np.ndarray, col: np.ndarray, n: int) -> None:
    if n == 1:
        row[0], col[0] = 0, 0
        row[1], col[1] = 1, 1
        row[2], col[2] = 0, 1
        row[3], col[3] = 1, 0
        return
    _pauli_ord(row, col, n - 1)
    pw = 1 << (2 * (n - 1))
    row[pw: 2 * pw] = row[:pw] + (1 << (n - 1))
    col[pw: 2 * pw] = col[:pw] + (1 << (n - 1))
    row[2 * pw: 3 * pw] = row[:pw]
    col[2 * pw: 3 * pw] = col[:pw] + (1 << (n - 1))
    row[3 * pw: 4 * pw] = row[:pw] + (1 << (n - 1))
    col[3 * pw: 4 * pw] = col[:pw]

def _mat_to_vec(A: np.ndarray, log2n: int) -> np.ndarray:
    row = np.zeros(4 ** log2n, dtype=np.int64)
    col = np.zeros(4 ** log2n, dtype=np.int64)
    _pauli_ord(row, col, log2n)
    flatindex = (1 << log2n) * row + col
    return A.reshape(-1)[flatindex].astype(np.complex128)

def _one_bit_indices(x: np.int64, n: int) -> list:
    x = int(x)
    one_bit_indexes = []
    while x:
        v = x & -x
        one_bit_indexes.append(n - v.bit_length())
        x ^= v
    return one_bit_indexes

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
    B = _mat_to_vec(matrix, log2n)
    h = 1
    while h < B.shape[0]:
        for i in range(0, B.shape[0], 4 * h):
            B[i : i + h], B[i + h : i + 2 * h] = (B[i : i + h] + B[i + h : i + 2 * h]) / 2, (B[i : i + h] - B[i + h : i + 2 * h]) / 2
            B[i + 2 * h : i + 3 * h], B[i + 3 * h : i + 4 * h] = (B[i + 2 * h : i + 3 * h] + B[i + 3 * h : i + 4 * h]) / 2, 1j * (B[i + 2 * h : i + 3 * h] - B[i + 3 * h : i + 4 * h]) / 2
        h *= 4
    pstr = get_identity(log2n)
    res = dict()
    # Find nonzero coefficient indices
    nnz_coeff_inds = np.where(np.abs(B) > tol)[0]
    if nnz_coeff_inds[0] == 0:
        res[pstr.copy()] = B[0]
        nnz_coeff_inds = nnz_coeff_inds[1:]
    npb = np.zeros(2 * log2n, dtype=np.bool)
    j = 0
    for i in nnz_coeff_inds:
        # Manually set bits for performance
        npb[_one_bit_indices(i ^ j, 2 * log2n)] ^= True
        pstr.bits = bitarray()
        pstr.bits.pack(npb.tobytes())
        pstr.bits_even = pstr.bits[::2]
        pstr.bits_odd = pstr.bits[1::2]
        res[pstr.copy()] = B[i]
        j = i
    return res

def matrix_decomposition_diagonal(diag: np.ndarray, tol: float=1e-8) -> dict[PauliString, np.complex128]:
    if diag.ndim != 1:
        raise ValueError("matrix must be a 1D ndarray")
    if diag.shape[0] == 1:
        raise ValueError(f"input cannot be scalar")
    if int(diag.shape[0]).bit_count() != 1:
        raise ValueError(f"expected 1D ndarray with power of two \
                         length but length is {diag.shape[0]}")
    n = diag.shape[0]
    log2n = int(n).bit_length() - 1
    B = diag.astype(np.complex128)
    h = 1
    while h < B.shape[0]:
        for i in range(0, B.shape[0], 2 * h):
            B[i : i + h], B[i + h : i + 2 * h] = (B[i : i + h] + B[i + h : i + 2 * h]) / 2, (B[i : i + h] - B[i + h : i + 2 * h]) / 2
        h *= 2
    pstr = get_identity(log2n)
    res = dict()
    # Find nonzero coefficient indices
    nnz_coeff_inds = np.where(np.abs(B) > tol)[0]
    if nnz_coeff_inds[0] == 0:
        res[pstr.copy()] = B[0]
        nnz_coeff_inds = nnz_coeff_inds[1:]
    npb = np.zeros(log2n, dtype=np.bool)
    j = 0
    for i in nnz_coeff_inds:
        # Manually set bits for performance
        npb[_one_bit_indices(i ^ j, log2n)] ^= True
        pstr.bits_odd = bitarray()
        pstr.bits_odd.pack(npb.tobytes())
        pstr.bits[1::2] = pstr.bits_odd
        res[pstr.copy()] = B[i]
        j = i
    return res
