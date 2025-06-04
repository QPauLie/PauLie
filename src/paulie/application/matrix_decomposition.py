"""
    Fast Pauli basis matrix decomposition algorithm.
"""
import numpy as np

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

def matrix_decomposition(matrix: np.ndarray, tol: float=1e-8) -> np.ndarray:
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
    return B

def matrix_decomposition_diagonal(diag: np.ndarray, tol: float=1e-8) -> np.ndarray:
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
    return B
