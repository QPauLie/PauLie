"""
    Fast Pauli basis matrix decomposition algorithm.
"""
import numpy as np

def _pauli_ord(row: np.ndarray, col: np.ndarray, n: int) -> None:
    r"""Recursively constructs the Pauli ordering."""
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

def _mat_to_vec(matrix: np.ndarray) -> np.ndarray:
    r"""
    Vectorizes a matrix in Pauli order.

    Given a matrix

    .. math::

        M =
        \begin{bmatrix}
            A & B \\
            C & D
        \end{bmatrix}

    we vectorize it as

    .. math::

        \text{vec}(M) =
        \begin{bmatrix}
            A \\
            B \\
            C \\
            D
        \end{bmatrix}

    Args:
        matrix (np.ndarray): Matrix to vectorize.
    Returns:
        np.ndarray: Vector of matrix.
    """
    log2n = int(matrix.shape[0]).bit_length() - 1
    row = np.zeros(4 ** log2n, dtype=np.int64)
    col = np.zeros(4 ** log2n, dtype=np.int64)
    _pauli_ord(row, col, log2n)
    flat_index = (1 << log2n) * row + col
    return matrix.reshape(-1)[flat_index].astype(np.complex128)

def matrix_decomposition(matrix: np.ndarray) -> np.ndarray:
    r"""
    Returns the weight vector corresponding to the Pauli basis decomposition of a matrix.

    Args:
        matrix (np.ndarray): Matrix to be decomposed.

    Returns:
        np.ndarray: Weight vector corresponding to the Pauli basis decomposition of a matrix.
    """
    if matrix.ndim != 2:
        raise ValueError("matrix must be a 2D ndarray")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"expected square matrix but matrix dimensions \
                         are ({matrix.shape[0]}, {matrix.shape[1]})")
    if matrix.shape[0] == 1:
        raise ValueError("input must be a matrix, not a scalar")
    if int(matrix.shape[0]).bit_count() != 1:
        raise ValueError(f"expected square matrix with power of two \
                         dimensions but matrix dimensions are \
                         ({matrix.shape[0]}, {matrix.shape[1]})")
    b = _mat_to_vec(matrix)
    h = 1
    while h < b.shape[0]:
        for i in range(0, b.shape[0], 4 * h):
            x, y = b[i : i + h], b[i + h : i + 2 * h]
            b[i : i + h], b[i + h : i + 2 * h] = (x + y) / 2, (x - y) / 2
            z, w = b[i + 2 * h : i + 3 * h], b[i + 3 * h : i + 4 * h]
            b[i + 2 * h : i + 3 * h], b[i + 3 * h : i + 4 * h] = (z + w) / 2, 1j * (z - w) / 2
        h *= 4
    return b

def matrix_decomposition_diagonal(diag: np.ndarray) -> np.ndarray:
    r"""
    Returns the weight vector corresponding to the Pauli basis decomposition of a diagonal matrix.

    Args:
        diag (np.ndarray): Main diagonal of the diagonal matrix to be decomposed.
    Returns:
        np.ndarray: Weight vector corresponding to the Pauli basis decomposition of a diagonal matrix.
    """
    if diag.ndim != 1:
        raise ValueError("matrix must be a 1D ndarray")
    if diag.shape[0] == 1:
        raise ValueError("input cannot be scalar")
    if int(diag.shape[0]).bit_count() != 1:
        raise ValueError(f"expected 1D ndarray with power of two \
                         length but length is {diag.shape[0]}")
    b = diag.astype(np.complex128)
    h = 1
    while h < b.shape[0]:
        for i in range(0, b.shape[0], 2 * h):
            x, y = b[i : i + h], b[i + h : i + 2 * h]
            b[i : i + h], b[i + h : i + 2 * h] = (x + y) / 2, (x - y) / 2
        h *= 2
    return b
