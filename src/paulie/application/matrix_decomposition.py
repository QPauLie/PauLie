"""
    Fast Pauli basis matrix decomposition algorithm (arXiV:2403.11644).
"""
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import get_pauli_string, get_identity
from collections import deque, Counter
from concurrent.futures import ProcessPoolExecutor

def matrix_decomposition(matrix: np.ndarray) -> dict[PauliString, complex]:
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
    n = matrix.shape[0]
    log2n = int(n).bit_length() - 1
    # Hardcoded parallel processing by splitting into 4 blocks
    if log2n >= 3:
        res = Counter()
        ops = [get_pauli_string('I'),
               get_pauli_string('X'),
               get_pauli_string('Y'),
               get_pauli_string('Z'),]
        Ms = [matrix[0 : n // 2, 0 : n // 2],   # -> 0.5 * (I + Z)
              matrix[0 : n // 2, n // 2 : n],   # -> 0.5 * (X + iY)
              matrix[n // 2 : n, 0 : n // 2],   # -> 0.5 * (X - iY)
              matrix[n // 2 : n, n // 2 : n],]  # -> 0.5 * (I - Z)
        # Break up matrix
        with ProcessPoolExecutor() as executor:
            for i, decomp in enumerate(executor.map(_decomp_general, Ms)):
                for pstr, coeff in decomp.items():
                    if i == 0:
                        res[ops[0] + pstr] += coeff / 2
                        res[ops[3] + pstr] += coeff / 2
                    elif i == 1:
                        res[ops[1] + pstr] += coeff / 2
                        res[ops[2] + pstr] += 1j * coeff / 2
                    elif i == 2:
                        res[ops[1] + pstr] += coeff / 2
                        res[ops[2] + pstr] -= 1j * coeff / 2
                    elif i == 3:
                        res[ops[0] + pstr] += coeff / 2
                        res[ops[3] + pstr] -= coeff / 2
        return res
    return dict(_decomp_general(matrix))

def _compute_coeff(k0: int, k: np.ndarray, m: np.ndarray, nY: int, matrix: np.ndarray) -> complex:
    n = matrix.shape[0]
    coeff = 0
    for j in range(n):
        coeff += m[j] * matrix[k0 + k[j]][j] / n
    return (-1j) ** (nY % 4) * coeff

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

def _decomp_general(matrix: np.ndarray) -> dict[PauliString, complex]:
    n = matrix.shape[0]
    log2n = int(n).bit_length() - 1
    k = np.zeros(n, dtype=np.int64)
    m = np.zeros(n, dtype=np.int64)
    nY = 0
    coeffs = Counter()
    # Initialize m[0] (k[0] is computed at leaves)
    m[0] = 1
    # Walk the tree
    dfs_stack = deque([('I', 0), ('X', 0), ('Y', 0), ('Z', 0)])
    pstr = get_identity(log2n)
    curlen = 0
    pos = log2n - 1
    while dfs_stack:
        node, depth = dfs_stack.popleft()
        while curlen > depth:
            rnode = pstr[pos + 1]
            pstr[pos + 1] = 'I'
            if str(rnode) == 'Y':
                nY -= 1
            curlen -= 1
            pos += 1
        if node == 'Y':
            nY += 1
        pstr[pos] = node
        curlen += 1
        pos -= 1
        _update_general(k, m, depth, node)
        if curlen == log2n:
            k0 = int(pstr.bits_even.to01(), base=2)
            coeff = _compute_coeff(k0, k, m, nY, matrix)
            if not np.isclose(np.abs(coeff), 0):
                coeffs[pstr.copy()] += coeff
        else:
            dfs_stack.appendleft(('Z', 1 + depth))
            dfs_stack.appendleft(('Y', 1 + depth))
            dfs_stack.appendleft(('X', 1 + depth))
            dfs_stack.appendleft(('I', 1 + depth))
    return dict(coeffs)
