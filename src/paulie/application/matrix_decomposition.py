"""
Decompose a complex matrix into the Pauli basis using 
PauliString bitarray representation.
"""

import numpy as np
from concurrent.futures import ThreadPoolExecutor
from paulie.common.pauli_string_bitarray import PauliString

def matrix_decomposition(matrix: np.ndarray, tol: float = 1e-12) -> dict[PauliString, complex]:
    """
    Decompose a 2^n × 2^n complex matrix A into the Pauli basis using the
    PTDR algorithm from Koska et al. (arXiv:2403.11644). Carries along index‐map k,
    sign‐map m arrays, and a Y‐count, doing one O(2^n) dot‐product per leaf:
        α_P = (−i)^{n_Y} · (1/2^n) · ∑_j m[j] · A[k[j], j].

    Args:
        matrix: square array of shape (2^n, 2^n)
        tol: drop any |α| ≤ tol

    Returns:
        Mapping PauliString → its complex coefficient α.
    """
    # — validate —
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    dim = matrix.shape[0]
    n = int(np.log2(dim))
    if 2**n != dim:
        raise ValueError("matrix dimension must be a power of 2")

    # precompute bit‐masks for each qubit position
    masks = [1 << (n - 1 - i) for i in range(n)]
    initial_k = np.arange(dim, dtype=np.int64)
    initial_m = np.ones(dim, dtype=np.int64)
    norm = 1 / (2**n)

    result: dict[PauliString, complex] = {}

    def recurse(
        depth: int,
        k: np.ndarray,
        m: np.ndarray,
        y_count: int,
        labels: list[str],
        acc: dict[PauliString, complex]
    ) -> None:
        if depth == n:
            # leaf: compute dot‐product ∑ m[j] A[k[j], j]
            idx = np.arange(dim)
            total = (m * matrix[k, idx]).sum()
            alpha = ((-1j) ** y_count) * norm * total
            if abs(alpha) > tol:
                ps = PauliString(pauli_str=''.join(labels))
                acc[ps] = alpha
            return

        mask = masks[depth]
        # I branch
        labels[depth] = 'I'
        recurse(depth + 1, k, m, y_count, labels, acc)

        # X branch
        labels[depth] = 'X'
        k ^= mask
        recurse(depth + 1, k, m, y_count, labels, acc)
        k ^= mask  # revert

        # Z branch
        labels[depth] = 'Z'
        bit = (k & mask) != 0
        m[bit] *= -1
        recurse(depth + 1, k, m, y_count, labels, acc)
        m[bit] *= -1  # revert

        # Y branch: apply Z then X to match Y = i·Z·X
        labels[depth] = 'Y'
        bit = (k & mask) != 0
        m[bit] *= -1
        k ^= mask
        recurse(depth + 1, k, m, y_count + 1, labels, acc)
        k ^= mask     # revert
        m[bit] *= -1  # revert

    def process_label(lbl0: str) -> dict[PauliString, complex]:
        k = initial_k.copy()
        m = initial_m.copy()
        y_c = 0
        labels = [''] * n
        mask0 = masks[0]

        if lbl0 == 'I':
            labels[0] = 'I'
        elif lbl0 == 'X':
            labels[0] = 'X'
            k ^= mask0
        elif lbl0 == 'Z':
            labels[0] = 'Z'
            bit0 = (k & mask0) != 0
            m[bit0] *= -1
        else:  # 'Y'
            labels[0] = 'Y'
            bit0 = (k & mask0) != 0
            m[bit0] *= -1
            k ^= mask0
            y_c = 1

        acc: dict[PauliString, complex] = {}
        recurse(1, k, m, y_c, labels, acc)
        return acc

    # parallel over I, X, Y, Z on qubit 0
    with ThreadPoolExecutor() as executor:
        for branch_map in executor.map(process_label, ['I', 'X', 'Y', 'Z']):
            result.update(branch_map)

    return result
