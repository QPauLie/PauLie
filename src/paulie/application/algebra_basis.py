"""
Defining-representation matrix bases for the four canonical DLA families.

Conventions are fixed for downstream use — basis ordering, sign choices,
and the choice of J for sp are documented per constructor.

  so(N)  real antisymmetric  {E_ij - E_ji : 1<=i<j<=N}  lex (i,j) order
         shape (N*(N-1)//2, N, N)   dtype float64

  su(N)  traceless anti-Hermitian (Gell-Mann generators x i)
         order: (a) sym off-diag i<j, (b) antisym off-diag i<j,
                (c) diagonal l=1..N-1
         shape (N**2-1, N, N)   dtype complex128

  sp(N)  N = half-dimension; matrices are 2N x 2N;
         J = [[0, I_N], [-I_N, 0]]  so  X in sp(N)  iff  X^T J + J X = 0
         order: (i)  E_ij - E_{j+N,i+N}       0<=i,j<N   (N^2 matrices)
                (ii) E_{i,j+N} + E_{j,i+N}    0<=i<=j<N  (N(N+1)/2)
                (iii)E_{i+N,j} + E_{j+N,i}    0<=i<=j<N  (N(N+1)/2)
         total dim = N(2N+1)  shape (N*(2N+1), 2N, 2N)  dtype float64

  u(1)   shape (1, 1, 1)  dtype complex128  single generator [[i]]
"""

from __future__ import annotations
import re
import numpy as np


def so_basis(N: int) -> np.ndarray:
    """Basis for so(N) in the N x N real anti-symmetric defining representation."""
    dim = N * (N - 1) // 2
    basis = np.zeros((dim, N, N), dtype=np.float64)
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            basis[k, i, j] = 1.0
            basis[k, j, i] = -1.0
            k += 1
    return basis


def su_basis(N: int) -> np.ndarray:
    """Basis for su(N): traceless anti-Hermitian matrices (Gell-Mann x i)."""
    dim = N * N - 1
    basis = np.zeros((dim, N, N), dtype=np.complex128)
    k = 0
    # (a) symmetric off-diagonal:  i*(|i><j| + |j><i|) / sqrt(2)
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), dtype=np.complex128)
            m[i, j] = m[j, i] = 1.0
            basis[k] = 1j * m / np.sqrt(2)
            k += 1
    # (b) antisymmetric off-diagonal:  (|i><j| - |j><i|) / sqrt(2)
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), dtype=np.complex128)
            m[i, j] = 1.0
            m[j, i] = -1.0
            basis[k] = m / np.sqrt(2)
            k += 1
    # (c) diagonal:  i * diag(1,...,1,-l,0,...) / sqrt(l*(l+1)/2)
    for ell in range(1, N):
        m = np.zeros((N, N), dtype=np.complex128)
        for i in range(ell):
            m[i, i] = 1.0
        m[ell, ell] = -ell
        basis[k] = 1j * m / np.sqrt(ell * (ell + 1) / 2)
        k += 1
    return basis


def sp_basis(N: int) -> np.ndarray:
    """
    Basis for sp(N) in the 2N x 2N defining representation.
    J = [[0, I_N], [-I_N, 0]];  X in sp(N)  iff  X^T J + J X = 0.
    N is the half-dimension: matrices are 2N x 2N.
    """
    size = 2 * N
    dim = N * (2 * N + 1)
    basis = np.zeros((dim, size, size), dtype=np.float64)
    k = 0
    # (i) E_ij - E_{j+N, i+N}
    for i in range(N):
        for j in range(N):
            basis[k, i, j] = 1.0
            basis[k, j + N, i + N] = -1.0
            k += 1
    # (ii) E_{i, j+N} + E_{j, i+N}  -- upper-right symmetric block
    for i in range(N):
        for j in range(i, N):
            basis[k, i, j + N] = 1.0
            basis[k, j, i + N] = 1.0
            k += 1
    # (iii) E_{i+N, j} + E_{j+N, i}  -- lower-left symmetric block
    for i in range(N):
        for j in range(i, N):
            basis[k, i + N, j] = 1.0
            basis[k, j + N, i] = 1.0
            k += 1
    return basis


def u1_basis() -> np.ndarray:
    """Basis for u(1): single generator i in the 1 x 1 defining rep."""
    return np.array([[[1j]]], dtype=np.complex128)


def algebra_basis_from_label(label: str) -> np.ndarray:
    """Return the defining-representation basis for the algebra named by *label*.

    Handles homogeneous direct sums (``"4*so(3)"``), heterogeneous direct sums
    (``"so(3)+u(1)"``), and single summands (``"sp(4)"``).  The result is a
    single ndarray with block-diagonal embedding; summand i occupies the i-th
    diagonal block.

    Parameters
    ----------
    label : str
        Label as returned by ``PauliStringCollection.get_algebra()``.

    Returns
    -------
    np.ndarray
        Shape ``(total_dim, total_M, total_M)`` where ``total_M`` is the sum
        of per-summand matrix sizes and ``total_dim`` is the sum of per-summand
        basis dimensions.

    Raises
    ------
    ValueError
        If a term cannot be parsed or uses an unknown family.
    NotImplementedError
        If ``u(N)`` with ``N > 1`` is requested.
    """

    def _primitive(family: str, N: int) -> np.ndarray:
        if family == "so":
            return so_basis(N)
        if family == "su":
            return su_basis(N)
        if family == "sp":
            return sp_basis(N)
        if family == "u":
            if N != 1:
                raise NotImplementedError(f"u({N}) is not implemented; only u(1) is supported.")
            return u1_basis()
        raise ValueError(f"Unknown Lie algebra family {family!r}")

    summands: list[np.ndarray] = []
    for term in label.split("+"):
        term = term.strip()
        m = re.match(r"^(so|sp|su|u)\((\d+)\)$", term)
        if m:
            summands.append(_primitive(m.group(1), int(m.group(2))))
            continue
        m = re.match(r"^(\d+)\*(so|sp|su|u)\((\d+)\)$", term)
        if m:
            k_rep, fam, N = int(m.group(1)), m.group(2), int(m.group(3))
            base = _primitive(fam, N)
            summands.extend([base] * k_rep)
            continue
        raise ValueError(f"Cannot parse term {term!r} in label {label!r}.")

    if not summands:
        raise ValueError(f"Empty label: {label!r}")

    if len(summands) == 1:
        return summands[0].copy()

    # Block-diagonal embedding — supports heterogeneous summand sizes.
    Ms = [b.shape[1] for b in summands]
    total_M = sum(Ms)
    dtype = np.result_type(*[b.dtype for b in summands])
    mats: list[np.ndarray] = []
    offset = 0
    for b, M in zip(summands, Ms):
        for mat in b:
            block = np.zeros((total_M, total_M), dtype=dtype)
            block[offset : offset + M, offset : offset + M] = mat
            mats.append(block)
        offset += M
    return np.array(mats, dtype=dtype)
