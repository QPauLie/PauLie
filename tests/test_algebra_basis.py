"""
Tests for the matrix-basis form of the classified Lie algebra.

Covers:
- The four single-algebra primitive bases (so, sp, su, u): defining
  property, dimension, real/complex dtype, bracket closure.
- The Classification / PauliStringCollection ``get_algebra_basis`` entry
  points: one canonical-type input per family plus one direct-sum input.
"""

from __future__ import annotations

import numpy as np
import pytest

from paulie import G_LIE, get_pauli_string as p
from paulie.application.algebra_basis import (
    so_basis,
    sp_basis,
    su_basis,
    u_basis,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
TOL = 1e-10


def _is_antisym_real(b: np.ndarray) -> bool:
    return bool(np.allclose(b + b.swapaxes(1, 2), 0, atol=TOL))


def _is_anti_hermitian(b: np.ndarray) -> bool:
    return bool(np.allclose(b + b.conj().swapaxes(1, 2), 0, atol=TOL))


def _is_traceless(b: np.ndarray) -> bool:
    return bool(np.allclose(np.trace(b, axis1=1, axis2=2), 0, atol=TOL))


def _is_symplectic_generator(b: np.ndarray) -> bool:
    """Check M^T J + J M = 0 for all generators in the stack."""
    n2 = b.shape[1]
    assert n2 % 2 == 0
    n = n2 // 2
    j_mat = np.zeros((n2, n2), dtype=b.dtype)
    j_mat[:n, n:] = np.eye(n)
    j_mat[n:, :n] = -np.eye(n)
    return bool(np.allclose(b.swapaxes(1, 2) @ j_mat + j_mat @ b, 0, atol=TOL))


def _linearly_independent(b: np.ndarray) -> bool:
    flat = b.reshape(b.shape[0], -1)
    return int(np.linalg.matrix_rank(flat)) == b.shape[0]


def _brackets_in_span(b: np.ndarray) -> bool:
    """Every pairwise commutator lies in span(b)."""
    flat = b.reshape(b.shape[0], -1)
    dim = b.shape[0]
    for i in range(dim):
        for j in range(i + 1, dim):
            comm = (b[i] @ b[j] - b[j] @ b[i]).reshape(-1)
            stacked = np.concatenate([flat, comm[np.newaxis]], axis=0)
            if np.linalg.matrix_rank(stacked) != dim:
                return False
    return True


# ---------------------------------------------------------------------------
# Primitive bases
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_so_basis_defining_property(n: int) -> None:
    """so(n) generators are real antisymmetric and have dim n(n-1)/2."""
    b = so_basis(n)
    assert b.shape == (n * (n - 1) // 2, n, n)
    assert b.dtype == np.float64
    assert _is_antisym_real(b)
    assert _linearly_independent(b)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_sp_basis_defining_property(n: int) -> None:
    """sp(n) generators (rank n, 2n*2n matrices) satisfy M^T J + J M = 0
    and have dim n(2n+1)."""
    b = sp_basis(n)
    assert b.shape == (n * (2 * n + 1), 2 * n, 2 * n)
    assert b.dtype == np.float64
    assert _is_symplectic_generator(b)
    assert _linearly_independent(b)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_su_basis_defining_property(n: int) -> None:
    """su(n) generators are traceless anti-Hermitian, dim n^2 - 1."""
    b = su_basis(n)
    assert b.shape == (n * n - 1, n, n)
    assert b.dtype == np.complex128
    assert _is_anti_hermitian(b)
    assert _is_traceless(b)
    assert _linearly_independent(b)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_u_basis_defining_property(n: int) -> None:
    """u(n) generators are anti-Hermitian (NOT traceless), dim n^2."""
    b = u_basis(n)
    assert b.shape == (n * n, n, n)
    assert b.dtype == np.complex128
    assert _is_anti_hermitian(b)
    assert _linearly_independent(b)


@pytest.mark.parametrize(
    "ctor,n",
    [
        (so_basis, 4),
        (sp_basis, 2),
        (su_basis, 3),
        (u_basis, 2),
    ],
)
def test_brackets_close(ctor, n: int) -> None:
    """For each primitive basis, every commutator lies in the span."""
    assert _brackets_in_span(ctor(n))


def test_su_1_is_trivial() -> None:
    """su(1) is the zero algebra."""
    assert su_basis(1).shape == (0, 1, 1)


def test_invalid_n_rejected() -> None:
    """All constructors reject n < 1."""
    with pytest.raises(ValueError):
        so_basis(0)
    with pytest.raises(ValueError):
        sp_basis(0)
    with pytest.raises(ValueError):
        su_basis(0)
    with pytest.raises(ValueError):
        u_basis(0)


# ---------------------------------------------------------------------------
# get_algebra_basis end-to-end on classified inputs
# ---------------------------------------------------------------------------
def test_get_algebra_basis_so() -> None:
    """A canonical so(3) input — one summand of dim 3, 3x3 antisymmetric."""
    collection = p(["XY"], n=3)
    assert collection.get_algebra() == "so(3)"
    basis = collection.get_algebra_basis()
    assert len(basis) == 1
    assert basis[0].shape == (3, 3, 3)
    assert _is_antisym_real(basis[0])


def test_get_algebra_basis_sp() -> None:
    """A canonical sp(4) input — one summand of dim 4*(2*4+1) = 36, 8x8."""
    collection = p(["XY", "XZ"], n=4)
    assert collection.get_algebra() == "sp(4)"
    basis = collection.get_algebra_basis()
    assert len(basis) == 1
    assert basis[0].shape == (36, 8, 8)
    assert _is_symplectic_generator(basis[0])


def test_get_algebra_basis_u_singleton() -> None:
    """A single Pauli string generates u(1)."""
    collection = p(["XY"])
    assert collection.get_algebra() == "u(1)"
    basis = collection.get_algebra_basis()
    assert len(basis) == 1
    assert basis[0].shape == (1, 1, 1)
    assert _is_anti_hermitian(basis[0])


def test_get_algebra_basis_direct_sum() -> None:
    """A canonical direct-sum input — paulie reports '4*so(5)', so we expect
    4 summands each of shape (10, 5, 5)."""
    collection = p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"])
    assert collection.get_algebra() == "4*so(5)"
    basis = collection.get_algebra_basis()
    assert len(basis) == 4
    for summand in basis:
        assert summand.shape == (10, 5, 5)
        assert _is_antisym_real(summand)


def test_get_algebra_basis_brackets_close_per_summand() -> None:
    """Bracket closure spot-check on a direct-sum input: each summand is a
    Lie algebra in its own right, so commutators must stay in the summand."""
    basis = p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]).get_algebra_basis()
    # Spot-check the first summand to keep the test fast.
    assert _brackets_in_span(basis[0])


def test_get_algebra_basis_su() -> None:
    """A canonical su(16) input from the G_LIE tabulated generators.

    a12 at n=4 produces a single su(16) summand. Its dim is
    :math:`16^2 - 1 = 255`, with 16x16 anti-Hermitian traceless matrices.
    """
    collection = p(G_LIE["a12"], n=4)
    assert collection.get_algebra() == "su(16)"
    basis = collection.get_algebra_basis()
    assert len(basis) == 1
    assert basis[0].shape == (255, 16, 16)
    assert _is_anti_hermitian(basis[0])
    assert _is_traceless(basis[0])
