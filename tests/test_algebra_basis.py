"""Tests for :meth:`get_algebra_basis` -- the defining-representation matrix
basis of the classified dynamical Lie algebra.

The suite covers one input per canonical type of Theorem 1 in Aguilar et al.
2024 (A / B1 / B2 / B3) plus direct-sum inputs with more than one summand, and
checks three things for every returned basis:

1. the matrices satisfy the *defining condition* of their family,
2. the matrices are linearly independent and close under the Lie bracket -- so
   their real span is a Lie subalgebra of the expected dimension,
3. that dimension agrees with an **independent** ground truth: the size of the
   brute-force Pauli Lie closure of the generators (the same closure used in
   ``docs/examples/benchmark_classification.py``).

Together (2) and (3) verify that the output is *actually a basis* of the
classified algebra, not merely a collection of matrices of the right shape.
"""
from itertools import product

import numpy as np
import pytest

from paulie import G_LIE
from paulie import get_pauli_string as p
from paulie.application.algebra_basis import symplectic_form

TOL = 1e-9


# --------------------------------------------------------------------------- #
# Linear-algebra helpers (a real Lie algebra is a real vector space)
# --------------------------------------------------------------------------- #
def _real_rows(stack: np.ndarray) -> np.ndarray:
    """Stack the real coordinates (Re and Im parts) of every matrix as rows."""
    flat = stack.reshape(stack.shape[0], -1)
    return np.hstack([flat.real, flat.imag])


def _dim_real_span(stack: np.ndarray) -> int:
    """Dimension of the real span of a stack of matrices."""
    if stack.shape[0] == 0:
        return 0
    return int(np.linalg.matrix_rank(_real_rows(stack), tol=TOL))


def _is_independent(stack: np.ndarray) -> bool:
    """The matrices are linearly independent over the reals."""
    return _dim_real_span(stack) == stack.shape[0]


def _is_closed_under_bracket(stack: np.ndarray) -> bool:
    """Every commutator of basis elements lies in the real span of the basis."""
    base_dim = _dim_real_span(stack)
    brackets = [
        stack[a] @ stack[b] - stack[b] @ stack[a]
        for a in range(stack.shape[0])
        for b in range(a + 1, stack.shape[0])
    ]
    if not brackets:
        return True
    augmented = np.concatenate([stack, np.stack(brackets, axis=0)], axis=0)
    return _dim_real_span(augmented) == base_dim


# --------------------------------------------------------------------------- #
# Independent ground truth: brute-force Pauli Lie closure
# --------------------------------------------------------------------------- #
def _lie_closure_dim(collection) -> int:
    """Dimension of the DLA via brute-force commutator enumeration.

    Mirrors ``lie_closure_brute_force`` from the benchmark example but is
    self-contained and deterministic.  The returned count is computed *without*
    consulting ``get_dla_dim``, so it is an independent check.
    """
    dla = list(collection.get())
    old_len, new_len = 0, len(dla)
    initial_len = new_len
    while new_len > old_len:
        for pw1, pw2 in product(dla[:initial_len], dla[old_len:]):
            if pw1.commutes_with(pw2):
                continue
            com = pw1 @ pw2
            if com not in dla:
                dla.append(com)
        old_len, new_len = new_len, len(dla)
    return len(dla)


# A representative generating set for each canonical type, small enough for the
# brute-force closure to terminate quickly.
#   A   -> so(N)        B1 -> sp(2^k)
#   B2  -> so(2^k)      B3 -> su(2^k)
CANONICAL_CASES = [
    ("A  so(3)", ["XY"], 3, "so(3)"),
    ("A  2*so(3)", ["XY"], 4, "2*so(3)"),
    ("B1 sp(4)", ["XY", "XZ"], 4, "sp(4)"),
    ("B2 so(16)", G_LIE["a11"], 4, "so(16)"),
    ("B3 su(16)", G_LIE["a12"], 4, "su(16)"),
    ("B3 2*su(8)", ["XX", "YY", "ZZ", "ZY"], 4, "2*su(8)"),
    ("u(1)", ["XY"], 2, "u(1)"),
]


# --------------------------------------------------------------------------- #
# Per-family defining-condition tests
# --------------------------------------------------------------------------- #
def test_so_summand_type_a():
    """Type A: real antisymmetric so(N)."""
    summands = p(["XY"], n=3).get_algebra_basis()
    assert len(summands) == 1
    basis = summands[0]
    assert basis.shape == (3, 3, 3)
    for mat in basis:
        assert np.allclose(mat.imag, 0.0, atol=TOL)
        assert np.allclose(mat + mat.T, 0.0, atol=TOL)
    assert _is_independent(basis)
    assert _is_closed_under_bracket(basis)


def test_sp_summand_type_b1():
    """Type B1: real symplectic sp(n), X^T J + J X = 0."""
    coll = p(["XY", "XZ"], n=4)
    assert coll.get_algebra() == "sp(4)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    basis = summands[0]
    assert basis.shape == (36, 8, 8)
    form = symplectic_form(4)
    for mat in basis:
        assert np.allclose(mat.imag, 0.0, atol=TOL)
        assert np.allclose(mat.T @ form + form @ mat, 0.0, atol=TOL)
    assert _is_independent(basis)
    assert _is_closed_under_bracket(basis)


def test_so_summand_type_b2():
    """Type B2: so(2^k); a11 at n=4 gives so(16)."""
    coll = p(G_LIE["a11"], n=4)
    assert coll.get_algebra() == "so(16)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    basis = summands[0]
    assert basis.shape == (120, 16, 16)
    for mat in basis:
        assert np.allclose(mat + mat.T, 0.0, atol=TOL)
    assert _is_independent(basis)
    assert _is_closed_under_bracket(basis)


def test_su_summand_type_b3():
    """Type B3: traceless anti-Hermitian su(2^k); a12 at n=4 gives su(16)."""
    coll = p(G_LIE["a12"], n=4)
    assert coll.get_algebra() == "su(16)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    basis = summands[0]
    assert basis.shape == (255, 16, 16)
    for mat in basis:
        assert np.allclose(mat + mat.conj().T, 0.0, atol=TOL)
        assert abs(np.trace(mat)) < TOL
    assert _is_independent(basis)
    assert _is_closed_under_bracket(basis)


def test_u1_summand():
    """A single Pauli string generates u(1) = span{i I}."""
    summands = p(["XY"], n=2).get_algebra_basis()
    assert len(summands) == 1
    basis = summands[0]
    assert basis.shape == (1, 1, 1)
    assert np.allclose(basis[0] + basis[0].conj().T, 0.0, atol=TOL)


def test_direct_sum_has_one_entry_per_summand():
    """B3 with n_1 >= 1: a direct sum 2*su(8) returns two su(8) summands."""
    coll = p(["XX", "YY", "ZZ", "ZY"], n=4)
    assert coll.get_algebra() == "2*su(8)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 2
    for basis in summands:
        assert basis.shape == (63, 8, 8)
        for mat in basis:
            assert np.allclose(mat + mat.conj().T, 0.0, atol=TOL)
            assert abs(np.trace(mat)) < TOL
        assert _is_independent(basis)
        assert _is_closed_under_bracket(basis)


# --------------------------------------------------------------------------- #
# Cross-summand brackets vanish (a direct sum acts on orthogonal blocks)
# --------------------------------------------------------------------------- #
def test_cross_summand_brackets_vanish():
    """Block-embedding two summands yields commuting blocks: [g_1, g_2] = 0."""
    summands = p(["XX", "YY", "ZZ", "ZY"], n=4).get_algebra_basis()
    assert len(summands) == 2
    first, second = summands[0], summands[1]
    size = first.shape[1]
    top = np.zeros((2 * size, 2 * size), dtype=np.complex128)
    top[:size, :size] = first[0]
    bottom = np.zeros((2 * size, 2 * size), dtype=np.complex128)
    bottom[size:, size:] = second[1]
    assert np.allclose(top @ bottom - bottom @ top, 0.0, atol=TOL)


# --------------------------------------------------------------------------- #
# The output is *actually a basis* of the classified algebra
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("label, gens, n, expected", CANONICAL_CASES)
def test_basis_spans_the_lie_closure(label, gens, n, expected):
    """Independent verification that the returned set is a genuine basis.

    The basis is genuine iff (a) its matrices are linearly independent, (b) their
    span closes under the bracket (so it is a Lie subalgebra), and (c) the total
    dimension equals the size of the brute-force Pauli Lie closure -- a ground
    truth computed without ``get_dla_dim``.
    """
    coll = p(gens, n=n)
    assert coll.get_algebra() == expected

    summands = coll.get_algebra_basis()
    total_dim = 0
    for basis in summands:
        assert _is_independent(basis)
        assert _is_closed_under_bracket(basis)
        total_dim += basis.shape[0]

    assert total_dim == _lie_closure_dim(coll)


@pytest.mark.parametrize("label, gens, n, expected", CANONICAL_CASES)
def test_total_dimension_matches_get_dla_dim(label, gens, n, expected):
    """The materialised basis has exactly ``get_dla_dim`` matrices in total."""
    coll = p(gens, n=n)
    total = sum(basis.shape[0] for basis in coll.get_algebra_basis())
    assert total == coll.get_dla_dim()


def test_summand_count_matches_label():
    """``len(basis)`` equals the number of direct summands in the label."""
    assert len(p(["XY"], n=3).get_algebra_basis()) == 1
    assert len(p(["XY"], n=4).get_algebra_basis()) == 2  # 2*so(3)
    assert len(p(["XX", "YY", "ZZ", "ZY"], n=4).get_algebra_basis()) == 2  # 2*su(8)
    assert len(p(G_LIE["a0"], n=4).get_algebra_basis()) == 3  # 3*u(1)


def test_ordering_is_deterministic():
    """Repeated calls return summands in the same order and values."""
    coll = p(["XX", "YY", "ZZ", "ZY"], n=4)
    first = coll.get_algebra_basis()
    second = p(["XX", "YY", "ZZ", "ZY"], n=4).get_algebra_basis()
    assert len(first) == len(second)
    for left, right in zip(first, second):
        assert left.shape == right.shape
        assert np.allclose(left, right, atol=TOL)


def test_returned_arrays_are_independent_copies():
    """Mutating a returned summand must not affect a fresh call."""
    basis = p(["XY"], n=3).get_algebra_basis()[0]
    basis[0] += 1.0
    fresh = p(["XY"], n=3).get_algebra_basis()[0]
    assert not np.allclose(basis[0], fresh[0], atol=TOL)
