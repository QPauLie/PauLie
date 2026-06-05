"""Tests for :meth:`get_algebra_basis` -- the defining-representation matrix
basis of the classified dynamical Lie algebra.

The suite covers one input per canonical type of Theorem 1 in Aguilar et al.
2024 (A / B1 / B2 / B3) plus direct-sum inputs, and proves that the returned set
is *actually a basis* of the named algebra -- not merely a set of matrices of the
right shape, and not merely a Lie algebra of the right dimension.

For each summand of family ``F`` and defining-rep size ``N`` we verify:

1. **Membership** -- every matrix satisfies the defining condition of ``F`` (real
   antisymmetric / traceless anti-Hermitian / ``X^T J + J X = 0`` / anti-Hermitian),
   so the set lies *inside* ``F``.
2. **Independence** -- the matrices are linearly independent over the reals.
3. **Completeness** -- the count equals the closed-form ``dim F``, and an
   independently generated, generic element of ``F`` lies in the real span of the
   basis.  (1)+(2)+(3) prove the set spans ``F`` exactly, i.e. *every* element of
   the algebra is reachable -- so it is a basis of ``F``.
4. **Closure** -- the span is closed under the Lie bracket.

Finally, the family is tied back to the *generated* algebra: ``dim F`` must equal
the size of the brute-force Pauli Lie closure of the generators (the closure used
in ``docs/examples/benchmark_classification.py``), computed independently of
``get_dla_dim``.  Matching dimension alone does not identify the algebra, which is
why it is only the *last* link after completeness has been established.
"""
from itertools import product

import numpy as np
import pytest

from paulie import G_LIE
from paulie import get_pauli_string as p

TOL = 1e-9


def _symplectic_J(rank: int) -> np.ndarray:
    """The symplectic form J = [[0, I], [-I, 0]], built independently of the
    module under test so the sp checks do not rely on its own ``symplectic_form``.
    """
    eye = np.eye(rank, dtype=np.complex128)
    zero = np.zeros((rank, rank), dtype=np.complex128)
    return np.block([[zero, eye], [-eye, zero]])


# --------------------------------------------------------------------------- #
# Linear algebra over the reals (a real Lie algebra is a real vector space)
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


def _in_real_span(stack: np.ndarray, mat: np.ndarray) -> bool:
    """Whether ``mat`` lies in the real span of the basis ``stack``."""
    base_dim = _dim_real_span(stack)
    augmented = np.concatenate([stack, mat[np.newaxis, :, :]], axis=0)
    return _dim_real_span(augmented) == base_dim


# --------------------------------------------------------------------------- #
# Independent constructions of each family (used to certify a *complete* span)
# --------------------------------------------------------------------------- #
def _family_dimension(family: str, size: int) -> int:
    """Closed-form dimension of the family in defining-rep size ``size``."""
    if family == "so":
        return size * (size - 1) // 2
    if family == "su":
        return size * size - 1
    if family == "u":
        return size * size
    if family == "sp":
        rank = size // 2
        return rank * (2 * rank + 1)
    raise ValueError(f"unknown family {family!r}")


def _satisfies_defining_condition(family: str, mat: np.ndarray) -> bool:
    """Whether ``mat`` satisfies the defining condition of ``family``."""
    size = mat.shape[0]
    if family == "so":
        return bool(np.allclose(mat.imag, 0.0, atol=TOL)
                    and np.allclose(mat + mat.T, 0.0, atol=TOL))
    if family == "u":
        return bool(np.allclose(mat + mat.conj().T, 0.0, atol=TOL))
    if family == "su":
        return bool(np.allclose(mat + mat.conj().T, 0.0, atol=TOL)
                    and abs(np.trace(mat)) < TOL)
    if family == "sp":
        form = _symplectic_J(size // 2)
        return bool(np.allclose(mat.imag, 0.0, atol=TOL)
                    and np.allclose(mat.T @ form + form @ mat, 0.0, atol=TOL))
    raise ValueError(f"unknown family {family!r}")


def _random_family_element(family: str, size: int, rng: np.random.Generator) -> np.ndarray:
    """A generic element of ``family``, built independently of the basis code."""
    if family == "so":
        mat = rng.standard_normal((size, size))
        return (mat - mat.T).astype(np.complex128)
    if family in ("u", "su"):
        mat = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))
        anti = mat - mat.conj().T
        if family == "su":
            anti -= np.trace(anti) / size * np.eye(size)
        return anti
    if family == "sp":
        rank = size // 2
        sym = rng.standard_normal((size, size))
        sym = sym + sym.T
        # X = -J S with S symmetric satisfies X^T J + J X = 0.
        return (-_symplectic_J(rank) @ sym).astype(np.complex128)
    raise ValueError(f"unknown family {family!r}")


def _assert_is_basis_of_family(stack: np.ndarray, family: str) -> None:
    """Assert that ``stack`` is a basis of the named family in its defining rep.

    Membership + independence + (count == dim F) already imply the span is *all*
    of F.  We additionally show, directly, that an independently constructed
    generic element of F is reachable -- the "reach every element" check.
    """
    size = stack.shape[1]
    # (1) every element lies inside the family algebra
    for mat in stack:
        assert _satisfies_defining_condition(family, mat)
    # (2) linearly independent and (3) of the full family dimension
    assert _is_independent(stack)
    assert stack.shape[0] == _family_dimension(family, size)
    # (3, direct) a generic element of the family is in the span -> span = F
    rng = np.random.default_rng(0)
    for _ in range(5):
        element = _random_family_element(family, size, rng)
        assert _satisfies_defining_condition(family, element)
        assert _in_real_span(stack, element)
    # (4) the span is a Lie subalgebra
    assert _is_closed_under_bracket(stack)


# --------------------------------------------------------------------------- #
# Independent ground truth: brute-force Pauli Lie closure
# --------------------------------------------------------------------------- #
def _lie_closure_dim(collection) -> int:
    """Dimension of the DLA via brute-force commutator enumeration.

    Mirrors ``lie_closure_brute_force`` from the benchmark example but is
    self-contained and deterministic.  Computed *without* ``get_dla_dim``.
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
# brute-force closure to terminate quickly.  ``family`` / ``size`` describe each
# summand of the resulting algebra.
#   A -> so(N)   B1 -> sp(2^k)   B2 -> so(2^k)   B3 -> su(2^k)
CANONICAL_CASES = [
    ("A  so(3)", ["XY"], 3, "so(3)", "so", 3),
    ("A  2*so(3)", ["XY"], 4, "2*so(3)", "so", 3),
    ("B1 sp(4)", ["XY", "XZ"], 4, "sp(4)", "sp", 8),
    ("B2 so(16)", G_LIE["a11"], 4, "so(16)", "so", 16),
    ("B3 su(16)", G_LIE["a12"], 4, "su(16)", "su", 16),
    ("B3 2*su(8)", ["XX", "YY", "ZZ", "ZY"], 4, "2*su(8)", "su", 8),
    ("u(1)", ["XY"], 2, "u(1)", "u", 1),
]


# --------------------------------------------------------------------------- #
# The returned set is actually a basis of the named algebra
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("label, gens, n, expected, family, size", CANONICAL_CASES)
def test_each_summand_is_a_basis_of_its_family(label, gens, n, expected, family, size):
    """Every summand spans -- exactly -- the family algebra named by the label."""
    coll = p(gens, n=n)
    assert coll.get_algebra() == expected
    summands = coll.get_algebra_basis()
    for basis in summands:
        assert basis.shape[1] == size and basis.shape[2] == size
        _assert_is_basis_of_family(basis, family)


@pytest.mark.parametrize("label, gens, n, expected, family, size", CANONICAL_CASES)
def test_basis_dimension_matches_brute_force_closure(label, gens, n, expected, family, size):
    """The (verified) basis has exactly as many elements as the DLA dimension.

    The basis is already shown to span its family by
    :func:`test_each_summand_is_a_basis_of_its_family`; here that family is tied
    to the *generated* algebra by matching its dimension against the brute-force
    Pauli Lie closure -- an independent ground truth.
    """
    coll = p(gens, n=n)
    total_dim = sum(basis.shape[0] for basis in coll.get_algebra_basis())
    assert total_dim == _lie_closure_dim(coll)


# --------------------------------------------------------------------------- #
# Per-canonical-type shape checks (readable documentation of the four cases)
# --------------------------------------------------------------------------- #
def test_so_summand_type_a():
    """Type A: real antisymmetric so(N)."""
    summands = p(["XY"], n=3).get_algebra_basis()
    assert len(summands) == 1
    assert summands[0].shape == (3, 3, 3)
    _assert_is_basis_of_family(summands[0], "so")


def test_sp_summand_type_b1():
    """Type B1: real symplectic sp(n), X^T J + J X = 0."""
    coll = p(["XY", "XZ"], n=4)
    assert coll.get_algebra() == "sp(4)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    assert summands[0].shape == (36, 8, 8)
    _assert_is_basis_of_family(summands[0], "sp")


def test_so_summand_type_b2():
    """Type B2: so(2^k); a11 at n=4 gives so(16)."""
    coll = p(G_LIE["a11"], n=4)
    assert coll.get_algebra() == "so(16)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    assert summands[0].shape == (120, 16, 16)
    _assert_is_basis_of_family(summands[0], "so")


def test_su_summand_type_b3():
    """Type B3: traceless anti-Hermitian su(2^k); a12 at n=4 gives su(16)."""
    coll = p(G_LIE["a12"], n=4)
    assert coll.get_algebra() == "su(16)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 1
    assert summands[0].shape == (255, 16, 16)
    _assert_is_basis_of_family(summands[0], "su")


def test_u1_summand():
    """A single Pauli string generates u(1) = span{i I}."""
    summands = p(["XY"], n=2).get_algebra_basis()
    assert len(summands) == 1
    assert summands[0].shape == (1, 1, 1)
    _assert_is_basis_of_family(summands[0], "u")


def test_direct_sum_has_one_entry_per_summand():
    """B3 with n_1 >= 1: a direct sum 2*su(8) returns two su(8) summands."""
    coll = p(["XX", "YY", "ZZ", "ZY"], n=4)
    assert coll.get_algebra() == "2*su(8)"
    summands = coll.get_algebra_basis()
    assert len(summands) == 2
    for basis in summands:
        assert basis.shape == (63, 8, 8)
        _assert_is_basis_of_family(basis, "su")


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
# Bookkeeping: summand count, total dimension, ordering, copies
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("label, gens, n, expected, family, size", CANONICAL_CASES)
def test_total_dimension_matches_get_dla_dim(label, gens, n, expected, family, size):
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
    first = p(["XX", "YY", "ZZ", "ZY"], n=4).get_algebra_basis()
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


def test_symplectic_form_matches_documented_convention():
    """The module's ``symplectic_form`` equals the documented J = [[0, I], [-I, 0]]."""
    from paulie.application.algebra_basis import symplectic_form
    for rank in (1, 2, 4):
        assert np.allclose(symplectic_form(rank), _symplectic_J(rank), atol=TOL)
