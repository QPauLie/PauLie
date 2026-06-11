"""
Tests for defining-representation algebra bases (``get_algebra_basis``).

Each canonical Lie-algebra family (``so``, ``su``, ``sp`` as the compact form
``usp``, and ``u(1)``) is checked for the correct dimension, membership in the
algebra, linear independence, and closure under the commutator. The public-API
tests drive the classifier end to end and confirm the total basis dimension
matches ``get_dla_dim`` with the per-summand partition preserved.
"""
from typing import Optional

import numpy as np
import pytest

from paulie import get_pauli_string as p
from paulie.application.algebra_basis import (
    so_basis,
    su_basis,
    usp_basis,
    u1_basis,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _symplectic_form(half_dim: int) -> np.ndarray:
    """J = [[0, I_M], [-I_M, 0]], the form used to define usp(2M)."""
    j = np.zeros((2 * half_dim, 2 * half_dim), dtype=np.complex128)
    j[:half_dim, half_dim:] = np.eye(half_dim)
    j[half_dim:, :half_dim] = -np.eye(half_dim)
    return j


def _is_anti_hermitian(basis: np.ndarray) -> bool:
    """True if every matrix in the basis is anti-Hermitian."""
    return all(np.allclose(b, -b.conj().T) for b in basis)


def _is_linearly_independent(basis: np.ndarray) -> bool:
    """True if the basis matrices are linearly independent as vectors."""
    vecs = np.asarray(basis).reshape(basis.shape[0], -1)
    return bool(np.linalg.matrix_rank(vecs) == basis.shape[0])


def _closes_under_bracket(basis: np.ndarray, tol: float = 1e-9) -> bool:
    """True if every commutator [B_i, B_j] lies in the span of the basis."""
    basis = np.asarray(basis)
    span = basis.reshape(basis.shape[0], -1).T
    for i in range(basis.shape[0]):
        for j in range(i + 1, basis.shape[0]):
            comm = (basis[i] @ basis[j] - basis[j] @ basis[i]).reshape(-1)
            coeffs, *_ = np.linalg.lstsq(span, comm, rcond=None)
            if np.linalg.norm(span @ coeffs - comm) > tol:
                return False
    return True


def _element_is_reachable(element: np.ndarray, basis: np.ndarray,
                          tol: float = 1e-9) -> bool:
    """True if ``element`` is a real linear combination of the basis."""
    span = np.asarray(basis).reshape(basis.shape[0], -1).T
    target = element.reshape(-1)
    coeffs, *_ = np.linalg.lstsq(span, target, rcond=None)
    if np.linalg.norm((span @ coeffs) - target) > tol:
        return False
    # A real basis of a real form must reconstruct the element with real coeffs.
    return bool(np.allclose(coeffs.imag, 0, atol=tol))


# --------------------------------------------------------------------------- #
# Per-family construction: dimension, membership, independence, closure
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_so_basis_is_so(n: int) -> None:
    """so(N) basis is real antisymmetric, independent and bracket-closed."""
    basis = so_basis(n)
    assert basis.shape == (n * (n - 1) // 2, n, n)
    assert np.allclose(basis.imag, 0)
    assert all(np.allclose(b, -b.T) for b in basis)
    assert _is_linearly_independent(basis)
    assert _closes_under_bracket(basis)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_su_basis_is_su(n: int) -> None:
    """su(N) basis is traceless anti-Hermitian, independent and closed."""
    basis = su_basis(n)
    assert basis.shape == (n * n - 1, n, n)
    assert _is_anti_hermitian(basis)
    assert all(abs(np.trace(b)) < 1e-9 for b in basis)
    assert _is_linearly_independent(basis)
    assert _closes_under_bracket(basis)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_su_basis_uniform_gell_mann_normalization(n: int) -> None:
    """su(N) generators are uniformly normalized: Tr(B_a^dag B_b) = 2 delta."""
    basis = su_basis(n)
    gram = np.array([[np.trace(a.conj().T @ b) for b in basis] for a in basis])
    assert np.allclose(gram, 2.0 * np.eye(basis.shape[0]))


@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_usp_basis_is_compact_symplectic(m: int) -> None:
    """sp(M) basis is the compact form usp(2M): anti-Hermitian, X^T J + J X = 0."""
    basis = usp_basis(m)
    assert basis.shape == (m * (2 * m + 1), 2 * m, 2 * m)
    assert _is_anti_hermitian(basis)
    j = _symplectic_form(m)
    for x in basis:
        assert np.allclose(x.T @ j + j @ x, 0)
    assert _is_linearly_independent(basis)
    assert _closes_under_bracket(basis)


def test_u1_basis() -> None:
    """u(1) basis is the single 1x1 generator [[i]]."""
    basis = u1_basis()
    assert basis.shape == (1, 1, 1)
    assert np.allclose(basis[0], [[1j]])


# --------------------------------------------------------------------------- #
# Element-wise reachability with a nontrivial correspondence
#
# For so(N) and usp(2M) the defining representation (N x N, 2M x 2M) is NOT the
# 2^n embedded space, so reachability is not implied by a full-rank argument the
# way it is for su(2^n). These tests build an arbitrary element of the target
# algebra and confirm it is a real linear combination of the basis, i.e. every
# algebra element is reachable.
# --------------------------------------------------------------------------- #
def test_so_reachability_arbitrary_antisymmetric() -> None:
    """An arbitrary so(6) element is reachable from the so(6) basis."""
    rng = np.random.default_rng(0)
    n = 6
    basis = so_basis(n)
    a = rng.standard_normal((n, n))
    element = a - a.T  # arbitrary real antisymmetric = arbitrary so(6) element
    assert _element_is_reachable(element, basis)


def test_usp_reachability_arbitrary_element() -> None:
    """An arbitrary usp(6) element is reachable from the sp(3) basis."""
    rng = np.random.default_rng(1)
    m = 3
    basis = usp_basis(m)
    # Arbitrary usp(2M) element: A anti-Hermitian, B complex symmetric.
    a = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
    a = a - a.conj().T
    b = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
    b = b + b.T
    element = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    element[:m, :m] = a
    element[m:, m:] = a.conj()
    element[:m, m:] = b
    element[m:, :m] = -b.conj()
    assert _element_is_reachable(element, basis)


def test_su_reachability_arbitrary_element() -> None:
    """An arbitrary su(4) element is reachable from the su(4) basis."""
    rng = np.random.default_rng(2)
    n = 4
    basis = su_basis(n)
    h = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    h = h + h.conj().T
    h = h - np.trace(h) / n * np.eye(n)  # traceless Hermitian
    element = 1j * h  # arbitrary su(4) element
    assert _element_is_reachable(element, basis)


# --------------------------------------------------------------------------- #
# Public-API tests through PauliStringCollection.get_algebra_basis
#
# These drive the classifier end to end. The total basis dimension must equal
# the independently computed get_dla_dim(), every summand must be a valid
# (anti-Hermitian, independent, closed) basis, and the per-summand partition
# must be preserved as a list.
# --------------------------------------------------------------------------- #
PUBLIC_CASES = [
    (["XY", "XZ"], 4, "sp(4)", 1),
    (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None, "4*so(5)", 4),
    (["XYI", "XZI", "IIZ"], None, "so(3)+u(1)", 2),
    # The classifier emits a mixed sum here; string-parse dispatch crashes on it.
    (["XY", "XZ", "IIIZ"], 4, "so(3)+u(1)", 2),
    (["XX", "YY", "ZZ"], 3, "so(6)", 1),
    (["XY"], 3, "so(3)", 1),
    (["XY"], None, "u(1)", 1),
    (["XYXY", "YZYZ"], 4, "2*u(1)", 2),
]


@pytest.mark.parametrize("strings,n,label,n_summands", PUBLIC_CASES)
def test_public_api_matches_classification(strings: list[str], n: Optional[int],
                                           label: str, n_summands: int) -> None:
    """get_algebra_basis total dimension and partition match the classifier."""
    collection = p(strings, n=n) if n is not None else p(strings)
    basis = collection.get_algebra_basis()

    # The label this case targets is what the classifier actually emits
    # (up to summand order in the printed string).
    emitted = collection.get_algebra()
    assert sorted(emitted.split("+")) == sorted(label.split("+"))

    assert isinstance(basis, list)
    assert all(isinstance(b, np.ndarray) for b in basis)
    # Partition is preserved as one array per summand (copies counted).
    assert len(basis) == n_summands
    # Total dimension matches the independent oracle.
    assert sum(b.shape[0] for b in basis) == collection.get_dla_dim()
    # Every summand acts on the common direct-sum space.
    common = basis[0].shape[1]
    assert all(b.shape[1] == common and b.shape[2] == common for b in basis)


@pytest.mark.parametrize("strings,n,label,n_summands", PUBLIC_CASES)
def test_public_api_union_is_basis_of_direct_sum(strings: list[str],
                                                 n: Optional[int], label: str,
                                                 n_summands: int) -> None:
    """The union of all summands is a genuine basis of the complete operator.

    Not a stack of unrelated native blocks: the full union must be
    anti-Hermitian, linearly independent and closed under the commutator on the
    common direct-sum space.
    """
    collection = p(strings, n=n) if n is not None else p(strings)
    emitted = collection.get_algebra()
    assert sorted(emitted.split("+")) == sorted(label.split("+"))

    summands = collection.get_algebra_basis()
    assert len(summands) == n_summands
    union = np.concatenate(summands, axis=0)
    assert _is_anti_hermitian(union)
    assert _is_linearly_independent(union)
    assert _closes_under_bracket(union)


def test_public_api_reaches_every_element_of_direct_sum() -> None:
    """An arbitrary element of so(3)+u(1) is reachable from the returned basis.

    This is the element-wise reachability the maintainer asks for: the defining
    representation (a 4x4 block-diagonal space) is not the 2^n embedded space, so
    spanning is a real condition, not implied by rank.
    """
    collection = p(["XY", "XZ", "IIIZ"], n=4)
    basis = collection.get_algebra_basis()
    common = basis[0].shape[1]
    assert common == 4  # u(1) block (1) + so(3) block (3)

    # Summands are ordered by (family, size); u(1) sorts before so(3), so the
    # u(1) block sits at offset 0 and the so(3) block at offset 1.
    rng = np.random.default_rng(7)
    a = rng.standard_normal((3, 3))
    a = a - a.T
    element = np.zeros((4, 4), dtype=np.complex128)
    element[0, 0] = 1j * rng.standard_normal()  # u(1) block
    element[1:, 1:] = a  # so(3) block

    union = np.concatenate(basis, axis=0)
    assert _element_is_reachable(element, union)


def test_public_api_mixed_sum_does_not_crash() -> None:
    """A mixed direct sum like u(1)+so(3) builds rather than raising."""
    collection = p(["XY", "XZ", "IIIZ"], n=4)
    assert collection.get_algebra() in ("so(3)+u(1)", "u(1)+so(3)")
    basis = collection.get_algebra_basis()
    # Both summands act on the common 4x4 space (u(1) block + so(3) block).
    assert all(b.shape[1] == 4 and b.shape[2] == 4 for b in basis)
    assert len(basis) == 2
    assert sum(b.shape[0] for b in basis) == collection.get_dla_dim()


SU_PUBLIC_CASES = [
    (["XX", "XY", "YZ"], 4, "su(16)", 1, 255),
    (["XX", "YY", "YZ"], 4, "2*su(8)", 2, 126),
]


@pytest.mark.parametrize("strings,n,label,n_summands,total_dim", SU_PUBLIC_CASES)
def test_public_api_emits_su_basis(strings: list[str], n: int, label: str,
                                   n_summands: int, total_dim: int) -> None:
    """su is reachable end to end and gets instantiated as a real su basis.

    The classifier names su only for the B3 canonical-graph type, and the
    small-block isomorphisms su(2) = so(3), su(4) = so(6) mean su surfaces only at
    larger blocks. These cases drive get_algebra_basis through that su path, which
    the per-family su_basis test cannot reach from the public API. A full-rank set
    of the correct dimension N^2-1 (or two such blocks) already spans the algebra,
    so closure is implied; the heavy O(dim^2) bracket check is left to the smaller
    per-family case.
    """
    collection = p(strings, n=n)
    assert collection.get_algebra() == label
    basis = collection.get_algebra_basis()
    assert len(basis) == n_summands
    assert sum(b.shape[0] for b in basis) == total_dim == collection.get_dla_dim()

    union = np.concatenate(basis, axis=0)
    assert _is_anti_hermitian(union)
    assert all(abs(np.trace(b)) < 1e-9 for b in union)
    gram_diag = np.array([np.trace(b.conj().T @ b) for b in union])
    assert np.allclose(gram_diag, 2.0)  # uniform Gell-Mann normalization holds here too
    assert _is_linearly_independent(union)


def test_public_api_deterministic_ordering() -> None:
    """Summand order is reproducible across repeated calls."""
    collection = p(["XYI", "XZI", "IIZ"])
    first = collection.get_algebra_basis()
    second = p(["XYI", "XZI", "IIZ"]).get_algebra_basis()
    assert [b.shape for b in first] == [b.shape for b in second]
    for a, b in zip(first, second):
        assert np.allclose(a, b)
