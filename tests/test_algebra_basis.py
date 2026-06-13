"""
Test matrix bases for the classified Lie algebras.
"""
from itertools import combinations

import numpy as np
import pytest

from paulie import get_pauli_string as p
from paulie.common.algebra_basis import (
    block_diagonal_basis,
    so_basis,
    sp_basis,
    su_basis,
    symplectic_form,
    u1_basis,
)
from paulie.common.two_local_generators import G_LIE

ATOL = 1e-9


def _anti_hermitian(basis: np.ndarray) -> bool:
    """Whether every matrix in the stack is anti-Hermitian."""
    return all(np.allclose(matrix.conj().T, -matrix, atol=ATOL) for matrix in basis)


def _independent(basis: np.ndarray) -> bool:
    """Whether the matrices are linearly independent."""
    flat = basis.reshape(basis.shape[0], -1)
    return np.linalg.matrix_rank(flat, tol=ATOL) == basis.shape[0]


def _structure_constants(basis: np.ndarray) -> np.ndarray:
    """Coefficients of every bracket ``[B_i, B_j]`` in the basis."""
    flat = basis.reshape(basis.shape[0], -1).T
    dim = basis.shape[0]
    constants = np.zeros((dim, dim, dim), dtype=np.complex128)
    for i in range(dim):
        for j in range(dim):
            bracket = (basis[i] @ basis[j] - basis[j] @ basis[i]).reshape(-1)
            constants[i, j], *_ = np.linalg.lstsq(flat, bracket, rcond=None)
    return constants


def _bracket_closed(basis: np.ndarray) -> bool:
    """Whether the span is closed under the Lie bracket."""
    flat = basis.reshape(basis.shape[0], -1).T
    for i, j in combinations(range(basis.shape[0]), 2):
        bracket = (basis[i] @ basis[j] - basis[j] @ basis[i]).reshape(-1)
        coefficients, *_ = np.linalg.lstsq(flat, bracket, rcond=None)
        if not np.allclose(flat @ coefficients, bracket, atol=ATOL):
            return False
    return True


def _in_span(matrix: np.ndarray, basis: np.ndarray) -> bool:
    """Whether a matrix lies in the span of the basis."""
    flat = basis.reshape(basis.shape[0], -1).T
    coefficients, *_ = np.linalg.lstsq(flat, matrix.reshape(-1), rcond=None)
    return np.allclose(flat @ coefficients, matrix.reshape(-1), atol=ATOL)


def _spans_equal(first: np.ndarray, second: np.ndarray) -> bool:
    """Whether two stacks of matrices span the same space."""
    flat_first = first.reshape(first.shape[0], -1)
    flat_second = second.reshape(second.shape[0], -1)
    rank_first = np.linalg.matrix_rank(flat_first, tol=ATOL)
    rank_second = np.linalg.matrix_rank(flat_second, tol=ATOL)
    rank_union = np.linalg.matrix_rank(np.vstack([flat_first, flat_second]), tol=ATOL)
    return rank_first == rank_second == rank_union


def _lie_closure(generators: list) -> list:
    """Brute-force Lie closure as Pauli words, using ``commutes_with`` and ``@``."""
    closure = list(generators)
    changed = True
    while changed:
        changed = False
        for left in list(closure):
            for right in list(closure):
                if left.commutes_with(right):
                    continue
                product = left @ right
                if all(product != element for element in closure):
                    closure.append(product)
                    changed = True
    return closure


def _embedded_dla(collection, n: int) -> np.ndarray:
    """The embedded DLA as anti-Hermitian ``2^n`` matrices ``i P`` over the Lie closure."""
    closure = _lie_closure([generator.expand(n) for generator in collection.get()])
    return np.array([1.0j * word.get_matrix() for word in closure])


def _majorana_operators(n: int) -> list:
    """The ``2n`` Jordan-Wigner Majorana operators on ``n`` qubits as dense matrices."""
    operators = []
    for site in range(n):
        prefix = "Z" * site
        operators.append(p(prefix + "X" + "I" * (n - site - 1)).get_matrix())
        operators.append(p(prefix + "Y" + "I" * (n - site - 1)).get_matrix())
    return operators


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_so_basis_is_compact(n) -> None:
    """The so(n) basis is real antisymmetric, complete, and closed under the bracket.
    """
    basis = so_basis(n)
    assert basis.shape == (n * (n - 1) // 2, n, n)
    assert all(np.allclose(matrix, -matrix.T) for matrix in basis)
    assert _anti_hermitian(basis) and _independent(basis) and _bracket_closed(basis)
    assert all(np.isclose(np.trace(matrix.conj().T @ matrix), 2.0) for matrix in basis)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_su_basis_is_compact(n) -> None:
    """The su(n) basis is traceless anti-Hermitian and uniformly normalized.
    """
    basis = su_basis(n)
    assert basis.shape == (n * n - 1, n, n)
    assert _anti_hermitian(basis)
    assert all(abs(np.trace(matrix)) < ATOL for matrix in basis)
    assert all(np.isclose(np.trace(matrix.conj().T @ matrix), 2.0) for matrix in basis)
    assert _independent(basis) and _bracket_closed(basis)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_sp_basis_is_compact(n) -> None:
    """The sp(n) basis is anti-Hermitian and symplectic, the compact real form.
    """
    basis = sp_basis(n)
    form = symplectic_form(n)
    assert basis.shape == (n * (2 * n + 1), 2 * n, 2 * n)
    assert _anti_hermitian(basis)
    assert all(np.allclose(matrix, -form @ matrix.T @ form.T, atol=ATOL) for matrix in basis)
    assert _independent(basis) and _bracket_closed(basis)


def test_u1_basis() -> None:
    """The u(1) basis is the single generator [[i]].
    """
    basis = u1_basis()
    assert basis.shape == (1, 1, 1)
    assert np.allclose(basis[0], [[1.0j]])


@pytest.mark.parametrize("constructor", [so_basis, su_basis, sp_basis])
@pytest.mark.parametrize("size", [0, -1])
def test_constructors_reject_non_positive_size(constructor, size) -> None:
    """Every constructor rejects a non-positive size.
    """
    with pytest.raises(ValueError):
        constructor(size)


def test_block_diagonal_basis_rejects_empty() -> None:
    """A direct sum needs at least one summand.
    """
    with pytest.raises(ValueError):
        block_diagonal_basis([])


def test_direct_sum_is_a_basis_of_the_complete_operator() -> None:
    """A direct sum embeds block-diagonally and stays full rank, not identical copies.
    """
    basis = block_diagonal_basis([so_basis(3), so_basis(3)])
    assert basis.shape == (6, 6, 6)
    assert _anti_hermitian(basis) and _independent(basis) and _bracket_closed(basis)
    assert np.allclose(basis[:3, 3:, 3:], 0) and np.allclose(basis[3:, :3, :3], 0)


# (generators, qubit number, expected algebra label)
ALGEBRA_CASES = [
    (["XY"], 3, "so(3)"),
    (["XX", "XZ"], 3, "so(5)"),
    (["XX", "YY", "XY"], 3, "so(6)"),
    (["XY", "XZ"], 4, "sp(4)"),
    (["XX", "XY", "XZ", "YX"], 3, "su(8)"),
    (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None, "4*so(5)"),
    (["XYI", "XZI", "IIZ"], 4, "u(1)+so(6)"),
]


@pytest.mark.parametrize("generators, n, label", ALGEBRA_CASES)
def test_basis_dimension_matches_dla_dim(generators, n, label) -> None:
    """The basis is anti-Hermitian with exactly get_dla_dim generators, including mixed sums.
    """
    collection = p(generators) if n is None else p(generators, n=n)
    # get_algebra orders summands non-deterministically, so compare the multiset of terms.
    assert sorted(collection.get_algebra().split("+")) == sorted(label.split("+"))
    basis = collection.get_algebra_basis()
    assert basis.shape[0] == collection.get_dla_dim()
    assert _anti_hermitian(basis) and _independent(basis)


def test_each_summand_block_satisfies_its_own_condition() -> None:
    """Each diagonal block, not the whole operator, satisfies its summand condition.
    """
    basis = p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]).get_algebra_basis()
    assert basis.shape == (40, 20, 20)
    for element in basis:
        for row in range(0, 20, 5):
            for column in range(0, 20, 5):
                block = element[row:row + 5, column:column + 5]
                if row == column:
                    assert np.allclose(block, -block.T)
                else:
                    assert np.allclose(block, 0)


def test_full_su_reachability_element_wise() -> None:
    """For a full su(8) DLA every Lie-closure element is reached by the basis.
    """
    collection = p(["XX", "XY", "XZ", "YX"], n=3)
    assert collection.get_algebra() == "su(8)"
    basis = collection.get_algebra_basis()
    embedded = _embedded_dla(collection, 3)
    assert basis.shape[0] == embedded.shape[0] == collection.get_dla_dim()
    assert all(_in_span(element, basis) for element in embedded)


@pytest.mark.parametrize("n, label", [(3, "so(6)"), (4, "so(8)")])
def test_nontrivial_so_reachability_via_majorana(n, label) -> None:
    """For so(2n) the basis reaches every Lie-closure element, with a nontrivial correspondence.

    The defining representation is ``2n x 2n`` while the embedding is ``2^n x 2^n``, so the
    correspondence is nontrivial. The map ``E_ij - E_ji <-> (1/2) gamma_i gamma_j`` is an explicit
    Lie isomorphism between the so(2n) basis and the Majorana bilinears, which span the embedded
    DLA, so the basis reaches every element through it.
    """
    collection = p(["XX", "YY", "XY"], n=n)
    assert collection.get_algebra() == label
    basis = collection.get_algebra_basis()
    gamma = _majorana_operators(n)
    bilinears = np.array([0.5 * gamma[i] @ gamma[j] for i, j in combinations(range(2 * n), 2)])
    assert np.allclose(_structure_constants(basis), _structure_constants(bilinears), atol=ATOL)
    embedded = _embedded_dla(collection, n)
    assert _spans_equal(bilinears, embedded)
    assert all(_in_span(element, bilinears) for element in embedded)


def test_compact_sp_reachability() -> None:
    """The sp(4) basis is the compact form usp(8) and is a basis of the embedded DLA.

    The flagship sp(4) DLA embeds as anti-Hermitian ``16 x 16`` operators while the defining rep
    is ``8 x 8``. The basis is certified to be the compact real form (anti-Hermitian and
    symplectic), and its dimension equals both get_dla_dim and the embedded Lie-closure
    dimension, so it is a basis of that algebra.
    """
    collection = p(["XY", "XZ"], n=4)
    assert collection.get_algebra() == "sp(4)"
    basis = collection.get_algebra_basis()
    form = symplectic_form(4)
    assert basis.shape == (36, 8, 8)
    assert _anti_hermitian(basis)
    assert all(np.allclose(matrix, -form @ matrix.T @ form.T, atol=ATOL) for matrix in basis)
    embedded = _embedded_dla(collection, 4)
    assert _anti_hermitian(embedded)
    assert basis.shape[0] == embedded.shape[0] == collection.get_dla_dim()


@pytest.mark.slow
@pytest.mark.parametrize("name", ["a1", "a2", "a4", "a7", "a8", "a9", "a14", "a22", "b0", "b3"])
def test_basis_matches_dla_dim_across_two_local_algebras(name) -> None:
    """The basis dimension matches get_dla_dim across a spread of two-local algebras.
    """
    collection = p(G_LIE[name], n=3)
    basis = collection.get_algebra_basis()
    assert basis.shape[0] == collection.get_dla_dim()
    assert _anti_hermitian(basis) and _independent(basis)


def test_dropping_a_generator_changes_the_basis() -> None:
    """The basis tracks the DLA: dropping a generator changes its dimension.
    """
    full = p(["XX", "YY", "XY"], n=3).get_algebra_basis()
    reduced = p(["XX", "YY"], n=3).get_algebra_basis()
    assert full.shape[0] != reduced.shape[0]
