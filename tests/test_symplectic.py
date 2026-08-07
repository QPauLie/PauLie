"""
    Tests for the GF(2) symplectic representation and the commutant it computes.
"""
import numpy as np
import pytest

from paulie import get_identity
from paulie import get_pauli_string as p
from paulie.common.symplectic import (
    from_symplectic,
    null_space,
    row_reduce,
    span,
    symplectic_gram,
    to_symplectic,
)

generator_sets = [
    (["XX", "Z"], 4),
    (["XX", "YY", "Z"], 4),
    (["XX", "YY", "ZZ"], 4),
    (["XY"], 5),
    (["XX"], 5),
    (["ZZ", "X"], 5),
    (["XZ", "YY"], 4),
    (["XX", "YY"], 6),
    (["Z"], 6),
]


def brute_force_commutants(collection) -> set[str]:
    """
    Find the commutant by filtering all 4**n Pauli strings.

    This is what ``get_commutants`` used to do, kept here as the reference the fast
    implementation is checked against.
    """
    candidates = list(get_identity(collection.get_size()).gen_all_pauli_strings())
    for g in collection.get():
        candidates = [q for q in candidates if g.commutes_with(q)]
    return {str(q) for q in candidates}


@pytest.mark.parametrize("pauli_strings", [
    ["XYZI", "IXXZ", "ZZII"],
    ["IIII"],
    ["XXXX", "YYYY", "ZZZZ"],
])
def test_symplectic_roundtrip(pauli_strings: list[str]) -> None:
    """
    Test that from_symplectic inverts to_symplectic.
    """
    original = [p(s) for s in pauli_strings]
    assert [str(q) for q in from_symplectic(to_symplectic(original))] == pauli_strings


def test_symplectic_gram_matches_commutes_with() -> None:
    """
    Test that the Gram matrix agrees with PauliString.commutes_with.
    """
    strings = [p(s) for s in ["XII", "ZII", "IXI", "YZX", "ZZZ", "IIY"]]
    gram = symplectic_gram(to_symplectic(strings))
    for i, a in enumerate(strings):
        for j, b in enumerate(strings):
            expected = 0 if i == j else int(not a.commutes_with(b))
            assert gram[i, j] == expected


def test_null_space_is_a_kernel_basis() -> None:
    """
    Test that null_space returns an independent basis of the kernel over GF(2).
    """
    rng = np.random.default_rng(0)
    for _ in range(200):
        rows, width = int(rng.integers(0, 6)), int(rng.integers(1, 9))
        matrix = rng.integers(0, 2, size=(rows, width)).astype(np.uint8)
        basis = null_space(matrix, width)
        _, pivots = row_reduce(matrix)

        assert basis.shape == (width - len(pivots), width)
        if rows and basis.size:
            assert not ((matrix @ basis.T) % 2).any()
            assert len(row_reduce(basis)[1]) == basis.shape[0]


def test_span_enumerates_the_subspace() -> None:
    """
    Test that span returns every distinct combination of the basis.
    """
    basis = np.array([[1, 0, 1], [0, 1, 1]], dtype=np.uint8)
    vectors = span(basis)
    assert vectors.shape == (4, 3)
    assert len({tuple(v) for v in vectors}) == 4
    assert (0, 0, 0) in {tuple(v) for v in vectors}


@pytest.mark.parametrize("generators,n", generator_sets)
def test_commutants_match_brute_force(generators: list[str], n: int) -> None:
    """
    Test that the null-space commutant is the same set as the filtered one.
    """
    collection = p(generators, n=n)
    assert {str(q) for q in collection.get_commutants()} == brute_force_commutants(collection)


@pytest.mark.parametrize("generators,n", generator_sets)
def test_commutant_basis_generates_the_commutant(generators: list[str], n: int) -> None:
    """
    Test that the basis is independent and spans exactly the commutant.
    """
    collection = p(generators, n=n)
    basis = collection.get_commutant_basis()
    commutants = {str(q) for q in collection.get_commutants()}

    assert len(commutants) == 2 ** len(basis)
    for b in basis:
        assert str(b) in commutants
        assert all(g.commutes_with(b) for g in collection.get())


@pytest.mark.parametrize("n", [2, 3, 4])
def test_identity_commutes_with_everything(n: int) -> None:
    """
    Test that the commutant of the identity is the whole Pauli group.
    """
    assert len({str(q) for q in get_identity(n).get_commutants()}) == 4**n


@pytest.mark.parametrize("n", [8, 12, 20])
def test_commutant_basis_scales_past_brute_force(n: int) -> None:
    """
    Test that the basis is available at sizes where filtering 4**n is hopeless.
    """
    basis = p(["XX", "Z"], n=n).get_commutant_basis()
    assert 0 < len(basis) <= 2 * n
