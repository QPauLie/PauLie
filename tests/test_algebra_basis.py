import numpy as np
import pytest
from paulie import get_pauli_string as p


def check_linear_independence(basis_tensor):
    """Verify linear independence of basis."""
    dim = basis_tensor.shape[0]
    # Flatten matrices [dim, N, N] do [dim, N*N]
    flat_shape = (dim, -1)
    flat_mats = basis_tensor.reshape(flat_shape)
    # Compute algebraic rank of vector set
    rank = np.linalg.matrix_rank(flat_mats)
    assert rank == dim


def test_canonical_type_so():
    gens = p(["IX", "XI"], n=2)
    basis = gens.get_algebra_basis()
    assert len(basis) >= 1
    for summand in basis:
        dim, N, _ = summand.shape
        # For u(1) dimension is 1, otherwise use so formula
        if N == 1:
            assert dim == 1
        else:
            assert dim == (N * (N - 1)) // 2

        # Sprawdzenie liniowej niezaleznosci
        check_linear_independence(summand)


def test_canonical_type_sp():
    gens = p(["XY", "XZ"], n=4)
    basis = gens.get_algebra_basis()
    assert len(basis) == 1
    summand = basis[0]
    dim, N, _ = summand.shape
    assert dim == (N // 2) * (N + 1)
    n_half = N // 2
    I = np.eye(n_half)
    J = np.block(
        [
            [np.zeros((n_half, n_half)), I],
            [-I, np.zeros((n_half, n_half))],
        ]
    )
    for mat in summand:
        cond = mat.T @ J + J @ mat
        assert np.allclose(cond, 0.0)

    # Sprawdzenie liniowej niezaleznosci
    check_linear_independence(summand)


def test_canonical_type_su():
    from paulie.helpers._lie_bases import (
        get_su_basis,
    )

    basis = get_su_basis(3)
    dim, N, _ = basis.shape
    assert dim == 8
    for mat in basis:
        assert np.allclose(
            mat, -mat.view(np.ndarray).conj().T
        )
        assert np.allclose(np.trace(mat), 0.0)

    # Sprawdzenie liniowej niezaleznosci
    check_linear_independence(basis)
