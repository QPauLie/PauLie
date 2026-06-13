"""
Tests for Lie algebra basis matrix generators.
"""

import pytest
import numpy as np

# Adjust this import to match your actual module structure
from paulie.classifier.classification import TypeAlgebra
from paulie.common.algebra_basis import (
    get_u_basis,
    get_n_so_basis,
    get_so_basis,
    get_n_su_basis,
    get_su_basis,
    get_n_sp_basis,
    get_sp_basis,
    get_group_basis,
    get_n_basis,
    get_algebras_basis,
)

# ---------------------------------------------------------------------------
# Mathematical Property Helpers
# ---------------------------------------------------------------------------


def is_skew_hermitian(mat: np.ndarray) -> bool:
    """Check if a matrix is skew-Hermitian (X^dagger = -X)."""
    return np.allclose(mat + mat.conj().T, 0)


def is_antisymmetric(mat: np.ndarray) -> bool:
    """Check if a matrix is antisymmetric (X^T = -X)."""
    return np.allclose(mat + mat.T, 0)


def is_traceless(mat: np.ndarray) -> bool:
    """Check if a matrix is traceless (Tr(X) = 0)."""
    return np.isclose(np.trace(mat), 0)


def is_symplectic(mat: np.ndarray, n: int) -> bool:
    """
    Check the symplectic condition (X^T J + J X = 0).
    J = [[0, I_n], [-I_n, 0]]
    """
    identity = np.eye(n)
    zeros = np.zeros((n, n))
    J = np.block([[zeros, identity], [-identity, zeros]])
    return np.allclose(mat.T @ J + J @ mat, 0)


# ---------------------------------------------------------------------------
# u(1) Tests
# ---------------------------------------------------------------------------


def test_u1_basis_properties() -> None:
    """Test dimension, shape, and properties of u(1) basis."""
    basis = get_u_basis(1)
    assert basis.shape == (1, 1, 1)
    assert np.allclose(basis[0], [[1j]])
    assert is_skew_hermitian(basis[0])


def test_u1_basis_errors() -> None:
    """Test that u(n) raises ValueError for n != 1."""
    with pytest.raises(ValueError, match="currently supported"):
        get_u_basis(2)


# ---------------------------------------------------------------------------
# so(n) Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n, expected_dim", [(2, 1), (3, 3), (4, 6)])
def test_so_basis_properties(n: int, expected_dim: int) -> None:
    """Test dimension calculation and matrix properties for so(n)."""
    assert get_n_so_basis(n) == expected_dim

    basis = get_so_basis(n)
    assert basis.shape == (expected_dim, n, n)

    for mat in basis:
        assert is_skew_hermitian(mat)
        assert is_antisymmetric(mat)
        assert np.all(np.isreal(mat))  # so(n) generators are strictly real


def test_so_basis_errors() -> None:
    """Test error handling for so(n) edge cases."""
    with pytest.raises(ValueError, match="positive integer"):
        get_n_so_basis(0)
    with pytest.raises(TypeError, match="positive integer"):
        get_n_so_basis(3.5)  # type: ignore


# ---------------------------------------------------------------------------
# su(n) Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n, expected_dim", [(5, 24), (6, 35)])
def test_su_basis_properties(n: int, expected_dim: int) -> None:
    """Test dimension calculation and matrix properties for su(n)."""
    assert get_n_su_basis(n) == expected_dim

    basis = get_su_basis(n)
    assert basis.shape == (expected_dim, n, n)

    for mat in basis:
        assert is_skew_hermitian(mat)
        assert is_traceless(mat)


def test_su_basis_errors() -> None:
    """Test error handling for su(n), bounded by n > 2^2 (4)."""
    with pytest.raises(ValueError, match=r"greater than 2\^2"):
        get_n_su_basis(4)
    with pytest.raises(TypeError):
        get_n_su_basis(5.5)  # type: ignore


# ---------------------------------------------------------------------------
# sp(n) Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n, expected_dim", [(1, 3), (2, 10), (4, 36)])
def test_sp_basis_properties(n: int, expected_dim: int) -> None:
    """Test dimension calculation and matrix properties for sp(n)."""
    assert get_n_sp_basis(n) == expected_dim

    basis = get_sp_basis(n)
    # sp(n) matrices are 2n x 2n
    assert basis.shape == (expected_dim, 2 * n, 2 * n)

    for mat in basis:
        assert is_skew_hermitian(mat)
        assert is_symplectic(mat, n)


def test_sp_basis_errors() -> None:
    """Test error handling for sp(n) power-of-2 requirements."""
    with pytest.raises(ValueError, match="power of 2 AND greater or equal to 2"):
        get_n_sp_basis(3)  # 3 is not a power of 2
    with pytest.raises(TypeError):
        get_n_sp_basis(2.0)  # type: ignore


# ---------------------------------------------------------------------------
# Dispatcher & Direct Sum Tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "algebra_type, n, expected_dim",
    [
        (TypeAlgebra.U, 1, 1),
        (TypeAlgebra.SO, 3, 3),  # n*(n-1)/2 -> 3*(2)/2 = 3
        (TypeAlgebra.SU, 5, 24),  # n^2 - 1 -> 25 - 1 = 24
        (TypeAlgebra.SP, 2, 10),  # n*(2n+1) -> 2*(5) = 10
    ],
)
def test_get_n_basis_dispatch(algebra_type: TypeAlgebra, n: int, expected_dim: int) -> None:
    """Test that get_n_basis correctly calculates total dimensions for each algebra type."""
    assert get_n_basis(algebra_type, n) == expected_dim


def test_get_n_basis_errors() -> None:
    """Test error handling for unsupported algebra types in get_n_basis."""
    with pytest.raises(ValueError, match="Unsupported algebra type"):
        get_n_basis("INVALID_TYPE", 2)  # type: ignore


def test_get_group_basis_dispatch() -> None:
    """Test that the dispatcher correctly routes to the right generator."""
    assert get_group_basis(TypeAlgebra.U, 1).shape == (1, 1, 1)
    assert get_group_basis(TypeAlgebra.SO, 3).shape == (3, 3, 3)
    assert get_group_basis(TypeAlgebra.SU, 5).shape == (24, 5, 5)
    assert get_group_basis(TypeAlgebra.SP, 2).shape == (10, 4, 4)

    # Use a dummy invalid string to test the default fallback
    with pytest.raises(ValueError, match="Unsupported algebra type"):
        get_group_basis("INVALID_TYPE", 2)  # type: ignore


def test_get_algebras_basis_direct_sum() -> None:
    """Test block-diagonal allocation for direct sums of Lie algebras."""
    multipliers = [2, 1]
    groups = [TypeAlgebra.U, TypeAlgebra.SO]
    sizes = [1, 3]

    # Total expected dimension: 2*(u1 dim 1) + 1*(so3 dim 3) = 5
    # Total expected matrix size: 2*(1) + 1*(3) = 5
    basis = get_algebras_basis(multipliers, groups, sizes)
    assert basis.shape == (5, 5, 5)

    # Validate the block diagonal isolation for the first matrix
    # The first U(1) generator should only occupy the [0, 0] cell.
    assert basis[0, 0, 0] != 0
    assert np.allclose(basis[0, 1:, :], 0)
    assert np.allclose(basis[0, :, 1:], 0)

    # Validate the last generator belongs to the SO(3) block
    # and leaves the upper blocks isolated (0 to 1).
    assert np.allclose(basis[-1, 0:2, 0:2], 0)


def test_get_algebras_basis_errors() -> None:
    """Test error handling for mismatched list inputs in direct sums."""
    with pytest.raises(ValueError, match="same length"):
        get_algebras_basis([1, 2], [TypeAlgebra.U], [1])
