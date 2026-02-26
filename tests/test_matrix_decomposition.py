"""
    Test Pauli decomposition routines.
"""
import pytest
import numpy as np
from paulie import (
    matrix_decomposition,
    matrix_decomposition_diagonal
)
from paulie.common.pauli_string_factory import get_pauli_string as p

def test_matrix_decomposition_errors() -> None:
    """
    Assert that matrix_decomposition raises errors for
    invalid input.
    """
    non_square = np.zeros((2, 3))
    not_power_of_two = np.zeros((3, 3))
    not_2d = np.zeros((6, 4, 3))
    scalar = np.zeros(1)
    with pytest.raises(ValueError):
        matrix_decomposition(non_square)
    with pytest.raises(ValueError):
        matrix_decomposition(not_power_of_two)
    with pytest.raises(ValueError):
        matrix_decomposition(not_2d)
    with pytest.raises(ValueError):
        matrix_decomposition(scalar)

@pytest.mark.parametrize("paulistr", [
    "I", "X", "Y", "Z",
    "II", "IX", "IY", "IZ",
    "XI", "XX", "XY", "XZ",
    "YI", "YX", "YY", "YZ",
    "ZI", "ZX", "ZY", "ZZ",
    "ZXIYZ", "XIYIZ",
    "YXZIXY", "XYXYZZ",
])
def test_pauli_string_matrices(paulistr: str) -> None:
    """
    Assert that matrix_decomposition recovers a Pauli string
    from its matrix.
    """
    pstr = p(paulistr)
    pmat = pstr.get_matrix()
    decomp = matrix_decomposition(pmat)
    # There should only be one entry
    assert np.count_nonzero(~np.isclose(np.abs(decomp), 0)) == 1
    # Non-zero entry must be equal to 1 and in correct position
    assert np.real(pstr.get_weight_in_matrix(decomp)) == pytest.approx(1.0)
    assert np.imag(pstr.get_weight_in_matrix(decomp)) == pytest.approx(0.0)

def test_matrix_decomposition_diagonal_errors() -> None:
    """
    Assert that matrix_decomposition_diagonal raises errors
    for invalid input.
    """
    not_1d = np.zeros((2, 3))
    not_power_of_two = np.zeros(3)
    scalar = np.zeros(1)
    with pytest.raises(ValueError):
        matrix_decomposition(not_1d)
    with pytest.raises(ValueError):
        matrix_decomposition(not_power_of_two)
    with pytest.raises(ValueError):
        matrix_decomposition(scalar)

@pytest.mark.parametrize("paulistr", [
    "I", "Z",
    "II", "IZ", "ZI", "ZZ",
    "ZIZIZ", "ZZZIZIZI",
    "IIIZZZIZIZZI",
])
def test_pauli_string_diagonal_matrices(paulistr: str) -> None:
    """
    Assert that matrix_decomposition_diagonal recovers
    a Pauli string of only I and Z from its matrix.
    """
    pstr = p(paulistr)
    pmat = np.diag(pstr.get_matrix())
    decomp = matrix_decomposition_diagonal(pmat)
    # There should only be one entry
    assert np.count_nonzero(~np.isclose(np.abs(decomp), 0)) == 1
    # Non-zero entry must be equal to 1 and in correct position
    assert np.real(pstr.get_weight_in_matrix(decomp)) == pytest.approx(1.0)
    assert np.imag(pstr.get_weight_in_matrix(decomp)) == pytest.approx(0.0)
