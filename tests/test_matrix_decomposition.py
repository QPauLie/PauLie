"""
    Test Pauli decomposition routines.
"""
import pytest
from paulie.application.matrix_decomposition import matrix_decomposition, matrix_decomposition_diagonal
from paulie.common.pauli_string_factory import get_pauli_string as p
import numpy as np

def test_matrix_decomposition_errors() -> None:
    """
    Assert that matrix_decomposition raises errors for
    invalid input.
    """
    A = np.zeros((2, 3))
    B = np.zeros((3, 3))
    C = np.zeros((6, 4, 3))
    D = np.zeros(1)
    with pytest.raises(ValueError):
        matrix_decomposition(A)
    with pytest.raises(ValueError):
        matrix_decomposition(B)
    with pytest.raises(ValueError):
        matrix_decomposition(C)
    with pytest.raises(ValueError):
        matrix_decomposition(D)

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
    assert len(decomp) == 1
    for ustr, coeff in decomp.items():
        # ustr must match pstr
        assert ustr == pstr
        # coeff must be real and equal to 1
        assert np.imag(coeff) == pytest.approx(0)
        assert np.real(coeff) == pytest.approx(1)

def test_matrix_decomposition_diagonal_errors() -> None:
    """
    Assert that matrix_decomposition_diagonal raises errors
    for invalid input.
    """
    A = np.zeros((2, 3))
    B = np.zeros(3)
    C = np.zeros(1)
    with pytest.raises(ValueError):
        matrix_decomposition(A)
    with pytest.raises(ValueError):
        matrix_decomposition(B)
    with pytest.raises(ValueError):
        matrix_decomposition(C)

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
    assert len(decomp) == 1
    for ustr, coeff in decomp.items():
        # ustr must match pstr
        assert ustr == pstr
        # coeff must be real and equal to 1
        assert np.imag(coeff) == pytest.approx(0)
        assert np.real(coeff) == pytest.approx(1)
