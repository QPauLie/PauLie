import numpy as np
import pytest
from paulie.common.pauli_string_factory import get_pauli_string
from paulie.application.matrix_decomposition import matrix_decomposition

# Single-qubit Pauli basis tests
@pytest.mark.parametrize("matrix,label", [
    (np.eye(2), "I"),
    (np.array([[0, 1], [1, 0]]), "X"),
    (np.array([[0, -1j], [1j, 0]]), "Y"),
    (np.array([[1, 0], [0, -1]]), "Z"),
])
def test_single_qubit_pauli_decomposition(matrix: np.ndarray, label: str) -> None:
    """
    Decomposing each single-qubit Pauli yields a coefficient of 1 on itself.
    """
    result = matrix_decomposition(matrix)
    assert len(result) == 1
    ps = get_pauli_string(label)
    assert ps in result
    assert result[ps] == pytest.approx(1.0)

# Two-qubit tensor-product tests
@pytest.mark.parametrize("label,matrix", [
    ("II", np.kron(np.eye(2), np.eye(2))),
    ("XZ", np.kron(np.array([[0, 1], [1, 0]]), np.array([[1, 0], [0, -1]]))),
])
def test_two_qubit_tensor_pauli(label: str, matrix: np.ndarray) -> None:
    """
    Decomposing a tensor product of single-qubit Paulis yields a single Pauli string.
    """
    result = matrix_decomposition(matrix)
    assert len(result) == 1
    ps = get_pauli_string(label)
    assert ps in result
    assert result[ps] == pytest.approx(1.0)

# Linear combination test
@pytest.mark.parametrize("coeff,label,matrix", [
    (0.5, "II", np.kron(np.eye(2), np.eye(2))),
    (0.25, "XY", np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]))),
])
def test_linear_combination_two_qubits(coeff: float, label: str, matrix: np.ndarray) -> None:
    """
    Decomposing a linear combination yields expected coefficients.
    """
    result = matrix_decomposition(matrix * coeff)
    ps = get_pauli_string(label)
    assert len(result) == 1
    assert ps in result
    assert result[ps] == pytest.approx(coeff)

# Invalid input tests
@pytest.mark.parametrize("matrix", [
    np.zeros((3, 3)),  # not power-of-2 dimension
    np.zeros((2, 4)),  # non-square
])
def test_invalid_dimensions(matrix: np.ndarray) -> None:
    """
    matrix_decomposition should raise ValueError on non-square or non-power-of-2 matrices.
    """
    with pytest.raises(ValueError):
        matrix_decomposition(matrix)
