"""Test the matrix_decomposition function"""
import numpy as np
import pytest
from paulie.common.pauli_string_factory import get_pauli_string
from paulie.application.matrix_decomposition import matrix_decomposition


def test_identity_decomposition_single_qubit() -> None:
    """Identity on 1 qubit decomposes to coefficient 1 on Pauli I"""
    Identity = np.eye(2)
    result = matrix_decomposition(Identity)
    assert len(result) == 1
    ps = get_pauli_string("I")
    assert ps in result
    assert result[ps] == pytest.approx(1.0)


def test_single_qubit_paulis() -> None:
    """Each single-qubit Pauli has coefficient 1 on itself"""
    paulis = {
        "X": np.array([[0, 1], [1, 0]]),
        "Y": np.array([[0, -1j], [1j, 0]]),
        "Z": np.array([[1, 0], [0, -1]])
    }
    for label, mat in paulis.items():
        result = matrix_decomposition(mat)
        assert len(result) == 1
        ps = get_pauli_string(label)
        assert ps in result
        assert result[ps] == pytest.approx(1.0)


def test_two_qubit_tensor_pauli() -> None:
    """Tensor product of Paulis decomposes to the corresponding Pauli string"""
    X = np.array([[0, 1], [1, 0]])
    Z = np.array([[1, 0], [0, -1]])
    XZ = np.kron(X, Z)
    result = matrix_decomposition(XZ)
    assert len(result) == 1
    ps = get_pauli_string("XZ")
    assert ps in result
    assert result[ps] == pytest.approx(1.0)


def test_linear_combination_two_qubits() -> None:
    """Check decomposition of a linear combination of II and XY"""
    Identity = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    # Construct 0.5*I⊗I + 0.25*X⊗Y
    combo = 0.5 * np.kron(Identity, Identity) + 0.25 * np.kron(X, Y)
    result = matrix_decomposition(combo)
    # Expect two non-zero coefficients: II and XY
    assert len(result) == 2
    ii = get_pauli_string("II")
    xy = get_pauli_string("XY")
    assert ii in result and xy in result
    assert result[ii] == pytest.approx(0.5)
    assert result[xy] == pytest.approx(0.25)
