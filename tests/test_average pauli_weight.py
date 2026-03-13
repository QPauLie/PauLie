"""
Test of average pauli weight
"""
import itertools
import pytest
import numpy as np
from paulie import (
    PauliString,
    quantum_fourier_entropy, average_pauli_weight
)

def get_O(n_qubit):
    """
    Create operator O = 2^(n/2) |0><0|^kron n
    """
    operator = np.zeros((2 ** n_qubit, 2 ** n_qubit))
    operator[0,0] = 2 ** (n_qubit / 2)

    return operator

test_case_pauli_str = [
    ("I", 0), ("X", 1), ("Y", 1), ("Z", 1),
    ("II", 0), ("IX", 1), ("IY", 1), ("IZ", 1),
    ("XI", 1), ("XX", 2), ("XY", 2), ("XZ", 2),
    ("YI", 1), ("YX", 2), ("YY", 2), ("YZ", 2),
    ("ZI", 1), ("ZX", 2), ("ZY", 2), ("ZZ", 2),
]

test_n_qubit = [1, 2, 3, 4]
test_case_weight_O = [
    (n_qubit1, n_qubit2)
    for (n_qubit1, n_qubit2)
    in itertools.product(test_n_qubit, repeat=2)
]

@pytest.mark.parametrize("P, expected_weight", test_case_pauli_str)
def test_avg_pauli_weight_pauli_str(P, expected_weight) -> None:
    """
    Test average Pauli weight from PauliString
    """

    pauli_string = PauliString(pauli_str=P)
    pauli_string_matrix = pauli_string.get_matrix()

    assert np.isclose(average_pauli_weight(pauli_string_matrix), expected_weight, atol=1e-10)

@pytest.mark.parametrize("n1, n2", test_case_weight_O)
def test_avg_pauli_weight_O(n1, n2) -> None:
    """
    Test average Pauli weight from operator O
    """

    operator1 = get_O(n1)
    operator2 = get_O(n2)
    operator12 = np.kron(operator1, operator2)

    weight1 = average_pauli_weight(operator1)
    weight2 = average_pauli_weight(operator2)
    weight12 = average_pauli_weight(operator12)

    assert np.isclose(weight1, n1 / 2, atol=1e-10)
    assert np.isclose(weight2, n2 / 2, atol=1e-10)
    assert np.isclose(weight12, weight1 + weight2, atol=1e-10)

@pytest.mark.parametrize("P, _", test_case_pauli_str)
def test_quantum_fourier_entropy_pauli_str(P, _) -> None:
    """
    Test quantum Fourier entropy from PauliString
    """

    pauli_string = PauliString(pauli_str=P)
    pauli_string_matrix = pauli_string.get_matrix()

    assert np.isclose(quantum_fourier_entropy(pauli_string_matrix), 0, atol=1e-10)

@pytest.mark.parametrize("n, _", test_case_weight_O)
def test_quantum_fourier_entropy_O(n, _) -> None:
    """
    Test quantum Fourier entropy from operator O
    """

    operator = get_O(n)
    entropy = quantum_fourier_entropy(operator)

    assert np.isclose(entropy, n, atol=1e-10)
