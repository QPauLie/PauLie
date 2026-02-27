"""
    Module to compute the average Pauli weight (influence) and quantum Fourier entropy of an
    operator O.
"""
import numpy as np
from paulie.application.matrix_decomposition import matrix_decomposition

def quantum_fourier_entropy(o: np.ndarray) -> float:
    r"""
    Finds the quantum Fourier entropy of an operator :math:`O`.

    The quantum Fourier entropy :math:`H` is given by

    .. math::

        O &= \sum_{P \in \text{all Pauli strings}} c_P P \\
        H(O) &= -\sum_{P} c_{P}^{2}\log_2 \left(c_{P}^{2}\right)

    Args:
        o (numpy.ndarray): The operator :math:`O` in matrix form.

    Returns:
        float: The quantum Fourier entropy of the operator.
    """
    # Get the coefficients c_P from the Pauli decomposition
    c_p = matrix_decomposition(o)
    # Calculate the probabilities p_P = c_P^2
    probs = np.abs(c_p)**2
    # Filter out zero probabilities to avoid log(0)
    non_zero_probs = probs[probs > 1e-12]
    # Calculate the Shannon entropy using base 2 for the logarithm
    entropy = -np.sum(non_zero_probs * np.log2(non_zero_probs))
    return entropy

def get_pauli_weights(num_qubits: int, identity_pos: int=0) -> np.ndarray:
    r"""
    Generates the weight :math:`|P|` for each of the :math:`4^{\text{number qubits}}` Pauli
    operators.

    .. currentmodule::
        paulie.application.matrix_decomposition

    The weight is the number of non-identity terms in the Pauli string. The ordering corresponds to
    the output of :func:`matrix_decomposition`, default is 'I' at position 0.

    Args:
        num_qubits (int): The number of qubits.
        identity_pos (int, optional): The position of the identity operator. Defaults to 0.

    Returns:
        numpy.ndarray: The weights for each of the :math:`4^{\text{number qubits}}` Pauli operators.
    """
    num_paulis = 4**num_qubits
    weights = np.zeros(num_paulis, dtype=int)
    for i in range(num_paulis):
        weight = 0
        temp_i = i
        # Convert index to base 4, number of non-zero digits is the weight
        for _ in range(num_qubits):
            if temp_i % 4 != identity_pos:  # Pauli I corresponds to digit "identity_pos"
                weight += 1
            temp_i //= 4
        weights[i] = weight
    return weights

def average_pauli_weight(o: np.ndarray, weights: np.ndarray) -> float:
    r"""
    Calculates the average Pauli weight (influence) for an operator :math:`O`.

    The influence :math:`I` is calculated as

    .. math::

        O &= \sum_{P \in \text{all Pauli strings}} c_P P \\
        I(O) &= \sum_{P} |P| * c_{P}^{2}
    
    where :math:`|P|` is the weight of the Pauli operator :math:`P`.

    .. currentmodule::
        paulie.application.matrix_decomposition

    Args:
        o (numpy.ndarray): The operator :math:`O` in matrix form.
        weights (numpy.ndarray): The weights of each Pauli operator in the ordering corresponding to
            the output of :func:`matrix_decomposition`.

    Returns:
        float: The average Pauli weight (influence) of the operator.
    """
    # Get the coefficients c_P from the Pauli decomposition
    coeffs = matrix_decomposition(o)
    # For a Hermitian operator O, the coefficients c_P are real.
    # The "probability" of a Pauli term P is c_P^2.
    # Note: sum(|c_P|^2) = 1 due to O^2=I.
    probs = np.abs(coeffs)**2
    # Calculate the influence I(O)
    influence = np.sum(weights * probs)
    return influence
