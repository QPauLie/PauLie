from paulie.application.matrix_decomposition import matrix_decomposition
from paulie.common.pauli_string_bitarray import PauliString
import numpy as np
from itertools import product

def quantum_fourier_entropy(O: np.ndarray) -> float:
    """
    Calculates the quantum Fourier entropy of an operator O.
    H(O) = -sum_P c_P**2 * log(c_P**2)
    """
    # Get the coefficients c_P from the Pauli decomposition
    c_P = matrix_decomposition(O)

    # Calculate the probabilities p_P = c_P^2
    probs = np.abs(c_P)**2

    # Filter out zero probabilities to avoid log(0)
    non_zero_probs = probs[probs > 1e-12]

    # Calculate the Shannon entropy using base 2 for the logarithm
    entropy = -np.sum(non_zero_probs * np.log2(non_zero_probs))
    return entropy


def avg_pauli_weights(O: np.ndarray) -> np.ndarray:
    """
    Calculate the average Pauli weights of an operator O.
    I(O) = sum_P |P| * c_P**2
    """
    # Get the coefficients c_P from the Pauli decomposition
    c_Ps = matrix_decomposition(O)
    # Get the number of qubits from the matrix decomposition
    # 4 Pauli matrices yield 4**n_qubits options
    dim = c_Ps.shape[0]
    n_qubits = int(np.emath.logn(4, dim))
    # get all options of Pauli Strings for n_qubits
    pauli_strings = product('IXYZ', repeat=n_qubits)
    I = 0
    for pauli_str in pauli_strings:
        pl = PauliString(pauli_str=pauli_str)
        c_P = pl.get_weight_in_matrix(c_Ps)
        abs_P = pl.get_count_non_trivially()
        I += abs(c_P) ** 2 * abs_P
    return I


def avg_pauli_weights_from_strings(O: np.ndarray, pauli_strings: list) -> np.ndarray:
    """
    Calculate the average Pauli weights of an operator O, given a list of Pauli strings.: This is useful to reduce the calculation, when testing.
    I(O) = sum_P |P| * c_P**2
    """
    # Get the coefficients c_P from the Pauli decomposition
    c_Ps = matrix_decomposition(O)
    # get all options of Pauli Strings for n_qubits
    I = 0
    for pauli_str in pauli_strings:
        pauli_str = ''.join(pauli_str)
        pl = PauliString(pauli_str=pauli_str)
        c_P = pl.get_weight_in_matrix(c_Ps)
        abs_P = pl.get_count_non_trivially()
        I += abs(c_P) ** 2 * abs_P
    return I



##### Alternate Implementation #####


def get_pauli_weights(num_qubits: int, identity_pos: int=0) -> np.ndarray:
    """
    Generates the weight |P| for each of the 4**num_qubits Pauli operators.
    The weight is the number of non-identity terms in the Pauli string.
    The ordering corresponds to the output of matrix_decomposition, default is 'I' at position 0.
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

def average_pauli_weight(O: np.ndarray, weights: np.ndarray) -> float:
    """
    Calculates the average Pauli weight (influence) for an operator O.
    I(O) = sum_P |P| * c_P**2
    """

    # Get the coefficients c_P from the Pauli decomposition
    coeffs = matrix_decomposition(O)
    # For a Hermitian operator O, the coefficients c_P are real.
    # The "probability" of a Pauli term P is c_P^2.
    # Note: sum(|c_P|^2) = 1 due to O^2=I.
    probs = np.abs(coeffs)**2
    
    # Calculate the influence I(O)
    influence = np.sum(weights * probs)
    return influence