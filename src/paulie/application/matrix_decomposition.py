"""
    Decompose a 2^n x 2^n complex matrix into the Pauli basis using 
    PauliString bitarray representation.
"""

import numpy as np
from typing import Dict
from paulie.common.pauli_string_bitarray import PauliString

def matrix_decomposition(matrix: np.ndarray) -> Dict[PauliString, complex]:
    """
    Decompose a 2^n x 2^n complex matrix into the Pauli basis using PauliString bitarray representation.

    Parameters
    ----------
    matrix : np.ndarray
        A square complex matrix of dimension 2^n x 2^n.

    Returns
    -------
    Dict[PauliString, complex]
        A mapping from each n-qubit PauliString to its coefficient in the decomposition.
    """
    # Validate input
    dim = matrix.shape[0]
    if matrix.shape[0] != matrix.shape[1] or (dim & (dim - 1)) != 0:
        raise ValueError("Matrix must be square with dimension 2^n.")

    # Number of qubits
    n = int(np.log2(dim))
    # Normalization factor for trace inner product
    normalization = float(dim)

    coeffs: Dict[PauliString, complex] = {}
    # Iterate over all Pauli strings of length n
    for ps in PauliString.all(n):
        # Convert bitarray representation to its matrix form
        P = ps.to_matrix()
        # Compute coefficient: Tr(P^† A) / 2^n
        coeff = np.trace(np.conjugate(P).T @ matrix) / normalization
        # Discard near-zero coefficients
        if not np.isclose(coeff, 0.0):
            coeffs[ps] = coeff

    return coeffs