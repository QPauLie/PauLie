import numpy as np
import itertools
from functools import reduce

def generate_soN_basis(N):
    """
    Generates the standard basis for the Lie algebra so(N).

    The Lie algebra so(N) consists of all N x N skew-symmetric matrices
    (matrices where X^T = -X).

    Convention Used:
    ----------------
    Mathematically, the standard basis is defined by the set:
        { E_ij - E_ji : 1 <= i < j <= N }
    where E_ij is an N x N matrix with a 1 in the i-th row and j-th column,
    and 0 everywhere else.

    Because Python uses 0-based indexing, this function translates the 
    mathematical bounds (1 to N) into Python bounds (0 to N-1). 
    The condition becomes:
        { E_ij - E_ji : 0 <= i < j < N }

    The algorithm strictly iterates through the upper-triangular indices 
    (where i < j). For each valid (i, j) pair, it:
      1. Initializes an N x N matrix of zeros.
      2. Sets the (i, j) element to 1.
      3. Sets the (j, i) element to -1 to satisfy skew-symmetry.

    Parameters:
    -----------
    N : int
        The dimension of the space (must be >= 2).

    Returns:
    --------
    list of numpy.ndarray
        A list containing the basis matrices. The total number of matrices
        returned will be exactly N * (N - 1) // 2.
    """
    if N < 2:
        raise ValueError("Dimension N must be at least 2.")

    basis = []
    
    # Iterate over rows (i)
    for i in range(N):
        # Iterate over columns strictly to the right of the diagonal (j)
        for j in range(i + 1, N):
            # Initialize an N x N zero matrix
            matrix = np.zeros((N, N), dtype=int)
            
            # Apply the E_ij - E_ji convention
            matrix[i, j] = 1
            matrix[j, i] = -1
            
            basis.append(matrix)
            
    return basis

# Example usage:


def generate_su2N_pauli_basis(N):
    """
    Generates the standard basis for the Lie algebra su(2^N) using the 
    Tensor Product (Pauli String) construction.

    This approach is foundational in quantum error correction and the 
    stabilizer formalism. Rather than building massive generalized 
    Gell-Mann matrices, it constructs the algebra's generators as 
    N-fold Kronecker products of the fundamental SU(2) Pauli matrices.

    Mathematical Convention:
    ------------------------
    1. Base Matrices:
       We use the standard Hermitian Pauli matrices, plus the Identity.
       I = [[1, 0], [0,  1]]
       X = [[0, 1], [1,  0]]
       Y = [[0, -1j], [1j, 0]]
       Z = [[1, 0], [0, -1]]

    2. Tensor Construction:
       A raw basis element is formed by taking the Kronecker product of N 
       Pauli matrices: P_1 (kron) P_2 (kron) ... (kron) P_N.
       We generate all 4^N possible strings using a Cartesian product.

    3. Lie Algebra Constraints for su(d):
       - Anti-Hermitian (X^dagger = -X): 
         The standard Pauli matrices (and their real tensor products) are 
         Hermitian. To map them to the compact Lie algebra su(2^N), we 
         multiply the final tensor product by the imaginary unit '1j'.
       - Traceless (Tr(X) = 0): 
         The trace of a tensor product is the product of the individual traces. 
         Because Tr(X) = Tr(Y) = Tr(Z) = 0, any string containing at least one 
         X, Y, or Z will have a trace of exactly 0. We must explicitly filter 
         out the all-Identity string (I (kron) I ... (kron) I), which has a 
         trace of 2^N.

    Parameters:
    -----------
    N : int
        The number of subsystems (e.g., qubits). The resulting matrices will 
        be of dimension 2^N x 2^N.

    Returns:
    --------
    list of numpy.ndarray
        A list containing the exactly (4^N - 1) basis matrices for su(2^N).
    """
    # 1. Generate the full u(2^N) basis
    full_basis = generate_u2N_pauli_basis(N)
    
    # 2. Slice off the first element (i * Identity) to enforce Tr(X) = 0
    return full_basis[1:]    

def generate_u2N_pauli_basis(N):
    """
    Generates the standard basis for the unitary Lie algebra u(2^N).
    
    Convention Used:
    ----------------
    This uses the same tensor product construction as su(2^N), mapping 
    to N-fold Pauli strings. 
    
    Because u(2^N) does NOT require generators to be traceless, we simply 
    include all 4^N combinations of {I, X, Y, Z}, including the all-Identity 
    string. Multiplying the all-Identity string by 1j provides the i*I 
    matrix, which mathematically satisfies the u(1) direct sum component.
    
    The dimension returned will be exactly 4^N.
    """
    if N < 1:
        raise ValueError("Dimension exponent N must be at least 1.")

    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    paulis = [I, X, Y, Z]
    basis = []

    # Iterate through ALL 4^N possible Pauli strings
    for string in itertools.product(paulis, repeat=N):
        
        tensor_product = reduce(np.kron, string)
        
        # Enforce Anti-Hermiticity
        basis_matrix = 1j * tensor_product
        
        basis.append(basis_matrix)
        
    return basis



def generate_spN_basis(N):
    """
    Generates the standard basis for the compact symplectic Lie algebra sp(N).
    
    The Lie algebra sp(N) consists of 2N x 2N complex matrices that satisfy 
    two independent constraints simultaneously:
    1. Symplectic: X^T * J + J * X = 0  (where J is the standard symplectic matrix)
    2. Anti-Hermitian: X^dagger + X = 0
    
    Convention Used:
    ----------------
    These constraints force the 2N x 2N matrix to have the rigid block structure:
        X = [[ A,  -conj(C) ],
             [ C,   conj(A) ]]
             
    Where A and C are N x N blocks:
    - A is an anti-Hermitian matrix.
    - C is a complex symmetric matrix.
    
    We generate the basis by isolating the independent real and imaginary 
    parameters of A and C into 5 distinct subsets:
    
    Part 1: Generators from the A block
      Subset 1: Imaginary diagonals (A = i * E_kk)
      Subset 2: Real skew-symmetric off-diagonals (A = E_jk - E_kj)
      Subset 3: Imaginary symmetric off-diagonals (A = i * (E_jk + E_kj))
      
    Part 2: Generators from the C block
      Subset 4: Real symmetric elements (C = E_jk + E_kj)
      Subset 5: Imaginary symmetric elements (C = i * (E_jk + E_kj))
      
    Parameters:
    -----------
    N : int
        The block dimension. The resulting matrices will be 2N x 2N.
        
    Returns:
    --------
    list of numpy.ndarray
        A list containing the basis matrices. The total number of matrices
        will be exactly N * (2N + 1).
    """
    if N < 1:
        raise ValueError("Block dimension N must be at least 1.")

    dim = 2 * N
    basis = []

    # Helper function to create an empty 2N x 2N complex matrix
    def empty_matrix():
        return np.zeros((dim, dim), dtype=complex)

    # ---------------------------------------------------------
    # PART I: Generators from the A block (anti-Hermitian)
    # A goes in top-left [0:N, 0:N]
    # conj(A) goes in bottom-right [N:2N, N:2N]
    # ---------------------------------------------------------

    # Subset 1: Imaginary Diagonals
    for k in range(N):
        M = empty_matrix()
        M[k, k] = 1j
        M[k + N, k + N] = -1j  # conj(1j) = -1j
        basis.append(M)

    # Subsets 2 & 3: Off-Diagonals
    for j in range(N):
        for k in range(j + 1, N):
            
            # Subset 2: Real Skew-Symmetric Off-Diagonals
            M2 = empty_matrix()
            M2[j, k] = 1
            M2[k, j] = -1
            M2[j + N, k + N] = 1   # conj(1) = 1
            M2[k + N, j + N] = -1
            basis.append(M2)
            
            # Subset 3: Imaginary Symmetric Off-Diagonals
            M3 = empty_matrix()
            M3[j, k] = 1j
            M3[k, j] = 1j
            M3[j + N, k + N] = -1j # conj(1j) = -1j
            M3[k + N, j + N] = -1j
            basis.append(M3)

    # ---------------------------------------------------------
    # PART II: Generators from the C block (complex symmetric)
    # C goes in bottom-left [N:2N, 0:N]
    # -conj(C) goes in top-right [0:N, N:2N]
    # ---------------------------------------------------------

    # Subsets 4 & 5: Complex Symmetric Elements
    for j in range(N):
        # We start k from j to include the diagonal elements of C
        for k in range(j, N):
            
            # Subset 4: Real Symmetric
            M4 = empty_matrix()
            M4[j + N, k] = 1
            M4[k + N, j] = 1
            M4[j, k + N] = -1      # -conj(1) = -1
            M4[k, j + N] = -1
            basis.append(M4)
            
            # Subset 5: Imaginary Symmetric
            M5 = empty_matrix()
            M5[j + N, k] = 1j
            M5[k + N, j] = 1j
            M5[j, k + N] = 1j      # -conj(1j) = -(-1j) = 1j
            M5[k, j + N] = 1j
            basis.append(M5)

    return basis
