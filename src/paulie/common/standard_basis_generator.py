import numpy as np

def generate_so_basis(N) -> list[np.ndarray]:
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

    A list containing the basis matrices each of size N * N. The total number of matrices
        returned will be exactly N * (N - 1) / 2.

    Parameters:
    -----------
    N : int
        The dimension of the space (must be >= 2).

    Returns:
        list[np.ndarray]
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

def generate_su_basis(N) -> list[np.ndarray]:
    """
    Generates the standard basis for the Lie algebra su(2^N) using the 
    Tensor Product (Pauli String) construction.

    Convention Used:
    he returned list of basis generators follows a strict lexicographical 
    ordering based on the binary symplectic representation of Pauli matrices.
    
    Each single-qubit Pauli operator is mapped to a 2-bit binary number:
        00 (0) -> Identity (I)
        01 (1) -> Pauli-Z (Z)
        10 (2) -> Pauli-X (X)
        11 (3) -> Pauli-Y (Y)

    An N-qubit Pauli string can therefore be interpreted as a 2N-bit integer. 
    Because we generate the basis by counting upwards in binary (from 1 to 
    4^N - 1), the basis is returned in the exact order of these growing 
    binary integers. (The integer 0 is exactly the all-Identity string, 
    which is sliced off to enforce the traceless condition of su).

    Example Sequence for N = 2:
    ---------------------------
    The 4^2 - 1 = 15 matrices will be returned in the following order 
    (where each string is implicitly multiplied by the imaginary unit 'i' 
    to enforce anti-Hermiticity):

        1.  IZ  (00 01)
        2.  IX  (00 10)
        3.  IY  (00 11)
        4.  ZI  (01 00)
        5.  ZZ  (01 01)
        6.  ZX  (01 10)
        7.  ZY  (01 11)
        8.  XI  (10 00)
        9.  XZ  (10 01)
        10. XX  (10 10)
        11. XY  (10 11)
        12. YI  (11 00)
        13. YZ  (11 01)
        14. YX  (11 10)
        15. YY  (11 11)
    Parameters:
        N (int): The number of subsyste
        The number of subsystems (e.g., qubits). The resulting matrices will 
        be of dimension 2^N x 2^N.

    Returns:
        list[np.ndarray]
        
    """
    # 1. Generate the full u(2^N) basis
    from paulie.common.pauli_string_factory  import gen_all_pauli_strings 
    all_pauli_strings = gen_all_pauli_strings(N)  # This generates all 4^N Pauli strings as PauliString objects
    basis = []
    for pauli_string in all_pauli_strings:
        basis.append(1j * pauli_string.get_matrix())  # Multiply by 1j to ensure anti-Hermiticity
    
    # 2. Slice off the first element (i * Identity) to enforce Tr(X) = 0
    return basis[1:]    



def generate_sp_basis(N) -> list[np.ndarray]:
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
