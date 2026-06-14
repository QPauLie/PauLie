"""
Basis matrix generators for u(1), so(n), su(n), sp(n) and their direct sums.
"""

import numpy as np

from paulie.classifier.types import TypeAlgebra


def get_u_basis(n: int) -> np.ndarray:
    """Returns the basis of u(1).

    Args:
        n (int): The dimension of the unitary algebra. Must be 1.

    Returns:
        np.ndarray: A 3D array of shape (1, 1, 1) representing the u(1) basis.

    Raises:
        ValueError: If n is not 1.
    """
    if n != 1:
        raise ValueError("u(n) basis: Only n=1 is currently supported.")
    return np.array([[[1j]]], dtype=np.complex128)


def get_n_so_basis(n: int) -> int:
    """Calculates the number of basis elements for so(n).

    Args:
        n (int): The dimension of the orthogonal algebra.

    Returns:
        int: Number of basis elements, calculated as n * (n - 1) / 2.

    Raises:
        TypeError: If n is not an integer.
        ValueError: If n <= 0.
    """
    err_msg = "so(n) basis: n must be a positive integer."
    if not isinstance(n, int):
        raise TypeError(err_msg)
    if n <= 0:
        raise ValueError(err_msg)
    return n * (n - 1) // 2


def get_so_basis(n: int) -> np.ndarray:
    """Generates the basis for so(n).

    Args:
        n (int): The dimension of the orthogonal algebra.

    Returns:
        np.ndarray: A 3D array of shape (dim, n, n) containing the basis matrices
            constructed from antisymmetric E_ij - E_ji generators.
    """
    dim: int = get_n_so_basis(n)
    basis = np.zeros((dim, n, n), dtype=np.complex128)
    rows, cols = np.triu_indices(n, k=1)
    k = np.arange(dim)
    basis[k, rows, cols] = 1.0
    basis[k, cols, rows] = -1.0

    return basis


def get_n_su_basis(n: int) -> int:
    """Calculates the number of basis elements for su(n).

    Args:
        n (int): The dimension of the special unitary algebra.

    Returns:
        int: Number of basis elements, calculated as n^2 - 1.

    Raises:
        TypeError: If n is not an integer.
        ValueError: If n <= 2^2.
    """
    err_msg = "su(n) basis: n must be integer(2^x) and greater than 2^2"
    if not isinstance(n, int):
        raise TypeError(err_msg)
    if n <= 2**2 or n.bit_count()!=1:
        raise ValueError(err_msg)
    return n**2 - 1


def get_su_basis(n: int) -> np.ndarray:
    """Generates the skew-Hermitian basis for su(n).

    Args:
        n (int): The dimension of the special unitary algebra.

    Returns:
        np.ndarray: A 3D array of shape (dim, n, n) containing the basis matrices:
            Symmetric Imaginary -> Anti-symmetric Real -> Diagonal Imaginary
    """
    dim: int = get_n_su_basis(n)
    basis = np.zeros((dim, n, n), dtype=np.complex128)

    # Number of symmetric basis is same with anti-symmetric basis
    n_symmetric: int = (n * (n - 1)) // 2
    n_diagonal: int = n - 1
    rows, cols = np.triu_indices(n, k=1)

    # Symmetric Imaginary Generators
    k = np.arange(n_symmetric)
    basis[k, rows, cols] = 1j
    basis[k, cols, rows] = 1j

    # Anti-symmetric Real Generator
    k = np.arange(n_symmetric, n_symmetric * 2)
    basis[k, rows, cols] = 1.0
    basis[k, cols, rows] = -1.0

    # Diagonal Imaginary Generator
    factors = np.asarray([np.sqrt(2 / (f * (f + 1))) 
                          for f in range(1, n_diagonal + 1)
                        ])

    i_inds, j_inds = np.tril_indices(n_diagonal)
    dim_inds = (n_symmetric*2) + i_inds
    basis[dim_inds, j_inds, j_inds] = 1j * factors[i_inds]


    i_arr = np.arange(n_diagonal)
    dim_inds = (n_symmetric *2) + i_arr
    diag_tails = i_arr + 1
    basis[dim_inds, diag_tails, diag_tails] = -1j * diag_tails * factors[i_arr]

    return basis


def get_n_sp_basis(n: int) -> int:
    """Calculates the number of basis elements for sp(n).

    Args:
        n (int): The dimension parameter of the symplectic algebra.

    Returns:
        int: Number of basis elements, calculated as n * (2n + 1).

    Raises:
        ValueError: If n is not a power of 2 or < 1.
    """
    err_msg = "sp(n) basis: n must be power of 2 AND greater or equal to 2 (at least 2^1)"
    if not isinstance(n, int):
        raise TypeError(err_msg)
    if (2 * n).bit_count() != 1 or n < 1:
        raise ValueError(err_msg)
    return n * (2 * n + 1)


def get_sp_basis(n: int) -> np.ndarray:
    """
    Generates the basis for the compact symplectic Lie algebra usp(n).
    Defining condition: X^T J + J X = 0
    Symplectic form: J = [[0, I_N], [-I_N, 0]]
    Matrix structure: X = [[A, B], [-conj(B), conj(A)]]
    where:
      - A is n x n skew-Hermitian (A^dagger = -A)
      - B is n x n complex symmetric (B^T = B)
      - conj(M) denotes the entry-wise complex conjugate

    Args:
        n (int): The dimension parameter of the symplectic algebra.

    Returns:
        np.ndarray: A 3D array of shape (dim, 2n, 2n) containing the basis matrices:
            k (Unitary Subgroup) {Block A: Real Anti-symmetric -> Block A: Imaginary Symmetric 
            -> Block A: Imaginary Diagonal}
            p or a a Off-diagonal Symplectic {Block B: Real Symmetric 
            -> Block B: Imaginary Symmetric 
            -> Block B: Imaginary Diagonal}
            a (Cartan Subgroup) {Block B: Real Diagonal}
    """
    dim: int = get_n_sp_basis(n)
    basis = np.zeros((dim, 2 * n, 2 * n), dtype=np.complex128)

    # nA_RealSymmetric = nA_ImaginarySymmetric = nB_RealSymmetric = nB_ImaginarySymmetric
    n_Sym = (n**2 - n) // 2
    # nB_ImaginaryDiagonal = nB_RealDiagonal = nA_DiagonalImaginary
    n_Diag = n

    # Prepare Sub Block
    zeros_block = np.zeros((n, n), dtype=np.complex128)
    sub_rows, sub_cols = np.triu_indices(n, 1)
    k_sym = np.arange(n_Sym)

    # E_{jk} - E_{kj}
    E1 = np.zeros((n_Sym, n, n), dtype=np.complex128)
    E1[k_sym, sub_rows, sub_cols] = 1.0
    E1[k_sym, sub_cols, sub_rows] = -1.0

    # E_{jk} + E_{kj}
    E2 = np.zeros((n_Sym, n, n), dtype=np.complex128)
    E2[k_sym, sub_rows, sub_cols] = 1.0
    E2[k_sym, sub_cols, sub_rows] = 1.0

    # E_{jj}
    k_diag = np.arange(n_Diag)
    E3 = np.zeros((n_Diag, n, n), dtype=np.complex128)
    di = np.arange(n)
    E3[k_diag, di, di] = 1.0

    # Assign Basis
    for i in range(n_Sym):
        # Block A: Real Anti-symmetric => k (Unitary Subgroup)
        basis[i] = np.block([[E1[i], zeros_block], [zeros_block, E1[i]]])

        # Block A: Imaginary Symmetric => k (Unitary Subgroup)
        basis[n_Sym + i] = np.block([[1j * E2[i], zeros_block], [zeros_block, -1j * E2[i]]])

        # Block B: Real Symmetric => p \ a Off-diagonal Symplectic
        basis[n_Sym * 2 + n_Diag + i] = np.block([
                    [zeros_block, E2[i]],
                    [-1 * E2[i], zeros_block]
            ])
        # Block B: Imaginary Symmetric => p \ a Off-diagonal Symplectic
        basis[n_Sym * 3 + n_Diag + i] = np.block([
                [zeros_block, 1j * E2[i]],
                [1j * E2[i], zeros_block]
            ])

    for i in range(n_Diag):
        # Block A: Imaginary Diagonal
        basis[n_Sym * 2 + i] = np.block([[1j * E3[i], zeros_block], [zeros_block, -1j * E3[i]]])

        # Block B: Imaginary Diagonal => p \ a Off-diagonal Symplectic
        basis[n_Sym * 4 + n_Diag + i] = np.block([
                [zeros_block, 1j * E3[i]],
                [1j * E3[i], zeros_block]
            ])

        # Block B: Real Diagonal => a (Cartan Subalgebra)
        basis[n_Sym * 4 + n_Diag * 2 + i] = np.block([
                [zeros_block, E3[i]],
                [-1 * E3[i], zeros_block]
            ])

    return basis


def get_group_basis(algebra_type: 'TypeAlgebra', n: int) -> np.ndarray:
    """Dispatcher to generate basis based on algebra type.

    Args:
        algebra_type (TypeAlgebra): The enum representing the Lie Algebra (U, SO, SU, SP).
        n (int): The dimension parameter for the algebra.

    Returns:
        np.ndarray: The generated basis matrices.

    Raises:
        ValueError: If algebra_type is unsupported.
    """
    match algebra_type:
        case TypeAlgebra.U:
            return get_u_basis(n)
        case TypeAlgebra.SO:
            return get_so_basis(n)
        case TypeAlgebra.SU:
            return get_su_basis(n)
        case TypeAlgebra.SP:
            return get_sp_basis(n)
    err_msg = f"Unsupported algebra type: {algebra_type}"
    raise ValueError(err_msg)


def get_n_basis(algebra_type: 'TypeAlgebra', n: int) -> int:
    """Calculates total dimension (number of generators) for a given algebra.

    Args:
        algebra_type: The enum representing the Lie Algebra.
        n: The dimension parameter.

    Returns:
        int: Total number of generators.
    """
    match algebra_type:
        case TypeAlgebra.U:
            return 1
        case TypeAlgebra.SO:
            return get_n_so_basis(n)
        case TypeAlgebra.SU:
            return get_n_su_basis(n)
        case TypeAlgebra.SP:
            return get_n_sp_basis(n)
    err_msg = f"Unsupported algebra type: {algebra_type}"
    raise ValueError(err_msg)


def get_algebras_basis(
        multipliers: list[int],
        groups: list['TypeAlgebra'],
        sizes: list[int]
    ) -> np.ndarray:
    """Constructs the block-diagonal basis for a direct sum of Lie algebras.

    This function pre-allocates a 3D array of shape (total_dim, total_n, total_n)
    and populates it by injecting the textbook basis matrices into the diagonal
    blocks.

    Args:
        multipliers (list[int]): A list of how many times each algebra is repeated.
        groups (list[TypeAlgebra]): A list of TypeAlgebra enums.
        sizes (list[int]): A list of integers defining the size (n) for each algebra.

    Returns:
        np.ndarray: The concatenated basis tensor of shape (total_dim, total_n, total_n).

    Raises:
        ValueError: If the input lists have inconsistent lengths.
    """
    n_pairs = len(multipliers)
    if not n_pairs == len(groups) == len(sizes):
        raise ValueError("Each input params should have same length")

    block_sizes = []
    n_basis_per_group = []
    total_dim = 0
    total_n = 0

    for i in range(n_pairs):
        actual_block_size = sizes[i] * 2 if groups[i] == TypeAlgebra.SP else sizes[i]
        block_sizes.append(actual_block_size)

        n_basis = get_n_basis(groups[i], sizes[i])
        n_basis_per_group.append(n_basis)

        total_dim += n_basis * multipliers[i]
        total_n += actual_block_size * multipliers[i]

    basis = np.zeros((total_dim, total_n, total_n), dtype=np.complex128)

    # Use rolling offsets to slice directly into memory
    gen_offset = 0  # Tracks depth in the 3D tensor
    diag_offset = 0  # Tracks slide down the 2D diagonal

    for i in range(n_pairs):
        base_matrices = get_group_basis(groups[i], sizes[i])

        # Stamp this basis out for however many copies (multipliers) we need
        for _ in range(multipliers[i]):
            n_gen = n_basis_per_group[i]
            b_size = block_sizes[i]

            # Vectorized assignment: Drop the entire chunk of matrices into the correct slice
            basis[
                gen_offset : gen_offset + n_gen,
                diag_offset : diag_offset + b_size,
                diag_offset : diag_offset + b_size
            ] = base_matrices

            gen_offset += n_gen
            diag_offset += b_size

    return basis
