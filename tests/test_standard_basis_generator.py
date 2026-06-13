"""
Unit tests for the basis generators in paulie.common.standard_basis_generator.py.
These tests verify that the generated bases have the correct size, element properties, and Lie closure properties for so(N), su(2^N), and sp(N) algebras.
"""

import numpy as np
import pytest
from paulie.common.standard_basis_generator import generate_so_basis, generate_su_basis, generate_sp_basis

@pytest.mark.parametrize("N", [2, 3, 4, 5, 10])
def test_soN_basis_size(N):
    """1. Ensure it returns the correct size of the basis (N*(N-1)/2)."""
    basis = generate_so_basis(N)
    expected_size = N * (N - 1) // 2
    
    assert len(basis) == expected_size, \
        f"Expected {expected_size} basis elements for N={N}, but got {len(basis)}"

@pytest.mark.parametrize("N", [2, 3, 4, 5])
def test_soN_element_size_and_symmetry(N):
    """2. Ensure it returns the correct element size and fundamental properties."""
    basis = generate_so_basis(N)
    
    for idx, mat in enumerate(basis):
        # Check that the matrix is N x N
        assert mat.shape == (N, N), \
            f"Matrix at index {idx} has wrong shape: {mat.shape}. Expected: {(N, N)}"
        
        # Check skew-symmetry (X^T = -X), which is the definition of so(N)
        assert np.array_equal(mat.T, -mat), \
            f"Matrix at index {idx} is not skew-symmetric."

@pytest.mark.parametrize("N", [2, 3, 4, 5])
def test_soN_lie_closure(N):
    """
    3. Ensure it spans the Lie closure.
    
    A basis spans the Lie closure if the Lie bracket (commutator) of ANY 
    two basis elements [A, B] = AB - BA can be expressed as a linear 
    combination of the basis elements themselves.
    """
    basis = generate_so_basis(N)
    num_elements = len(basis)
    
    # Flatten basis matrices into column vectors to form a linear algebra system
    # Matrix shape will be: (N^2, num_basis_elements)
    basis_vectors = np.array([mat.flatten() for mat in basis]).T
    
    # Test every pair of basis elements
    for i in range(num_elements):
        for j in range(i, num_elements):
            A = basis[i]
            B = basis[j]
            
            # Compute the commutator [A, B] = AB - BA
            # Note: We use np.dot or '@' for matrix multiplication
            commutator = A @ B - B @ A
            
            # To prove the commutator is in the span of the basis, we solve the
            # linear equation: basis_vectors * coefficients = flattened_commutator
            flat_comm = commutator.flatten()
            
            # np.linalg.lstsq finds the best-fit coefficients. 
            coeffs, residuals, rank, s = np.linalg.lstsq(basis_vectors, flat_comm, rcond=None)
            
            # Reconstruct the matrix using the found coefficients
            reconstructed_flat = basis_vectors @ coeffs
            reconstructed_matrix = reconstructed_flat.reshape(N, N)
            
            # If the basis spans the closure, the reconstructed matrix will 
            # perfectly match the commutator.
            assert np.allclose(commutator, reconstructed_matrix, atol=1e-10), \
                f"Commutator of basis[{i}] and basis[{j}] is not in the span of the basis!"
            
@pytest.mark.parametrize("N", [1, 2, 3])
def test_su2N_basis_size(N):
    """1. Ensure it returns the correct size of the basis (4^N - 1)."""
    basis = generate_su_basis(N)
    
    # Mathematical dimension of su(2^N) is (2^N)^2 - 1 = 4^N - 1
    expected_size = 4**N - 1
    
    assert len(basis) == expected_size, \
        f"Expected {expected_size} basis elements for N={N}, but got {len(basis)}"

@pytest.mark.parametrize("N", [1, 2, 3])
def test_su2N_element_properties(N):
    """2. Ensure it returns the correct element size and fundamental properties."""
    basis = generate_su_basis(N)
    dim = 2**N
    
    for idx, mat in enumerate(basis):
        # Check that the matrix is 2^N x 2^N
        assert mat.shape == (dim, dim), \
            f"Matrix at index {idx} has wrong shape: {mat.shape}. Expected: {(dim, dim)}"
        
        # Check Anti-Hermitian (X^dagger = -X)
        assert np.allclose(mat.conj().T, -mat), \
            f"Matrix at index {idx} is not anti-Hermitian."
            
        # Check Traceless (Tr(X) = 0)
        assert np.isclose(np.trace(mat), 0.0, atol=1e-10), \
            f"Matrix at index {idx} is not traceless. Trace: {np.trace(mat)}"

@pytest.mark.parametrize("N", [1, 2])
def test_su2N_lie_closure(N):
    """
    3. Ensure it spans the Lie closure.
    
    Because su(2^N) is a Lie algebra over the REAL numbers, the commutator 
    [A, B] must be a REAL linear combination of the basis elements. 
    To test this robustly with numpy, we split the complex matrices into 
    stacked real and imaginary vectors to force a real-coefficient solution.
    """
    basis = generate_su_basis(N)
    num_elements = len(basis)
    
    # Flatten the basis matrices and split into real and imaginary parts
    real_parts = np.array([mat.real.flatten() for mat in basis])
    imag_parts = np.array([mat.imag.flatten() for mat in basis])
    
    # FIX: Use hstack so the arrays are glued side-by-side, then transpose
    # Shape becomes (2 * 4^N, num_elements)
    basis_vectors = np.hstack([real_parts, imag_parts]).T
    
    # Test pairs of basis elements
    for i in range(num_elements):
        for j in range(i + 1, num_elements):
            A = basis[i]
            B = basis[j]
            
            # Compute the commutator [A, B] = AB - BA
            commutator = A @ B - B @ A
            
            # Split the commutator exactly as we split the basis
            flat_comm = np.hstack([commutator.real.flatten(), commutator.imag.flatten()])
            
            # If they commute (commutator is 0), it trivially spans. Skip solver.
            if np.allclose(flat_comm, 0):
                continue
            
            # Solve for the linear coefficients
            coeffs, residuals, rank, s = np.linalg.lstsq(basis_vectors, flat_comm, rcond=None)
            
            # Reconstruct the vector using the found coefficients
            reconstructed_flat = basis_vectors @ coeffs
            
            # If the basis spans the closure, the reconstruction will perfectly match
            assert np.allclose(flat_comm, reconstructed_flat, atol=1e-10), \
                f"Commutator of basis[{i}] and basis[{j}] is not in the REAL span of the basis!"            


@pytest.mark.parametrize("N", [1, 2, 3, 4])
def test_spN_basis_size(N):
    """1. Ensure it returns the correct size of the basis: N(2N + 1)."""
    basis = generate_sp_basis(N)
    
    expected_size = N * (2 * N + 1)
    
    assert len(basis) == expected_size, \
        f"Expected {expected_size} basis elements for N={N}, but got {len(basis)}"

@pytest.mark.parametrize("N", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
def test_spN_element_properties(N):
    """
    2. Ensure it returns the correct element size and fundamental properties.
    Specifically, testing Anti-Hermiticity and the Symplectic condition.
    """
    basis = generate_sp_basis(N)
    dim = 2 * N
    
    # Construct the standard symplectic matrix J = [[0, I], [-I, 0]]
    I_N = np.eye(N)
    J = np.block([
        [np.zeros((N, N)), I_N],
        [-I_N, np.zeros((N, N))]
    ])
    
    for idx, mat in enumerate(basis):
        # Check that the matrix is 2N x 2N
        assert mat.shape == (dim, dim), \
            f"Matrix at index {idx} has wrong shape: {mat.shape}. Expected: {(dim, dim)}"
        
        # Check Anti-Hermitian (X^dagger = -X)
        assert np.allclose(mat.conj().T, -mat), \
            f"Matrix at index {idx} is not anti-Hermitian."
            
        # Check Symplectic (X^T J + J X = 0)
        symplectic_condition = mat.T @ J + J @ mat
        assert np.allclose(symplectic_condition, 0), \
            f"Matrix at index {idx} does not satisfy the symplectic condition."

@pytest.mark.parametrize("N", [1, 2, 3, 4, ])
def test_spN_lie_closure(N):
    """
    3. Ensure it spans the Lie closure.
    
    Like su(N), sp(N) is a Lie algebra over the REAL numbers. We split the 
    complex matrices into real and imaginary vectors and stack them side-by-side 
    (hstack) to force np.linalg.lstsq to find a strictly REAL coefficient solution.
    """
    basis = generate_sp_basis(N)
    num_elements = len(basis)
    
    # Flatten the basis matrices and split into real and imaginary parts
    real_parts = np.array([mat.real.flatten() for mat in basis])
    imag_parts = np.array([mat.imag.flatten() for mat in basis])
    
    # Glue them side-by-side and transpose to create the coefficient matrix
    basis_vectors = np.hstack([real_parts, imag_parts]).T
    
    # Test pairs of basis elements
    for i in range(num_elements):
        for j in range(i + 1, num_elements):
            A = basis[i]
            B = basis[j]
            
            # Compute the commutator [A, B] = AB - BA
            commutator = A @ B - B @ A
            
            # Split the commutator exactly as we split the basis
            flat_comm = np.hstack([commutator.real.flatten(), commutator.imag.flatten()])
            
            # If they commute (commutator is 0), it trivially spans. Skip solver.
            if np.allclose(flat_comm, 0):
                continue
            
            # Solve for the real linear coefficients
            coeffs, residuals, rank, s = np.linalg.lstsq(basis_vectors, flat_comm, rcond=None)
            
            # Reconstruct the vector using the found coefficients
            reconstructed_flat = basis_vectors @ coeffs
            
            # If the basis spans the closure, the reconstruction will perfectly match
            assert np.allclose(flat_comm, reconstructed_flat, atol=1e-10), \
                f"Commutator of basis[{i}] and basis[{j}] is not in the REAL span of the basis!"

@pytest.mark.parametrize("N", [1, 2, 3])
def test_su2N_linear_independence(N):
    """Ensure the generated basis elements are linearly independent over the reals."""
    basis = generate_su_basis(N)
    
    real_parts = np.array([mat.real.flatten() for mat in basis])
    imag_parts = np.array([mat.imag.flatten() for mat in basis])
    
    # Glue real and imaginary parts side-by-side to form strictly real vectors
    basis_vectors = np.hstack([real_parts, imag_parts])
    
    rank = np.linalg.matrix_rank(basis_vectors)
    
    assert rank == len(basis), \
        f"su(2^{N}) basis is not linearly independent over the reals! Expected rank {len(basis)}, got {rank}"

@pytest.mark.parametrize("N", [1, 2, 3, 4])
def test_spN_linear_independence(N):
    """Ensure the generated basis elements are linearly independent over the reals."""
    basis = generate_sp_basis(N)
    
    real_parts = np.array([mat.real.flatten() for mat in basis])
    imag_parts = np.array([mat.imag.flatten() for mat in basis])
    
    basis_vectors = np.hstack([real_parts, imag_parts])
    
    rank = np.linalg.matrix_rank(basis_vectors)
    
    assert rank == len(basis), \
        f"sp({N}) basis is not linearly independent over the reals! Expected rank {len(basis)}, got {rank}"


def get_structure_constants(basis, is_complex=False):
    """
    Computes the 3D tensor of structure constants C_{ij}^k for a given basis.
    [B_i, B_j] = sum_k (C_{ij}^k * B_k)
    """
    n = len(basis)
    
    # Flatten the basis to use a linear solver
    if is_complex:
        # For complex real-algebras like su(N), split real and imag
        real_parts = np.array([mat.real.flatten() for mat in basis])
        imag_parts = np.array([mat.imag.flatten() for mat in basis])
        basis_flat = np.hstack([real_parts, imag_parts]).T
    else:
        # For purely real algebras like so(N)
        basis_flat = np.array([mat.flatten() for mat in basis]).T

    # Initialize the 3D tensor for structure constants
    C = np.zeros((n, n, n))
    
    for i in range(n):
        for j in range(n):
            A = basis[i]
            B = basis[j]
            
            # Compute commutator
            comm = A @ B - B @ A
            
            if is_complex:
                flat_comm = np.hstack([comm.real.flatten(), comm.imag.flatten()])
            else:
                flat_comm = comm.flatten()
                
            # Solve for the coefficients: basis_flat * coeffs = flat_comm
            coeffs, _, _, _ = np.linalg.lstsq(basis_flat, flat_comm, rcond=None)
            
            # Store the coefficients in the tensor
            C[i, j, :] = coeffs
            
    return C

def test_su2_so3_exact_equivalence():
    """
    Proves the exact isomorphism su(2) == so(3) by showing their 
    structure constants are identical.
    """
    # 1. Generate the bases
    # su(2) corresponds to N=1 (1 qubit)
    # so(3) corresponds to N=3 (3D space)
    su2_basis = generate_su_basis(1) 
    so3_basis = generate_so_basis(3)
    
    # 2. Normalize the su(2) basis
    # Pauli matrices have an inherent factor of 2 in their commutation 
    # relations ([X, Y] = 2iZ). Standard elementary matrices do not. 
    # We scale the su(2) basis by 1/2 to align the geometries perfectly.
    su2_basis_normalized = [mat / 2.0 for mat in su2_basis]
    
    # 3. Compute structure constants
    C_su2 = get_structure_constants(su2_basis_normalized, is_complex=True)
    C_so3 = get_structure_constants(so3_basis, is_complex=False)
    
    # 4. Assert exact algebraic equivalence
    # Remarkably, because of the binary counting order of your su() generator 
    # and the upper-triangular loop of your so() generator, the bases are 
    # generated in the exact same matching sequence. No sorting is required!
    assert np.allclose(C_su2, C_so3, atol=1e-10), \
        "The structure constants of su(2) and so(3) do not match! They are not equivalent."
    
def test_su2_sp1_exact_equivalence():
    """
    Proves the exact isomorphism su(2) == sp(1) by showing their 
    structure constants are identical.
    """
    # 1. Generate the bases
    # su(2) corresponds to N=1
    # sp(1) corresponds to N=1
    su2_basis = generate_su_basis(1) 
    sp1_basis = generate_sp_basis(1)
    
    # 2. Compute structure constants (both contain complex matrices)
    C_su2 = get_structure_constants(su2_basis, is_complex=True)
    C_sp1 = get_structure_constants(sp1_basis, is_complex=True)
    
    # 3. Assert exact algebraic equivalence
    assert np.allclose(C_su2, C_sp1, atol=1e-10), \
        "The structure constants of su(2) and sp(1) do not match! They are not equivalent."