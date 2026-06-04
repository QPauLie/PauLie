r"""Defining-representation bases for the classical Lie algebras.

The classifier (:meth:`paulie.Classification.get_algebra`) returns only the
*name* of the dynamical Lie algebra, e.g. ``"sp(4)"`` or ``"2*su(8)"``.  This
module turns such a name into a concrete, usable mathematical object: an
explicit basis of matrices in the *defining* (fundamental) representation of
each simple summand.

Every routine here is purely table-driven -- the matrices depend only on the
family and its size parameter, never on the Pauli strings that happened to
generate the algebra.  Because downstream consumers (e.g. the recursive Cartan
decomposition pipeline of Wierichs et al., arXiv:2503.19014) need a *stable*
basis, the ordering and sign conventions below are fixed and documented.

Conventions
-----------
All matrices are returned with ``complex128`` dtype so that the four families
share a single return type, even though the orthogonal and symplectic families
are real-valued.  Each algebra is realised as a *real* Lie algebra: its basis
spans a real vector space and is closed under the matrix commutator
:math:`[A, B] = AB - BA` with real structure constants.

``so(N)`` -- real antisymmetric :math:`N \times N` matrices.
    Basis :math:`\{E_{ij} - E_{ji} : 1 \le i < j \le N\}` ordered
    lexicographically in ``(i, j)``.  Dimension :math:`N(N-1)/2`.

``su(N)`` -- traceless anti-Hermitian :math:`N \times N` matrices.
    An (unnormalised) anti-Hermitian generalized Gell-Mann basis.  For each
    pair ``(i, j)`` with :math:`i < j`, taken lexicographically, two
    off-diagonal generators are emitted consecutively,

    * the real generator   :math:`E_{ij} - E_{ji}`,
    * the imaginary generator :math:`i\,(E_{ij} + E_{ji})`,

    followed by the :math:`N - 1` diagonal generators
    :math:`i\,\mathrm{diag}(1, \dots, 1, -k, 0, \dots, 0)` (``k`` leading ones)
    for :math:`k = 1, \dots, N - 1`.  Dimension :math:`N^2 - 1`.

``sp(n)`` -- the rank-``n`` symplectic algebra :math:`C_n`, realised as real
    :math:`2n \times 2n` matrices :math:`X` obeying :math:`X^{T} J + J X = 0`
    with :math:`J = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix}`.
    Writing :math:`X = \begin{bmatrix} A & B \\ C & -A^{T} \end{bmatrix}` with
    ``B`` and ``C`` symmetric, the basis is ordered as

    1. the ``A`` block, :math:`E_{ij}` (top-left) paired with :math:`-E_{ji}`
       (bottom-right), for every ``(i, j)`` lexicographically,
    2. the symmetric ``B`` block (top-right), for :math:`i \le j`,
    3. the symmetric ``C`` block (bottom-left), for :math:`i \le j`.

    Dimension :math:`n(2n + 1)`.  Note that the *rank* ``n`` is the argument;
    the matrices are :math:`2n \times 2n`.

``u(N)`` -- anti-Hermitian :math:`N \times N` matrices.
    The ``su(N)`` basis with the trace direction :math:`i\,I_N` appended last.
    Dimension :math:`N^2`.
"""

import numpy as np
import numpy.typing as npt

# Single return dtype shared by every family; keeps direct sums homogeneous.
Basis = npt.NDArray[np.complex128]


def symplectic_form(rank: int) -> Basis:
    r"""Return the symplectic form used by :func:`symplectic_basis`.

    Args:
        rank: Rank ``n`` of the symplectic algebra; ``J`` is ``2n x 2n``.

    Returns:
        The block matrix :math:`J = \begin{bmatrix} 0 & I_n \\ -I_n & 0
        \end{bmatrix}`.
    """
    eye = np.eye(rank, dtype=np.complex128)
    zero = np.zeros((rank, rank), dtype=np.complex128)
    return np.block([[zero, eye], [-eye, zero]])


def orthogonal_basis(size: int) -> Basis:
    r"""Return a basis of ``so(size)`` in the defining representation.

    Args:
        size: Matrix dimension ``N``.

    Returns:
        Array of shape ``(N(N-1)/2, N, N)`` whose entries are the real
        antisymmetric matrices :math:`E_{ij} - E_{ji}` for :math:`i < j`,
        ordered lexicographically in ``(i, j)``.
    """
    basis = []
    for i in range(size):
        for j in range(i + 1, size):
            mat = np.zeros((size, size), dtype=np.complex128)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            basis.append(mat)
    return _stack(basis, size)


def special_unitary_basis(size: int) -> Basis:
    r"""Return a basis of ``su(size)`` in the defining representation.

    The basis is the anti-Hermitian generalized Gell-Mann basis described in
    the module docstring: paired off-diagonal generators followed by the
    diagonal Cartan generators.

    Args:
        size: Matrix dimension ``N``.

    Returns:
        Array of shape ``(N^2 - 1, N, N)`` of traceless anti-Hermitian
        matrices.
    """
    basis = []
    for i in range(size):
        for j in range(i + 1, size):
            real = np.zeros((size, size), dtype=np.complex128)
            real[i, j] = 1.0
            real[j, i] = -1.0
            basis.append(real)

            imag = np.zeros((size, size), dtype=np.complex128)
            imag[i, j] = 1j
            imag[j, i] = 1j
            basis.append(imag)
    for k in range(1, size):
        diagonal = np.zeros(size, dtype=np.complex128)
        diagonal[:k] = 1.0
        diagonal[k] = -k
        basis.append(np.diag(1j * diagonal))
    return _stack(basis, size)


def symplectic_basis(rank: int) -> Basis:
    r"""Return a basis of ``sp(rank)`` in the defining representation.

    The matrices are real :math:`2\,\text{rank} \times 2\,\text{rank}` and obey
    :math:`X^{T} J + J X = 0` with ``J`` from :func:`symplectic_form`.  See the
    module docstring for the block ordering.

    Args:
        rank: Rank ``n`` of the symplectic algebra; matrices are ``2n x 2n``.

    Returns:
        Array of shape ``(n(2n+1), 2n, 2n)``.
    """
    dim = 2 * rank
    basis = []
    # A block: top-left E_ij paired with bottom-right -E_ji = -(E_ij)^T.
    for i in range(rank):
        for j in range(rank):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[i, j] = 1.0
            mat[rank + j, rank + i] = -1.0
            basis.append(mat)
    # Symmetric B block (top-right).
    for i in range(rank):
        for j in range(i, rank):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[i, rank + j] = 1.0
            mat[j, rank + i] = 1.0
            basis.append(mat)
    # Symmetric C block (bottom-left).
    for i in range(rank):
        for j in range(i, rank):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[rank + i, j] = 1.0
            mat[rank + j, i] = 1.0
            basis.append(mat)
    return _stack(basis, dim)


def unitary_basis(size: int) -> Basis:
    r"""Return a basis of ``u(size)`` in the defining representation.

    Args:
        size: Matrix dimension ``N``.

    Returns:
        Array of shape ``(N^2, N, N)`` of anti-Hermitian matrices: the
        ``su(N)`` basis followed by the trace direction :math:`i I_N`.
    """
    trace_direction = 1j * np.eye(size, dtype=np.complex128)
    if size == 1:
        return _stack([trace_direction], size)
    return np.concatenate(
        [special_unitary_basis(size), trace_direction[np.newaxis, :, :]], axis=0
    )


def _stack(matrices: list[Basis], size: int) -> Basis:
    """Stack basis matrices into one array, preserving an empty leading axis.

    Args:
        matrices: List of equally shaped ``(size, size)`` matrices.
        size: The matrix dimension, used to shape an empty result.

    Returns:
        Array of shape ``(len(matrices), size, size)``.
    """
    if not matrices:
        return np.empty((0, size, size), dtype=np.complex128)
    return np.stack(matrices, axis=0)
