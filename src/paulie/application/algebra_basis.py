r"""
Defining-representation bases for the canonical Lie algebras returned by the
classifier (:meth:`paulie.Classification.get_algebra`).

Every function here is purely table-driven: the basis of an algebra is fully
determined by its name, with no dependence on the Pauli strings that generated
it.  The conventions below are **fixed** so that downstream consumers (e.g. a
Cartan-decomposition pipeline) get a stable basis ordering.

Conventions
-----------
``so(N)`` : real :math:`N\times N` antisymmetric matrices.
    Basis :math:`\{E_{ij} - E_{ji} : 1 \le i < j \le N\}`, ordered by
    ``(i, j)`` lexicographically.  Dimension :math:`N(N-1)/2`.

``su(N)`` : traceless anti-Hermitian :math:`N\times N` matrices.
    Generalized Gell-Mann basis multiplied by :math:`i`, ordered as
    (1) symmetric off-diagonal :math:`i(E_{ij} + E_{ji})`,
    (2) antisymmetric off-diagonal :math:`E_{ij} - E_{ji}`, both for
    :math:`1 \le i < j \le N` lexicographically in ``(i, j)``, then
    (3) diagonal :math:`i\,\mathrm{diag}(1, \dots, 1, -k, 0, \dots, 0)`
    (``k`` leading ones) for :math:`k = 1, \dots, N-1`.
    Dimension :math:`N^2 - 1`.

``sp(n)`` : the rank-``n`` symplectic algebra :math:`C_n`, realized as real
    :math:`2n \times 2n` matrices :math:`X` with :math:`X^{T} J + J X = 0`,
    where :math:`J = \begin{bmatrix} 0 & I_n \\ -I_n & 0 \end{bmatrix}`.
    Writing :math:`X = \begin{bmatrix} A & B \\ C & -A^{T} \end{bmatrix}` with
    ``B``, ``C`` symmetric, the basis is ordered as
    (1) the ``A`` block :math:`E_{ij}` (paired with :math:`-E_{ji}` in the
    lower-right block) for all ``(i, j)``,
    (2) the symmetric ``B`` block for :math:`i \le j`,
    (3) the symmetric ``C`` block for :math:`i \le j`.
    Dimension :math:`n(2n+1)`.  Note the label argument ``n`` is the *rank*;
    the matrices are :math:`2n \times 2n`.

``u(N)`` : anti-Hermitian :math:`N\times N` matrices.
    Basis is the ``su(N)`` basis together with :math:`i\,I_N`, appended last.
    Dimension :math:`N^2`.
"""
import numpy as np


def so_basis(size: int) -> np.ndarray:
    r"""
    Basis of ``so(size)`` in the defining representation.

    Args:
        size (int): Matrix dimension ``N``.
    Returns:
        np.ndarray: Array of shape ``(N(N-1)/2, N, N)`` of real antisymmetric
        matrices :math:`E_{ij} - E_{ji}` ordered by ``(i, j)``, ``i < j``.
    """
    basis = []
    for i in range(size):
        for j in range(i + 1, size):
            mat = np.zeros((size, size), dtype=np.complex128)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            basis.append(mat)
    return np.array(basis, dtype=np.complex128)


def su_basis(size: int) -> np.ndarray:
    r"""
    Basis of ``su(size)`` in the defining representation.

    The generalized Gell-Mann (Hermitian, traceless) matrices multiplied by
    :math:`i`, giving traceless anti-Hermitian matrices.

    Args:
        size (int): Matrix dimension ``N``.
    Returns:
        np.ndarray: Array of shape ``(N^2 - 1, N, N)`` of traceless
        anti-Hermitian matrices.  See module docstring for the ordering.
    """
    basis = []
    # Symmetric off-diagonal: i (E_ij + E_ji)
    for i in range(size):
        for j in range(i + 1, size):
            mat = np.zeros((size, size), dtype=np.complex128)
            mat[i, j] = 1j
            mat[j, i] = 1j
            basis.append(mat)
    # Antisymmetric off-diagonal: E_ij - E_ji
    for i in range(size):
        for j in range(i + 1, size):
            mat = np.zeros((size, size), dtype=np.complex128)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            basis.append(mat)
    # Diagonal: i diag(1, ..., 1, -k, 0, ..., 0), k leading ones
    for k in range(1, size):
        diag = np.zeros(size, dtype=np.complex128)
        diag[:k] = 1.0
        diag[k] = -k
        basis.append(np.diag(1j * diag))
    return np.array(basis, dtype=np.complex128)


def sp_basis(rank: int) -> np.ndarray:
    r"""
    Basis of ``sp(rank)`` (the algebra :math:`C_{\text{rank}}`) in the defining
    representation as real :math:`2\,\text{rank} \times 2\,\text{rank}`
    matrices :math:`X` satisfying :math:`X^{T} J + J X = 0`.

    Args:
        rank (int): Rank ``n`` of the symplectic algebra.  Matrices are
            ``2n x 2n``.
    Returns:
        np.ndarray: Array of shape ``(n(2n+1), 2n, 2n)``.  See module docstring
        for the ordering and the choice of ``J``.
    """
    n = rank
    dim = 2 * n
    basis = []
    # A block: [[E_ij, 0], [0, -E_ji]]
    for i in range(n):
        for j in range(n):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[i, j] = 1.0
            mat[n + j, n + i] = -1.0
            basis.append(mat)
    # B block (symmetric, upper-right): [[0, B], [0, 0]]
    for i in range(n):
        for j in range(i, n):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[i, n + j] = 1.0
            mat[j, n + i] = 1.0
            basis.append(mat)
    # C block (symmetric, lower-left): [[0, 0], [C, 0]]
    for i in range(n):
        for j in range(i, n):
            mat = np.zeros((dim, dim), dtype=np.complex128)
            mat[n + i, j] = 1.0
            mat[n + j, i] = 1.0
            basis.append(mat)
    return np.array(basis, dtype=np.complex128)


def u_basis(size: int) -> np.ndarray:
    r"""
    Basis of ``u(size)`` in the defining representation.

    The ``su(size)`` basis together with the trace direction :math:`i I`,
    appended last.

    Args:
        size (int): Matrix dimension ``N``.
    Returns:
        np.ndarray: Array of shape ``(N^2, N, N)`` of anti-Hermitian matrices.
    """
    trace_dir = 1j * np.eye(size, dtype=np.complex128)
    if size == 1:
        return np.array([trace_dir], dtype=np.complex128)
    return np.concatenate([su_basis(size), trace_dir[np.newaxis, :, :]], axis=0)


def symplectic_form(rank: int) -> np.ndarray:
    r"""
    The symplectic form :math:`J = \begin{bmatrix} 0 & I_n \\ -I_n & 0
    \end{bmatrix}` used by :func:`sp_basis`.

    Args:
        rank (int): Rank ``n``; ``J`` is ``2n x 2n``.
    Returns:
        np.ndarray: The ``2n x 2n`` symplectic form.
    """
    n = rank
    eye = np.eye(n, dtype=np.complex128)
    zero = np.zeros((n, n), dtype=np.complex128)
    return np.block([[zero, eye], [-eye, zero]])
