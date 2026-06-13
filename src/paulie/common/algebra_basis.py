"""
Defining-representation matrix bases for the classified Lie algebras.

The bases are the compact real forms used by the recursive Cartan decomposition pipeline of
Wierichs et al. (arXiv:2503.19014, Eq. (6)-(9)), so every generator is anti-Hermitian and
``exp`` of any real combination is unitary:

* :math:`\\mathfrak{so}(m)`: real antisymmetric, :math:`x = -x^{T}`.
* :math:`\\mathfrak{su}(m)`: traceless anti-Hermitian, :math:`x = -x^{\\dagger}`.
* :math:`\\mathfrak{sp}(m) = \\{x \\in \\mathfrak{su}(2m) : x = -J_m x^{T} J_m^{T}\\}`, the compact
  symplectic algebra, with :math:`J_m = [[0, \\mathbb{1}_m], [-\\mathbb{1}_m, 0]]`.

Each constructor returns a stack of shape ``(dim, d, d)`` with a fixed, documented ordering so
that downstream consumers get a stable basis. A direct sum is assembled block-diagonally by
:func:`block_diagonal_basis` into a single basis of the complete operator. The ordering follows
the conventions of Wierichs et al. (Tab. II) so the basis feeds directly into the recursive
Cartan-decomposition (KAK) pipeline that motivates the issue.
"""
import numpy as np


def _require_positive(n: int, name: str) -> None:
    """
    Validate that a size parameter is a positive integer.

    Args:
        n (int): Size parameter.
        name (str): Algebra name, used in the error message.
    Returns:
        None
    Raises:
        ValueError: If ``n`` is not a positive integer.
    """
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError(f"{name} requires a positive integer size, got {n!r}")


def symplectic_form(n: int) -> np.ndarray:
    """
    Standard symplectic form :math:`J_n = [[0, I_n], [-I_n, 0]]`.

    Args:
        n (int): Half the matrix dimension.
    Returns:
        numpy.ndarray: The ``2n x 2n`` symplectic form.
    """
    _require_positive(n, "symplectic_form")
    j = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    j[:n, n:] = np.eye(n)
    j[n:, :n] = -np.eye(n)
    return j


def so_basis(n: int) -> np.ndarray:
    """
    Basis of :math:`\\mathfrak{so}(n)`: real antisymmetric matrices.

    The generators are :math:`E_{ij} - E_{ji}` for :math:`i < j`, ordered by ``i`` then ``j``.
    Each is real antisymmetric, hence anti-Hermitian, with :math:`\\mathrm{tr}(B^{\\dagger}B) = 2`.

    Args:
        n (int): Matrix dimension.
    Returns:
        numpy.ndarray: Stack of shape ``(n(n-1)/2, n, n)``.
    Raises:
        ValueError: If ``n`` is not a positive integer.
    """
    _require_positive(n, "so(n)")
    generators = []
    for i in range(n):
        for j in range(i + 1, n):
            mat = np.zeros((n, n), dtype=np.float64)
            mat[i, j] = 1.0
            mat[j, i] = -1.0
            generators.append(mat)
    return np.array(generators, dtype=np.float64).reshape((-1, n, n))


def su_basis(n: int) -> np.ndarray:
    """
    Basis of :math:`\\mathfrak{su}(n)`: traceless anti-Hermitian matrices.

    The generators are the anti-Hermitian generalized Gell-Mann matrices, ordered as
    antisymmetric :math:`E_{ij} - E_{ji}`, then symmetric :math:`i(E_{ij} + E_{ji})`
    (both for :math:`i < j`), then the :math:`n - 1` diagonal generators. Every generator is
    uniformly normalized to :math:`\\mathrm{tr}(B^{\\dagger}B) = 2`.

    Args:
        n (int): Matrix dimension.
    Returns:
        numpy.ndarray: Stack of shape ``(n^2 - 1, n, n)``.
    Raises:
        ValueError: If ``n`` is not a positive integer.
    """
    _require_positive(n, "su(n)")
    generators = []
    for i in range(n):
        for j in range(i + 1, n):
            antisymmetric = np.zeros((n, n), dtype=np.complex128)
            antisymmetric[i, j] = 1.0
            antisymmetric[j, i] = -1.0
            generators.append(antisymmetric)
            symmetric = np.zeros((n, n), dtype=np.complex128)
            symmetric[i, j] = 1.0j
            symmetric[j, i] = 1.0j
            generators.append(symmetric)
    for level in range(1, n):
        diagonal = np.zeros((n, n), dtype=np.complex128)
        scale = np.sqrt(2.0 / (level * (level + 1)))
        for k in range(level):
            diagonal[k, k] = 1.0j * scale
        diagonal[level, level] = -1.0j * scale * level
        generators.append(diagonal)
    return np.array(generators, dtype=np.complex128).reshape((-1, n, n))


def sp_basis(n: int) -> np.ndarray:
    """
    Basis of the compact symplectic algebra :math:`\\mathfrak{sp}(n) = \\mathfrak{usp}(2n)`.

    Generators are the ``2n x 2n`` anti-Hermitian matrices :math:`x = -J_n x^{T} J_n^{T}`,
    parametrized as :math:`[[A, B], [-\\bar{B}, \\bar{A}]]`
    with ``A`` anti-Hermitian and ``B`` complex symmetric. The ordering is the ``A`` block
    (diagonal :math:`iE_{kk}`, then antisymmetric, then symmetric off-diagonal generators of
    :math:`\\mathfrak{u}(n)`) followed by the ``B`` block (real symmetric, then imaginary
    symmetric).

    Args:
        n (int): Rank; the matrices act on :math:`\\mathbb{C}^{2n}`.
    Returns:
        numpy.ndarray: Stack of shape ``(n(2n+1), 2n, 2n)``.
    Raises:
        ValueError: If ``n`` is not a positive integer.
    """
    _require_positive(n, "sp(n)")
    generators = []
    for a_block in _anti_hermitian_basis(n):
        mat = np.zeros((2 * n, 2 * n), dtype=np.complex128)
        mat[:n, :n] = a_block
        mat[n:, n:] = np.conjugate(a_block)
        generators.append(mat)
    for b_block in _symmetric_basis(n):
        mat = np.zeros((2 * n, 2 * n), dtype=np.complex128)
        mat[:n, n:] = b_block
        mat[n:, :n] = -np.conjugate(b_block)
        generators.append(mat)
    return np.array(generators, dtype=np.complex128).reshape((-1, 2 * n, 2 * n))


def u1_basis() -> np.ndarray:
    """
    Basis of :math:`\\mathfrak{u}(1)`: the single generator :math:`[[i]]`.

    Returns:
        numpy.ndarray: Stack of shape ``(1, 1, 1)``.
    """
    return np.array([[[1.0j]]], dtype=np.complex128)


def block_diagonal_basis(summands: list[np.ndarray]) -> np.ndarray:
    """
    Assemble a basis of the direct sum from per-summand bases by block-diagonal embedding.

    The direct sum :math:`\\bigoplus_k \\mathfrak{g}_k` is realized in the common space
    :math:`\\mathbb{C}^{D \\times D}` with :math:`D = \\sum_k d_k`. Summand ``k`` occupies the
    ``k``-th diagonal block, so the result is a basis of the complete operator (its dimension is
    :math:`\\sum_k \\dim \\mathfrak{g}_k`) rather than a list of disjoint block bases.

    Args:
        summands (list[numpy.ndarray]): Per-summand bases, each of shape ``(dim_k, d_k, d_k)``.
    Returns:
        numpy.ndarray: Stack of shape ``(sum dim_k, D, D)``.
    Raises:
        ValueError: If no summands are given.
    """
    if not summands:
        raise ValueError("a direct sum needs at least one summand")
    block_dims = [summand.shape[1] for summand in summands]
    total_dim = sum(summand.shape[0] for summand in summands)
    size = sum(block_dims)
    basis = np.zeros((total_dim, size, size), dtype=np.complex128)
    row = 0
    offset = 0
    for summand, dim in zip(summands, block_dims):
        count = summand.shape[0]
        basis[row:row + count, offset:offset + dim, offset:offset + dim] = summand
        row += count
        offset += dim
    return basis


def _anti_hermitian_basis(n: int) -> list[np.ndarray]:
    """
    Basis of :math:`\\mathfrak{u}(n)`: anti-Hermitian matrices (the ``A`` block of sp).

    Args:
        n (int): Matrix dimension.
    Returns:
        list[numpy.ndarray]: The ``n^2`` anti-Hermitian generators.
    """
    generators = []
    for k in range(n):
        diagonal = np.zeros((n, n), dtype=np.complex128)
        diagonal[k, k] = 1.0j
        generators.append(diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            antisymmetric = np.zeros((n, n), dtype=np.complex128)
            antisymmetric[i, j] = 1.0
            antisymmetric[j, i] = -1.0
            generators.append(antisymmetric)
            symmetric = np.zeros((n, n), dtype=np.complex128)
            symmetric[i, j] = 1.0j
            symmetric[j, i] = 1.0j
            generators.append(symmetric)
    return generators


def _symmetric_basis(n: int) -> list[np.ndarray]:
    """
    Basis of complex symmetric matrices (the ``B`` block of sp).

    Args:
        n (int): Matrix dimension.
    Returns:
        list[numpy.ndarray]: The ``n(n+1)`` symmetric generators (real then imaginary).
    """
    generators = []
    for scale in (1.0, 1.0j):
        for i in range(n):
            for j in range(i, n):
                symmetric = np.zeros((n, n), dtype=np.complex128)
                symmetric[i, j] = scale
                symmetric[j, i] = scale
                generators.append(symmetric)
    return generators
