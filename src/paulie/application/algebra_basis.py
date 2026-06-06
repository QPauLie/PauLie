"""
Matrix bases for the classical Lie algebras :math:`\\mathfrak{u}(N)`,
:math:`\\mathfrak{su}(N)`, :math:`\\mathfrak{so}(N)` and :math:`\\mathfrak{sp}(N)`
in their defining representations.

The functions exposed here produce a list of generators that span the named
Lie algebra. Each generator is a square NumPy array; the full basis is
returned as a stack of shape ``(dim, M, M)`` where ``dim`` is the algebra's
dimension and ``M`` is the defining-representation matrix size.

These primitives are consumed by :meth:`Classification.get_algebra_basis`
and :meth:`PauliStringCollection.get_algebra_basis`, which turn an
isomorphism label like ``"4*so(5)"`` into the concrete matrices a user
needs for Cartan decompositions, structure-constant computations, or any
downstream task that requires the algebra as a mathematical object rather
than a name.
"""

from __future__ import annotations

import numpy as np


def so_basis(n: int) -> np.ndarray:
    r"""
    Basis of :math:`\mathfrak{so}(n)` in its defining representation.

    :math:`\mathfrak{so}(n)` consists of the real antisymmetric
    :math:`n \times n` matrices. The basis returned here is the
    standard one made from the matrix units
    :math:`E_{ij}` (one at row :math:`i`, column :math:`j`, zero
    elsewhere):

    .. math::

        T_{(i,j)} = E_{ij} - E_{ji}, \qquad 0 \le i < j < n.

    There are :math:`\binom{n}{2} = n(n-1)/2` such generators, which is
    the dimension of :math:`\mathfrak{so}(n)`. They satisfy
    :math:`T^\top = -T`, span the real antisymmetric matrices, and the
    bracket :math:`[T_{(i,j)}, T_{(k,l)}]` is again real antisymmetric so
    the algebra closes.

    Args:
        n: Matrix size of the defining representation. ``n >= 1``.
           ``n == 1`` returns an empty basis (the trivial algebra).

    Returns:
        np.ndarray: Real array of shape ``(n*(n-1)//2, n, n)``.

    Examples:
        >>> from paulie.application.algebra_basis import so_basis
        >>> b = so_basis(3)
        >>> b.shape
        (3, 3, 3)
        >>> # Antisymmetry of every generator
        >>> import numpy as np
        >>> np.allclose(b + b.swapaxes(1, 2), 0)
        True

    References:
        Hall, B. C. *Lie Groups, Lie Algebras, and Representations*,
        2nd ed., Springer (2015), Section 3.4.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    dim = n * (n - 1) // 2
    out = np.zeros((dim, n, n), dtype=np.float64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            out[idx, i, j] = 1.0
            out[idx, j, i] = -1.0
            idx += 1
    return out


def sp_basis(n: int) -> np.ndarray:
    r"""
    Basis of :math:`\mathfrak{sp}(n)` in its defining representation.

    Following the convention used elsewhere in PauLie (consistent with
    :cite:t:`Aguilar_2024`), :math:`\mathfrak{sp}(n)` denotes the
    rank-:math:`n` real symplectic algebra. Its defining representation
    acts on :math:`\mathbb{R}^{2n}` by :math:`2n \times 2n` matrices
    :math:`M` satisfying

    .. math::

        M^\top J + J M = 0, \qquad
        J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}.

    Writing :math:`M` in block form
    :math:`\begin{pmatrix} A & B \\ C & -A^\top \end{pmatrix}` with
    :math:`A` arbitrary and :math:`B, C` symmetric :math:`n \times n`
    blocks gives :math:`n^2 + n(n+1) = n(2n+1)` real parameters, which
    is the dimension of :math:`\mathfrak{sp}(n)`. The basis returned
    here exposes each parameter as a separate generator: the :math:`n^2`
    "A-type" generators come first, followed by the :math:`n(n+1)/2`
    "B-type" generators and the :math:`n(n+1)/2` "C-type" generators.

    Args:
        n: Rank of the symplectic algebra. ``n >= 1``. Matrices have
           shape ``(2n, 2n)``.

    Returns:
        np.ndarray: Real array of shape ``(n*(2n+1), 2n, 2n)``.

    Examples:
        >>> from paulie.application.algebra_basis import sp_basis
        >>> b = sp_basis(1)
        >>> b.shape           # dim sp(1) = 3, matrices 2x2
        (3, 2, 2)
        >>> # Defining condition: M^T J + J M = 0
        >>> import numpy as np
        >>> J = np.array([[0, 1], [-1, 0]], dtype=float)
        >>> np.allclose(b.swapaxes(1, 2) @ J + J @ b, 0)
        True

    References:
        Hall, B. C. *Lie Groups, Lie Algebras, and Representations*,
        2nd ed., Springer (2015), Section 3.4 (Example 4).
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    matrix_size = 2 * n
    dim = n * (2 * n + 1)
    out = np.zeros((dim, matrix_size, matrix_size), dtype=np.float64)
    idx = 0

    # A-block: n*n arbitrary generators. M placed at (i, j) of upper-left
    # block forces -1 at (n+j, n+i) of lower-right to keep -A^T relation.
    for i in range(n):
        for j in range(n):
            out[idx, i, j] = 1.0
            out[idx, n + j, n + i] = -1.0
            idx += 1

    # B-block: symmetric n*n in the upper-right corner.
    for i in range(n):
        for j in range(i, n):
            out[idx, i, n + j] = 1.0
            if i != j:
                out[idx, j, n + i] = 1.0
            idx += 1

    # C-block: symmetric n*n in the lower-left corner.
    for i in range(n):
        for j in range(i, n):
            out[idx, n + i, j] = 1.0
            if i != j:
                out[idx, n + j, i] = 1.0
            idx += 1

    return out


def su_basis(n: int) -> np.ndarray:
    r"""
    Basis of :math:`\mathfrak{su}(n)` in its defining representation.

    :math:`\mathfrak{su}(n)` consists of the traceless anti-Hermitian
    :math:`n \times n` complex matrices. The basis returned here is
    constructed from generalised Gell-Mann matrices, each multiplied by
    :math:`i` so the generators are anti-Hermitian (this matches the
    physics convention :math:`U = \exp(iH) \in \mathrm{SU}(n)` with
    Hermitian :math:`H`, viewed as the Lie algebra element :math:`iH`):

    - **Symmetric off-diagonal**, for :math:`j < k`:
      :math:`i(\lvert j\rangle\!\langle k\rvert + \lvert k\rangle\!\langle j\rvert)`.
    - **Antisymmetric off-diagonal**, for :math:`j < k`:
      :math:`\lvert j\rangle\!\langle k\rvert - \lvert k\rangle\!\langle j\rvert`.
    - **Diagonal traceless**, for :math:`0 \le m < n-1`:
      :math:`i\bigl(\lvert m\rangle\!\langle m\rvert - \lvert m+1\rangle\!\langle m+1\rvert\bigr)`.

    Total: :math:`n^2 - 1` generators, matching
    :math:`\dim \mathfrak{su}(n)`. The diagonal generators are NOT
    orthonormalised; they are linearly independent, which is enough for
    a basis.

    Args:
        n: Matrix size of the defining representation. ``n >= 1``.
           ``n == 1`` returns an empty basis (since :math:`\mathfrak{su}(1) = 0`).

    Returns:
        np.ndarray: Complex array of shape ``(n*n - 1, n, n)``.

    Examples:
        >>> from paulie.application.algebra_basis import su_basis
        >>> b = su_basis(2)
        >>> b.shape           # 3 Pauli-like generators
        (3, 2, 2)
        >>> # Anti-Hermitian and traceless
        >>> import numpy as np
        >>> np.allclose(b + b.conj().swapaxes(1, 2), 0)
        True
        >>> np.allclose(np.trace(b, axis1=1, axis2=2), 0)
        True

    References:
        Georgi, H. *Lie Algebras in Particle Physics*, 2nd ed.,
        Westview Press (1999), Chapter 1; generalised Gell-Mann matrices
        also discussed in Bertlmann & Krammer, J. Phys. A 41, 235303
        (2008).
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    dim = n * n - 1
    out = np.zeros((dim, n, n), dtype=np.complex128)
    idx = 0

    # Symmetric off-diagonal: i*(|j><k| + |k><j|)
    for j in range(n):
        for k in range(j + 1, n):
            out[idx, j, k] = 1j
            out[idx, k, j] = 1j
            idx += 1

    # Antisymmetric off-diagonal: |j><k| - |k><j|
    for j in range(n):
        for k in range(j + 1, n):
            out[idx, j, k] = 1.0
            out[idx, k, j] = -1.0
            idx += 1

    # Diagonal traceless: i*(|m><m| - |m+1><m+1|)
    for m in range(n - 1):
        out[idx, m, m] = 1j
        out[idx, m + 1, m + 1] = -1j
        idx += 1

    return out


def u_basis(n: int) -> np.ndarray:
    r"""
    Basis of :math:`\mathfrak{u}(n)` in its defining representation.

    :math:`\mathfrak{u}(n) = \mathfrak{su}(n) \oplus \mathfrak{u}(1)`,
    i.e. the anti-Hermitian :math:`n \times n` complex matrices. The
    basis returned here is :func:`su_basis` extended by the additional
    central generator :math:`i I_n`, giving :math:`n^2` matrices total.

    Args:
        n: Matrix size of the defining representation. ``n >= 1``.

    Returns:
        np.ndarray: Complex array of shape ``(n*n, n, n)``.

    Examples:
        >>> from paulie.application.algebra_basis import u_basis
        >>> b = u_basis(1)            # u(1) = iR
        >>> b.shape
        (1, 1, 1)
        >>> b[0, 0, 0]                # generator is i
        1j
        >>> import numpy as np
        >>> np.allclose(b + b.conj().swapaxes(1, 2), 0)  # anti-Hermitian
        True

    References:
        Hall, B. C. *Lie Groups, Lie Algebras, and Representations*,
        2nd ed., Springer (2015), Section 1.2.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    su = su_basis(n)
    central = (1j * np.eye(n, dtype=np.complex128))[np.newaxis, ...]
    return np.concatenate([su, central], axis=0)
