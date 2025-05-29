Matrix Decomposition
====================

This tutorial demonstrates how to use :code:`paulie` to decompose an arbitrary \(2^n \times 2^n\) complex matrix into the Pauli basis.  Any operator \(A\) on \(n\) qubits can be written as:

.. math::
    A = \sum_{P \in \mathcal{P}_n} \alpha_P P,\quad
    \alpha_P = \frac{1}{2^n} \mathrm{Tr}(P^\dagger A)

where :math:`\mathcal{P}_n` is the set of all :math:`4^n` Pauli strings.

The function :func:`paulie.application.matrix_decomposition.matrix_decomposition` computes this mapping and returns a dictionary from each :class:`~paulie.common.pauli_string_bitarray.PauliString` to its complex coefficient.

Usage Example
-------------

.. code-block:: python

    import numpy as np
    from paulie.common.pauli_string_factory import get_pauli_string as p
    from paulie.application.matrix_decomposition import matrix_decomposition

    # 2-qubit example: A = 0.5 * I\u00d7I + 0.25 * X\u00d7Y
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    A = 0.5 * np.kron(I, I) + 0.25 * np.kron(X, Y)

    coeffs = matrix_decomposition(A)
    for ps, coeff in coeffs.items():
        print(f"{ps}: {coeff}")

Expected Output
---------------

.. code-block:: text

    II: (0.5+0j)
    XY: (0.25+0j)