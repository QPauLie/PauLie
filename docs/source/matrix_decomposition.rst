Matrix decomposition
====================

This tutorial demonstrates how to use :code:`paulie` to decompose a matrix in :math:`\mathbb{C}^{2^n \times 2^n}`
as a sum of Pauli strings of length :math:`n`.

The set of Pauli strings of length :math:`n` given by :math:`\mathcal{P}=\{\bigotimes_{i=1}^{n}U_{i} \mid U_{i} \in (I, X, Y, Z)\}`
forms a basis of the vector space of complex matrices of dimensions :math:`2^n \times 2^n` denoted by :math:`\mathbb{C}^{2^n \times 2^n}`.

:code:`paulie` provides the method :code:`matrix_decomposition` for decomposing general matrices which returns a weight vector:

.. code-block:: python

    import numpy as np
    from paulie.application.matrix_decomposition import matrix_decomposition

    A = np.array([[0. +0.j, -2.9-0.9j, 0.5+0.j, 0. +0.j]
                  [3.1+1.5j, 0. +0.j, 0. +0.j, -0.5+0.j]
                  [0.5+0.j, 0. +0.j, 0. +0.j, 3.1+1.5j]
                  [0. +0.j, -0.5+0.j, -2.9-0.9j, 0. +0.j]])
    decomp = matrix_decomposition(A)
    print("Decomposition weight vector:")
    print(decomp)

which outputs:

.. code-block:: bash

    Decomposition weight vector:
    [0. +0.j  0. +0.j  0.1+0.3j 0. +0.j  0. +0.j  0. +0.j  0. +0.j  1.2-3.j
     0. +0.j  0.5+0.j  0. +0.j  0. +0.j  0. +0.j  0. +0.j  0. +0.j  0. +0.j ]

From this weight vector, we can obtain the coefficient of a given Pauli string in the decomposition:

.. code-block:: python

    from paulie.common.pauli_string_factory import get_pauli_string as p

    p1 = p("IX")
    p2 = p("ZY")
    p3 = p("XZ")

    print("Coefficient of IX: ", np.round(p1.get_weight_in_matrix(decomp), 5))
    print("Coefficient of ZY: ", np.round(p2.get_weight_in_matrix(decomp), 5))
    print("Coefficient of XZ: ", np.round(p3.get_weight_in_matrix(decomp), 5))

which outputs:

.. code-block:: bash

    Coefficient of IX:  (0.1+0.3j)
    Coefficient of ZY:  (1.2-3j)
    Coefficient of XZ:  (0.5+0j)

There is also a specialized method :code:`matrix_decomposition_diagonal` for decomposing diagonal matrices:

.. code-block:: python

    from paulie.application.matrix_decomposition import matrix_decomposition_diagonal

    A = np.array([0.41+0.16j, 0.1 +0.72j, 0.31+0.12j, 0.28+0.67j])
    decomp = matrix_decomposition_diagonal(A)
    pauli_strs = ["II", "IZ", "ZI", "ZZ"]
    for pauli_str in pauli_strs:
        pstr = p(pauli_str)
        coeff = pstr.get_weight_in_matrix(decomp)
        print(f"Coefficient of {pauli_str}: {np.round(coeff, 5)}")

which outputs:

.. code-block:: bash

    Coefficient of II: (0.275+0.4175j)
    Coefficient of IZ: (0.085-0.2775j)
    Coefficient of ZI: (-0.02+0.0225j)
    Coefficient of ZZ: (0.07-0.0025j)

