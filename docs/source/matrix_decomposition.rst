Matrix decomposition
====================

This tutorial demonstrates how to use :code:`paulie` to decompose a matrix in :math:`\mathbb{C}^{2^n \times 2^n}`
as a sum of Pauli strings of length :math:`n`.

The set of Pauli strings of length :math:`n` given by :math:`\mathcal{P}=\{\bigotimes_{i=1}^{n}U_{i} \mid U_{i} \in (I, X, Y, Z)\}`
forms a basis of the vector space of complex matrices of dimensions :math:`2^n \times 2^n` denoted by :math:`\mathbb{C}^{2^n \times 2^n}`.

:code:`paulie` provides the method :code:`matrix_decomposition` for decomposing general matrices which returns a dictionary of Pauli strings
with the corresponding coefficients:

.. code-block:: python

    from paulie.application.matrix_decomposition import matrix_decomposition
    import numpy as np

    A = np.array([[0.48+0.9j , 0.64+0.3j , 0.47+0.97j, 0.96+0.11j],
                  [0.52+0.05j, 0.96+0.82j, 1.  +0.27j, 0.89+0.22j],
                  [0.01+0.2j , 0.22+0.87j, 0.37+0.25j, 0.4 +0.68j],
                  [0.04+0.2j , 0.4 +0.83j, 0.22+0.25j, 0.69+0.19j]])
    decomp = matrix_decomposition(A)
    print('Pauli decomposition: ' + ' + '.join([f'{np.round(coeff, 5)}*{str(p)}' for p, coeff in decomp.items()]))
    print('\nObtained matrix:')
    print(sum(coeff * p.get_matrix() for p, coeff in decomp.items()))

which outputs:

.. code-block:: bash

    Pauli decomposition:
    (0.625+0.54j)*II + (0.4425+0.555j)*XI + (-0.04+0.2375j)*YI + (0.095+0.32j)*ZI + (0.445+0.32j)*IX + (0.555+0.3625j)*XX + (0.1725+0.425j)*YX + (0.135-0.145j)*ZX + (-0.17+0.075j)*IY + (-0.1275+0.035j)*XY + (0.055+0.2075j)*YY + (0.045-0.015j)*ZY + (-0.2+0.035j)*IZ + (-0.2025+0.03j)*XZ + (-0.345-0.0075j)*YZ + (-0.04+0.005j)*ZZ

    Obtained matrix:
    [[0.48+0.9j  0.64+0.3j  0.47+0.97j 0.96+0.11j]
    [0.52+0.05j 0.96+0.82j 1.  +0.27j 0.89+0.22j]
    [0.01+0.2j  0.22+0.87j 0.37+0.25j 0.4 +0.68j]
    [0.04+0.2j  0.4 +0.83j 0.22+0.25j 0.69+0.19j]]
    
There is also a specialized method :code:`matrix_decomposition_diagonal` for decomposing diagonal matrices:

.. code-block:: python

    from paulie.application.matrix_decomposition import matrix_decomposition_diagonal
    import numpy as np

    A = np.array([0.41+0.16j, 0.1 +0.72j, 0.31+0.12j, 0.28+0.67j])
    decomp = matrix_decomposition_diagonal(A)
    print('Pauli decomposition: ' + ' + '.join([f'{np.round(coeff, 5)}*{str(p)}' for p, coeff in decomp.items()]))
    print('\nObtained matrix:')
    print(sum(coeff * p.get_matrix() for p, coeff in decomp.items()))

which outputs:

.. code-block:: bash

    Pauli decomposition: (0.275+0.4175j)*II + (-0.02+0.0225j)*ZI + (0.085-0.2775j)*IZ + (0.07-0.0025j)*ZZ

    Obtained matrix:
    [[0.41+0.16j 0.  +0.j   0.  +0.j   0.  +0.j  ]
    [0.  +0.j   0.1 +0.72j 0.  +0.j   0.  +0.j  ]
    [0.  +0.j   0.  +0.j   0.31+0.12j 0.  +0.j  ]
    [0.  +0.j   0.  +0.j   0.  +0.j   0.28+0.67j]]
