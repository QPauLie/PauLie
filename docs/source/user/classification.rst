Classification of Pauli DLAs
============================

This tutorial will illustrate how to use :code:`paulie` to classify the dynamical Lie algebra of a circuit given
the generators consisting of Pauli strings.
A Pauli string is a tensor product of Pauli matrices

.. math::
    P = \bigotimes_i  \sigma_i , \, \sigma_i \in \{I,X,Y,Z\}

and is represented as a string indicating the Pauli matrices successively.
Given a set of Pauli strings, the closure under the commutator defines a Lie algebra.

In :cite:t:`Aguilar_2024`, an efficient algorithm for classifying which Lie algebra is generated is given. PauLie implements
a modified version of this algorithm.
The function :code:`get_algebra` returns exactly which algebra is generated when
given the generator set which can be extended with identities to arbitrary qubit numbers
specified.
We can reproduce Example I.5 in :cite:t:`Wiersema_2024`:

.. code-block:: python

    from paulie import get_pauli_string as p

    generators = p(["XY"])
    algebra = generators.get_algebra()
    print(f"algebra = {algebra}")

outputs

.. code-block:: bash

    algebra = u(1)

whereas changing to a three qubit system, results in another algebra:

.. code-block:: python

    size = 3
    generators = p(["XY"], n=size)
    algebra = generators.get_algebra()
    print(f"algebra = {algebra}")

outputs

.. code-block:: bash

    algebra = so(3)

The algorithm is based on the concept of an anticommutation graph. Given a set of n-qubit Pauli strings
:math:`\mathcal{G} = \{P_1,\dots ,P_{n_G}\}`, the anticommutation graph has as a vertex set :math:`\mathcal{G}`
and edges between all vertices that do not commute. The fundamental operation we use is a *contraction* between two Pauli strings :math:`P_i` and :math:`P_j` which maps :math:`P_i \mapsto \pm \frac{1}{2} i [P_i,P_j] = P_i^\star`. Crucially, it can be shown that this operation leaves the Lie algebra invariant. Now if :math:`P_i^\star` is already in :math:`\mathcal{G}`, the size of the generator set has been reduced. Otherwise, the operation results the complementation of the edge set between :math:`P_i^\star` and vertices in the neighbourhood of :math:`P_j`. Simply put, we can use such contractions to manipulates the edges of the anticommutation graph so as to bring it to a canonical form.

For any generator set consisting of Pauli strings, the anticommutation graph can be transformed to four canonical types (Theorem 1 in :cite:t:`Aguilar_2024`). Each canonical type corresponds to a particular algebra. There is an exception if :math:`\mathcal{P}` only has one Pauli string. A single Pauli string generates the algebra :math:`\mathfrak{u}(1)`.

.. table:: Canonical types and associated Lie algebras for a starlike graph with :math:`n_1 + 1` :math:`\left(n_1 \geq 0\right)` legs of length :math:`1` and :math:`n_2` :math:`\left(n_2 \geq 1\right)` legs of length :math:`2`. Note that some graphs may belong to multiple canonical types, which reflects the existence of certain *exceptional isomorphisms* between Lie algebras.

   +----------------+------------------------------------------+-------------------------------------------------------------+
   | Canonical type | Structure                                | Lie algebra                                                 |
   +================+==========================================+=============================================================+
   | A              | :math:`\left|N\right| = n_1`,            | :math:`\bigoplus_{i=1}^{2^{n_1}}\mathfrak{so}(n + 3)`       |
   |                | :math:`\left|T\right| = 0`,              |                                                             |
   |                | :math:`\left|L^{\mathcal{B}}\right| = n` |                                                             |
   +----------------+------------------------------------------+-------------------------------------------------------------+
   | B1             | :math:`\left|N\right|=n_1`,              | :math:`\bigoplus_{i=1}^{2^{n_1}}\mathfrak{sp}(2^{n_2})`     |
   |                | :math:`\left|T\right| = n_2`,            |                                                             |
   |                | :math:`\left|L^{\mathcal{B}}\right|=0`   |                                                             |
   +----------------+------------------------------------------+-------------------------------------------------------------+
   | B2             | :math:`\left|N\right|=n_1`,              | :math:`\bigoplus_{i=1}^{2^{n_1}}\mathfrak{so}(2^{n_2 + 3})` |
   |                | :math:`\left|T\right| = n_2`,            |                                                             |
   |                | :math:`\left|L^{\mathcal{B}}\right|=4`   |                                                             |
   +----------------+------------------------------------------+-------------------------------------------------------------+
   | B3             | :math:`\left|N\right|=n_1`,              | :math:`\bigoplus_{i=1}^{2^{n_1}}\mathfrak{su}(2^{n_2 + 2})` |
   |                | :math:`\left|T\right| = n_2`,            |                                                             |
   |                | :math:`\left|L^{\mathcal{B}}\right|=3`   |                                                             |
   +----------------+------------------------------------------+-------------------------------------------------------------+


Classification of A-type canonical graph
----------------------------------------
Let's try to classify a generator set that corresponds to an A-type canonical graph. This algebra is generated by :math:`\mathcal{P}=\{IYZI,IIXX,IIYZ,IXXI,XXII,YZII\}`.

.. code-block:: python

    generators = p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"])
    algebra = generators.get_algebra()
    print(f"algebra = {algebra}")

outputs

.. code-block:: bash

    algebra = 4*so(5)

.. Here we should add some definitions like lit, dependent etc.

According to the table, the resultant graph corresponds to :math:`\mathfrak{so}(5)\oplus \mathfrak{so}(5)\oplus \mathfrak{so}(5)\oplus \mathfrak{so}(5)`. But it is worth noting that it also corresponds to :math:`\mathfrak{sp}(2)\oplus \mathfrak{sp}(2)\oplus \mathfrak{sp}(2)\oplus \mathfrak{sp}(2)`. This shows that there is an exceptional isomorphism between :math:`\mathfrak{so}(5)` and :math:`\mathfrak{sp}(2)`.

Classification of B-type canonical graph
----------------------------------------
Let's try to classify a generator set that corresponds to a B-type canonical graph, that is a anticommutation graph that is a star graph. We demonstrate it by the algebra :math:`\mathfrak{a}_9` [1], generated by :math:`XY` and :math:`XZ`.

.. code-block:: python

    n_qubits = 4
    generators = p(["XY", "XZ"], n=n_qubits)
    algebra = generators.get_algebra()
    print(f"algebra = {algebra}")

outputs

.. code-block:: bash

    algebra = sp(4)

From the algebra to concrete matrices
-------------------------------------
The label returned by :code:`get_algebra` fixes the algebra up to isomorphism, but it is just a name.
To do anything algebraic with the result — dimension checks, structure constants, or feeding a
Cartan-decomposition pipeline — the algebra has to be instantiated as concrete operators.
:code:`get_algebra_basis` returns the named algebra as a list of matrices in its defining
representation, one stack of basis matrices per direct summand. The output is table-driven: it
depends only on the label, not on the generating Pauli strings.

Reusing the two examples above, :math:`\mathfrak{so}(3)` is returned as a single summand of
:math:`3\times 3` real antisymmetric matrices:

.. code-block:: python

    generators = p(["XY"], n=3)
    basis = generators.get_algebra_basis()
    print(len(basis), basis[0].shape)   # number of summands, (dim, N, N)

outputs

.. code-block:: bash

    1 (3, 3, 3)

while :math:`\mathfrak{sp}(4)` is a single summand of :math:`8\times 8` symplectic matrices
(:math:`X^{T} J + J X = 0`):

.. code-block:: python

    generators = p(["XY", "XZ"], n=4)
    basis = generators.get_algebra_basis()
    print(len(basis), basis[0].shape)

outputs

.. code-block:: bash

    1 (36, 8, 8)

For a direct sum such as :math:`2\,\mathfrak{su}(8)` the list has one entry per summand.
The basis ordering and sign conventions are fixed and documented in the docstring of
:code:`get_algebra_basis`, so downstream consumers get a stable basis.

The Lie algebra plays a pivotal role in quantum control theory to understand the reachability of states.
Also measures of operator spread complexity rely on this concept.
Furthermore, determining moments of circuits can be significantly simplified when the Lie algebra is known.
All these applications are functionalities of :code:`paulie`.







