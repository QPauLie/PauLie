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
and edges between all vertices that do not commute. Starting with an empty generator set, we incrementally add Pauli strings to our set. At each step, we must keep the anticommutation graph in a *canonical form* (explained later). The fundamental operation we use is a *contraction* between two Pauli strings :math:`P_i` and :math:`P_j` which maps :math:`P_i \mapsto \pm \frac{1}{2} i [P_i,P_j] = P_i^\star`. Crucially, it can be shown that this operation leaves the Lie algebra invariant. Now if :math:`P_i^\star` is already in :math:`\mathcal{G}`, the size of the generator set has been reduced. Otherwise, the operation results the complementation of the edge set between :math:`P_i^\star` and vertices in the neighbourhood of :math:`P_j`. Simply put, we can use such contractions to manipulates the edges of the anticommutation graph so as to bring it to a canonical form.

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

All four canonical types are starlike graphs. The algorithm essentially proceeds by first creating a "core" of vertices (the central vertex of the star and at most two other vertices), and then appending vertices to legs or splitting legs and connecting them to the center, as required. There are essentially four key steps:

- Step 1: Build the core.
- Step 2: If the vertex to be added is only connected to a subset of the isolated (only connected to center excluding the vertex to be added) vertices,
  use contractions until it is only connected to one isolated vertex and exit. Otherwise use contractions to disconnect the isolated vertices from the
  vertex to be added and go to Step 3.
- Step 3: Use contractions until the vertex is only connected to vertices on the longest leg of the graph.
- Step 4: Use contractions until we can attach the vertex to the end of the longest leg or to the center. This step cannot always be achieved
  using just contractions: we may have to remove a vertex to complete this step. But it is fine since we can add it again afterwards without affecting the algebra.

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

We can also animate the transformation to a canonical graph.
We use the following color convention.

.. Here we should add some definitions like lit, dependent etc.

.. role:: raw-html(raw)
   :format: html

.. list-table::
   :header-rows: 1
   :widths: 40 20

   * - Node
     - Color
   * - Lighting
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:red; border:1px solid #666;"></span>`
   * - Dependent
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#2F4F4F; border:1px solid #666;"></span>`
   * - Unlit vertex in a leg of length one
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#FF00FF; border:1px solid #666;"></span>`
   * - Removing
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:black; border:1px solid #666;"></span>`
   * - Unlit
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#cccccc; border:1px solid #666;"></span>`
   * - Appending
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#00FF00; border:1px solid #666;"></span>`
   * - Contracting
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#008080; border:1px solid #666;"></span>`
   * - Lit vertex in a leg of length one
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#6A5ACD; border:1px solid #666;"></span>`
   * - Replacing
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:#8B008B; border:1px solid #666;"></span>`
   * - Lit
     - :raw-html:`<span style="display:inline-block; width:1.4em; height:1.4em; background:cyan; border:1px solid #666;"></span>`

.. raw:: html
   :file: ../media/example_c.html

This is what happens:

1. :math:`IYZI` is added.
2. :math:`IIXX` is added as one leg.
3. :math:`IIYZ` is added as the second leg. Now the core is complete.
4. :math:`IXXI` is added. It is connected to :math:`IIYZ`.
5. :math:`XXII` is added. It is connected to :math:`IYZI`.
6. :math:`YZII` is added. It is connected to :math:`IYZI` and :math:`IXXI`. We are now in Step 4. But the configuration of attached vertices is
   such that we must remove :math:`IXXI` so that it is valid for us to attach :math:`YZII` to :math:`IYZI`.
7. Now we must add :math:`IXXI` again. It is connected to :math:`YZII` and :math:`IIYZ`. We are now in Step 2. In this case, we pick
   two isolated vertices :math:`P` and :math:`Q` such that :math:`P` is connected to the vertex to be added and :math:`Q` is not.
   Then for every remaining vertex :math:`W` which is connected to the vertex to be added, we perform :math:`W\mapsto PQW` and now :math:`PQW`
   commutes with the vertex to be added, so it is no longer connected. It is guaranteed that such a transformation is always possible by
   Theorem 4 of :cite:t:`Aguilar_2024`. In our case, we can pick :math:`P=YZII` and :math:`Q=XXII`, which means we have :math:`IIYZ\mapsto ZYYZ`
   and then we can connect :math:`IXXI` to :math:`YZII`.

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

We can also animate the transformation to a star graph:

.. raw:: html
    :file: ../media/example_d.html

This is what happens:

1. :math:`IXZI` is added.
2. :math:`IIXZ` is added as one leg.
3. :math:`IIXY` is added. It is connected to :math:`IXZI` and :math:`IIXZ`, so we contract with :math:`IXZI` to get :math:`[IXZI, IIXY]\propto IXYY`.
   We add it as the second leg. Now the core is complete.
4. :math:`IXYI` is added. It is connected to :math:`IXZI` and :math:`IIXZ`. We are now in Step 2. We pick :math:`P=IIXZ` and :math:`Q=IXYY` and we get
   :math:`IXZI\mapsto IIIX`. Then we can attach :math:`IXYI` to :math:`IIXZ`.
5. :math:`XZII` is added. It is connected to :math:`IXYY` and :math:`IXYI`. We are now in Step 2. We can contract with :math:`IXYY` and :math:`IIIX`
   in that order to get :math:`[IIIX,[IXYY,XZII]]\propto XYYZ`. Now it is connected to :math:`IIIX`, :math:`IIXZ`, and :math:`IXYI`. Now we are in
   Step 4 of the algorithm. We can contract with :math:`IXYI` to get :math:`[IXYI,XYYZ]\propto XZIZ`. Now it is no longer connected to :math:`IIXZ`.
   Next we contract with :math:`IIIX` to get :math:`[IIIX,XZIZ]\propto XZIY`. Now :math:`XZIY` is connected to all vertices. Then we contract with
   :math:`IXYY` to get :math:`[IXYY,XZIY]\propto XYYI`. Now it is connected to all except :math:`IIIX`. Then we contract with :math:`IIXZ` to get
   :math:`[IIXZ,XYYI]\propto XYZZ`. Now it is connected to all except :math:`IXYI`. Finally, we contract with :math:`IIIX` to get :math:`[IIIX,XYZZ]\propto XYZY`.
   Now it is attached to only :math:`IIIX`. We can consider :math:`IIIX` as the center.
6. :math:`XYII` is added. It is connected to :math:`IXYY` and :math:`IXYI`. We are now in Step 2. We pick :math:`P=IXYY` and :math:`Q=XYZY`
   and we get :math:`IXYI\mapsto XYZI`. Then we can attach :math:`XYII` to :math:`IXYY`.

According to the table, the resultant graph corresponds to :math:`\mathfrak{sp}(4)`.

The Lie algebra plays a pivotal role in quantum control theory to understand the reachability of states.
Also measures of operator spread complexity rely on this concept.
Furthermore, determining moments of circuits can be significantly simplified when the Lie algebra is known.
All these applications are functionalities of :code:`paulie`.







