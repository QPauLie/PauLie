Multi-angle MaxCut QAOA
=======================

This tutorial illustrates how to use :code:`paulie` to classify the dynamical Lie algebra (DLA)
of the *multi-angle* MaxCut Quantum Approximate Optimization Algorithm (QAOA) ansatz, reproducing the
six-family result of :cite:t:`Kazi_2025` and the classification shown in Figure 8 of
:cite:t:`Shaya_2025`.

The ansatz and its Lie algebra
------------------------------

For a connected graph :math:`G = (V, E)` on :math:`n = |V|` qubits, the multi-angle (also called
*free*) MaxCut QAOA ansatz assigns an independent angle to every generator in

.. math::
    \mathcal{G}_\mathrm{free} = \{\, X_v : v \in V \,\} \cup \{\, Z_u Z_v : \{u,v\} \in E \,\}.

Each :math:`X_v` is the Pauli string carrying :math:`X` on qubit :math:`v` (e.g. ``XII``) and each
:math:`Z_u Z_v` carries :math:`Z` on the two endpoints of an edge (e.g. ``ZZI``). The object we
classify is the DLA, the real span of these generators together with all of their nested
commutators,

.. math::
    \mathfrak{g}_\mathrm{free} = \langle\, i\,\mathcal{G}_\mathrm{free} \,\rangle_\mathrm{Lie}
    \subseteq \mathfrak{su}(2^n).

The six connected graph families and their Lie algebras
-------------------------------------------------------

The six families are distinguished by standard properties of the graph. A graph is *bipartite* if its
vertices split into two sets, its two *parts*, so that every edge joins one part to the other
(equivalently, it contains no odd-length cycle). Among the graphs used below, a *path* :math:`P_n`
joins :math:`n` vertices in a line, a *cycle* :math:`C_n` joins them in a closed loop, the *complete
bipartite* graph :math:`K_{a,b}` joins every vertex of one part of size :math:`a` to every vertex of
the other of size :math:`b`, and the *complete* graph :math:`K_n` joins every pair of its :math:`n`
vertices.

:cite:t:`Kazi_2025` prove that, for *any* connected graph, the multi-angle DLA falls into exactly one
of six families, determined by these properties. A graph is called *archetypal* when it is connected
but neither bipartite nor a cycle graph; the even--even, odd--odd and even--odd labels refer to the
parities of the sizes of the two parts of a connected bipartite graph.

.. table:: The six multi-angle (free) DLA families (:cite:t:`Kazi_2025`, Table II). The path and cycle rows grow polynomially in :math:`n`; the remaining four grow exponentially.

   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Graph family | Bipartite?        | :math:`\dim\,\mathfrak{g}_\mathrm{free}` | Lie algebra :math:`\mathfrak{g}_\mathrm{free}`               |
   +==============+===================+==========================================+==============================================================+
   | Path         | yes               | :math:`2n^2 - n`                         | :math:`\mathfrak{so}(2n)`                                    |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Cycle        | if :math:`n` even | :math:`4n^2 - 2n`                        | :math:`\mathfrak{so}(2n) \oplus \mathfrak{so}(2n)`           |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Even--even   | yes               | :math:`2^{2n-2} - 2^{n-1}`               | :math:`\mathfrak{so}(2^{n-1}) \oplus \mathfrak{so}(2^{n-1})` |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Odd--odd     | yes               | :math:`2^{2n-2} + 2^{n-1}`               | :math:`\mathfrak{sp}(2^{n-1}) \oplus \mathfrak{sp}(2^{n-1})` |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Even--odd    | yes               | :math:`2^{2n-2} - 1`                     | :math:`\mathfrak{su}(2^{n-1})`                               |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+
   | Archetypal   | no                | :math:`2^{2n-1} - 2`                     | :math:`\mathfrak{su}(2^{n-1}) \oplus \mathfrak{su}(2^{n-1})` |
   +--------------+-------------------+------------------------------------------+--------------------------------------------------------------+

.. note::
   The six-family result is specific to the multi-angle (free) ansatz. Thus, the classification should
   not be extended to QAOA in general.

Diagnosing barren plateaus
--------------------------

The dimension of the DLA controls how the cost-function gradients scale, and therefore whether a
variational circuit can be trained. When this dimension grows exponentially with the number of
qubits :math:`n`, the gradient variance vanishes exponentially, which is the barren-plateau
condition (:cite:t:`Kazi_2025`, Corollary 2; :cite:t:`Shaya_2025`). For the multi-angle MaxCut
ansatz, every family except the path and cycle has an exponentially large DLA (see the table above),
so those circuits are extremely prone to barren plateaus even at a single layer, whereas the
polynomially small path and cycle families may instead be classically simulable. Classifying the DLA
therefore diagnoses in advance whether the variational circuit is trainable.

Reproducing the classification
------------------------------

The six-family result can be reproduced with :code:`paulie` alone, by giving it the generators of a
small representative of each family as Pauli strings:

.. code-block:: python

    from paulie import get_pauli_string as p

    families = {
        "path P3":        ["XII", "IXI", "IIX", "ZZI", "IZZ"],
        "cycle C4":       ["XIII", "IXII", "IIXI", "IIIX", "ZZII", "IZZI", "IIZZ", "ZIIZ"],
        "even-odd K1,4":  ["XIIII", "IXIII", "IIXII", "IIIXI", "IIIIX",
                           "ZZIII", "ZIZII", "ZIIZI", "ZIIIZ"],
        "odd-odd K1,3":   ["XIII", "IXII", "IIXI", "IIIX", "ZZII", "ZIZI", "ZIIZ"],
        "archetypal K4":  ["XIII", "IXII", "IIXI", "IIIX",
                           "ZZII", "ZIZI", "ZIIZ", "IZZI", "IZIZ", "IIZZ"],
        "even-even K2,4": ["XIIIII", "IXIIII", "IIXIII", "IIIXII", "IIIIXI", "IIIIIX",
                           "ZIZIII", "ZIIZII", "ZIIIZI", "ZIIIIZ", "IZZIII", "IZIZII", "IZIIZI", "IZIIIZ"],
    }
    for name, gens in families.items():
        g = p(gens)
        print(f"{name}: {g.get_algebra()}, dim {g.get_dla_dim()}")

outputs

.. code-block:: bash

    path P3: so(6), dim 15
    cycle C4: 2*so(8), dim 56
    even-odd K1,4: su(16), dim 255
    odd-odd K1,3: 2*sp(4), dim 72
    archetypal K4: 2*su(8), dim 126
    even-even K2,4: 2*so(32), dim 992

Each line reproduces the corresponding row of the table above. The path gives a single
:math:`\mathfrak{so}` algebra and the cycle a sum of two, both of polynomial dimension; the even-odd
graph gives :math:`\mathfrak{su}`, the odd-odd graph :math:`\mathfrak{sp}`, and the even-even and
archetypal graphs sums of :math:`\mathfrak{so}` and :math:`\mathfrak{su}`, all of exponential
dimension. The smallest case is checkable by hand: the path :math:`P_3` gives :math:`\mathfrak{so}(6)`
of dimension :math:`2n^2 - n = 15` at :math:`n = 3`. (:code:`paulie` writes a direct sum of two
identical factors with a numeric prefix, so ``2*so(8)`` is :math:`\mathfrak{so}(8) \oplus
\mathfrak{so}(8)`.) For each family, :code:`paulie` returns the algebra and dimension predicted by
:cite:t:`Kazi_2025`; Figure 8 of :cite:t:`Shaya_2025` shows this agreement across :math:`n = 4` to
:math:`24`.

.. seealso::
   For the algorithm behind :code:`get_algebra`, see :doc:`classification`.