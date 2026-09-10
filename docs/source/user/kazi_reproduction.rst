Multi-angle MaxCut QAOA
=======================

This tutorial illustrates how to use :code:`paulie` to classify the dynamical Lie algebra (DLA)
of the *multi-angle* MaxCut Quantum Approximate Optimization Algorithm (QAOA) ansatz, and how the
resulting classification reproduces the six-family result of :cite:t:`Kazi_2025`. The accompanying
example (:code:`docs/examples/kazi_reproduction.py`) runs the full check across graph families and
regenerates the scaling figure reported in :cite:t:`Shaya_2025`.

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

Building :math:`\mathfrak{g}_\mathrm{free}` by brute-force Lie closure costs time and memory that
are exponential in :math:`n`. :code:`paulie` avoids this by contracting the anticommutation graph of
the generators (which preserves the Lie algebra) into a canonical form. The algebra and its
dimension are then determined from that form. It does not require generating the exponentially
large set of elements.

Example: classifying the DLA of a graph
---------------------------------------

Take the path graph :math:`P_3` on vertices :math:`0\text{--}1\text{--}2`. Its generators are
:math:`X_0, X_1, X_2` and :math:`Z_0 Z_1, Z_1 Z_2`:

.. code-block:: python

    from paulie import get_pauli_string as p

    generators = p([
        "XII", "IXI", "IIX",   # X_v on each vertex
        "ZZI", "IZZ",          # Z_u Z_v on each edge
    ])
    print(f"algebra   = {generators.get_algebra()}")
    print(f"dimension = {generators.get_dla_dim()}")

outputs

.. code-block:: bash

    algebra   = so(6)
    dimension = 15

This matches the closed form for path graphs, :math:`\mathfrak{g}_\mathrm{free} \cong
\mathfrak{so}(2n)` of dimension :math:`2n^2 - n = 15` at :math:`n = 3`.

.. seealso::
   For the algorithm behind :code:`get_algebra`, see :doc:`classification`.

The six connected graph families and their Lie algebras
-------------------------------------------------------

:cite:t:`Kazi_2025` prove that, for *any* connected graph, the multi-angle DLA falls into exactly
one of six families, fixed by simple combinatorial properties of the graph. A graph is called
*archetypal* when it is connected but neither bipartite nor a cycle graph; the even--even,
odd--odd and even--odd labels refer to the parities of the two parts of a connected bipartite
graph's bipartition.

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

Running the reproduction
------------------------

The reproduction is provided as the example script :code:`docs/examples/kazi_reproduction.py`. Run it
with :code:`paulie` installed (Python 3.12 or newer):

.. code-block:: bash

    pip install paulie          # networkx and matplotlib are also used by this example

    python docs/examples/kazi_reproduction.py   # sweeps n = 4,6,...,24

It checks one representative of each family (:math:`P_8`, :math:`C_8`, :math:`K_{2,4}`,
:math:`K_{3,3}`, :math:`K_{2,3}`, :math:`K_5`), sweeps :math:`n` to show how the dimension scales, and writes
:code:`kazi_reproduction.csv` together with the log-scale figure :code:`kazi_reproduction.pdf`
(the plot reported in :cite:t:`Shaya_2025`).