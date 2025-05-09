Commutator graph
======================

This tutorial will illustrate how to use :code:`paulie` to analyse short term dynamics of quantum systems.
The Hamiltonian of the system is described by a linear combination of Pauli strings :math:`H = \sum_{P \in \mathcal{G}} c_P P` where
:math:`\mathcal{G}` is the generator set. The commutator graph of an :math:`n`-qubit system has the :math:`4^n` Pauli strings as
vertices and there exists and edge :math:`(P,Q)` if there exists a generator :math:`G \in \mathcal{G}` s.t.
:math:`[P,G] \prop Q`.
While the anticommutation graph is isomorphic for generator sets with isomorphic DLA's and only capture the
long term dynamics, the commutator graph captures short term dynamics.

.. code-block:: python
    from paulie.common.pauli_string_factory import get_pauli_string as p
    generators = p(["XI", "ZI", "IX", "IZ", "ZZ"], n=3)
    vertices, edges = generators.get_commutator_graph()




