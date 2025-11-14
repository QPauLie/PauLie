Commutator graph
======================

This tutorial will illustrate how to use :code:`paulie` to analyse short term dynamics of quantum systems.
The Hamiltonian of the system is described by a linear combination of Pauli strings :math:`H = \sum_{P \in \mathcal{G}} c_P P` where
:math:`\mathcal{G}` is the generator set. The commutator graph of an :math:`n`-qubit system has the :math:`4^n` Pauli strings as
vertices and there exists and edge :math:`(P,Q)` if there exists a generator :math:`G \in \mathcal{G}` s.t.
:math:`[P,G] \propto Q`.
While the anticommutation graph is isomorphic for generator sets with isomorphic DLA's and only capture the
long term dynamics, the commutator graph captures short term dynamics.
We can plot the commutator graph

.. code-block:: python

    from paulie.common.pauli_string_factory import get_pauli_string as p
    generators = p(["XI", "ZI", "IX", "IZ", "ZZ"], n=2)
    vertices, edges = generators.get_commutator_graph()
    plot_graph(vertices, edges)


.. image:: media/commutator_graph.png

With the help of the commutator graph we can analyse the short-term dynamics of the quantum system.
In particular, we can determine how chaotic the system behaves. A quantifier of quantum chaos is the
the out-of-time-order correlator (OTOC) between two unitary operators :math:`W` and :math:`V` defined as

.. math::

    F(W, V) = \frac{1}{d}\text{tr}\left[W^{\dagger}V^{\dagger}WV\right]

where :math:`d` is the total Hilbert space dimension. For two initially commuting operators, the OTOC is equal to one.

If we Heisenberg-evolve one of the operators under some dynamics as :math:`V_t = U^{\dagger}VU` where :math:`U=e^{-iHt}`, the operators :math:`W` and :math:`V_t` begin to fail to commute as the operator :math:`V_t` becomes global under the dynamics.
The exponential decay of the OTOC :math:`F(W, V_t)` with time is a probe of quantum choas.

Given a Pauli string dynamical Lie algebra (DLA), we can compute the average OTOC between two Pauli strings where the Heisenberg-evolution is averaged over the dynamics using the function :code:`average_otoc`:

.. code-block:: python

    import numpy as np
    from paulie.common.pauli_string_factory import get_pauli_string as p
    from paulie.application.otoc import average_otoc

    # These are the generators of n=3 "matchgate" dynamics
    generators = p(["ZII", "IZI", "IIZ", "XXI", "IXX"])
    v = p("XYZ")
    w = p("YZZ")
    print(f"Do {v} and {w} initially commute?", v.commutes_with(w))
    print("Average OTOC of {v} and {w} under matchgate dynamics:", np.round(average_otoc(generators, v, w), 3))

which outputs:

.. code-block:: bash

    Do XYZ and YZZ initially commute? True
    Average OTOC of {v} and {w} under matchgate dynamics: -0.333

As we can see, two initially commuting operators may fail to commute after Heisenberg-evolution under some dynamics.

