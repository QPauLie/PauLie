Average Out-of-Time-Order Correlator
====================================

The out-of-time-order correlator (OTOC) is a four-point correlator between two unitary operators :math:`W` and :math:`V` defined as

.. math::

    F(W, V) = \frac{1}{d}\text{tr}\left[W^{\dagger}V^{\dagger}WV\right]

where :math:`d` is the total Hilbert space dimension. For two initially commuting operators, the OTOC is equal to one.

If we Heisenberg-evolve one of the operators under some dynamics as :math:`V_t = U^{\dagger}VU` where :math:`U=e^{-iHt}`, the operators :math:`W` and :math:`V_t` begin to fail to commute as the operator :math:`V` gets scrambled under the dynamics. The decay of the OTOC :math:`F(W, V_t)` is a well-studied probe of operator scrambling in quantum systems.

Given a Pauli string dynamical Lie algebra (DLA), we can compute the average OTOC between two Pauli strings where the Heisenberg-evolution is averaged over the dynamics using the function `average_otoc`:

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

As we can see, after Heisenberg-evolution, two initially commuting operators may fail to commute after Heisenberg-evolution under some dynamics.
