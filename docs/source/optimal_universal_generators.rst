Optimal universal generator sets
======================

Here we demonstarte the feature `get_optimal_su_2_n_generators` with which we obtain generating
sets of :math:`\mathfrak{su}(2^n)` that have optimal generation rate.
To this end, we start with some universal generator set and find the graph with an anticommutation
fraction nearest to :math:`0.706`.

.. code-block:: python

    from paulie.application.get_optimal_su2_n import get_optimal_su_2_n_generators
    from paulie.common.two_local_generators import G_LIE
    from paulie.common.pauli_string_factory import get_pauli_string as p

    n = 4  # dimension of the system
    initial_generators = p(G_LIE["a12"], n=n)
    optimal_generators = get_optimal_su_2_n_generators(initial_generators)
    print(f" {optimal_generators} fraction={optimal_generators.get_anticommutation_fraction()}")


which outputs:

.. code-block:: bash

     XZYX,ZZXZ,IYYY,XZYY,IXZX,YYXY,YYYY,ZZIX,XXZY fraction=0.6944444444444444
