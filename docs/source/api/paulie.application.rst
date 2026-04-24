.. module:: paulie.application

Application (:mod:`paulie.application`)
=======================================

Visualization
-------------
Utilities for visualizing the anti-commutation graph.

.. autosummary::
   :toctree: generated/

   plot.plot_anti_commutation_graph

Average Pauli weight
--------------------
Utilities related to average Pauli weights.

.. autosummary::
   :toctree: generated/

   average_pauli_weight.average_pauli_weight
   average_pauli_weight.get_pauli_weights
   average_pauli_weight.quantum_fourier_entropy


Optimal :math:`\mathfrak{su}(2^n)` generators
---------------------------------------------
Utilities related to :math:`\mathfrak{su}(2^n)` generators.

.. autosummary::
   :toctree: generated/

   get_optimal_su2_n.get_optimal_edges_su_2_n
   get_optimal_su2_n.get_optimal_su_2_n_generators

Chaos metrics
-------------
Utilities related to different DLA properties.

.. autosummary::
   :toctree: generated/

   average_graph_complexity.average_graph_complexity
   otoc.average_otoc
   otoc.fourpoint
   second_moment.second_moment

Pauli instability and fixed-unitary OTOC
-----------------------------------------
Expectations over uniform Pauli pairs for a given unitary (contrast with Haar-averaged
:func:`~paulie.application.otoc.average_otoc` over a DLA).

.. autosummary::
   :toctree: generated/

   otoc.otoc_fixed_unitary
   otoc.mean_abs_otoc_uniform
   otoc.pauli_instability

Matrix decomposition
--------------------
Utilities related to matrix decompositions.

.. autosummary::
   :toctree: generated/

   matrix_decomposition.matrix_decomposition
   matrix_decomposition.matrix_decomposition_diagonal
