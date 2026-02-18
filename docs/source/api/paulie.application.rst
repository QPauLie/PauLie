.. module:: paulie.application

Application (:mod:`paulie.application`)
=======================================

Visualization
-------------
Utilities for visualizing the anti-commutation graph.

.. autosummary::
   :toctree: generated/

   animation.animation_anti_commutation_graph
   plot.plot_anti_commutation_graph

Average Pauli weight
--------------------
Utilities related to average Pauli weights.

.. autosummary::
   :toctree: generated/

   average_pauli_weight.average_pauli_weight
   average_pauli_weight.avg_pauli_weights
   average_pauli_weight.avg_pauli_weights_from_strings
   average_pauli_weight.get_pauli_weights
   average_pauli_weight.quantum_fourier_entropy

Charges
-------
Utilities for finding charges.

.. autosummary::
   :toctree: generated/

   charges.non_commuting_charges

Four-point correlation functions
--------------------------------
Utilities for finding four-point correlation functions.

.. autosummary::
   :toctree: generated/

   fourpoint.fourpoint

Optimal :math:`\mathfrak{su}(2^n)` generators
---------------------------------------------
Utilities related to :math:`\mathfrak{su}(2^n)` generators.

.. autosummary::
   :toctree: generated/

   get_optimal_su2_n.get_optimal_connections_su_2_n
   get_optimal_su2_n.get_optimal_su_2_n_generators

DLA metrics
-----------
Utilities related to different DLA properties.

.. autosummary::
   :toctree: generated/

   graph_complexity.average_graph_complexity
   otoc.average_otoc
   second_moment.second_moment

Matrix decomposition
--------------------
Utilities related to matrix decompositions.

.. autosummary::
   :toctree: generated/

   matrix_decomposition.matrix_decomposition
   matrix_decomposition.matrix_decomposition_diagonal

Pauli compiler
--------------
Utilities related to Pauli compilation.

.. autosummary::
   :toctree: generated/
   
   pauli_compiler.OptimalPauliCompiler
   pauli_compiler.PauliCompilerConfig
   pauli_compiler.SubsystemCompiler
   pauli_compiler.SubsystemCompilerConfig
   pauli_compiler.choose_u_for_b
   pauli_compiler.compile_target
   pauli_compiler.construct_universal_set
   pauli_compiler.left_a_minimal
   pauli_compiler.left_map_over_a
   

