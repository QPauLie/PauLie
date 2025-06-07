Average Out-of-Time-Order Correlator (OTOC)
###########################################

This page explains how to compute the average Out-of-Time-Order Correlator (OTOC),
a key diagnostic for quantum chaos and information scrambling.

Theoretical Background
======================

The OTOC is a measure of how much two initially commuting operators fail to commute
after one of them has been evolved in time. Its decay is often taken as a signature
of quantum chaos.

While the OTOC for a specific evolution can be complex, its *average* over an entire
family of dynamics (an ensemble) can be calculated with graph theory. As shown in
`arXiv:2502.16404 <https://arxiv.org/abs/2502.16404>`_ (Corollary 3), this average is
determined by the structure of the **commutator graph** for the system.

The formula is given by:

.. math::

   \mathbb{E}_{U \sim G} [F(V, U^\dagger W U)] = 1 - \frac{2|\{T \in C(V) : \{W, T\} = 0\}|}{|C(V)|}

where :math:`C(V)` is the connected component of the operator :math:`V` in the commutator
graph, and :math:`\{W, T\}` denotes the anticommutator.

``PauLie`` provides a function to compute this value directly.

Usage Example
=============

To compute the average OTOC, you need to define your system's generators and the
two Pauli string operators you wish to compare.

.. code-block:: python

   from paulie.common.pauli_string_factory import get_pauli_string as p
   from paulie.application.OTOC import OTOC

   # Define the system generators for the "matchgate" model on 3 qubits
   system = p(["Z", "XX"], n=3)
   
   # Define two Pauli operators to compare
   V = p("XXI")
   W = p("XYI")

   # Compute the average OTOC between V and W for the given system
   average_otoc = OTOC(V, W, system)

   print(f"The average OTOC between {V} and {W} is: {average_otoc}")

   # The OTOC with the Identity operator should always be 1
   I = p("III")
   otoc_with_identity = OTOC(V, I, system)
   print(f"The average OTOC between {V} and {I} is: {otoc_with_identity}")

**Expected Output:**

.. code-block:: text

   The average OTOC between XXI and XYI is: 0.8
   The average OTOC between XXI and III is: 1.0

*(Note: The value 0.8 is an example; the actual calculated value will depend on the specific graph structure for n=3.)*
