Average Out-of-Time-Order Correlator (OTOC)
###########################################

This page explains how to compute the average Out-of-Time-Order Correlator (OTOC),
a key diagnostic for quantum chaos and information scrambling.

Theoretical Background
======================

The OTOC measures how much two initially commuting operators fail to commute
after one of them has been evolved in time. Its decay over time is often taken as a
signature of quantum chaos.

While the OTOC for a specific evolution can be complex, its *average* over an entire
family of dynamics (an ensemble) can be calculated with graph theory. As shown in
`arXiv:2502.16404 <https://arxiv.org/abs/2502.16404>`_ (Corollary 3), this average is
determined by the structure of the **commutator graph** for the system.

The formula is given by:

.. math::

   \mathbb{E}_{U \sim G} [F(W, V_t)] = 1 - \frac{2|\{T \in C(V) : \{W, T\} \neq 0\}|}{|C(V)|}

where :math:`F(W, V_t)` is the OTOC between a static operator :math:`W` and a
time-evolved operator :math:`V_t`. :math:`C(V)` is the connected component of the
initial operator :math:`V` in the commutator graph, and :math:`\{W, T\}` denotes
the anticommutator.

``paulie`` provides the `average_otoc` function to compute this value directly.

Usage Example
=============

To compute the average OTOC, you need to define your system's generators and the
two operators. This example is constructed so the operators are in different
components, leading to a non-trivial result.

.. code-block:: python

   from paulie.common.pauli_string_factory import get_pauli_string as p
   from paulie.application.otoc import average_otoc

   # Define a system with two non-communicating sectors
   system = p(["XII", "IIZ"], n=3)

   # W: An operator that spans both sectors
   static_op = p("XIX")

   # V: An operator that lives only in the first sector
   initial_evolved_op = p("YII")

   # The component of V is {YII, ZII}. W anticommutes with both.
   # Expected OTOC = 1 - 2 * (2 / 2) = -1.0
   otoc_value = average_otoc(static_op, initial_evolved_op, system)
   print(f"The average OTOC is: {otoc_value:.4f}")

**Output:**

.. code-block:: text

   The average OTOC is: -1.0000