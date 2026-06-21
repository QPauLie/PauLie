Pauli compiler
==============

In quantum computing, direct access to a target operation is usually not
available due to a fixed set of native gates determined by the physical
architecture of the device. This requires decomposing the target operation
in terms of those available operators before it can be executed.

Constructing such decompositions manually is often non-trivial, especially for systems with many qubits.
The Pauli compiler addresses this problem by automatically finding a **constructive decomposition**
of a target Pauli string into a sequence of available generators. This removes the need for manual
algebraic derivations and provides a direct route from high-level targets to executable quantum circuits.

Core idea
---------

A target operation :math:`P` is considered reachable if it can be generated from the available system
operators :math:`A_m` via nested commutators. The compiler searches for a sequence

.. math::

    (A_{i_1}, A_{i_2}, \dots, A_{i_k})

such that :math:`P` can be expressed (up to proportionality) as

.. math::

    [A_{i_1}, [A_{i_2}, [\dots, [A_{i_{k-1}}, A_{i_k}] \dots ]]] \propto P

This provides a constructive synthesis of the target operator from the available algebra.

Algorithm
---------

To compile a target Pauli string :math:`P`, the string is partitioned into left and right subsystem components.
Depending on this structure, the compiler applies one of three cases:

1. **Left-only** (:math:`V \neq I, W = I`): The target has support only on the left subsystem.
   A BFS over adjoint maps finds a path from a seed generator to the target.

2. **Both sides** (:math:`V \neq I, W \neq I`): When the target has support on both subsystems,
   the right subsystem component is first compiled using a subsystem compiler, yielding undesired
   operator :math:`V'` on the left subsystem. This operator
   is then corrected via BFS over adjoint mapping.

3. **Right-only** (:math:`V = I, W \neq I`): The target acts only on the right subsystem and
   is decomposed as :math:`W \propto [W_1, W_2]`. Each :math:`W_1` and :math:`W_2` is compiled
   separately using a subsystem compiler, yielding operators :math:`V_1' \otimes W_1` and
   :math:`V_2' \otimes W_2`. A sequence mapping :math:`V_1'` and :math:`V_2'` via adjoint maps is
   then found via BFS. Finally, a permutation :math:`\sigma` of this generator sequence that
   yields a nonzero adjoint chain is selected, giving the desired compilation.

In all cases, the resulting sequence has :math:`\mathcal{O}(N)` length.

Quick start
-----------

For 4-qubits system, if the target Pauli string is IZXI, the simplest way to use the compiler is through the
:func:`~paulie.application.pauli_compiler.compile_target` function:

.. code-block:: python

    from paulie import compile_target, get_pauli_string as p

    target = p("IZXI")
    sequence = compile_target(target, k_left=2)
    print(f"Sequence length: {len(sequence)}")

In the example, the target Pauli string is split into a left partition on ``k_left`` qubits
and a right part on the remaining qubits.
The algorithm requires that the number of qubits in the left partition be greater than or equal to 2.
The :func:`~paulie.application.pauli_compiler.compile_target` returns a list of Pauli operators forming
a sequence that reconstructs the target operator.

To verify the result:

.. code-block:: python

    from paulie import PauliStringCollection

    result = PauliStringCollection(sequence[:-1]).nested_adjoint(sequence[-1])
    print(f"Target:  {target}")
    print(f"Result:  {result}")
    assert result == target

