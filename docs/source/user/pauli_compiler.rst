Pauli compiler
==============

The Pauli compiler implements the algorithm from :cite:t:`Smith_2025` to compile a target
Pauli string into a sequence of generators that produces the target via nested commutators.
The resulting sequence has length :math:`\mathcal{O}(N)` for :math:`N` qubits.

Overview
--------
In a closed quantum system, the evolution of the state :math:`|\psi\rangle` is governed by the time dependent
Schrödinger equation:

.. math::
   i \hbar \frac{d}{dt} |\psi\rangle = H |\psi\rangle

where :math:`H` is the Hamiltonian, which can be discretized into a sequence of operators :math:`H_m`,
each associated with a discrete time step :math:`t_m`. Each :math:`H_m` can be expressed as
a sum of Pauli strings acting on the qubits in the system.
By solving the Schrödinger equation and assuming each :math:`H_m` is time-independent,
we obtain the target state :math:`|\psi_1\rangle` from the initial state :math:`|\psi_0\rangle` as

.. math::
    |\psi_1\rangle = e^{-iH_m t_m} \;\cdots\; e^{-iH_1 t_1} |\psi_0\rangle

Using Baker-Campbell-Hausdorff (BCH) formula :cite:t:`Hall2013LieGroups` and defining :math:`A_m = -iH_m t_m`,
this product of exponentials can be expressed as a single exponential whose generator is built from
the system operators :math:`A_m` through commutators.

In quantum computing, a target Pauli string :math:`P` is considered reachable if it can be generated from the
available system operators :math:`A_m`. The Pauli compiler searches for a sequence
:math:`(A_{i_1}, A_{i_2}, \dots, A_{i_k})` such that :math:`P` can be expressed
(up to proportionality) as a nested commutator,

.. math::

    [A_{i_1}, [A_{i_2}, [\dots, [A_{i_{k-1}}, A_{i_k}] \dots ]]] \propto P

This provides a sequence of Pauli operations in a circuit that implement target operator using a fixed generator set
through nested commutators (adjoint maps).

Algorithm
---------

To compile a target Pauli string :math:`P`, the string is partitioned into left and right subsystem components.
Depending on this structure, the compiler applies one of three cases:

1. **Left-only** (:math:`V \neq I, W = I`): The target has support only on the left subsystem.
   A BFS over adjoint maps finds a path from a seed generator to the target.

2. **Both sides** (:math:`V \neq I, W \neq I`): When the target has support on both subsystems,
   the right subsystem component is first compiled using a subsystem compiler. This introduces an undesired operator
   on the left subsystem, which is then corrected via adjoint mapping.

3. **Right-only** (:math:`V = I, W \neq I`): The target acts only on the right subsystem and
   is decomposed as :math:`W = W_1 W_2` with :math:`[W_1, W_2] \neq 0`.
   Each factor is compiled separately and the sequences are interleaved to cancel the left residue.
   If deterministic methods fail, a bounded BFS is used as fallback.

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

Using the class-based API
-------------------------

For more control, use the :class:`~paulie.application.pauli_compiler.OptimalPauliCompiler`
class directly:

.. code-block:: python

    from paulie import (
        OptimalPauliCompiler,
        PauliCompilerConfig,
        get_pauli_string as p,
        PauliStringCollection,
    )

    config = PauliCompilerConfig(k_left=2, n_total=4)
    compiler = OptimalPauliCompiler(config)

    v_left = p("IZ")
    w_right = p("XI")
    sequence = compiler.compile(v_left, w_right)

    target = p("IZXI")
    result = PauliStringCollection(sequence[:-1]).nested_adjoint(sequence[-1])
    assert result == target

The :class:`~paulie.application.pauli_compiler.PauliCompilerConfig` accepts additional
parameters for the fallback search:

- ``fallback_depth`` (default 8): Maximum BFS depth for the fallback search.
- ``fallback_nodes`` (default 200000): Maximum number of nodes explored in the fallback BFS.
