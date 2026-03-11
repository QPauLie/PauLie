Pauli compiler
==============

The Pauli compiler implements the algorithm from :cite:t:`Bittel_2024` to compile a target
Pauli string into a sequence of generators that produces the target via nested commutators.
The resulting sequence has length :math:`\mathcal{O}(N)` for :math:`N` qubits.

Overview
--------

Given a set of universal generators and a target Pauli string :math:`P`, the compiler finds a
sequence :math:`(A_1, A_2, \dots, A_m)` of generators such that

.. math::

    [A_1, [A_2, [\dots, [A_{m-1}, A_m] \dots ]]] \propto P

This is useful for constructing circuits that realize arbitrary Pauli strings using only
a fixed generator set through nested commutators (adjoint maps).

Quick start
-----------

The simplest way to use the compiler is through the :func:`~paulie.application.pauli_compiler.compile_target`
function:

.. code-block:: python

    from paulie import compile_target, get_pauli_string as p

    target = p("IZXI")
    sequence = compile_target(target, k_left=2)
    print(f"Sequence length: {len(sequence)}")

The ``k_left`` parameter specifies the left-right partition of the qubit system. The target
is split into a left part on ``k_left`` qubits and a right part on the remaining qubits.
The parameter must satisfy ``k_left >= 2``.

To verify the result:

.. code-block:: python

    from paulie import PauliStringCollection

    result = PauliStringCollection(sequence[:-1]).nested_adjoint(sequence[-1])
    print(f"Target:  {target}")
    print(f"Result:  {result}")
    assert result == target

Universal generator set
-----------------------

The compiler works with a universal generator set consisting of:

- **Left generators:** :math:`\{X_i, Z_i\}_{i=1}^{k} \cup \{Z^{\otimes k}\}`, extended
  with identities on the right subsystem.
- **Right generators:** Single-qubit :math:`X_i` and :math:`Z_i` on the right subsystem,
  tagged with a fixed left Pauli string.

You can construct this set explicitly:

.. code-block:: python

    from paulie import construct_universal_set

    n_total = 5
    k_left = 2
    universal_set = construct_universal_set(n_total, k_left)
    print(f"Universal set size: {len(universal_set)}")
    # Output: Universal set size: 11  (= 2 * n_total + 1)

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

Algorithm
---------

The compiler handles three cases depending on the structure of the target:

1. **Left-only** (:math:`W = I`): The target is supported only on the left subsystem.
   A BFS over adjoint maps finds a path from a seed generator to the target.

2. **Both sides** (:math:`V \neq I, W \neq I`): The right part is compiled using a
   subsystem compiler, and the left factor is corrected via adjoint mapping.

3. **Right-only** (:math:`V = I, W \neq I`): The target is decomposed as :math:`W = W_1 W_2`
   with :math:`[W_1, W_2] \neq 0`. Each factor is compiled separately and the sequences
   are interleaved to cancel the left residue. If deterministic methods fail, a bounded
   BFS is used as fallback.

In all cases, the resulting sequence has :math:`\mathcal{O}(N)` length.
