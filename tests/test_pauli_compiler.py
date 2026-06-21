"""
Test pauli compiler
"""
from itertools import product, combinations
import pytest
from paulie import PauliString, PauliStringCollection
from paulie.application.pauli_compiler import (
    compile_target, left_a_minimal, left_map_over_a,
    SubsystemCompilerConfig, SubsystemCompiler
)

gates = ["I", "X", "Y", "Z"]
pauli_string_length = 5
test_cases = [
    "".join(p)
    for n in range(1, pauli_string_length + 1)
    for p in product(gates, repeat=n)
]

@pytest.mark.parametrize("input_arg", test_cases)
def test_compile_target(input_arg) -> None:
    """
    Test compile target for all possible combination
    """
    target = PauliString(pauli_str=input_arg)
    n = len(input_arg)

    for k in range(1, n):

        if k < 2 or k == n:
            with pytest.raises(ValueError):
                compile_target(target, k_left=k)
        else:
            sequence = compile_target(target, k_left=k)
            result = PauliStringCollection(sequence).evaluate_commutator_sequence()
            assert result == target

@pytest.mark.parametrize("input_arg", test_cases)
def test_subcompiler(input_arg) -> None:
    """
    Test subcompiler (internal) for all possible combination
    """
    target = PauliString(pauli_str=input_arg)
    n = len(input_arg)

    for k in range(1, n):

        w = target.get_substring(k, n)

        if k < 2 or k == n:
            with pytest.raises(ValueError):
                compile_target(target, k_left=k)
        else:
            sub_config = SubsystemCompilerConfig(k_left=k, n_total=n)
            sub_compiler = SubsystemCompiler(sub_config)

            g = sub_compiler.subsystem_compiler(w)
            result = PauliStringCollection(g).evaluate_commutator_sequence()

            if w.is_identity():
                assert result is None
            else:
                assert result.get_substring(k, n) == w

left_length = 3
test_adjoint_maps_strings = [
    "".join(p)
    for p in product(gates, repeat=left_length)
    if "".join(p) != left_length*"I"
]

test_adjoint_maps_cases = list(
    combinations(test_adjoint_maps_strings, 2)
)

@pytest.mark.parametrize("v_from, v_to", test_adjoint_maps_cases)
def test_left_adjoint_map(v_from, v_to) -> None:
    """
    Test left adjoint map (internal) for all possible combination
    """
    pauli_str_v_from = PauliString(pauli_str=v_from)
    pauli_str_v_to = PauliString(pauli_str=v_to)
    n = len(v_to)
    generators = left_a_minimal(n)

    sequence = left_map_over_a(pauli_str_v_from, pauli_str_v_to, generators)
    sequence.append(pauli_str_v_from)
    result = PauliStringCollection(sequence).evaluate_commutator_sequence()

    assert result == pauli_str_v_to
