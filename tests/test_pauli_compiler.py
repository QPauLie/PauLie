"""Tests for the Pauli Compiler implemented in PauLie.

Validates that the returned sequence G satisfies
    nested_adjoint(G[:-1], G[-1]) == target
for a selection of small k, N, and target strings.
"""
from paulie.application.pauli_compiler import (
    PauliCompilerConfig,
    OptimalPauliCompiler,
    _evaluate_paulie_orientation,
    compile_target,
    construct_universal_set,
)
from paulie.common.pauli_string_factory import get_pauli_string


def _assert_compiles(target_str: str, k_left: int) -> None:
    target = get_pauli_string(target_str)
    sequence = compile_target(target, k_left=k_left)

    assert sequence, f"Empty sequence for target={target_str}, k_left={k_left}"

    result = _evaluate_paulie_orientation(sequence)
    assert result is not None
    assert result == target, f"Compilation failed for target={target_str}, got={result}"


def test_compile_target_smoke_cases() -> None:
    cases = [
        (2, 3, "XII"),
        (2, 3, "IIY"),
        (2, 3, "YII"),
        (2, 4, "IZXI"),
        (2, 4, "IIYZ"),
        (3, 5, "IIZXI"),
        (3, 5, "IIIYZ"),
    ]

    for k_left, n_total, target_str in cases:
        assert len(target_str) == n_total
        _assert_compiles(target_str, k_left)


def test_class_compile_api() -> None:
    k_left, n_total = 2, 4
    compiler = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n_total))

    v_left = get_pauli_string("IZ")
    w_right = get_pauli_string("XI")
    sequence = compiler.compile(v_left, w_right)

    target = get_pauli_string("IZXI")
    result = _evaluate_paulie_orientation(sequence)

    assert result is not None
    assert result == target


def test_universal_set_size_minimal() -> None:
    for k_left, n_total in [(2, 3), (2, 4), (3, 5)]:
        universal_set = construct_universal_set(n_total, k_left)
        assert len(universal_set) == 2 * n_total + 1