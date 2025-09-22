"""Tests for the Optimal Pauli Compiler implemented in PauLie.

Validates that the returned sequence G satisfies
    nested_adjoint(G[:-1], G[-1]) == target
for a selection of small k, N, and target strings.
"""
from paulie.common.pauli_string_factory import get_pauli_string
from paulie.application.pauli_compiler import nested_adjoint
from paulie.application.pauli_compiler_optimal import (
    compile_target,
    OptimalPauliCompiler,
    PauliCompilerConfig,
    construct_universal_set,
)


def _assert_compiles(target_str: str, k: int) -> None:
    target = get_pauli_string(target_str)
    seq = compile_target(target, k_left=k)
    assert seq, f"Empty sequence for target={target_str}, k={k}"
    res = nested_adjoint(seq[:-1], seq[-1])
    assert res is not None and res == target, (
        f"Compilation failed: target={target_str}, k={k}, got={res}"
    )


def test_compile_target_smoke_cases() -> None:
    # (k, N, target)
    cases = [
        (2, 3, "XII"),    # V=X1, W=I
        (2, 3, "IIY"),    # V=I,  W=Y
        (2, 3, "YII"),    # V=Y1, W=I
        (2, 4, "IZXI"),   # V=Z2, W=XI
        (2, 4, "IIYZ"),   # V=I,  W=YZ
        (3, 5, "IIZXI"),  # V=Z3, W=XI
        (3, 5, "IIIYZ"),  # V=I,  W=YZ
    ]
    for (k, N, t) in cases:
        assert len(t) == N
        _assert_compiles(t, k)


def test_class_compile_api() -> None:
    # One class-based call to ensure V/W API works
    k, N = 2, 4
    cfg = PauliCompilerConfig(k_left=k, n_total=N)
    opc = OptimalPauliCompiler(cfg)
    # Target: V=Z2, W=XI => target "IZXI"
    V = get_pauli_string("IZ").get_substring(0, k)
    W = get_pauli_string("XI").get_substring(0, N - k)
    seq = opc.compile(V, W)
    target = get_pauli_string("IZXI")
    res = nested_adjoint(seq[:-1], seq[-1])
    assert res is not None and res == target


def test_universal_set_size_minimal() -> None:
    for (k, N) in [(2, 3), (2, 4), (3, 5)]:
        U = construct_universal_set(N, k)
        assert len(U) == 2 * N + 1

