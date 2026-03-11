"""Tests for the Pauli Compiler implemented in PauLie.

Validates that the returned sequence G satisfies
    nested_adjoint(G[:-1], G[-1]) == target
for a selection of small k, N, and target strings.
"""

import pytest

from paulie.application.pauli_compiler import (
    PauliCompilerConfig,
    OptimalPauliCompiler,
    SubsystemCompiler,
    SubsystemCompilerConfig,
    _evaluate_paulie_orientation,
    _evaluate_sequence,
    _sequence_to_paulie_orientation,
    compile_target,
    construct_universal_set,
    left_a_minimal,
    left_map_over_a,
    choose_u_for_b,
)
from paulie.common.pauli_string_factory import get_pauli_string, get_identity, get_single


def _assert_compiles(target_str: str, k_left: int) -> None:
    """Assert that the compiler produces a valid sequence for the given target."""
    target = get_pauli_string(target_str)
    sequence = compile_target(target, k_left=k_left)

    assert sequence, f"Empty sequence for target={target_str}, k_left={k_left}"

    result = _evaluate_paulie_orientation(sequence)
    assert result is not None
    assert result == target, f"Compilation failed for target={target_str}, got={result}"


# ---------------------------------------------------------------------------
# compile_target smoke tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "k_left, n_total, target_str",
    [
        # Case 1: left-only targets (W = I)
        (2, 3, "XII"),
        (2, 3, "YII"),
        (2, 4, "XZII"),
        # Case 2: both sides (V != I, W != I)
        (2, 4, "IZXI"),
        (3, 5, "IIZXI"),
        # Case 3: right-only targets (V = I, W != I)
        (2, 3, "IIY"),
        (2, 4, "IIYZ"),
        (3, 5, "IIIYZ"),
    ],
)
def test_compile_target(k_left: int, n_total: int, target_str: str) -> None:
    """Test compilation across all three cases."""
    assert len(target_str) == n_total
    _assert_compiles(target_str, k_left)


@pytest.mark.parametrize("target_str", ["IIYZ", "IIXY"])
def test_compile_target_multi_pauli_right(target_str: str) -> None:
    """Multi-Pauli targets on the right subsystem (n=4, k=2)."""
    _assert_compiles(target_str, k_left=2)


@pytest.mark.parametrize("target_str", ["XII", "ZII", "YII", "XZI", "ZXI"])
def test_compile_target_left_only(target_str: str) -> None:
    """Targets supported only on the left partition."""
    _assert_compiles(target_str, k_left=2)


# ---------------------------------------------------------------------------
# Class-based API
# ---------------------------------------------------------------------------


def test_compiler_class_api() -> None:
    """Test OptimalPauliCompiler with explicit left/right factors."""
    k_left, n_total = 2, 4
    compiler = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n_total))

    v_left = get_pauli_string("IZ")
    w_right = get_pauli_string("XI")
    sequence = compiler.compile(v_left, w_right)

    target = get_pauli_string("IZXI")
    result = _evaluate_paulie_orientation(sequence)

    assert result is not None
    assert result == target


def test_compiler_left_only() -> None:
    """Compile a target with identity on the right."""
    k_left, n_total = 2, 3
    compiler = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n_total))

    v_left = get_pauli_string("XZ")
    w_right = get_identity(1)
    sequence = compiler.compile(v_left, w_right)

    result = _evaluate_paulie_orientation(sequence)
    assert result is not None
    assert result == get_pauli_string("XZI")


def test_compiler_right_only() -> None:
    """Compile a target with identity on the left."""
    k_left, n_total = 2, 4
    compiler = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n_total))

    v_left = get_identity(k_left)
    w_right = get_pauli_string("YZ")
    sequence = compiler.compile(v_left, w_right)

    result = _evaluate_paulie_orientation(sequence)
    assert result is not None
    assert result == get_pauli_string("IIYZ")


def test_compiler_k_left_too_small() -> None:
    """k_left < 2 must raise ValueError."""
    with pytest.raises(ValueError, match="k_left must be >= 2"):
        OptimalPauliCompiler(PauliCompilerConfig(k_left=1, n_total=3))


def test_compiler_custom_fallback_config() -> None:
    """Custom fallback parameters are respected."""
    cfg = PauliCompilerConfig(k_left=2, n_total=3, fallback_depth=4, fallback_nodes=1000)
    compiler = OptimalPauliCompiler(cfg)
    assert compiler.fallback_depth == 4
    assert compiler.fallback_nodes == 1000


# ---------------------------------------------------------------------------
# Universal generator set
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "k_left, n_total",
    [(2, 3), (2, 4), (3, 5), (2, 5), (3, 6)],
)
def test_universal_set_size(k_left: int, n_total: int) -> None:
    """Universal set size should be 2 * n_total + 1."""
    universal_set = construct_universal_set(n_total, k_left)
    assert len(universal_set) == 2 * n_total + 1


def test_universal_set_same_length() -> None:
    """All elements should have the same qubit count."""
    n_total, k = 5, 2
    universal_set = construct_universal_set(n_total, k)
    for ps in universal_set:
        assert len(ps) == n_total


def test_universal_set_no_duplicates() -> None:
    """Universal set should contain no duplicate elements."""
    universal_set = construct_universal_set(5, 2)
    assert len(set(str(ps) for ps in universal_set)) == len(universal_set)


def test_universal_set_invalid_k() -> None:
    """Invalid k should raise ValueError."""
    with pytest.raises(ValueError, match="Require 1 <= k < N"):
        construct_universal_set(3, 0)


def test_universal_set_k_equals_n() -> None:
    """k equal to n should raise ValueError."""
    with pytest.raises(ValueError, match="Require 1 <= k < N"):
        construct_universal_set(3, 3)


# ---------------------------------------------------------------------------
# compile_target input validation
# ---------------------------------------------------------------------------


def test_compile_target_k_left_zero() -> None:
    """k_left=0 should raise ValueError."""
    with pytest.raises(ValueError):
        compile_target(get_pauli_string("XYZ"), k_left=0)


def test_compile_target_k_left_one() -> None:
    """k_left=1 should raise ValueError (minimum is 2)."""
    with pytest.raises(ValueError):
        compile_target(get_pauli_string("XYZ"), k_left=1)


def test_compile_target_k_left_equals_n() -> None:
    """k_left equal to target length should raise ValueError."""
    with pytest.raises(ValueError):
        compile_target(get_pauli_string("XYZ"), k_left=3)


def test_compile_target_k_left_exceeds_n() -> None:
    """k_left exceeding target length should raise ValueError."""
    with pytest.raises(ValueError):
        compile_target(get_pauli_string("XY"), k_left=5)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [2, 3, 4])
def test_left_a_minimal_size(k: int) -> None:
    """Should return 2k + 1 elements."""
    result = left_a_minimal(k)
    assert len(result) == 2 * k + 1


def test_left_a_minimal_contains_all_x_z() -> None:
    """Should contain X_i and Z_i for each qubit."""
    k = 3
    result = left_a_minimal(k)
    result_strs = [str(ps) for ps in result]
    for i in range(k):
        x_str = "I" * i + "X" + "I" * (k - i - 1)
        z_str = "I" * i + "Z" + "I" * (k - i - 1)
        assert x_str in result_strs
        assert z_str in result_strs


def test_left_a_minimal_contains_z_tensor() -> None:
    """Should contain Z^{otimes k}."""
    k = 3
    result = left_a_minimal(k)
    assert get_pauli_string("Z" * k) in result


@pytest.mark.parametrize("k", [2, 3, 4])
def test_choose_u_for_b(k: int) -> None:
    """Should return X on the first qubit."""
    u = choose_u_for_b(k)
    assert len(u) == k
    assert u == get_single(k, 0, "X")


def test_left_map_identity_path() -> None:
    """Same source and target should return empty path."""
    v = get_pauli_string("XZ")
    path = left_map_over_a(v, v, left_a_minimal(2))
    assert not path


def test_left_map_finds_path() -> None:
    """Should find a non-trivial path between different Pauli strings."""
    generators = left_a_minimal(2)
    v_from = get_single(2, 0, "X")
    v_to = get_single(2, 0, "Z")
    path = left_map_over_a(v_from, v_to, generators)
    assert len(path) > 0


# ---------------------------------------------------------------------------
# Sequence orientation
# ---------------------------------------------------------------------------


def test_sequence_orientation_empty() -> None:
    """Empty input should return empty output."""
    assert not _sequence_to_paulie_orientation([])


def test_sequence_orientation_single() -> None:
    """Single element should be preserved."""
    ps = get_pauli_string("XY")
    result = _sequence_to_paulie_orientation([ps])
    assert result == [ps]


def test_sequence_orientation_reversal() -> None:
    """[base, A1, A2] -> [A2, A1, base]."""
    base = get_pauli_string("XI")
    a1 = get_pauli_string("ZI")
    a2 = get_pauli_string("IX")
    result = _sequence_to_paulie_orientation([base, a1, a2])
    assert result == [a2, a1, base]


def test_evaluate_orientations_agree() -> None:
    """Evaluating in both orientations should give the same result."""
    base = get_pauli_string("XII")
    a1 = get_pauli_string("ZII")
    sequence = [base, a1]
    result_seq = _evaluate_sequence(sequence)
    paulie_seq = _sequence_to_paulie_orientation(sequence)
    result_paulie = _evaluate_paulie_orientation(paulie_seq)
    assert result_seq == result_paulie


def test_evaluate_paulie_orientation_empty() -> None:
    """Empty sequence should return None."""
    assert _evaluate_paulie_orientation([]) is None


# ---------------------------------------------------------------------------
# SubsystemCompiler
# ---------------------------------------------------------------------------


def test_subsystem_compiler_k_left_too_small() -> None:
    """k_left < 2 should raise ValueError."""
    with pytest.raises(ValueError, match="k_left must be >= 2"):
        SubsystemCompiler(SubsystemCompilerConfig(k_left=1, n_total=3))


def test_subsystem_compiler_extend_left() -> None:
    """Extending a left Pauli string should pad with identities."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=4))
    a = get_pauli_string("XZ")
    extended = sub.extend_left(a)
    assert len(extended) == 4
    assert str(extended) == "XZII"


def test_subsystem_compiler_extend_pair() -> None:
    """Combining left and right parts should concatenate them."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=4))
    u = get_pauli_string("XI")
    b = get_pauli_string("ZY")
    combined = sub.extend_pair(u, b)
    assert str(combined) == "XIZY"


def test_factor_w_orders_single_x() -> None:
    """An X-only target should have exactly one factoring order."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=3))
    w = get_pauli_string("X")
    orders = sub.factor_w_orders(w)
    assert len(orders) == 1


def test_factor_w_orders_y_has_two() -> None:
    """A Y target should have two factoring orders (XZ and ZX)."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=3))
    w = get_pauli_string("Y")
    orders = sub.factor_w_orders(w)
    assert len(orders) == 2


def test_factor_w_orders_identity() -> None:
    """An identity target should produce empty factor lists."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=3))
    w = get_pauli_string("I")
    orders = sub.factor_w_orders(w)
    assert len(orders) == 1
    assert len(orders[0]) == 0


def test_subsystem_compiler_single_x() -> None:
    """Compile a single X on the right subsystem."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=3))
    w = get_pauli_string("X")
    result = sub.subsystem_compiler(w)
    assert len(result) > 0


def test_subsystem_compiler_identity() -> None:
    """Identity on the right subsystem should return empty sequence."""
    sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=2, n_total=3))
    w = get_pauli_string("I")
    result = sub.subsystem_compiler(w)
    assert not result
