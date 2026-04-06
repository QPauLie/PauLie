"""
Tests for Pauli indexing and fixed-unitary OTOC.
"""
import numpy as np
import pytest

from paulie import PauliString, get_pauli_string as p
from paulie.application.pauli_instability import (
    mean_abs_otoc_uniform,
    otoc_fixed_unitary,
    pauli_instability,
)
from paulie.common.pauli_string_factory import pauli_string_from_index


def _haar_unitary(rng: np.random.Generator, d: int) -> np.ndarray:
    """Haar-random ``d`` x ``d`` unitary (QR with phase fix)."""
    z = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    phases = diag / np.abs(diag)
    return (q * phases).astype(np.complex128)


def test_pauli_string_from_index_round_trip_n1() -> None:
    for i in range(4):
        s = pauli_string_from_index(i, 1)
        assert s.get_index() == i
        assert len(s) == 1


def test_pauli_string_from_index_round_trip_n2() -> None:
    for i in range(16):
        s = pauli_string_from_index(i, 2)
        assert s.get_index() == i
        assert len(s) == 2


def test_pauli_string_from_index_matches_string_literals() -> None:
    assert str(pauli_string_from_index(0, 2)) == "II"
    assert str(pauli_string_from_index(2, 2)) == "IX"
    assert str(pauli_string_from_index(8, 2)) == "XI"


def test_pauli_string_from_index_errors() -> None:
    with pytest.raises(ValueError):
        pauli_string_from_index(-1, 1)
    with pytest.raises(ValueError):
        pauli_string_from_index(4, 1)
    with pytest.raises(ValueError):
        pauli_string_from_index(1, 0)


def test_otoc_fixed_unitary_identity_n1_all_pairs_abs_one() -> None:
    d = 2
    eye = np.eye(d, dtype=np.complex128)
    for i in range(4):
        for j in range(4):
            p1 = pauli_string_from_index(i, 1)
            p2 = pauli_string_from_index(j, 1)
            val = otoc_fixed_unitary(eye, p1, p2)
            assert abs(val - round(val.real)) < 1e-10
            assert abs(val) == pytest.approx(1.0)


def test_otoc_fixed_unitary_hadamard_n1_explicit() -> None:
    """Compare to explicit 2x2 multiply for one (P1, P2) pair."""
    h = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
    p_x = p("X")
    p_z = p("Z")
    p1_m = p_x.get_matrix().astype(np.complex128)
    p2_m = p_z.get_matrix().astype(np.complex128)
    ud = h.conj().T
    m = ud @ p1_m @ h @ p2_m @ ud @ p1_m @ h @ p2_m
    expected = np.trace(m) / 2
    got = otoc_fixed_unitary(h, p_x, p_z)
    assert got == pytest.approx(expected)


def test_otoc_fixed_unitary_rejects_non_unitary() -> None:
    bad = np.array([[1.0, 2.0], [0.0, 1.0]], dtype=np.complex128)
    with pytest.raises(ValueError, match="unitary"):
        otoc_fixed_unitary(bad, p("X"), p("I"))


def test_otoc_fixed_unitary_dimension_mismatch() -> None:
    eye2 = np.eye(2, dtype=np.complex128)
    with pytest.raises(ValueError, match="length"):
        otoc_fixed_unitary(eye2, PauliString(pauli_str="II"), PauliString(pauli_str="II"))


def test_otoc_fixed_unitary_n0() -> None:
    u = np.array([[1.0]], dtype=np.complex128)
    p0 = PauliString(n=0)
    assert otoc_fixed_unitary(u, p0, p0) == pytest.approx(1.0)


def test_mean_abs_otoc_uniform_identity_exact_n1() -> None:
    eye = np.eye(2, dtype=np.complex128)
    assert mean_abs_otoc_uniform(eye, method="exact") == pytest.approx(1.0)


def test_mean_abs_otoc_uniform_identity_exact_n2() -> None:
    eye = np.eye(4, dtype=np.complex128)
    assert mean_abs_otoc_uniform(eye, method="exact") == pytest.approx(1.0)


def test_mean_abs_otoc_uniform_identity_exact_n0() -> None:
    u = np.array([[1.0]], dtype=np.complex128)
    assert mean_abs_otoc_uniform(u, method="exact") == pytest.approx(1.0)


def test_mean_abs_otoc_uniform_identity_monte_carlo_is_one() -> None:
    """For ``U = I``, every pair gives ``|OTOC| = 1``, so MC is exact for any sample count."""
    eye = np.eye(4, dtype=np.complex128)
    m = mean_abs_otoc_uniform(
        eye, method="monte_carlo", num_samples=17, check_unitary=True
    )
    assert m == pytest.approx(1.0)


def test_mean_abs_otoc_uniform_invalid_method() -> None:
    eye = np.eye(2, dtype=np.complex128)
    with pytest.raises(ValueError, match="Invalid method"):
        mean_abs_otoc_uniform(eye, method="not_a_method")  # type: ignore[arg-type]


def test_mean_abs_otoc_uniform_exact_rejects_more_than_four_qubits() -> None:
    eye32 = np.eye(32, dtype=np.complex128)  # n=5
    with pytest.raises(ValueError, match="at most 4 qubits"):
        mean_abs_otoc_uniform(eye32, method="exact")


def test_mean_abs_otoc_uniform_monte_carlo_allows_five_qubits() -> None:
    eye32 = np.eye(32, dtype=np.complex128)
    m = mean_abs_otoc_uniform(
        eye32, method="monte_carlo", num_samples=20, check_unitary=False
    )
    assert m == pytest.approx(1.0)


def test_pauli_instability_exact_rejects_more_than_four_qubits() -> None:
    eye32 = np.eye(32, dtype=np.complex128)
    with pytest.raises(ValueError, match="at most 4 qubits"):
        pauli_instability(eye32, method="exact")


def test_mean_abs_otoc_uniform_matches_bruteforce_random_unitary_n2() -> None:
    rng = np.random.default_rng(0)
    u = _haar_unitary(rng, 4)
    n = 2
    total = 0.0
    for i in range(4**n):
        p1 = pauli_string_from_index(i, n)
        for j in range(4**n):
            p2 = pauli_string_from_index(j, n)
            o = otoc_fixed_unitary(u, p1, p2, check_unitary=False)
            total += abs(complex(o))
    brute = total / (16**n)
    got = mean_abs_otoc_uniform(u, method="exact", check_unitary=False)
    assert got == pytest.approx(brute, rel=0, abs=1e-9)


def test_pauli_instability_identity_exact_is_zero() -> None:
    eye = np.eye(2, dtype=np.complex128)
    assert pauli_instability(eye, method="exact") == pytest.approx(0.0, abs=1e-12)


def test_pauli_instability_identity_log_base_2() -> None:
    eye = np.eye(4, dtype=np.complex128)
    assert pauli_instability(eye, method="exact", base=2.0) == pytest.approx(0.0, abs=1e-12)


def test_pauli_instability_hadamard_near_zero() -> None:
    """Single-qubit Hadamard: uniform mean of ``|OTOC|`` is 1, so ``I(H)`` vanishes."""
    h = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
    assert mean_abs_otoc_uniform(h, method="exact") == pytest.approx(1.0, abs=1e-12)
    assert pauli_instability(h, method="exact") == pytest.approx(0.0, abs=1e-10)


def test_mean_abs_otoc_uniform_monte_carlo_seed_reproducible_large_n() -> None:
    """``n > 4`` uses the per-sample path; seed still fixes the estimate."""
    eye32 = np.eye(32, dtype=np.complex128)
    a = mean_abs_otoc_uniform(
        eye32,
        method="monte_carlo",
        num_samples=30,
        seed=7,
        check_unitary=False,
    )
    b = mean_abs_otoc_uniform(
        eye32,
        method="monte_carlo",
        num_samples=30,
        seed=7,
        check_unitary=False,
    )
    assert a == pytest.approx(b, rel=0, abs=0)


def test_mean_abs_otoc_uniform_monte_carlo_seed_reproducible() -> None:
    rng = np.random.default_rng(0)
    u = _haar_unitary(rng, 4)
    a = mean_abs_otoc_uniform(
        u,
        method="monte_carlo",
        num_samples=200,
        seed=12345,
        check_unitary=False,
    )
    b = mean_abs_otoc_uniform(
        u,
        method="monte_carlo",
        num_samples=200,
        seed=12345,
        check_unitary=False,
    )
    assert a == pytest.approx(b, rel=0, abs=0)


def test_pauli_instability_monte_carlo_seed_reproducible() -> None:
    rng = np.random.default_rng(1)
    u = _haar_unitary(rng, 4)
    a = pauli_instability(
        u, method="monte_carlo", num_samples=150, seed=99, check_unitary=False
    )
    b = pauli_instability(
        u, method="monte_carlo", num_samples=150, seed=99, check_unitary=False
    )
    assert a == pytest.approx(b, rel=0, abs=0)


def test_mean_abs_otoc_uniform_exact_matches_bruteforce_random_unitary_n3() -> None:
    """Optimized exact path agrees with naive ``otoc_fixed_unitary`` double sum."""
    rng = np.random.default_rng(2)
    u = _haar_unitary(rng, 8)
    n = 3
    total = 0.0
    for i in range(4**n):
        p1 = pauli_string_from_index(i, n)
        for j in range(4**n):
            p2 = pauli_string_from_index(j, n)
            o = otoc_fixed_unitary(u, p1, p2, check_unitary=False)
            total += abs(complex(o))
    brute = total / (16**n)
    got = mean_abs_otoc_uniform(u, method="exact", check_unitary=False)
    assert got == pytest.approx(brute, rel=0, abs=1e-9)
