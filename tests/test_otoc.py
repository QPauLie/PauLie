"""
Tests for OTOC: Haar-averaged ``average_otoc``, Pauli indexing helpers, and fixed-unitary
quantities (``otoc_fixed_unitary``, ``mean_abs_otoc_uniform``, ``pauli_instability``).
"""

import networkx as nx
import numpy as np
import pytest
from scipy.stats import unitary_group

from paulie import (
    PauliString,
    PauliStringCollection,
    average_otoc,
    get_identity,
    get_pauli_string as p,
    mean_abs_otoc_uniform,
    otoc_fixed_unitary,
    pauli_instability,
)
from collections.abc import Callable

from paulie.common.pauli_string_factory import pauli_string_from_index

generators_list = [
    ["X"],
    ["Y"],
    ["Z"],
    ["XX", "YY", "ZZ"],
    ["ZI", "IZ", "XX"],
    ["XI", "IX", "XX", "YY"],
    ["XY", "YX", "YZ", "ZY"],
    ["XI", "IX", "YI", "IY", "ZZ"],
    ["XI", "ZZ", "YI", "IY", "XY", "YX"],
    ["XI", "IX", "YI", "IY", "ZI", "IZ", "XX"],
]


def naive_otoc(
    generators: PauliStringCollection, v: PauliString, w: PauliString
) -> float:
    """
    Computes the Haar averaged out-of-time-order correlator between two Pauli strings
    using NetworkX.

    We can compute this as
    1 - 2 * |{W, P} = 0 : P in connected component of V| / |connected component of V|
    where we take the commutator graph. (arXiV:2502.16404)

    Args:
        generators: Generating set of the Pauli string DLA.
        v: Pauli string V
        w: Pauli string W
    """
    # Generate commutator graph
    vertices, edges = generators.get_commutator_graph()
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    graph.add_edges_from(edges)
    # Get connected component of V
    v_connected_component = nx.node_connected_component(graph, str(v))
    # Count the number of elements t in the connected component of V
    # that anticommute with W
    anticommute_count = sum(not w | p(t) for t in v_connected_component)
    return 1 - 2 * anticommute_count / len(v_connected_component)


@pytest.mark.parametrize("generators", generators_list)
def test_average_otoc_matches_naive(generators: list[str]) -> None:
    """
    Test that average_otoc(g, v, w) == naive_otoc(g, v, w).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    for v in all_paulis:
        for w in all_paulis:
            assert average_otoc(g, v, w) == pytest.approx(naive_otoc(g, v, w))


@pytest.mark.parametrize("generators", generators_list)
def test_average_otoc_is_symmetric(generators: list[str]) -> None:
    """
    Test that average_otoc(g, v, w) == average_otoc(g, w, v).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    for v in all_paulis:
        for w in all_paulis:
            assert average_otoc(g, v, w) == pytest.approx(average_otoc(g, w, v))


@pytest.mark.parametrize("generators", generators_list)
def test_average_eq_initial_otoc_for_commutants(generators: list[str]) -> None:
    """
    Test that the average_otoc(g, v, w) == tr[W @ V @ W @ V] / d if V is a
    commutant of the DLA (since it would commute with the evolution unitary by
    the Baker-Campbell-Haussdorf formula).
    """
    g = p(generators)
    i = get_identity(len(generators[0]))
    all_paulis = i.get_commutants()
    commutants = g.get_commutants()
    d = 2 ** len(generators[0])
    for v in commutants:
        vmat = v.get_matrix()
        for w in all_paulis:
            wmat = w.get_matrix()
            analytical_value = np.trace(wmat @ vmat @ wmat @ vmat) / d
            assert average_otoc(g, v, w) == pytest.approx(analytical_value)


def test_su_otoc():
    """
    for n in range(4, 12):
    """
    n = 4
    g_su = p(["XX", "XY", "YZ"], n=n)
    d = 2**n
    min_val = -1 / (d**2 - 1)
    i = get_identity(n)
    all_paulis = i.get_commutants()
    del all_paulis[0]
    for v, w in zip(all_paulis, all_paulis):
        assert pytest.approx(min_val) == average_otoc(g_su, v, w)


@pytest.mark.parametrize("n", [1, 2])
def test_pauli_string_from_index_round_trip(n: int) -> None:
    num = 4**n
    for i in range(num):
        s = pauli_string_from_index(i, n)
        assert s.get_index() == i
        assert len(s) == n


@pytest.mark.parametrize(
    ("index", "n", "expected"),
    [
        (0, 2, "II"),
        (2, 2, "IX"),
        (8, 2, "XI"),
    ],
)
def test_pauli_string_from_index_string_form(index: int, n: int, expected: str) -> None:
    assert str(pauli_string_from_index(index, n)) == expected


@pytest.mark.parametrize(
    ("index", "n"),
    [
        (-1, 1),
        (4, 1),
        (1, 0),
    ],
)
def test_pauli_string_from_index_errors(index: int, n: int) -> None:
    with pytest.raises(ValueError):
        pauli_string_from_index(index, n)


@pytest.mark.parametrize("i", list(range(4)))
@pytest.mark.parametrize("j", list(range(4)))
def test_otoc_fixed_unitary_identity_n1_all_pairs_abs_one(i: int, j: int) -> None:
    eye = np.eye(2, dtype=np.complex128)
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


@pytest.mark.parametrize(
    ("u", "pair", "match"),
    [
        pytest.param(
            np.array([[1.0, 2.0], [0.0, 1.0]], dtype=np.complex128),
            (lambda: p("X"), lambda: p("I")),
            "unitary",
            id="non_unitary",
        ),
        pytest.param(
            np.eye(2, dtype=np.complex128),
            (
                lambda: PauliString(pauli_str="II"),
                lambda: PauliString(pauli_str="II"),
            ),
            "length",
            id="dimension_mismatch",
        ),
    ],
)
def test_otoc_fixed_unitary_validation_errors(
    u: np.ndarray,
    pair: tuple,
    match: str,
) -> None:
    f1, f2 = pair
    with pytest.raises(ValueError, match=match):
        otoc_fixed_unitary(u, f1(), f2())


def test_otoc_fixed_unitary_n0() -> None:
    u = np.array([[1.0]], dtype=np.complex128)
    p0 = PauliString(n=0)
    assert otoc_fixed_unitary(u, p0, p0) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "u",
    [
        pytest.param(np.array([[1.0]], dtype=np.complex128), id="n0"),
        pytest.param(np.eye(2, dtype=np.complex128), id="n1"),
        pytest.param(np.eye(4, dtype=np.complex128), id="n2"),
    ],
)
def test_mean_abs_otoc_uniform_identity_exact(u: np.ndarray) -> None:
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


@pytest.mark.parametrize(
    "fn",
    [mean_abs_otoc_uniform, pauli_instability],
    ids=["mean_abs_otoc_uniform", "pauli_instability"],
)
def test_exact_method_rejects_large_n(fn) -> None:
    eye64 = np.eye(64, dtype=np.complex128)  # n=6, above exact-mode cap
    with pytest.raises(ValueError, match=r"method='exact' supports at most \d+ qubits"):
        fn(eye64, method="exact")


def test_mean_abs_otoc_uniform_monte_carlo_allows_five_qubits() -> None:
    eye32 = np.eye(32, dtype=np.complex128)
    m = mean_abs_otoc_uniform(
        eye32, method="monte_carlo", num_samples=20, check_unitary=False
    )
    assert m == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("u", "base", "abs_tol"),
    [
        (np.eye(2, dtype=np.complex128), None, 1e-12),
        (np.eye(4, dtype=np.complex128), 2.0, 1e-12),
    ],
    ids=["natural_log", "log_base_2"],
)
def test_pauli_instability_identity_vanishes(
    u: np.ndarray,
    base: float | None,
    abs_tol: float,
) -> None:
    if base is None:
        out = pauli_instability(u, method="exact")
    else:
        out = pauli_instability(u, method="exact", base=base)
    assert out == pytest.approx(0.0, abs=abs_tol)


@pytest.mark.parametrize(
    ("n", "rng_seed"),
    [(2, 0), (3, 2)],
    ids=["n2", "n3"],
)
def test_mean_abs_otoc_uniform_exact_matches_bruteforce_random_unitary(
    n: int,
    rng_seed: int,
) -> None:
    rng = np.random.default_rng(rng_seed)
    d = 2**n
    u = np.asarray(unitary_group.rvs(d, random_state=rng), dtype=np.complex128)
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


def test_pauli_instability_hadamard_near_zero() -> None:
    """Single-qubit Hadamard: uniform mean of ``|OTOC|`` is 1, so ``I(H)`` vanishes."""
    h = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2)
    assert mean_abs_otoc_uniform(h, method="exact") == pytest.approx(1.0, abs=1e-12)
    assert pauli_instability(h, method="exact") == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize(
    ("u_factory", "fn", "num_samples", "seed"),
    [
        pytest.param(
            lambda: np.eye(32, dtype=np.complex128),
            mean_abs_otoc_uniform,
            30,
            7,
            id="mean_abs_identity_n5",
        ),
        pytest.param(
            lambda: np.asarray(
                unitary_group.rvs(4, random_state=np.random.default_rng(0)),
                dtype=np.complex128,
            ),
            mean_abs_otoc_uniform,
            200,
            12345,
            id="mean_abs_haar_n2",
        ),
        pytest.param(
            lambda: np.asarray(
                unitary_group.rvs(4, random_state=np.random.default_rng(1)),
                dtype=np.complex128,
            ),
            pauli_instability,
            150,
            99,
            id="pauli_instability_haar_n2",
        ),
    ],
)
def test_monte_carlo_seed_reproducible(
    u_factory: Callable[[], np.ndarray],
    fn,
    num_samples: int,
    seed: int,
) -> None:
    """Repeated calls with the same ``seed`` return the same Monte Carlo estimate."""
    u = u_factory()
    kw = dict(
        method="monte_carlo",
        num_samples=num_samples,
        seed=seed,
        check_unitary=False,
    )
    a = fn(u, **kw)
    b = fn(u, **kw)
    assert a == pytest.approx(b, rel=0, abs=0)
