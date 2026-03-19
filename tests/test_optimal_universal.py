from __future__ import annotations

import math
import random

import pytest

from paulie import G_LIE, get_optimal_su_2_n_generators, get_pauli_string as p
from paulie.application.get_optimal_su2_n import get_optimal_edges_su_2_n


def su_dim(n: int) -> int:
    return 4**n - 1


def growth_score_5_rounds(generators, rounds: int = 5) -> int:
    """
    Early-growth score used for the slow benchmark:
    A_0 = G,
    A_{r+1} = A_r union {[g, a] != 0 : g in G, a in A_r}.

    We count distinct non-identity Paulis reached up to the given depth.
    """
    generators = list(generators)
    discovered = {str(op): op for op in generators}
    frontier = list(discovered.values())

    for _ in range(rounds):
        new_frontier = {}
        for g in generators:
            for a in frontier:
                comm = g ^ a
                if comm is None or comm.is_identity():
                    continue
                key = str(comm)
                if key not in discovered:
                    discovered[key] = comm
                    new_frontier[key] = comm
        if not new_frontier:
            break
        frontier = list(new_frontier.values())

    return len(discovered)


@pytest.fixture(params=[4, 5])
def universal_a12(request):
    """
    Universal fixture from the public API:
        p(G_LIE["a12"], n=n)
    """
    n = request.param
    g = p(G_LIE["a12"], n=n)
    assert g.get_dla_dim() == su_dim(n)
    return n, g


@pytest.mark.parametrize(
    ("ng", "expected"),
    [
        (0, -1),
        (1, -1),
        (2, 0),
        (5, math.floor(0.706 * 10)),
        (7, math.floor(0.706 * 21)),
        (9, math.floor(0.706 * 36)),
    ],
)
def test_get_optimal_edges_su_2_n_formula(ng: int, expected: int) -> None:
    assert get_optimal_edges_su_2_n(ng) == expected


def test_get_optimal_su_2_n_generators_hits_target_edge_count(universal_a12) -> None:
    n, g = universal_a12

    g_ind = g.copy().get_independents()
    target_edges = get_optimal_edges_su_2_n(len(g_ind))

    random.seed(target_edges)
    g_opt = get_optimal_su_2_n_generators(g)

    assert g_opt is not None
    assert len(g_opt) == len(g_ind)
    assert g_opt.get_dla_dim() == su_dim(n)
    assert g_opt.get_anticommutation_pair() == target_edges


@pytest.mark.slow
def test_get_optimal_su_2_n_generators_is_locally_best_for_5_round_growth(
    universal_a12,
) -> None:
    """
    Paper-style regression test for 'fast early growth': compare the 5-round
    commutator-growth score in a small window around the target edge count.
    """
    n, g = universal_a12

    g_ind = g.copy().get_independents()
    target_edges = get_optimal_edges_su_2_n(len(g_ind))

    random.seed(target_edges)
    g_opt = get_optimal_su_2_n_generators(g)

    assert g_opt is not None
    assert g_opt.get_dla_dim() == su_dim(n)

    target_score = growth_score_5_rounds(g_opt, rounds=5)

    total_pairs = g_ind.get_pair()
    window = min(4, total_pairs)
    scores_by_realized_edges = {}

    for desired_edges in range(
        max(0, target_edges - window),
        min(total_pairs, target_edges + window) + 1,
    ):
        cand = g_ind.find_generators_with_connection(desired_edges)
        realized_edges = cand.get_anticommutation_pair()
        score = growth_score_5_rounds(cand, rounds=5)
        scores_by_realized_edges[realized_edges] = max(
            score,
            scores_by_realized_edges.get(realized_edges, -1),
        )

    realized_target_edges = g_opt.get_anticommutation_pair()
    assert realized_target_edges in scores_by_realized_edges

    best_score = max(scores_by_realized_edges.values())
    tol = max(1, math.ceil(0.1 * best_score))
    assert target_score >= best_score - tol
