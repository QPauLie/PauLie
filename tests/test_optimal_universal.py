"""Tests for optimal universal generator set construction."""

from __future__ import annotations

import math
import random

import pytest

from paulie import G_LIE,  get_optimal_universal_generators, get_pauli_string as p
from paulie.application.get_optimal_su2_n import get_optimal_edges_su_2_n


def su_dim(n: int) -> int:
    """Return the dimension of su(2^n)."""
    return 4**n - 1


@pytest.fixture(params=list(range(2, 16)))
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
    """Verify the edge-count formula for small inputs."""
    assert get_optimal_edges_su_2_n(ng) == expected


def test_get_optimal_universal_generators_hits_target_edge_count(universal_a12) -> None:
    """Optimal generators should achieve the target anticommutation edge count."""
    n, g = universal_a12

    g_ind = g.copy().get_independents()
    target_edges = get_optimal_edges_su_2_n(len(g_ind))

    random.seed(target_edges)
    g_opt = get_optimal_universal_generators(n)

    assert g_opt is not None
    assert len(g_opt) == len(g_ind)
    assert g_opt.get_dla_dim() == su_dim(n)
    assert g_opt.get_anticommutation_pair() == target_edges