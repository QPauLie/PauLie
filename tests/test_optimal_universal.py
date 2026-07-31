"""Tests for optimal universal generator set construction."""

from __future__ import annotations

import math

import pytest

from paulie import get_optimal_universal_generators
from paulie.application.get_optimal_su2_n import get_optimal_edges_su_2_n


def su_dim(n: int) -> int:
    """Return the dimension of su(2^n)."""
    return 4**n - 1


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


@pytest.mark.parametrize("n", list(range(2, 16)))
def test_get_optimal_universal_generators(n: int) -> None:
    """The generated set is minimal, universal, and hits the target fraction.

    These three are the function's contract, checked directly rather than
    against the a12 construction (which is only universal for n >= 4).
    """
    g = get_optimal_universal_generators(n)

    assert g is not None
    # minimal: 2n+1 generators (n=1 is the special {X, Z} case, not tested here)
    assert len(g) == 2 * n + 1
    # universal: the dynamical Lie algebra is all of su(2^n)
    assert g.get_dla_dim() == su_dim(n)
    # optimal fraction: hits the target anticommuting-pair count
    assert g.get_anticommutation_pair() == get_optimal_edges_su_2_n(len(g))
