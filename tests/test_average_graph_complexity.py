"""
    Tests for average graph complexity
"""
from math import comb
from itertools import combinations
import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_factory import get_identity
from paulie.application.graph_complexity import average_graph_complexity
from paulie.common.pauli_string_collection import PauliString

def matchgate_generators(n: int) -> list[str]:
    """
    Construct the "matchgate" generators (arXiV:2502.16404 eq.(44)).
    """
    c = []
    for i in range(n):
        c.append("I" * i + "Z" + "I" * (n - i - 1))
    for i in range(n - 1):
        c.append("I" * i + "XX" + "I" * (n - i - 2))
    return c

def majorana_operators(n: int) -> list[PauliString]:
    """
    Construct the Majorana operators (arXiV:2502.16404 eq.(44)).
    """
    c = []
    for i in range(n):
        c.append(p("Z" * i + "X" + "I" * (n - i - 1)))
        c.append(p("Z" * i + "Y" + "I" * (n - i - 1)))
    return c

@pytest.mark.parametrize("n", [2, 3])
def test_matchgate_average_graph_complexity(n: int) -> None:
    """
    Test the average graph complexities for the matchgate set (arXiV:2502.16404 eq.(C3)).
    """
    gens = matchgate_generators(n)
    ops = majorana_operators(n)
    mg = p(gens, n=n)
    for k in range(1, 2 * n + 1):
        connected_component_size = comb(2 * n, k)
        for i_a in combinations(range(1, 2 * n + 1), k):
            # The Pauli string is given by c[i_a[0]] * c[i_a[1]] * ...
            c = get_identity(n=n)
            for i in i_a:
                c @= ops[i - 1]
            # Get computed value
            computed_val = average_graph_complexity(mg, c)
            # Get analytical value
            analytical_val = 0
            for j_a in combinations(range(1, 2 * n + 1), k):
                analytical_val += sum(abs(i - j) for i, j in zip(i_a, j_a))
            analytical_val /= connected_component_size
            assert computed_val == pytest.approx(analytical_val)


@pytest.mark.parametrize("generators", [
    ["X"],
    ["ZZ", "XX", "YY"],
    matchgate_generators(2),
    matchgate_generators(3),
])
def test_dla_commutants_have_zero_graph_complexity(generators: list[str]) -> None:
    """
    Test that the commutants of a DLA have zero average graph complexity.
    """
    gens = p(generators)
    commutants = gens.get_commutants()
    for c in commutants:
        assert average_graph_complexity(gens, c) == pytest.approx(0.0)
