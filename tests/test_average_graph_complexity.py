"""
    Tests for average graph complexity
"""
import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_factory import get_identity
from paulie.application.graph_complexity import average_graph_complexity
from paulie.common.pauli_string_collection import PauliString, PauliStringCollection
from math import comb
from itertools import combinations

def matchgate_generators(n: int) -> list[str]:
    """
    Construct the "matchgate" generators (arXiV:2502.16404 eq.(44)).
    """
    c = []
    for i in range(n):
        c.append('I' * i + 'Z' + 'I' * (n - i - 1))
    for i in range(n - 1):
        c.append('I' * i + 'XX' + 'I' * (n - i - 2))
    return c

def majorana_operators(n: int) -> list[PauliString]:
    """
    Construct the Majorana operators (arXiV:2502.16404 eq.(44)).
    """
    c = []
    for i in range(n):
        c.append(p('Z' * i + 'X' + 'I' * (n - i - 1)))
        c.append(p('Z' * i + 'Y' + 'I' * (n - i - 1)))
    return c

testdata = [2, 3]

@pytest.mark.parametrize("n", testdata)
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
                analytical_val += sum([abs(i - j) for i, j in zip(i_a, j_a)])
            analytical_val /= connected_component_size
            assert computed_val == pytest.approx(analytical_val)