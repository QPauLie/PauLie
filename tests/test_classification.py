"""
Test classification
"""
from operator import itemgetter
from typing import Any
import pytest
from paulie import (
    G_LIE,
    two_local_algebras,
    get_pauli_string as p
)


@pytest.mark.parametrize("generators_to_compare", [
    ([["XY", "XZ"], ["IX", "XY"]]), #Example III.1 Wie+24
    ([["XX", "YZ"], ["YY", "ZX"]]), #Example III.4
    ([["ZZ", "YX", "XY"],["XX", "YZ", "ZY"],["YY", "ZX", "XZ"]]), #Example III.5
    (itemgetter("a2", "a4")(G_LIE)), #so(n)+so(n)
    (itemgetter("a11", "a16")(G_LIE)), #so(2^n)
    (itemgetter("a6", "a7", "a10")(G_LIE)), #su(2^(n-1))
    (itemgetter("a13", "a20", "a15")(G_LIE)), #su(2^(n-1))+su(2^(n-1))
    (itemgetter("a12", "a17", "a18", "a19", "a21", "a22")(G_LIE)), #su(2^n)

])
def test_multiple_algebra_equivalences(generators_to_compare:Any) -> None:
    """Test algebra equivalences"""
    algebra = []
    for generator_set in generators_to_compare:
        algebra.append(p(generator_set, n=4).get_algebra()) # equivalences hold for n>=3

    assert all(x==algebra[0] for x in algebra), "Expected to represent the same algebra"

#Consider isomorphisms and improve is_algebra
def test_explicit_algebras() -> None:
    """
    Test explicit algebras
    """
    assert p(["XX", "YY", "ZZ", "ZY"]).is_algebra("u(1)+2*su(2)") # Example III.8
    assert p(["XY"]).is_algebra("u(1)") # Example I.4
    assert p(["XY"], n = 3).is_algebra("so(3)")
    # n>=3 for n = 3 eg. a6 should be su(4)
    # but Theorem 2 does not possibly lead to this result
    for n in [6,10]:
        algs = two_local_algebras(n)
        for name, alg in algs.items():
            assert p(G_LIE[name], n = n).is_algebra(alg)
