"""
    Tests for the average out-of-time-order correlator.
"""
import pytest
import numpy as np
from paulie.application.otoc import average_otoc
from paulie.common.pauli_string_factory import (
    get_pauli_string as p,
    get_identity
)

generators_list = [
    ["I"], ["X"], ["Y"], ["Z"],
    ["XX", "YY", "ZZ"],
    ["ZI", "IZ", "XX"],
    ["XI", "IX", "XX", "YY"],
    ["XY", "YX", "YZ", "ZY"],
    ["XI", "IX", "YI", "IY", "ZZ"],
    ["XI", "ZZ", "YI", "IY", "XY", "YX"],
    ["XI", "IX", "YI", "IY", "ZI", "IZ", "XX"],
]

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
            assert average_otoc(g, v, w) == pytest.approx(average_otoc(g, v, w))

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
