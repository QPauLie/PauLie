"""
    Test second-order twirl routine.
"""
import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.application.second_moment import second_moment

generator_list = [
    ["I"], ["X"], ["Y"], ["Z"],
    ["ZI", "IZ", "XX"],
    ["XI", "IX", "XY"],
    ["XX", "YY", "ZZ"],
]

@pytest.mark.parametrize("generators", generator_list)
def test_twirl_in_quadsym_basis_is_equal(generators: list[str]):
    """
    Test that the twirl of a quadratic symmetry is equal to the quadratic symmetry.

    second_moment(quadratic_symmetry) == quadratic_symmetry
    """
    gens = p(generators)
    quadsyms = gens.get_quadratic_symmetries()
    for quadsym in quadsyms:
        twirl = second_moment(gens, quadsym)
        print(quadsym)
        print(twirl)
        assert twirl == quadsym
