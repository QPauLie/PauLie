import pytest
from paulie.application.classify import get_algebra
from paulie.common.pauli_string_factory import get_pauli_string as p

@pytest.mark.parametrize("generators1, generators2, should_equal", [
    (["XY", "XZ"], ["IX", "XY"], True), #Example III.1 Wie+24
    (["XX", "YZ"],["YY", "ZX"], True), #Example III.4
    (["ZZ", "YX", "XY"],["XX", "YZ", "ZY"], True), #Example III.5
    (["ZZ", "YX", "XY"], ["YY", "ZX", "XZ"], True),
    # Add more test cases here as needed
])
def test_multiple_algebra_equivalences(generators1, generators2, should_equal):
    """Test multiple cases of algebra classifications using parametrization."""
    algebra1 = get_algebra(p(generators1))
    algebra2 = get_algebra(p(generators2))
    
    if should_equal:
        assert algebra1 == algebra2, f"Expected {generators1} and {generators2} to represent the same algebra"
    else:
        assert algebra1 != algebra2, f"Expected {generators1} and {generators2} to represent different algebras"

def test_explicit_algebras():
    assert p(["XX", "YY", "ZZ", "ZY"]).get_class().is_algebra("u(1)+2*so(2)") # Example III.8
    assert p(["XY"]).get_class().is_algebra("u(1)") # Example I.4
    assert p(["XY"], n = 3).get_class().is_algebra("so(3)")





