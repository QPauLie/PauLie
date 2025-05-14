"""
Test of charges
"""
import pytest
from paulie.common.two_local_generators import get_lie_algebra
from paulie.application.charges import non_commuting_charges
from paulie.common.pauli_string_factory import get_pauli_string as p

@pytest.mark.parametrize("algebra_name, should_have_charges", [
    ("a8", True),
    ("a9", True),
    # Add more algebras to test as needed
])
def test_algebras_with_non_commuting_charges(algebra_name:str, should_have_charges:bool) -> None:
    """
    Test that specific algebras possess non-commuting charges.
    
    Args:
        algebra_name: The name of the Lie algebra to test
        should_have_charges: Whether the algebra should have non-commuting charges
    """
    # Get the generators for the specified algebra
    generators = p(get_lie_algebra(algebra_name))
    # Find non-commuting charges
    charges = non_commuting_charges(generators)
    # Check if non-commuting charges exist as expected
    if should_have_charges:
        assert len(charges) > 0, (f"Expected {algebra_name} to "
               "have non-commuting charges, but found none")
    else:
        assert len(charges) == 0, (f"Expected {algebra_name} to "
               f"have no non-commuting charges, but found {len(charges)}")
