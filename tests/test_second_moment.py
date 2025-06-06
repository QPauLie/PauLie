"""
Unit tests for the second_moment function in paulie.application.second_moment.
"""
import unittest

# Import the classes we need to construct objects directly
from paulie.common.pauli_string_collection import PauliStringCollection

# We still need the factory for simple generators
from paulie.common.pauli_string_factory import get_pauli_string as p

# Import the function we are testing
from paulie.application.second_moment import second_moment


class TestSecondMoment(unittest.TestCase):
    """
    Unit tests for the second_moment function.
    """
    def test_projection_to_zero(self):
        """
        Tests that an operator with no overlap with the quadratic symmetries
        is correctly projected to the zero operator.
        """
        system = p(["Z"])
        # Construct operator_m directly to ensure it has the correct type
        operator_m = p([(1.0, "XI")])
        twirled_m = second_moment(operator_m, system)
        self.assertTrue(twirled_m.is_zero(), f"Expected a zero operator, but got {twirled_m}")

    def test_projection_of_symmetry_is_identity_mapping(self):
        """
        Tests that twirling a quadratic symmetry basis vector returns the
        same vector, confirming the projection logic.
        """
        system = p(["Z"])
        # Construct operator_m directly
        operator_m = p([(1.0, "IZ")])
        twirled_m = second_moment(operator_m, system)
        self.assertEqual(twirled_m, operator_m, 
                         f"Expected {operator_m}, but twirling returned {twirled_m}")

        # Test a more complex symmetry, again using the factory
        operator_m_complex = p([(1.0, "XX"), (1.0, "YY")])
        twirled_m_complex = second_moment(operator_m_complex, system)
        self.assertEqual(twirled_m_complex, operator_m_complex)

    
    
