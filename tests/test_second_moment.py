import unittest
# Import the classes we need to construct objects directly
from paulie.common.pauli_string_linear import PauliStringLinear
from paulie.common.pauli_string_collection import PauliStringCollection
# We still need the factory for simple generators
from paulie.common.pauli_string_factory import get_pauli_string as p
# Import the function we are testing
from paulie.application.second_moment import second_moment

class TestSecondMoment(unittest.TestCase):

    def test_projection_to_zero(self):
        """
        Tests that an operator with no overlap with the quadratic symmetries
        is correctly projected to the zero operator.
        """
        system = PauliStringCollection(p(["Z"]))
        # Construct M directly to ensure it has the correct type
        M = PauliStringLinear([(1.0, "XI")])
        twirled_M = second_moment(M, system)
        self.assertTrue(twirled_M.is_zero(), 
                        f"Expected a zero operator, but got {twirled_M}")

    def test_projection_of_symmetry_is_identity_mapping(self):
        """
        Tests that twirling a quadratic symmetry basis vector returns the
        same vector, confirming the projection logic.
        """
        system = PauliStringCollection(p(["Z"]))
        # Construct M directly
        M = PauliStringLinear([(1.0, "IZ")])
        twirled_M = second_moment(M, system)
        self.assertEqual(twirled_M, M,
                         f"Expected {M}, but twirling returned {twirled_M}")

        # Test a more complex symmetry
        M_complex = PauliStringLinear([(1.0, "XX"), (1.0, "YY")])
        twirled_M_complex = second_moment(M_complex, system)
        self.assertEqual(twirled_M_complex, M_complex)

    def test_second_moment_skips_on_zero_denominator(self):
        """
        Tests that the function gracefully handles a quadratic symmetry
        that results in a zero-valued trace denominator, triggering the
        `continue` statement and preventing a ZeroDivisionError.
        """
        system = PauliStringCollection(p(["Z"]))
        M = PauliStringLinear([(1.0, "XI")])
        
        zero_term_q = PauliStringLinear([])
        overlapping_q = PauliStringLinear([(1.0, "XI")])
        
        # Our handcrafted basis for the test
        fake_q_basis = [zero_term_q, overlapping_q]
        
        # Temporarily replace the real get_quadratic_symmetries with our fake one
        system.get_quadratic_symmetries = lambda: fake_q_basis
        twirled_M = second_moment(M, system)
        
        # The result should be `overlapping_q` itself, as we are projecting it
        # onto a basis that contains only it and a term that should be skipped.
        self.assertEqual(twirled_M, overlapping_q,
                         "Function did not correctly skip the zero-denominator term.")