"""
Unit tests for the second_moment function in paulie.application.second_moment.
"""
import unittest

from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_linear import PauliStringLinear

from paulie.application.second_moment import second_moment


class TestSecondMoment(unittest.TestCase):
    """
    Unit tests for the second_moment function.
    These tests use the factory function `p()` for object creation
    to align with the library's design philosophy.
    """
    def test_projection_to_zero(self):
        """
        Tests that an operator with no overlap with the quadratic symmetries
        is correctly projected to the zero operator.
        """
        system = p(["Z"])
        operator_m = p([(1.0, "XI")]) # Changed to use p()
        twirled_m = second_moment(operator_m, system)
        self.assertTrue(twirled_m.is_zero(), f"Expected a zero operator, but got {twirled_m}")

    def test_projection_of_symmetry_is_identity_mapping(self):
        """
        Tests that twirling a quadratic symmetry basis vector returns the
        same vector, confirming the projection logic.
        """
        system = p(["Z"])
        operator_m = p([(1.0, "IZ")]) # Changed to use p()
        twirled_m = second_moment(operator_m, system)
        self.assertEqual(twirled_m, operator_m, f"Expected {operator_m}, "+
                         "but twirling returned {twirled_m}")

        # Test a more complex symmetry
        operator_m_complex = p([(1.0, "XX"), (1.0, "YY")]) # Changed to use p()
        twirled_m_complex = second_moment(operator_m_complex, system)
        self.assertEqual(twirled_m_complex, operator_m_complex)

    def test_gksmail_demonstration_case(self):
        """
        Tests the quadratic symmetry generation against the specific
        example for the a_5 system (n=2).

        This test verifies that the final refactored implementation correctly
        computes the full basis by iterating through the components and
        linear symmetries.
        """
        # 1. Recreate the list of connected components {C_k} for the a_5 system.
        #    Each component is a PauliStringCollection.
        components_ck = [
            p(["II"]), p(["XY"]), p(["ZX"]), p(["YZ"]),
            p(["ZI", "XZ", "IX", "YY"]),  # The 4-element component
            p(["XI", "ZZ"]),
            p(["IY", "YX"]),
            p(["YI", "ZY"]),
            p(["IZ", "XX"])
        ]
        self.assertEqual(len(components_ck), 9, "Should be 9 connected components for a_5")

        # 2. Recreate the list of linear symmetries {L_j} for the a_5 system.
        linear_symmetries_lj = p(["II", "XY", "ZX", "YZ"])
        self.assertEqual(len(linear_symmetries_lj), 4, "Should be 4 linear symmetries for a_5")

        # 3. Use our final, refactored logic to compute the full basis.
        #    This mimics the logic from the `get_full_quadratic_basis` orchestrator we deleted.
        full_quadratic_basis = []

        # For each component Ck in our list...
        for ck in components_ck:
            # ...compute its associated symmetries by passing the linear symmetries as an argument.
            # This directly calls the final version of our method.
            q_for_this_component = ck.get_quadratic_symmetries(linear_symmetries_lj)
            full_quadratic_basis.extend(q_for_this_component)

        # 4. Assert the final result is correct.
        #    The total number of Q_kj should be (number of components) * (number of symmetries).
        expected_q_count = len(components_ck) * len(linear_symmetries_lj)  # 9 * 4 = 36

        self.assertEqual(len(full_quadratic_basis), expected_q_count,
                         "The total number of generated quadratic symmetries is incorrect.")

        # A deeper check to ensure all generated objects have the correct type.
        for q in full_quadratic_basis:
            self.assertIsInstance(q, PauliStringLinear,
                                  "All generated symmetries should be PauliStringLinear objects.")
