"""
Unit tests for the average_otoc application function.
"""
import unittest

from paulie.application.otoc import average_otoc
from paulie.common.pauli_string_factory import get_pauli_string as p


class TestAverageOTOC(unittest.TestCase):
    """
    Test suite for the average Out-of-Time-Order Correlator function.
    """

    def setUp(self):
        """Set up a non-trivial system for multiple tests."""
        self.system = p(["XXI", "IYY", "ZIZ"])
        self.op_v = p("IZI")
        self.op_w = p("YXI")

    def test_mismatched_qubit_count_raises_value_error(self):
        """Tests that a ValueError is raised for inconsistent operator sizes."""
        system = p(["ZI"])  # n=2
        op_v = p("X")       # n=1
        op_w = p("II")      # n=2
        with self.assertRaises(ValueError):
            average_otoc(op_v, op_w, system)

    def test_otoc_with_identity_is_one(self):
        """
        Tests that the average OTOC with the Identity operator is always 1.
        """
        identity_op = p("III")
        # Test both as static and evolved operator
        self.assertAlmostEqual(average_otoc(self.op_v, identity_op, self.system), 1.0)
        self.assertAlmostEqual(average_otoc(identity_op, self.op_v, self.system), 1.0)

    def test_operators_in_same_component_is_one(self):
        """
        Tests that the average OTOC is 1 if operators are in the same component.
        """
        system = p(["Z"], n=3)
        op_v = p("XII")
        op_w = p("YII")
        self.assertAlmostEqual(average_otoc(op_v, op_w, system), 1.0)

    def test_known_value_for_simple_system(self):
        """
        Tests the core counting logic against a simple, manually verifiable case.
        This provides a robust check of the implementation's correctness.
        """
        # System: n=2, Generator G = {ZI}
        system = p(["ZI"])
        static_op = p("IX")
        initial_evolved_op = p("XI")

        # Manual Derivation:
        # 1. Component of evolved_op (XI) is C(V) = {XI, YI}. Size = 2.
        # 2. Check anticommutations with static_op (IX):
        #    - IX and XI commute (0 non-identity overlaps).
        #    - IX and YI commute (0 non-identity overlaps).
        # 3. Anticommuting count = 0.
        # 4. OTOC = 1 - 2 * (0 / 2) = 1.0.
        expected_otoc = 1.0
        actual_otoc = average_otoc(static_op, initial_evolved_op, system)
        self.assertAlmostEqual(actual_otoc, expected_otoc)

    def test_corollary_4_symmetry_property(self):
        """
        Tests the fundamental property E[F(V, W_t)] == E[F(W, V_t)], as proven
        in Corollary 4. This is a robust check of the logic.
        """
        # Use the Matchgate system again with different operators
        mg_system = p(["Z", "XX"], n=3)
        op_1 = p("XII")   # In C_1
        op_2 = p("ZZI")   # In C_2
        otoc_12 = average_otoc(op_1, op_2, mg_system)
        otoc_21 = average_otoc(op_2, op_1, mg_system)

        # Assert that the two values are equal, confirming the symmetry
        self.assertAlmostEqual(otoc_12, otoc_21, delta=1e-9)
        self.assertLess(otoc_12, 1.0) # Ensure it's not a trivial case
