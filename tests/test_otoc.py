"""
Unit tests for the average_otoc application function.
"""
import unittest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.application.otoc import average_otoc


class TestOTOC(unittest.TestCase):
    """
    Test suite for the average Out-of-Time-Order Correlator function.
    """

    def test_corollary_4_symmetry(self):
        """Tests that OTOC(V, W) == OTOC(W, V), a key symmetry."""
        # A non-trivial system from the paper's Table I
        system = p(["XX", "YY", "X"], n=3)
        op_v = p("IZI")
        op_w = p("YXI")
        otoc_vw = average_otoc(op_v, op_w, system)
        otoc_wv = average_otoc(op_w, op_v, system)
        self.assertAlmostEqual(otoc_vw, otoc_wv, delta=1e-9)

    def test_otoc_with_identity_is_one(self):
        """Tests that the average OTOC with the Identity operator is 1."""
        system = p(["Z", "XX"], n=3)
        op_v = p("XXI")
        identity_op = p("III")
        self.assertAlmostEqual(average_otoc(op_v, identity_op, system), 1.0, delta=1e-9)
        self.assertAlmostEqual(average_otoc(identity_op, op_v, system), 1.0, delta=1e-9)

    def test_operators_in_same_component_is_one(self):
        """Tests that OTOC is 1 if operators are in the same component."""
        system = p(["Z", "XX"], n=3)
        op_v = p("XII")
        op_w = p("YII")
        otoc_val = average_otoc(op_v, op_w, system)
        self.assertAlmostEqual(otoc_val, 1.0, delta=1e-9)

    def test_mismatched_qubit_count_raises_value_error(self):
        """
        Covers: raise ValueError(...)
        Tests that a ValueError is raised if operators have inconsistent sizes.
        """
        system = p(["ZI"])  # n=2
        op_v = p("X")       # n=1
        op_w = p("II")      # n=2
        # Use assertRaises to confirm that the expected error is thrown.
        with self.assertRaises(ValueError):
            average_otoc(op_v, op_w, system)
