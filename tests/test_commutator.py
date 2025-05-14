"""
Test commutator
"""
import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_bitarray import PauliString

@pytest.fixture(scope="module")
def pauli_setup() -> dict[str,PauliString]:
    """Set up Pauli operators for testing."""
    return {
        "I": p("I"),
        "X": p("X"),
        "Y": p("Y"),
        "Z": p("Z"),
        "XI": p("XI"),
        "IX": p("IX"),
        "XY": p("XY"),
        "YX": p("YX"),
        "XIIII": p("XIIII"),
        "ZIIII": p("ZIIII"),
        "IYXII": p("IYXII"),
        "IXIXI": p("IXIXI"),
        "IIIZI": p("IIIZI"),
        "IZIXZ": p("IZIXZ"),
        "IIIIX": p("IIIIX"),
        "IZXXI": p("IZXXI")
    }

def test_single_qubit_commutation(pauli_setup:dict[str,PauliString]) -> None:
    """Test commutation relations for single-qubit Pauli operators."""
    # Identity commutes with everything
    assert pauli_setup["I"] | pauli_setup["I"]
    assert pauli_setup["I"] | pauli_setup["X"]
    assert pauli_setup["I"] | pauli_setup["Y"]
    assert pauli_setup["I"] | pauli_setup["Z"]
    # X commutation relations
    assert pauli_setup["X"] | pauli_setup["I"]
    assert pauli_setup["X"] | pauli_setup["X"]
    assert not pauli_setup["X"] | pauli_setup["Y"]
    assert not pauli_setup["X"] | pauli_setup["Z"]
    # Y commutation relations
    assert pauli_setup["Y"] | pauli_setup["I"]
    assert not pauli_setup["Y"] | pauli_setup["X"]
    assert pauli_setup["Y"] | pauli_setup["Y"]
    assert not pauli_setup["Y"] | pauli_setup["Z"]
    # Z commutation relations
    assert pauli_setup["Z"] | pauli_setup["I"]
    assert not pauli_setup["Z"] | pauli_setup["X"]
    assert not pauli_setup["Z"] | pauli_setup["Y"]
    assert pauli_setup["Z"] | pauli_setup["Z"]

def test_single_qubit_products(pauli_setup:dict[str,PauliString]) -> None:
    """Test multiplication of single-qubit Pauli operators."""
    # Extract operators for cleaner code
    i, x, y, z = pauli_setup["I"], pauli_setup["X"], pauli_setup["Y"], pauli_setup["Z"]
    # Identity products
    assert i @ i == i
    assert i @ x == x
    assert i @ y == y
    assert i @ z == z
    # X products
    assert x @ i == x
    assert x @ x == i
    assert x @ y == z
    assert x @ z == y
    # Y products
    assert y @ i == y
    assert y @ x == z
    assert y @ y == i
    assert y @ z == x
    # Z products
    assert z @ i == z
    assert z @ x == y
    assert z @ y == x
    assert z @ z == i

def test_operator_chaining(pauli_setup:dict[str,PauliString]) -> None:
    """Test chaining of Pauli operators."""
    x, y, z = pauli_setup["X"], pauli_setup["Y"], pauli_setup["Z"]
    assert z ^ y == x
    assert x @ (y @ z) == pauli_setup["I"]
    assert x @ y @ z == pauli_setup["I"]

def test_multi_qubit_commutation(pauli_setup:dict[str,PauliString]) -> None:
    """Test commutation relations for multi-qubit Pauli strings."""
    # Same operator commutes with itself
    assert pauli_setup["XI"] | pauli_setup["XI"]
    # Different operators on different qubits
    assert pauli_setup["XY"] | pauli_setup["YX"]
    # Non-commuting operators on same qubit
    assert not pauli_setup["XIIII"] | pauli_setup["ZIIII"]
    # Test commutation with the first pattern
    assert pauli_setup["XIIII"] | pauli_setup["IYXII"]
    assert pauli_setup["XIIII"] | pauli_setup["IXIXI"]
    assert pauli_setup["XIIII"] | pauli_setup["IIIZI"]
    assert pauli_setup["XIIII"] | pauli_setup["IZIXZ"]
    assert pauli_setup["XIIII"] | pauli_setup["IIIIX"]

def test_complex_commutation_matrix(pauli_setup:dict[str,PauliString]) -> None:
    """Test commutation relations between multiple multi-qubit operators."""
    # Define the operators to test
    operators = [
        pauli_setup["ZIIII"], pauli_setup["IYXII"], pauli_setup["IXIXI"],
        pauli_setup["IIIZI"], pauli_setup["IZIXZ"], pauli_setup["IIIIX"]
    ]
    # Expected commutation results (True if commutes, False otherwise)
    expected_results = [
        # ZIIII  IYXII  IXIXI  IIIZI  IZIXZ  IIIIX
        [True,  True,  True,  True,  True,  True],   # ZIIII
        [True,  True,  False, True,  False, True],   # IYXII
        [True,  False, True,  False, False, True],   # IXIXI
        [True,  True,  False, True,  False, True],   # IIIZI
        [True,  False, False, False, True,  False],  # IZIXZ
        [True,  True,  True,  True,  False, True]    # IIIIX
    ]
    # Test all combinations
    for i, op1 in enumerate(operators):
        for j, op2 in enumerate(operators):
            assert (op1 | op2) == expected_results[i][j], (f"Commutation failed "
                   f"for operators at position {i},{j}")

def test_multi_qubit_products(pauli_setup:dict[str,PauliString]) -> None:
    """Test products of multi-qubit Pauli strings."""
    assert pauli_setup["IXIXI"] ^ pauli_setup["IYXII"] == pauli_setup["IZXXI"]
