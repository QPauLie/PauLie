import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p

@pytest.fixture(scope="module")
def pauli_setup():
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

def test_single_qubit_commutation(pauli_setup):
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

def test_single_qubit_products(pauli_setup):
    """Test multiplication of single-qubit Pauli operators."""
    # Extract operators for cleaner code
    I, X, Y, Z = pauli_setup["I"], pauli_setup["X"], pauli_setup["Y"], pauli_setup["Z"]
    
    # Identity products
    assert I @ I == I
    assert I @ X == X
    assert I @ Y == Y
    assert I @ Z == Z
    
    # X products
    assert X @ I == X
    assert X @ X == I
    assert X @ Y == Z
    assert X @ Z == Y
    
    # Y products
    assert Y @ I == Y
    assert Y @ X == Z
    assert Y @ Y == I
    assert Y @ Z == X
    
    # Z products
    assert Z @ I == Z
    assert Z @ X == Y
    assert Z @ Y == X
    assert Z @ Z == I

def test_operator_chaining(pauli_setup):
    """Test chaining of Pauli operators."""
    X, Y, Z = pauli_setup["X"], pauli_setup["Y"], pauli_setup["Z"]
    
    assert Z ^ Y == X
    assert X @ (Y @ Z) == pauli_setup["I"]
    assert X @ Y @ Z == pauli_setup["I"]

def test_multi_qubit_commutation(pauli_setup):
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

def test_complex_commutation_matrix(pauli_setup):
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
            assert (op1 | op2) == expected_results[i][j], f"Commutation failed for operators at position {i},{j}"

def test_multi_qubit_products(pauli_setup):
    """Test products of multi-qubit Pauli strings."""
    assert pauli_setup["IXIXI"] ^ pauli_setup["IYXII"] == pauli_setup["IZXXI"]
