"""
Test nested
"""
import pytest
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_collection import PauliStringCollection

# Define expected nested Pauli string pairs for each node
NESTED_PAIRS = {
    "IX": [["ZY", "ZZ"], ["XY", "XZ"], ["YZ", "YY"], ["IZ", "IY"]],
    "IY": [["IX", "IZ"], ["XX", "XZ"], ["YX", "YZ"], ["ZX", "ZZ"]],
    "IZ": [["IX", "IY"], ["XX", "XY"], ["ZX", "ZY"], ["YX", "YY"]],
    "XI": [["ZY", "YY"], ["ZX", "YX"], ["ZZ", "YZ"], ["ZI", "YI"]],
    "XX": [["XZ", "IY"], ["IZ", "XY"], ["ZI", "YX"], ["YI", "ZX"]],
    "XY": [["IX", "XZ"], ["IZ", "XX"], ["YI", "ZY"], ["YY", "ZI"]],
    "XZ": [["IX", "XY"], ["IY", "XX"], ["YI", "ZZ"], ["YZ", "ZI"]],
    "YI": [["XX", "ZX"], ["XY", "ZY"], ["XZ", "ZZ"], ["ZI", "XI"]],
    "YX": [["IY", "YZ"], ["IZ", "YY"], ["ZI", "XX"], ["ZX", "XI"]],
    "YY": [["YX", "IZ"], ["ZI", "XY"], ["ZY", "XI"], ["IX", "YZ"]],
    "YZ": [["IY", "YX"], ["ZZ", "XI"], ["ZI", "XZ"], ["YY", "IX"]],
    "ZX": [["ZY", "IZ"], ["YI", "XX"], ["IY", "ZZ"], ["YX", "XI"]],
    "ZY": [["ZX", "IZ"], ["YI", "XY"], ["XI", "YY"], ["IX", "ZZ"]],
    "ZZ": [["XZ", "YI"], ["IY", "ZX"], ["ZY", "IX"], ["XI", "YZ"]]
}

def check_nested(node:PauliStringCollection) -> bool:
    """
    Check if a node's nested pairs match the expected pairs.
    
    Args:
        node: A Pauli string node to test
        
    Returns:
        bool: True if the nested pairs match expectations, False otherwise
    """
    # Get the actual nested pairs from the node
    nested = node.get_nested()
    # Get the expected nested pairs
    nested_source = NESTED_PAIRS[str(node)]
    # Check that we have the expected number of pairs
    if len(nested) != len(nested_source):
        return False
    # Check that all actual pairs are in the expected set
    for pair in nested:
        found = False
        for source_pair in nested_source:
            if str(pair[0]) in source_pair and str(pair[1]) in source_pair:
                found = True
                break
        if not found:
            return False
    # Check that all expected pairs are in the actual set
    for source_pair in nested_source:
        found = False
        for pair in nested:
            if str(pair[0]) in source_pair and str(pair[1]) in source_pair:
                found = True
                break
        if not found:
            return False
    return True

@pytest.mark.parametrize("pauli_string", list(NESTED_PAIRS.keys()))
def test_nested(pauli_string:str):
    """
    Test that all Pauli strings have the correct nested pairs.
    
    Args:
        pauli_string: The string representation of a Pauli operator
    """
    pauli_operator = p(pauli_string)
    assert check_nested(pauli_operator), (f"Nested pairs for {pauli_string} "
           "do not match expected values")
