"""
Test k-local
"""
from paulie.common.pauli_string_factory import get_pauli_string as p

def test_k_local() -> None:
    """
    Test k-local
    """
    generators = p(["XY", "XZ"], n=5)
    assert len(generators) == 8
    assert p("XYIII") in generators
    assert p("IXYII") in generators
    assert p("IIXYI") in generators
    assert p("IIIXY") in generators
    assert p("XZIII") in generators
    assert p("IXZII") in generators
    assert p("IIXZI") in generators
    assert p("IIIXZ") in generators
