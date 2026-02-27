"""Test subgraph"""
from paulie import get_pauli_string as p

def test_subgraphs() -> None:
    """Test subgraph"""
    generators = p(["XIIII", "ZIIII", "IYXII", "IXIXI", "IIIZI", "IZIXZ", "IIIIX"])
    subgraphs = generators.get_subgraphs()
    assert len(subgraphs) == 2
