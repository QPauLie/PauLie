"""Test subgraph"""
import pytest

from paulie import get_pauli_string as p


def test_subgraphs() -> None:
    """Test subgraph"""
    generators = p(["XIIII", "ZIIII", "IYXII", "IXIXI", "IIIZI", "IZIXZ", "IIIIX"])
    subgraphs = generators.get_subgraphs()
    assert len(subgraphs) == 2


@pytest.mark.parametrize("generators", [
    ["XXII", "IXXI", "ZIII", "IZII"],
    ["XXII", "IIZZ"],
    ["XXII", "IXXI", "IIXX"],
    ["XY", "YX"],
    ["XXI", "IXX", "ZII", "IZI", "IIZ"],
])
def test_graph_components_agree_with_subgraphs(generators: list[str]) -> None:
    """The two ways of splitting on the anticommutation graph must agree.

    get_graph_components used to pass the collection to get_graph, which keeps only
    those edges whose commutator is itself a generator. That is a strictly smaller
    graph: for a transverse-field Ising chain [XX, Z] ~ YX is not a generator, so no
    edges survived and every term looked like its own component.
    """
    collection = p(generators)
    components = {frozenset(str(s) for s in c)
                  for c in collection.get_graph_components("anticommutator")}
    subgraphs = {frozenset(str(s) for s in c) for c in collection.get_subgraphs()}
    assert components == subgraphs


def test_anticommuting_terms_land_in_one_component() -> None:
    """A transverse-field Ising chain is connected, so its evolution does not factorize."""
    components = p(["XXII", "IXXI", "ZIII", "IZII"]).get_graph_components("anticommutator")
    assert len(components) == 1


def test_disjoint_support_gives_separate_components() -> None:
    """Terms that commute with everything else split off, exactly."""
    components = p(["XXII", "IIZZ"]).get_graph_components("anticommutator")
    assert len(components) == 2
