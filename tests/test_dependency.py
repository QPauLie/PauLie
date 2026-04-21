"""
Tests for linear dependency checks during canonicalization.
"""
from paulie import get_pauli_string as p
from paulie.classifier.canonicalizer import Canonicalizer
from paulie.common.pauli_string_bitarray import PauliString

def test_central_dependency_dropping():
    """
    Test that a generator dependent on the central vertex and other legs is dropped.

    Bug: 'Incomplete Dependency Check'
    Previously, _dependency_check only performed Gaussian elimination on length-1 legs,
    ignoring the central_vertex. This allowed generators dependent on {central, leg1, ...}
    to be incorrectly retained.
    """
    # Central = XI
    # L1 = ZX
    # L2 = IZ
    # L3 = YY (Product of XI, ZX, IZ)
    # symp(XI, ZX) = 1
    # symp(XI, IZ) = 1
    # symp(XI, YY) = 1
    # All are connected to Central=XI.
    # L3 = Central @ L1 @ L2 (ignoring phases).
    gens = ["XI", "ZX", "IZ", "YY"]

    cc = Canonicalizer()
    v_stack = [PauliString(pauli_str=g) for g in gens]

    # Stack is LIFO. To ensure XI is Central, put it last so it's popped first.
    # v_stack[0] is XI
    ordered_stack = [v_stack[1], v_stack[2], v_stack[3], v_stack[0]]

    cc.build_canonical_graph(ordered_stack)

    # After fix, YY should be dropped because it is dependent on {XI, ZX, IZ}.
    # We expect 2 legs: [ZX] and [IZ].
    assert len(cc.legs) == 2
    assert cc.central_vertex == v_stack[0]

    # Total vertices in canonical graph should be 3: XI, ZX, IZ.
    total_verts = 1 + sum(len(leg) for leg in cc.legs)
    assert total_verts == 3

def test_dla_dimension_with_redundant_generators():
    """
    Verify that redundant generators do not affect the calculated DLA dimension.
    """
    # {XI, ZX, IZ} generates so(4) ~ su(2) + su(2), dimension 6.
    # Adding YY (which is XI*ZX*IZ) should not change the DLA.
    gens_independent = ["XI", "ZX", "IZ"]
    gens_redundant = ["XI", "ZX", "IZ", "YY"]

    coll_indep = p(gens_independent)
    coll_redun = p(gens_redundant)

    dim_indep = coll_indep.get_dla_dim()
    dim_redun = coll_redun.get_dla_dim()

    assert dim_indep == 6
    assert dim_redun == 6
    assert coll_redun.get_algebra() == coll_indep.get_algebra()
