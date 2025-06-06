"""
    Test second-order twirl routine.
"""
import pytest
import numpy as np
from paulie.common.pauli_string_factory import get_pauli_string as p
from paulie.common.pauli_string_bitarray import PauliString
from paulie.application.second_moment import second_moment

def create_swap_matrix(n: int) -> np.ndarray:
    """
    Construct SWAP matrix on n qubits.

    The SWAP matrix is given by sum_{i} sum_{j} |i, j> <j, i|.
    """
    swap = np.zeros((2 ** n, 2 ** n))
    for i in range(2 ** n):
        j = int(np.binary_repr(i, width=n)[::-1], base=2)
        swap[i, j] = 1
    return swap

@pytest.mark.parametrize("gens", [
    ["I"], ["X"], ["Y"], ["Z"],
    ["I", "X"], ["X", "Y"], ["Y", "Z"],
    ["X", "Y", "Z"], ["I", "X", "Y", "Z"]
])
def test_frame_potential_from_second_order_twirl(gens: list[str]) -> None:
    """
    Reference: arxiv:2502.16404

    There is a method get_frame_potential to get the frame
    potential for any generating set. We can also compute
    this quantity using second-order twirls (Appendix A).
    """
    n = len(gens[0])
    gens = p(gens)
    # Frame potential from commutator graph
    fp_from_graph = gens.get_frame_potential()
    # Frame potential from second-order twirls
    fp_from_twirl = 0
    i = PauliString(n=n)
    swap = create_swap_matrix(2 * n)
    all_pauli_strings = list(i.gen_all_pauli_strings())
    for a1 in all_pauli_strings:
        for a2 in all_pauli_strings:
            a1krona2 = a1 + a2
            for b1 in all_pauli_strings:
                for b2 in all_pauli_strings:
                    b1kronb2 = p([(1, b1 + b2)])
                    twirlb1b2 = second_moment(gens, b1kronb2)
                    qq = a1krona2 @ twirlb1b2
                    qqmat = qq.get_matrix()
                    tr = np.trace(qqmat @ swap)
                    fp_from_twirl += np.abs(tr) ** 2
    fp_from_twirl /= 2 ** (4 * n)
    assert fp_from_twirl == pytest.approx(fp_from_graph)
