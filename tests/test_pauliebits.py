"""
Test pauliebits methods
"""
from pauliebits import pauliebits
from pauliebits.util import count_and, count_or, ba2int
from paulie.common.random_pauli_strings import get_random_list


CODEC = {
    "I": pauliebits([0, 0]),
    "X": pauliebits([1, 0]),
    "Y": pauliebits([1, 1]),
    "Z": pauliebits([0, 1]),
}


def test_pauliebits_new_methods() -> None:
    """Test pauliebits new methods of pauliebits"""
    paulies = get_random_list(1000, 4000)
    for paulie in paulies:
        temp_bits = pauliebits()
        temp_bits.encode_ixyz(paulie)
        assert paulie == temp_bits.decode_ixyz()

        temp_bits2 = pauliebits()
        temp_bits2.encode(CODEC, paulie)
        assert temp_bits == temp_bits2

        bits_even = temp_bits[::2]
        bits_odd = temp_bits[1::2]
        assert temp_bits.count_non_trivially() == count_or(bits_even, bits_odd)
        diagonal_index = ba2int(bits_odd) if ba2int(bits_even) == 0 else -1
        assert temp_bits.diagonal_index() == diagonal_index
        ys = count_and(bits_odd, bits_even)
        assert temp_bits.complex_conjugate() == ys
        assert temp_bits.not_identity_mask() == bits_even | bits_odd
        others = get_random_list(1000, 2)
        for other in others:
            temp_other_bits = pauliebits()
            temp_other_bits.encode_ixyz(other)
            other_bits_even = temp_other_bits[::2]
            other_bits_odd = temp_other_bits[1::2]
            f = 2 * count_and(bits_even, other_bits_odd) + \
                count_and(bits_odd, bits_even) + \
                count_and(other_bits_odd, other_bits_even) - \
                count_and(bits_even ^ other_bits_even,
                bits_odd ^ other_bits_odd)
            assert temp_bits.phase(temp_other_bits) == f
            old_commute_with = ((count_and(bits_even, other_bits_odd)
                + count_and(bits_odd, other_bits_even)) %2 == 0)
            assert temp_bits.commutes_with(temp_other_bits) == old_commute_with
