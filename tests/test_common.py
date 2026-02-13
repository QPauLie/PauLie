"""
Test common
"""
import itertools
import pytest
from bitarray import bitarray
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection

test_cases = [
    ({"n": 0}, "", "PauliString()"),
    ({"pauli_str": ""}, "", "PauliString()"),
    ({"n": 4}, "IIII", "PauliString(IIII)"),
    ({"pauli_str": "IXYZ"}, "IXYZ", "PauliString(IXYZ)"),
    ({"bits": bitarray("00101101")}, "IXYZ", "PauliString(IXYZ)"),
    ({"n": 5, "pauli_str": "IXYZ"}, "IXYZI", "PauliString(IXYZI)"),
]

test_case_collection = [
    (input_arg1, expected_str1, _expected_repr1, input_arg2, expected_str2, _expected_repr2)
    for (input_arg1, expected_str1, _expected_repr1), (input_arg2, expected_str2, _expected_repr2)
    in itertools.combinations(test_cases, 2)
]

test_coeff = [-1j, -1, 0, 1, 1j]
test_case_linear = [
    (
        coeff1, input_arg1, expected_str1, expected_repr1,
        coeff2, input_arg2, expected_str2, expected_repr2,
    )
    for (coeff1, coeff2) in itertools.product(test_coeff, repeat=2)
    for (
        (input_arg1, expected_str1, expected_repr1),
        (input_arg2, expected_str2, expected_repr2),
    ) in itertools.combinations(test_cases, 2)
]

@pytest.mark.parametrize("input_arg, expected_str, expected_repr", test_cases)
def test_PauliString_str_repr(input_arg, expected_str, expected_repr) -> None:
    """
    Test PauliString for consistent __repr__ and __str__
    """

    pauli_string = PauliString(**input_arg)

    assert str(pauli_string) == expected_str
    assert repr(pauli_string) == expected_repr

@pytest.mark.parametrize(
    "input_arg1, expected_str1, expected_repr1,"
    "input_arg2, expected_str2, expected_repr2",
    test_case_collection
)
def test_PauliStringCollection_str_repr(
    input_arg1, expected_str1, expected_repr1,
    input_arg2, expected_str2, expected_repr2
    ) -> None:
    """
    Test PauliStringCollection for consistent __repr__ and __str__
    """

    pauli_string1 = PauliString(**input_arg1)
    pauli_string2 = PauliString(**input_arg2)

    pauli_string_collection = PauliStringCollection([pauli_string1, pauli_string2])

    max_len = max(len(pauli_string1),len(pauli_string2))

    expected_str1 = expected_str1.ljust(max_len, "I")
    expected_str2 = expected_str2.ljust(max_len, "I")

    expected_str = f"[{expected_str1}, {expected_str2}]"
    expected_repr = f"PauliStringCollection({expected_str})"

    assert str(pauli_string_collection) == expected_str
    assert repr(pauli_string_collection) == expected_repr
