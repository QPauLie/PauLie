

from paulie.common.pauli import (
    commutator,
    get_pauli_array,
    is_commutative,
    multiply_pauli_arrays,
)


def test_commmutator():
    I = get_pauli_array("I")   # noqa
    X = get_pauli_array("X")
    Y = get_pauli_array("Y")
    Z = get_pauli_array("Z")
    assert is_commutative(I, I)
    assert is_commutative(I, X)
    assert is_commutative(I, Y)
    assert is_commutative(I, Z)
    assert is_commutative(X, I)
    assert is_commutative(X, X)
    assert is_commutative(X, Y) is False
    assert is_commutative(X, Z) is False
    assert is_commutative(Y, I)
    assert is_commutative(Y, X) is False
    assert is_commutative(Y, Y)
    assert is_commutative(Y, Z) is False
    assert is_commutative(Z, I)
    assert is_commutative(Z, X) is False
    assert is_commutative(Z, Y) is False
    assert is_commutative(Z, Z)
    assert multiply_pauli_arrays(I, I) == I
    assert multiply_pauli_arrays(I, X) == X
    assert multiply_pauli_arrays(I, Y) == Y
    assert multiply_pauli_arrays(I, Z) == Z
    assert commutator(X, Y) == Z
    assert commutator(X, Z) == Y
    assert commutator(Y, X) == Z
    assert commutator(Y, Z) == X
    assert commutator(Z, X) == Y
    assert commutator(Z, Y) == X
    assert(commutator(X, commutator(Y, Z))) == I

    XI = get_pauli_array("XI")
    get_pauli_array("IX")
    assert is_commutative(XI, XI)

    XY = get_pauli_array("XY")
    YX = get_pauli_array("YX")
    assert is_commutative(XY, YX)
    XIIII = get_pauli_array("XIIII")
    ZIIII = get_pauli_array("ZIIII")
    assert is_commutative(XIIII, ZIIII) is False

    IYXII = get_pauli_array("IYXII")
    IXIXI = get_pauli_array("IXIXI")
    IIIZI = get_pauli_array("IIIZI")
    IZIXZ = get_pauli_array("IZIXZ")
    IIIIX = get_pauli_array("IIIIX")

    assert is_commutative(XIIII, IYXII)
    assert is_commutative(XIIII, IXIXI)
    assert is_commutative(XIIII, IIIZI)
    assert is_commutative(XIIII, IZIXZ)
    assert is_commutative(XIIII, IIIIX)

    assert is_commutative(ZIIII, IYXII)
    assert is_commutative(ZIIII, IXIXI)
    assert is_commutative(ZIIII, IIIZI)
    assert is_commutative(ZIIII, IZIXZ)
    assert is_commutative(ZIIII, IIIIX)

    assert is_commutative(IYXII, IYXII)
    assert is_commutative(IYXII, IXIXI) is False
    assert is_commutative(IYXII, IIIZI)
    assert is_commutative(IYXII, IZIXZ) is False
    assert is_commutative(IYXII, IIIIX)

    assert is_commutative(IXIXI, IYXII) is False
    assert is_commutative(IXIXI, IXIXI)
    assert is_commutative(IXIXI, IIIZI) is False
    assert is_commutative(IXIXI, IZIXZ) is False
    assert is_commutative(IXIXI, IIIIX)

    assert is_commutative(IIIZI, IYXII)
    assert is_commutative(IIIZI, IXIXI) is False
    assert is_commutative(IIIZI, IIIZI)
    assert is_commutative(IIIZI, IZIXZ) is False
    assert is_commutative(IIIZI, IIIIX)

    assert is_commutative(IZIXZ, IYXII) is False
    assert is_commutative(IZIXZ, IXIXI) is False
    assert is_commutative(IZIXZ, IIIZI) is False
    assert is_commutative(IZIXZ, IZIXZ)
    assert is_commutative(IZIXZ, IIIIX) is False

    assert is_commutative(IIIIX, IYXII)
    assert is_commutative(IIIIX, IXIXI)
    assert is_commutative(IIIIX, IIIZI)
    assert is_commutative(IIIIX, IZIXZ) is False
    assert is_commutative(IIIIX, IIIIX)

    IZXXI = get_pauli_array("IZXXI")
    assert commutator(IXIXI, IYXII) == IZXXI
