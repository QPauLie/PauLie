

from paulie.common.pauli import (
    commutator,
    get_pauli_array,
    is_commutate,
    multi_pauli_arrays,
)


def test_commmutator():
    I = get_pauli_array("I")   # noqa
    X = get_pauli_array("X")
    Y = get_pauli_array("Y")
    Z = get_pauli_array("Z")
    assert is_commutate(I, I)
    assert is_commutate(I, X)
    assert is_commutate(I, Y)
    assert is_commutate(I, Z)
    assert is_commutate(X, I)
    assert is_commutate(X, X)
    assert is_commutate(X, Y) is False
    assert is_commutate(X, Z) is False
    assert is_commutate(Y, I)
    assert is_commutate(Y, X) is False
    assert is_commutate(Y, Y)
    assert is_commutate(Y, Z) is False
    assert is_commutate(Z, I)
    assert is_commutate(Z, X) is False
    assert is_commutate(Z, Y) is False
    assert is_commutate(Z, Z)
    assert multi_pauli_arrays(I, I) == I
    assert multi_pauli_arrays(I, X) == X
    assert multi_pauli_arrays(I, Y) == Y
    assert multi_pauli_arrays(I, Z) == Z
    assert commutator(X, Y) == Z
    assert commutator(X, Z) == Y
    assert commutator(Y, X) == Z
    assert commutator(Y, Z) == X
    assert commutator(Z, X) == Y
    assert commutator(Z, Y) == X
    assert(commutator(X, commutator(Y, Z))) == I

    XI = get_pauli_array("XI")
    get_pauli_array("IX")
    assert is_commutate(XI, XI)

    XY = get_pauli_array("XY")
    YX = get_pauli_array("YX")
    assert is_commutate(XY, YX)
    XIIII = get_pauli_array("XIIII")
    ZIIII = get_pauli_array("ZIIII")
    assert is_commutate(XIIII, ZIIII) is False

    IYXII = get_pauli_array("IYXII")
    IXIXI = get_pauli_array("IXIXI")
    IIIZI = get_pauli_array("IIIZI")
    IZIXZ = get_pauli_array("IZIXZ")
    IIIIX = get_pauli_array("IIIIX")

    assert is_commutate(XIIII, IYXII)
    assert is_commutate(XIIII, IXIXI)
    assert is_commutate(XIIII, IIIZI)
    assert is_commutate(XIIII, IZIXZ)
    assert is_commutate(XIIII, IIIIX)

    assert is_commutate(ZIIII, IYXII)
    assert is_commutate(ZIIII, IXIXI)
    assert is_commutate(ZIIII, IIIZI)
    assert is_commutate(ZIIII, IZIXZ)
    assert is_commutate(ZIIII, IIIIX)

    assert is_commutate(IYXII, IYXII)
    assert is_commutate(IYXII, IXIXI) is False
    assert is_commutate(IYXII, IIIZI)
    assert is_commutate(IYXII, IZIXZ) is False
    assert is_commutate(IYXII, IIIIX)

    assert is_commutate(IXIXI, IYXII) is False
    assert is_commutate(IXIXI, IXIXI)
    assert is_commutate(IXIXI, IIIZI) is False
    assert is_commutate(IXIXI, IZIXZ) is False
    assert is_commutate(IXIXI, IIIIX)

    assert is_commutate(IIIZI, IYXII)
    assert is_commutate(IIIZI, IXIXI) is False
    assert is_commutate(IIIZI, IIIZI)
    assert is_commutate(IIIZI, IZIXZ) is False
    assert is_commutate(IIIZI, IIIIX)

    assert is_commutate(IZIXZ, IYXII) is False
    assert is_commutate(IZIXZ, IXIXI) is False
    assert is_commutate(IZIXZ, IIIZI) is False
    assert is_commutate(IZIXZ, IZIXZ)
    assert is_commutate(IZIXZ, IIIIX) is False

    assert is_commutate(IIIIX, IYXII)
    assert is_commutate(IIIIX, IXIXI)
    assert is_commutate(IIIIX, IIIZI)
    assert is_commutate(IIIIX, IZIXZ) is False
    assert is_commutate(IIIIX, IIIIX)

    IZXXI = get_pauli_array("IZXXI")
    assert commutator(IXIXI, IYXII) == IZXXI
