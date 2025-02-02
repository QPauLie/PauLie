from bitarray import bitarray

# Object for handling Pauli strings that relies on the binary symplectic form
# See Section 2 of https://quantum-journal.org/papers/q-2020-06-04-278/
# The binary symplectic form works as follows. For N = 1 we have
#I = (00)
#X = (10)
#Y = (11)
#Z = (01)
# This extends obviously for N>1, for example XYZ = (1,1,0|0,1,1)
# By performing modular arithmetic on this array we can implement the Pauli algebra.
codec = {
    "I": bitarray([0, 0]), 
    "X": bitarray([1, 0]), 
    "Y": bitarray([1, 1]), 
    "Z": bitarray([0, 1]),
}


def get_I(n):
    if n == 0:
        return bitarray()
    return bitarray(2*n)


def get_Y(n):
    if n == 0:
        return bitarray()
    onePauliArray = get_I(n)
    onePauliArray.setall(1)
    return onePauliArray


def get_Z(n):
    if n == 0:
        return bitarray()
    zetPauliArray = get_I(n)
    i = 1
    while(i < 2*n):
        zetPauliArray[i] = 1
        i += 2
    return zetPauliArray


def get_pauli_string(bitstring):
    return "".join(bitstring.decode(codec))


def get_pauli_array(pauli_string):
    pauli_array = bitarray()
    pauli_array.encode(codec, pauli_string)
    return pauli_array


def inc_pauli_array(pauli_array):
    n = len(pauli_array) - 1
    stop = False
    while stop is not True:
        if pauli_array[n] == 0:
           pauli_array[n] = 1
           break
        if pauli_array[n] == 1:
           pauli_array[n] = 0
           n = n - 1
    return pauli_array

   
def inc_IZ_pauli_array(pauli_array):
    n = len(pauli_array) - 1
    stop = False
    while stop is not True:
        if pauli_array[n] == 0:
           pauli_array[n] = 1
           break
        if pauli_array[n] == 1:
           pauli_array[n] = 0
           n = n - 2
    return pauli_array


def is_commutate(a, b):
    if len(a) != len(b):
        raise ValueError("gates must have the same length")
    a_dot_b = 0
    b_dot_a = 0
    i = 0

    while i < len(a):
       a_dot_b += 1 if a[i] and  b[i + 1] else 0
       b_dot_a += 1 if a[i + 1] and b[i] else 0
       i = i + 2

    a_dot_b %= 2
    b_dot_a %= 2
    return a_dot_b == b_dot_a


def is_commutate_by_string(a, b):
     return is_commutate(get_pauli_array(a), get_pauli_array(b))


def multi_pauli_arrays(a, b):
    if len(a) != len(b):
        raise ValueError("gates must have the same length")
    c = bitarray(len(a))
    i = 0
    while i < len(c):
       c[i] = (a[i] + b[i]) % 2
       i = i + 1
    return c


def multi_pauli_string_to_array(a, b):
    aArray = get_pauli_array(a)
    bArray = get_pauli_array(b)
    cArray = multi_pauli_arrays(aArray, bArray)
    return cArray


def multi_pauli_string(a, b):
    return get_pauli_string(multi_pauli_string_to_array(a, b))


def commutator(a, b):
    if is_commutate(a, b):
        return bitarray(len(a))
    return multi_pauli_arrays(a, b)


def commutator_pauli_string(a, b):
    aArray = get_pauli_array(a)
    bArray = get_pauli_array(b)
    cArray = commutator(aArray, bArray)
    return get_pauli_string(cArray)


def is_IZ_string(a):
    i = 0
    size = len(a)
    while(i < size):
        if a[i] != 0:
            return False
        i += 2
    return True


def is_sub_in_array(sub, a, pos=0):
    index = a.find(sub, pos)
    if index == -1:
        return False
    if index % 2 == 0:
        return True
    return is_sub_in_array(sub, a, index+1)


def gen_all_nodes(n):
    a = get_I(n)
    yield a
    last = get_Y(n)
    while True:
        a = inc_pauli_array(a)
        yield a
        if a == last:
            break


def gen_all_IZ(n):
    a = get_I(n)
    yield a
    last = get_Z(n)
    while True:
        a = inc_IZ_pauli_array(a)
        yield a
        if a == last:
            break


def get_array_pauli_strings(bit_arrays):
    return list(map(get_pauli_string, bit_arrays))


def get_array_pauli_arrays(pauli_strings):
    return list(map(get_pauli_array, pauli_strings))

