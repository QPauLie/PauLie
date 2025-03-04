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
CODEC = {
    "I": bitarray([0, 0]), 
    "X": bitarray([1, 0]), 
    "Y": bitarray([1, 1]), 
    "Z": bitarray([0, 1]),
}


def get_I(n: int) -> bitarray:
    """Returns the identity Pauli bitarray of size `2 * n`."""
    return bitarray(2 * n)


def get_Y(n: int) -> bitarray:
    """Returns a Pauli bitarray filled with `1`."""
    return bitarray([1] * (2 * n))


def get_Z(n: int) -> bitarray:
    """Returns a Pauli bitarray where every second bit is set to 1."""
    z_array = get_I(n)
    z_array[1::2] = 1  # Set every odd index to 1
    return z_array


def get_pauli_string(bitstring: bitarray) -> str:
    """Decodes a Pauli bitarray into its string representation."""
    return "".join(bitstring.decode(CODEC))


def get_pauli_array(pauli_string: str) -> bitarray:
    """Encodes a Pauli string into a bitarray representation."""
    pauli_array = bitarray()
    pauli_array.encode(CODEC, pauli_string)
    return pauli_array


def inc_pauli_array(pauli_array: bitarray) -> bitarray:
    """Increments the binary representation of a Pauli array."""
    for i in reversed(range(len(pauli_array))):
        if pauli_array[i] == 0:
            pauli_array[i] = 1
            break
        pauli_array[i] = 0
    return pauli_array

   
def inc_IZ_pauli_array(pauli_array: bitarray) -> bitarray:
    """Increments a binary Pauli array, skipping every second bit."""
    for i in reversed(range(0, len(pauli_array), 2)):
        if pauli_array[i] == 0:
            pauli_array[i] = 1
            break
        pauli_array[i] = 0
    return pauli_array


def is_commutative(a: bitarray, b: bitarray) -> bool:
    """Checks if two Pauli operators commute using the symplectic inner product."""
    if len(a) != len(b):
        raise ValueError("Pauli arrays must be of equal length")
    return sum(a[i] & b[i + 1] for i in range(0, len(a), 2)) % 2 == sum(a[i + 1] & b[i] for i in range(0, len(a), 2)) % 2


def is_commutative_by_string(a: str, b: str) -> bool:
    """Checks commutativity by converting Pauli strings to bitarrays first."""
    return is_commutative(get_pauli_array(a), get_pauli_array(b))


def multiply_pauli_arrays(a: bitarray, b: bitarray) -> bitarray:
    """Performs element-wise modulo-2 addition (XOR) for two Pauli bitarrays."""
    if len(a) != len(b):
        raise ValueError("Pauli arrays must have the same length")
    return a ^ b  # Bitwise XOR is equivalent to mod-2 addition


def multiply_pauli_strings(a: str, b: str) -> str:
    """Multiplies two Pauli strings and returns the result as a string."""
    return get_pauli_string(multiply_pauli_arrays(get_pauli_array(a), get_pauli_array(b)))


def commutator(a: bitarray, b: bitarray) -> bitarray:
    """Computes the commutator [A, B] = AB - BA."""
    return bitarray(len(a)) if is_commutative(a, b) else multiply_pauli_arrays(a, b)


def commutator_pauli_string(a: str, b: str) -> str:
    """Computes the commutator of two Pauli strings."""
    return get_pauli_string(commutator(get_pauli_array(a), get_pauli_array(b)))


def commutant(generators: list[str]) -> list[str]:
    """Returns the commutant of a set of Pauli generators."""
    size = len(generators[0])
    generators_array = get_array_pauli_arrays(generators)
    commutant_list = [
        node.copy() for node in gen_all_nodes(size) if all(is_commutative(node, gen) for gen in generators_array)
    ]
    return get_array_pauli_strings(commutant_list)


def is_IZ_string(a: bitarray) -> bool:
    """Checks if a Pauli bitarray consists only of 'I' and 'Z'."""
    return all(a[i] == 0 for i in range(0, len(a), 2))


def is_sub_in_array(sub: bitarray, a: bitarray, pos: int = 0) -> bool:
    """Checks if `sub` exists in `a` starting from `pos`, ensuring it aligns with even indices."""
    while (index := a.find(sub, pos)) != -1:
        if index % 2 == 0:
            return True
        pos = index + 1
    return False


def gen_all_nodes(n: int):
    """Generator yielding all possible Pauli bitarrays of size `2*n`."""
    a = get_I(n)
    yield a.copy()
    last = get_Y(n)
    while True:
        a = inc_pauli_array(a)
        yield a.copy()
        if a == last:
            break


def gen_all_IZ(n: int):
    """Generator yielding all Pauli bitarrays consisting of 'I' and 'Z'."""
    node = get_I(n)
    last = get_Z(n)
    while True:
        yield node
        node = inc_IZ_pauli_array(node)
        if node == last:
            break


def get_array_pauli_strings(bit_arrays: list[bitarray]) -> list[str]:
    """Converts a list of Pauli bitarrays to a list of Pauli strings."""
    return list(map(get_pauli_string, bit_arrays))


def get_array_pauli_arrays(pauli_strings: list[str]) -> list[bitarray]:
    """Converts a list of Pauli strings to a list of Pauli bitarrays."""
    return list(map(get_pauli_array, pauli_strings))


def complex_conj(a: list[bitarray]) -> list[bitarray]:
    """Computes the complex conjugate of a Pauli array."""
    pauli_strings = get_array_pauli_strings(a)
    return pauli_strings.count("Y") * get_array_pauli_arrays(pauli_strings)