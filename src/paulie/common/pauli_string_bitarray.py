"""Representation of a Pauli string as a bitarray."""
from typing import Self, Generator
from six.moves import reduce
import numpy as np
from bitarray import bitarray
from bitarray.util import count_and, ba2int
from paulie.common.pauli_string_parser import pauli_string_parser

CODEC = {
    "I": bitarray([0, 0]), 
    "X": bitarray([1, 0]), 
    "Y": bitarray([1, 1]), 
    "Z": bitarray([0, 1]),
}

Si = np.array([[1,0],[0,1]])
Sx = np.array([[0,1],[1,0]])
Sy = np.array([[0,-1j],[1j,0]])
Sz = np.array([[1,0],[0,-1]])

class PauliString:
    """Representation of a Pauli string as a bitarray."""

    def __init__(self, n: int = None, pauli_str: str = None, bits: bitarray = None) -> None:
        """Initialize a Pauli string.
        
        Args:
            n: length of the Pauli string
            pauli_str: String representation of a Pauli string
            bits: Bit representation of a Pauli string
        """
        super().__init__()
        self.nextpos = 0
        if bits is not None:
            self.bits = bits.copy()
        elif n is not None and pauli_str is None:
            self.bits = bitarray(2 * n)
        elif pauli_str is not None:
            pauli_str = pauli_string_parser(pauli_str)
            self.bits = bitarray()
            self.bits.encode(CODEC, pauli_str)
            if n is not None and n > len(self):
                o = self + PauliString(n = n - len(self))
                self.bits = o.bits.copy()
        self.bits_even = self.bits[::2]
        self.bits_odd  = self.bits[1::2]

    def get_index(self) -> int:
        """
         Return index in matrix decomposition vector
        """
        return ba2int(self.bits)

    def get_diagonal_index(self) -> int:
        """
         Return index in diagonal matrix decomposition vector
        """
        if ba2int(self.bits_even) == 0:
            return ba2int(self.bits_odd)
        return -1

    def get_weight_in_matrix(self, b_matrix: np.ndarray) -> np.complex128:
        """
         Return weight in diagonal matrix decomposition vector
        """
        len_matrix = len(b_matrix)
        len_string = len(self)
        if len_matrix not in (2**len_string, 4**len_string):
            raise ValueError("Incorrect matrix size")
        if len_matrix == 2**len_string:
            index = self.get_diagonal_index()
            if index > -1:
                return b_matrix[index]
            return 0.0
        return b_matrix[self.get_index()]

    def create_instance(self, n: int = None, pauli_str: str = None):
        """Create a Pauli string instance.
           Args:
                n: legnth of the Pauli string
                pauli_str: String representation of a Pauli string
           Returns the intensity of a Pauli string
        """
        return PauliString(n=n, pauli_str=pauli_str)

    def __str__(self) -> str:
        """Convert PauliString to readable string (e.g., "XYZI")."""
        return "".join(self.bits.decode(CODEC))

    def _ensure_pauli_string(self, other:str|Self):
        return other if isinstance(other, PauliString) else PauliString(pauli_str=str(other))

    def __eq__(self, other:str|Self) -> bool:
        """Overloading the equality operator relating two Pauli strings.
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits == other.bits

    def __lt__(self, other:str|Self) -> bool:
        """
        Overloading < operator for two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits < other.bits

    def __le__(self, other:str|Self) -> bool:
        """
        Overloading <= operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits <= other.bits

    def __gt__(self, other:str|Self) -> bool:
        """
        Overloading > operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits > other.bits

    def __ge__(self, other:str|Self) -> bool:
        """
        Overloading >= operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits >= other.bits

    def __ne__(self, other:str|Self) -> bool:
        """
        Overloading != operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        other = self._ensure_pauli_string(other)
        return self.bits != other.bits

    def __hash__(self) -> int:
        """Make PauliString hashable so it can be used in sets"""
        return hash(str(self.bits))

    def __len__(self) -> int:
        """
        Returns the lenght of the Pauli string
        """
        return len(self.bits) // 2

    def __iter__(self) -> Self:
        """
        Pauli String Iterator
        """
        self.nextpos = 0
        return self

    def __next__(self) -> Self:
        """
        The value of the next position of the Pauli string
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = PauliString(bits=self.bits[2*self.nextpos:2*self.nextpos+2])
        self.nextpos += 1
        return value

    def __setitem__(self, position: int, pauli_string: str|Self):
        """
        Sets a specified Pauli at a given position in the Paulistring
        """
        self.set_substring(position, pauli_string)

    def __getitem__(self, position: int) -> Self:
        """
        Gets the Pauli at specified position
        """
        return self.get_substring(position)

    def __copy__(self) -> Self:
        """
        Pauli string copy operator
        """
        return PauliString(bits=self.bits)

    def copy(self) -> Self:
        """ Copy Pauli string """
        return PauliString(bits=self.bits)

    def __add__(self, other:str|Self):
        """
        Pauli string addition operator
        """
        other = self._ensure_pauli_string(other)
        return self.tensor(other)

    def __or__(self, other:str|Self)->bool:
        """
        Overloading | operator of two Pauli strings like commutes_with
        """
        return self.commutes_with(other)

    def __xor__(self, other:str|Self):
        """
        Overloading ^ operator of two Pauli strings like adjoint_map
        """
        return self.adjoint_map(other)

    def __matmul__(self, other:str|Self):
        """
        Overloading @ operator of two Pauli strings like multiply
        """
        return self.multiply(other)

    def sign(self, other:Self):
        """
        Sign of multiplication of two pauli string
        return +- 1
        """
        if self|other:
            return 1
        s = 1
        for i in range(0, len(self)):
            if ((self[i] == 'I' and other[i] == "I") or
                (self[i] == 'I' and other[i] == "X") or
                (self[i] == 'I' and other[i] == "Y") or
                (self[i] == 'I' and other[i] == "Z") or
                (self[i] == 'X' and other[i] == "I") or
                (self[i] == 'Y' and other[i] == "I") or
                (self[i] == 'Z' and other[i] == "I")):
                continue

            if ((self[i] == 'X' and other[i] == "Y") or
                (self[i] == 'Y' and other[i] == "Z") or
                (self[i] == 'Z' and other[i] == "X")):
                s *= complex(0,1)
                continue
            s *= complex(0,-1)
        return s

    def commutes_with(self, other:str|Self) -> bool:
        """
        Check if this Pauli string commutes with another
        Returns True if they commute, False if they anticommute
        """
        # Compute symplectic product mod 2
        # Paulis commute iff the symplectic product is 0
        other = self._ensure_pauli_string(other)

        if len(self) != len(other):
            raise ValueError("Pauli arrays must be of equal length")
        return (count_and(self.bits_even, other.bits_odd) % 2 ==
               count_and(other.bits_even, self.bits_odd) % 2)

    def get_substring(self, start: int, length: int = 1) -> Self:
        """
        Get a substring of Paulis inside the Pauli string

        Args:
            start: Index to begin extracting the string.
            length: Length of each substring.

        Returns:
            substring of the Pauli string.
        """
        return PauliString(bits=self.bits[2*start:2*start+2*length])

    def set_substring(self, start: int, pauli_string:str|Self) -> None:
        """
        Set substring starting at position `start`
        """
        pauli_string = self._ensure_pauli_string(pauli_string)

        for i in range(0, len(pauli_string)):
            self.bits[2*start + 2*i] = pauli_string.bits[2*i]
            self.bits[2*start + 2*i + 1] = pauli_string.bits[2*i + 1]
            self.bits_even[start  + i] = pauli_string.bits_even[i]
            self.bits_odd[start + i] = pauli_string.bits_odd[i]

    def is_identity(self) -> bool:
        """Check if this Pauli string is the identity"""
        return bitarray(len(self.bits)) == self.bits

    def tensor(self, other: Self) -> Self:
        """Tensor product of this Pauli string with another"""
        new_bits = bitarray(len(self.bits) + len(other.bits))
        for i in range(len(new_bits)):
            s = self.bits if i < len(self.bits) else other.bits
            j = i if i < len(self.bits) else i - len(self.bits)
            new_bits[i] = s[j]

        return PauliString(bits=new_bits)

    def multiply(self, other:str|Self) -> Self:
        """
        Proportional multiplication operator of two Pauli strings
        Returns a PauliString proportional to the multiplication 
        """
        other = self._ensure_pauli_string(other)

        if len(self.bits) != len(other.bits):
            raise ValueError("Pauli arrays must have the same length")
        # Bitwise XOR is equivalent to mod-2 addition
        return PauliString(bits = self.bits ^ other.bits)

    def adjoint_map(self, other:str|Self) -> Self:
        """
        Compute the adjoint map ad_A(B) = [A,B]
        Returns None if the commutator is zero (i.e., if A and B commute)
        Otherwise returns a PauliString proportional to the commutator
        """
        other = self._ensure_pauli_string(other)

        if self.commutes_with(other):
            return None
        # For Pauli strings, if they anticommute, [A,B] = 2AB
        # In the context of generating a Lie algebra, we only care about
        # the result up to a constant factor
        # For anticommuting Paulis, the product gives a new Pauli
        # XOR of the bit vectors gives the non-phase part of the product
        if len(self.bits) != len(other.bits):
            raise ValueError("Pauli arrays must have the same length")
        # Bitwise XOR is equivalent to mod-2 addition
        return PauliString(bits = self.bits ^ other.bits)

    def inc(self) -> None:
        """
        Pauli string increment operator
        """
        for i in reversed(range(len(self.bits))):
            if self.bits[i] == 0:
                self.bits[i] = 1
                break
            self.bits[i] = 0
        self.bits_even = self.bits[::2]
        self.bits_odd  = self.bits[1::2]

    def expand(self, n: int) -> Self:
        """
        Increasing the size of the Pauli string by taking the tensor product
        with identities in the end
        Args:
            n (int): New Pauli string length
        Returns the Pauli string of extend length
        """
        return self + PauliString(n = n - len(self))

    def gen_all_pauli_strings(self) -> Generator[list[Self], None, None]:
        """
        Generate a list of Pauli strings that commute with this string
        Yields the commutant of the Pauli string
        """
        n = len(self)
        pauli_string = PauliString(n=n)

        last = PauliString(bits = bitarray([1] * (2 * n)))

        while pauli_string !=last:
            yield pauli_string.copy()
            pauli_string.inc()
        yield pauli_string.copy()

    def get_commutants(self, generators:list[Self] = None) -> list[Self]:
        """
        Get a list of Pauli strings that commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if self|g]

    def get_anti_commutants(self, generators:list[Self] = None) -> list[Self]:
        """
        Get a list of Pauli strings that no-commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if not self|g]



    def get_nested(self, generators:list[Self] = None) ->list[tuple[Self, Self]]:
        """
        Get nested of Pauli string
        Args:
            generators: Collection of Pauli strings on which nested is searched
                        If not specified, then the search area is all Pauli strings of the same size
        """

        # Retrieve the Pauli strings that anticommute with self.
        anti_commuting = self.get_anti_commutants(generators=generators)
        nested_pairs = set()

        # Iterate through all anti-commuting Pauli strings
        for g in anti_commuting:
            # Compute the adjoint map (or the product) once
            adj = g @ self
            # Use canonical ordering to ensure the pair is unique: store the pair as (min, max).
            canonical_pair = (g, adj) if g < adj else (adj, g)
            nested_pairs.add(canonical_pair)

        return list(nested_pairs)

    def _match_matrix(self, v:str) -> np.array:
        """
         Matching matrix for the string item
         Args: v - a item of PauliString
         Returns: Matrix representation for the string item 
        """
        match v:
            case "I":
                return Si
            case "X":
                return Sx
            case "Y":
                return Sy
            case "Z":
                return Sz

    def get_matrix(self) -> np.array:
        """
        Get matrix representation for Pauli string
        Returns: Matrix representation for the Pauli string
        """
        return reduce(lambda matrix, v: np.kron(matrix, self._match_matrix(v))
                      if matrix is not None else self._match_matrix(v), str(self), None)
