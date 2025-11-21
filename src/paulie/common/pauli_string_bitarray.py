"""Representation of a Pauli string as a bitarray."""
from typing import Self, Generator
from six.moves import reduce
import numpy as np
from bitarray import bitarray
from bitarray.util import count_and, count_or, ba2int

from paulie.common.pauli_string_parser import pauli_string_parser

CODEC = {
    "I": bitarray([0, 0]), 
    "X": bitarray([1, 0]), 
    "Y": bitarray([1, 1]), 
    "Z": bitarray([0, 1]),
}

SI = np.array([[1,0],[0,1]])
SX = np.array([[0,1],[1,0]])
SY = np.array([[0,-1j],[1j,0]])
SZ = np.array([[1,0],[0,-1]])

class PauliString:
    """Representation of a Pauli string as a bitarray."""

    def __init__(self, n: int = None, pauli_str: str = None, bits: bitarray = None) -> None:
        """Initialize a Pauli string.
        
        Args:
            n: Length of the Pauli string.
            pauli_str: String representation of a Pauli string.
            bits: Bits representation of a Pauli string.
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
        Get index in matrix decomposition vector.
        Returns:
            Index in matrix decomposition vector.
        """
        return ba2int(self.bits)

    def get_diagonal_index(self) -> int:
        """
        Get index in diagonal matrix decomposition vector.
        Returns:
            Index in matrix decomposition vector.
        """

        if ba2int(self.bits_even) == 0:
            return ba2int(self.bits_odd)
        return -1

    def get_weight_in_matrix(self, b_matrix: np.ndarray) -> np.complex128:
        """
        Get weight in diagonal matrix decomposition vector.
        Args: 
            b_matrix: Matrix decomposition vector.
        Returns:
            Weight in diagonal matrix decomposition vector.
        Raises:
            ValueError:
                Incorrect matrix size.
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
                n: Length of the Pauli string.
                pauli_str: String representation of a Pauli string.
           Returns:
               Intensity of a Pauli string.
        """
        return PauliString(n=n, pauli_str=pauli_str)

    def __str__(self) -> str:
        """
        Convert PauliString to readable string (e.g., "XYZI").
        Returns:
            String representation.
        """
        return "".join(self.bits.decode(CODEC))

    def _ensure_pauli_string(self, other:str|Self):
        """
        Get self Pauli string representation.
        Args: 
            other: Pauli string representation.
        Returns:
            Self representation.
        """
        return other if isinstance(other, PauliString) else PauliString(pauli_str=str(other))

    def __eq__(self, other:str|Self) -> bool:
        """
        Overloading the equality operator relating two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits == other.bits

    def __lt__(self, other:str|Self) -> bool:
        """
        Overloading < operator for two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits < other.bits

    def __le__(self, other:str|Self) -> bool:
        """
        Overloading <= operator of two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits <= other.bits

    def __gt__(self, other:str|Self) -> bool:
        """
        Overloading > operator of two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits > other.bits

    def __ge__(self, other:str|Self) -> bool:
        """
        Overloading >= operator of two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits >= other.bits

    def __ne__(self, other:str|Self) -> bool:
        """
        Overloading != operator of two Pauli strings.
        Args:
             other: Pauli string to compare with.
        Returns:
            Result of the comparison.
        """
        other = self._ensure_pauli_string(other)
        return self.bits != other.bits

    def __hash__(self) -> int:
        """
        Make Pauli string hashable so it can be used in sets.
        Returns:
            Hash of Pauli string.
        """
        return hash(str(self.bits))

    def __len__(self) -> int:
        """
        Get the length of the Pauli string.
        Returns:
            Length of the Pauli string.
        """
        return len(self.bits) // 2

    def __iter__(self) -> Self:
        """
        Get Pauli String Iterator.
        Returns:
            Iterator of the Pauli string.
        """
        self.nextpos = 0
        return self

    def __next__(self) -> Self:
        """
        Get the value of the next position of the Pauli string.
        Returns:
            Value of the next position of the Pauli string.
        Raises:
            StopIteration: End of Pauli string.
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = PauliString(bits=self.bits[2*self.nextpos:2*self.nextpos+2])
        self.nextpos += 1
        return value

    def __setitem__(self, position: int, pauli_string: str|Self):
        """
        Sets a specified Pauli at a given position in the Pauli string.
        Args:
            pauli_string: Pauli string to set in position.
        Returns:
            None
        """
        self.set_substring(position, pauli_string)

    def __getitem__(self, position: int) -> Self:
        """
        Gets the PauliString at specified position.
        Args:
            position: Position in Pauli string.
        Returns:
            Pauli string at specified position.
        """
        return self.get_substring(position)

    def __copy__(self) -> Self:
        """
        Pauli string copy operator.
        Returns:
            Copy of self.
        """
        return PauliString(bits=self.bits)

    def copy(self) -> Self:
        """ 
        Copy Pauli string.
        Returns:
            Copy of self.
        """
        return PauliString(bits=self.bits)

    def __add__(self, other:str|Self):
        """
        Pauli string addition operator.
        Args:
             other: Pauli string to add.
        Returns:
            Result of adding.
        """
        other = self._ensure_pauli_string(other)
        return self.tensor(other)

    def __or__(self, other:str|Self)->bool:
        """
        Overloading | operator of two Pauli strings like commutes_with.
        Args:
             other: Pauli string to commutes_with.
        Returns:
            Result of commutes_with.
        """
        return self.commutes_with(other)

    def __xor__(self, other:str|Self):
        """
        Overloading ^ operator of two Pauli strings like adjoint_map.
        Args:
             other: Pauli string to adjoint_map.
        Returns:
            Result of adjoint_map.
        """
        return self.adjoint_map(other)

    def __matmul__(self, other:str|Self):
        """
        Overloading @ operator of two Pauli strings like multiply.
        Args:
             other: Pauli string to multiply.
        Returns:
            Result of multiply.
        """
        return self.multiply(other)

    def sign(self, other: 'PauliString') -> complex:
        """
        Calculates the phase of the product of two Pauli strings: self * other.
        The product is defined as P1 * P2 = phase * P3. This method returns the phase.

        This implementation uses the correct symplectic product formalism, which can
        be found in various literature, including the supplemental material of
        the paper referenced in the related GitHub issue.
        See also: arxiv.org:2405.19287

        Args:
            other (PauliString): Pauli string to multiply with.
        Returns:
            Complex phase of the product (1, -1, 1j, or -1j) .
        Raises:
            ValueError:
                len(self) != len(other)
        """
        other = self._ensure_pauli_string(other)
        if len(self) != len(other):
            raise ValueError("Pauli arrays must have the same length for multiplication.")

        # This is the full, correct formula for the exponent f in phase = i^f.
        # It is based on the bit-array representations of the two Pauli strings.
        # self.bits_even corresponds to the X part, self.bits_odd to the Z part.
        f = 2 * count_and(self.bits_even, other.bits_odd) + \
            count_and(self.bits_odd, self.bits_even) + \
            count_and(other.bits_odd, other.bits_even) - \
            count_and(self.bits_even ^ other.bits_even,
                      self.bits_odd ^ other.bits_odd)

        # The final phase is (-1j)^f mod 4.
        return (-1j) ** (f % 4)

    def complex_conj(self: Self) -> tuple[complex, str | Self]:
        """
        Get the complex conjugate of the Pauli string.
        Returns:
            Complex conjugate of the Pauli string.
        """
        ys = count_and(self.bits_odd, self.bits_even)
        return ((-1)**(ys), self)

    def commutes_with(self, other:str|Self) -> bool:
        """
        Check if this Pauli string commutes with another.
        Args:
            other: Pauli string to commutes.
        Returns:
            True if they commute, False if they anticommute.
        Raises:
            ValueError:
                len(self) != len(other)
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
        Get a substring of Paulis inside the Pauli string.
        Args:
            start: Index to begin extracting the string.
            length: Length of each substring.
        Returns:
            Substring of the Pauli string.
        """
        return PauliString(bits=self.bits[2*start:2*start+2*length])

    def set_substring(self, start: int, pauli_string:str|Self) -> None:
        """
        Set substring starting at position `start`.
        Args:
            start: Index to begin the string.
            pauli_string: Substring of PauliString.
        Returns:
            None
        """
        pauli_string = self._ensure_pauli_string(pauli_string)

        for i in range(0, len(pauli_string)):
            self.bits[2*start + 2*i] = pauli_string.bits[2*i]
            self.bits[2*start + 2*i + 1] = pauli_string.bits[2*i + 1]
            self.bits_even[start  + i] = pauli_string.bits_even[i]
            self.bits_odd[start + i] = pauli_string.bits_odd[i]

    def is_identity(self) -> bool:
        """
        Check if this Pauli string is the identity.
        Returns:
            True if self is the identity.
        """
        return bitarray(len(self.bits)) == self.bits

    def tensor(self, other: Self) -> Self:
        """
        Tensor product of this Pauli string with another.
        Args:
            other: Pauli string to tensor product.
        Returns:
            Result of the tensor product self on another.
        """
        new_bits = bitarray(len(self.bits) + len(other.bits))
        for i in range(len(new_bits)):
            s = self.bits if i < len(self.bits) else other.bits
            j = i if i < len(self.bits) else i - len(self.bits)
            new_bits[i] = s[j]

        return PauliString(bits=new_bits)

    def multiply(self, other:str|Self) -> Self:
        """
        Proportional multiplication operator of two Pauli strings.
        Args:
            other: Pauli string to multiplication.
        Returns:
            PauliString proportional to the multiplication.
        Raises:
            ValueError:
                 len(self.bits) != len(other.bits)
        """
        other = self._ensure_pauli_string(other)

        if len(self.bits) != len(other.bits):
            raise ValueError("Pauli arrays must have the same length")
        # Bitwise XOR is equivalent to mod-2 addition
        return PauliString(bits = self.bits ^ other.bits)

    def adjoint_map(self, other:str|Self) -> Self:
        """
        Compute the adjoint map ad_A(B) = [A,B].
        Args:
            other: Pauli string to adjoint map with self.
        Returns: 
            None if the commutator is zero (i.e., if A and B commute). Otherwise, returns a PauliString proportional to the commutator.
        Raises:
            ValueError:
                len(self.bits) != len(other.bits)
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

    def inc(self) -> Self:
        """
        Pauli string increment operator.
        Returns:
            Pauli string whose bit representation is greater than 1.
        """
        for i in reversed(range(len(self.bits))):
            if self.bits[i] == 0:
                self.bits[i] = 1
                break
            self.bits[i] = 0
        self.bits_even = self.bits[::2]
        self.bits_odd  = self.bits[1::2]
        return self

    def expand(self, n: int) -> Self:
        """
        Increasing the size of the Pauli string by taking the tensor product with identities in the end.
        Args:
            n (int): New Pauli string length.
        Returns:
            Pauli string of extend length.
        """
        return self + PauliString(n = n - len(self))

    def gen_all_pauli_strings(self) -> Generator[list[Self], None, None]:
        """
        Generate a list of Pauli strings that commute with this string.
        Yields:
            Commutant of the Pauli string.
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
        Get a list of Pauli strings that commute with this string.
        Args:
            generators: Collection of Pauli strings on which commutant is searched. If not specified, then the search area is all Pauli strings of the same size.
        Returns:
            List of Pauli strings that commute with this string.
        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if self|g]

    def get_anti_commutants(self, generators:list[Self] = None) -> list[Self]:
        """
        Get a list of Pauli strings that no-commute with this string.
        Args:
            generators: Collection of Pauli strings on which commutant is searched. If not specified, then the search area is all Pauli strings of the same size.
        Returns:
            List of Pauli strings that no-commute with this string.
        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if not self|g]

    def get_nested(self, generators:list[Self] = None) ->list[tuple[Self, Self]]:
        """
        Get nested of Pauli string.
        Args:
            generators: Collection of Pauli strings on which commutant is searched. If not specified, then the search area is all Pauli strings of the same size.
        Returns:
            Nested of Pauli string.
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
         Matching matrix for the string item.
         Args:
             v: Item of Pauli string.
         Returns:
             Matrix representation for the string item.
        """
        match v:
            case "I":
                return SI
            case "X":
                return SX
            case "Y":
                return SY
            case "Z":
                return SZ

    def get_matrix(self) -> np.array:
        """
        Get matrix representation for Pauli string.
        Returns:
            Matrix representation for the Pauli string.
        """
        return reduce(lambda matrix, v: np.kron(matrix, self._match_matrix(v))
                      if matrix is not None else self._match_matrix(v), str(self), None)

    def get_count_non_trivially(self) -> int:
        """ 
        Get count non-trivially.
        Returns:
            Count non-trivially.
        """
        return count_or(self.bits_even, self.bits_odd)
