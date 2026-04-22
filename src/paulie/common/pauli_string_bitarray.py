"""Representation of a Pauli string as a bitarray."""
from __future__ import annotations

from collections.abc import Generator
from typing import Self
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

DECODEC = {
    (0, 0): "I",
    (1, 0): "X",
    (1, 1): "Y",
    (0, 1): "Z",
}

SI = np.array([[1, 0], [0, 1]])
SX = np.array([[0, 1], [1, 0]])
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]])


class PauliString:
    """Representation of a Pauli string as a bitarray."""

    def __init__(self, n: int | None = None, pauli_str: str | None = None,
                 bits: bitarray | None = None,
                 bits_even: bitarray | None = None,
                 bits_odd: bitarray | None = None) -> None:
        """Initialize a Pauli string.

        Args:
            n (int, optional): Length of the Pauli string.
            pauli_str (str, optional): String representation of a Pauli string.
            bits (bitarray, optional): Bits representation of a Pauli string.
            bits_even (bitarray, optional): X-bits representation.
            bits_odd (bitarray, optional): Z-bits representation.
        """
        self.nextpos = 0
        self._bits = None
        if bits_even is not None and bits_odd is not None:
            self.bits_even = bits_even
            self.bits_odd = bits_odd
        elif bits is not None:
            self._bits = bits.copy()
            self.bits_even = self._bits[::2]
            self.bits_odd = self._bits[1::2]
        elif n is not None and pauli_str is None:
            self.bits_even = bitarray(n)
            self.bits_even.setall(0)
            self.bits_odd = bitarray(n)
            self.bits_odd.setall(0)
        elif pauli_str is not None:
            pauli_str = pauli_string_parser(pauli_str)
            temp_bits = bitarray()
            temp_bits.encode(CODEC, pauli_str)
            if n is not None and n > len(temp_bits) // 2:
                # Padding with identity
                pad_n = n - len(temp_bits) // 2
                temp_bits.extend(bitarray([0, 0] * pad_n))
            self._bits = temp_bits
            self.bits_even = self._bits[::2]
            self.bits_odd = self._bits[1::2]
        else:
            self.bits_even = bitarray(0)
            self.bits_odd = bitarray(0)

    @property
    def bits(self) -> bitarray:
        """
        Return the interleaved bit representation.

        The bitarray is constructed by interleaving ``bits_even`` and
        ``bits_odd`` such that even indices correspond to ``bits_even``
        and odd indices correspond to ``bits_odd``. The result is cached
        after the first construction.
        """
        if self._bits is None:
            n = len(self.bits_even)
            self._bits = bitarray(2 * n)
            self._bits[::2] = self.bits_even
            self._bits[1::2] = self.bits_odd
        return self._bits

    @bits.setter
    def bits(self, value: bitarray) -> None:
        """
        Set the interleaved bit representation.

        Updates the cached bitarray and reconstructs ``bits_even`` and
        ``bits_odd`` by splitting the input into even and odd indices.
        """
        self._bits = value
        self.bits_even = value[::2]
        self.bits_odd = value[1::2]

    def get_index(self) -> int:
        """
        Get index in matrix decomposition vector.

        Returns:
            int: Index in matrix decomposition vector.
        """
        return ba2int(self.bits)

    def get_diagonal_index(self) -> int:
        """
        Get index in diagonal matrix decomposition vector.

        Returns:
            int: Index in diagonal matrix decomposition vector.
        """

        if ba2int(self.bits_even) == 0:
            return ba2int(self.bits_odd)
        return -1

    def get_weight_in_matrix(self, b_matrix: np.ndarray) -> np.complex128:
        """
        Get weight in matrix decomposition vector.

        Args:
            b_matrix (numpy.ndarray): Matrix decomposition vector.
        Returns:
            complex: Weight in matrix decomposition vector.

        Raises:
            ValueError: Incorrect matrix size.
        """
        len_matrix = len(b_matrix)
        len_string = len(self)
        if len_matrix not in (2 ** len_string, 4 ** len_string):
            raise ValueError("Incorrect matrix size")
        if len_matrix == 2 ** len_string:
            index = self.get_diagonal_index()
            if index > -1:
                return np.complex128(b_matrix[index])
            return np.complex128(0.0)
        return np.complex128(b_matrix[self.get_index()])

    def create_instance(self, n: int | None = None, pauli_str: str | None = None) -> PauliString:
        """Create a Pauli string instance.

           Args:
                n (int, optional): Length of the Pauli string.
                pauli_str (str, optional): String representation of a Pauli string.
           Returns:
               PauliString: Instance of a Pauli string.
        """
        return PauliString(n=n, pauli_str=pauli_str)

    def __repr__(self) -> str:
        """Convert PauliString to readable string.

        Returns:
            str: String representation.
        """

        return f"PauliString({str(self)})"

    def __str__(self) -> str:
        """Convert PauliString to readable string.

        Returns:
            str: String representation.
        """
        return "".join(DECODEC[(self.bits_even[i], self.bits_odd[i])] for i in range(len(self)))

    def _ensure_pauli_string(self, other: object) -> PauliString:
        """
        Get self Pauli string representation.

        Args:
            other: Pauli string representation.
        Returns:
            PauliString: Self representation.
        """
        return other if isinstance(other, PauliString) else PauliString(pauli_str=str(other))

    def __eq__(self, other: object) -> bool:
        """
        Overloading the equality operator relating two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits_even == other.bits_even and self.bits_odd == other.bits_odd

    def __lt__(self, other: object) -> bool:
        """
        Overloading < operator for two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits < other.bits

    def __le__(self, other: object) -> bool:
        """
        Overloading <= operator of two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits <= other.bits

    def __gt__(self, other: object) -> bool:
        """
        Overloading > operator of two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits > other.bits

    def __ge__(self, other: object) -> bool:
        """
        Overloading >= operator of two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits >= other.bits

    def __ne__(self, other: object) -> bool:
        """
        Overloading != operator of two Pauli strings.

        Args:
             other: Pauli string to compare with.
        Returns:
            bool: Result of the comparison.
        """
        if not isinstance(other, PauliString):
            other = self._ensure_pauli_string(other)
        return self.bits_even != other.bits_even or self.bits_odd != other.bits_odd

    def __hash__(self) -> int:
        """
        Make Pauli string hashable so it can be used in sets.

        Returns:
            int: Hash of Pauli string.
        """
        return hash((self.bits_even.tobytes(), self.bits_odd.tobytes()))

    def __len__(self) -> int:
        """
        Get the length of the Pauli string.

        Returns:
            int: Length of the Pauli string.
        """
        return len(self.bits_even)

    def __iter__(self) -> Self:
        """
        Get Pauli string iterator.

        Returns:
            Iterator of the Pauli string.
        """
        self.nextpos = 0
        return self

    def __next__(self) -> PauliString:
        """
        Get the value of the next position of the Pauli string.

        Returns:
            PauliString: Value of the next position of the Pauli string.

        Raises:
            StopIteration: End of Pauli string.
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = PauliString(bits_even=self.bits_even[self.nextpos:self.nextpos+1],
                          bits_odd=self.bits_odd[self.nextpos:self.nextpos+1])
        self.nextpos += 1
        return value

    def __setitem__(self, position: int, pauli_string: str | Self) -> None:
        """
        Sets a specified Pauli at a given position in the Pauli string.

        Args:
            position (int): Position in the Pauli string to modify.
            pauli_string (str|Self): Pauli to set at `position`.
        Returns:
            None
        """
        self.set_substring(position, pauli_string)

    def __getitem__(self, position: int) -> PauliString:
        """
        Gets the Pauli at the specified position.

        Args:
            position (int): Position in Pauli string.
        Returns:
            PauliString: Pauli at the specified position.
        """
        return self.get_substring(position)

    def __copy__(self) -> PauliString:
        """
        Pauli string copy operator.

        Returns:
            PauliString: Copy of self.
        """
        return self.copy()

    def copy(self) -> PauliString:
        """
        Copy Pauli string.

        Returns:
            PauliString: Copy of self.
        """
        return PauliString(bits_even=self.bits_even.copy(),
                          bits_odd=self.bits_odd.copy())

    def __add__(self, other: object) -> PauliString:
        """
        Pauli string addition operator.

        Args:
             other: Pauli string to add.
        Returns:
            PauliString: Result of adding.
        """
        other = self._ensure_pauli_string(other)
        return self.tensor(other)

    def __or__(self, other: object) -> bool:
        """
        Overloading | operator of two Pauli strings like commutes_with.

        Args:
             other: Pauli string to commutes_with.
        Returns:
            bool: Result of commutes_with.
        """
        return self.commutes_with(other)

    def __xor__(self, other: object) -> PauliString | None:
        """
        Overloading ^ operator of two Pauli strings like adjoint_map.

        Args:
             other: Pauli string to adjoint_map.
        Returns:
            PauliString|None: Result of adjoint_map.
        """
        return self.adjoint_map(other)

    def __matmul__(self, other: object) -> PauliString:
        """
        Overloading @ operator of two Pauli strings like multiply.

        Args:
             other: Pauli string to multiply.
        Returns:
            PauliString: Result of multiply.
        """
        return self.multiply(other)

    def sign(self, other: PauliString) -> complex:
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
            complex: Complex phase of the product (1, -1, 1j, or -1j).

        Raises:
            ValueError: If the length of `self` is not equal to the length of `other`.
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

    def complex_conj(self: Self) -> tuple[complex, Self]:
        """
        Get the complex conjugate of the Pauli string.

        Returns:
            tuple[complex, Self]: Complex conjugate of the Pauli string.
        """
        ys = count_and(self.bits_odd, self.bits_even)
        return ((-1)**(ys), self)

    def commutes_with(self, other:object) -> bool:
        """
        Check if this Pauli string commutes with another.

        Args:
            other: Pauli string to commutes.
        Returns:
            bool: True if they commute, False if they anticommute.
        Raises:
            ValueError: If the length of `self` is not equal to the length of `other`.
        """
        # Compute symplectic product mod 2
        # Paulis commute iff the symplectic product is 0
        other = self._ensure_pauli_string(other)

        if len(self) != len(other):
            raise ValueError("Pauli arrays must be of equal length")
        return (count_and(self.bits_even, other.bits_odd) % 2 ==
               count_and(other.bits_even, self.bits_odd) % 2)

    def get_substring(self, start: int, length: int = 1) -> PauliString:
        """
        Get a substring of Paulis inside the Pauli string.

        Args:
            start (int): Starting index of the substring.
            length (int, optional): Length of the substring. Default is 1.
        Returns:
            PauliString: Substring of the Pauli string.
        """
        return PauliString(bits_even=self.bits_even[start:start+length],
                          bits_odd=self.bits_odd[start:start+length])

    def set_substring(self, start: int, pauli_string:str|PauliString) -> None:
        """
        Set substring starting at position `start`.

        Args:
            start (int): Starting index of the substring.
            pauli_string (str|PauliString): Pauli string to substitute.
        Returns:
            None
        """
        pauli_string = self._ensure_pauli_string(pauli_string)

        for i in range(0, len(pauli_string)):
            self.bits_even[start  + i] = pauli_string.bits_even[i]
            self.bits_odd[start + i] = pauli_string.bits_odd[i]
        # Invalidate cached bits
        self._bits = None

    def is_identity(self) -> bool:
        """
        Check if this Pauli string is the identity.

        Returns:
            bool: True if `self` is the identity.
        """
        return self.bits_even.count(1) == 0 and self.bits_odd.count(1) == 0

    def tensor(self, other: Self) -> PauliString:
        """
        Tensor product of this Pauli string with another.

        Args:
            other: Pauli string to tensor product.
        Returns:
            PauliString: Result of the tensor product of `self` with `other`.
        """
        new_even = self.bits_even.copy()
        new_even.extend(other.bits_even)
        new_odd = self.bits_odd.copy()
        new_odd.extend(other.bits_odd)
        return PauliString(bits_even=new_even, bits_odd=new_odd)

    def multiply(self, other:object) -> PauliString:
        """
        Proportional multiplication operator of two Pauli strings.

        Args:
            other: Pauli string to multiply.
        Returns:
            PauliString: Pauli string proportional to the multiplication of `self` with `other`.

        Raises:
            ValueError: If the length of `self` is not equal to the length of `other`.
        """
        other = self._ensure_pauli_string(other)

        if len(self) != len(other):
            raise ValueError("Pauli arrays must have the same length")
        return PauliString(bits_even=self.bits_even ^ other.bits_even,
                          bits_odd=self.bits_odd ^ other.bits_odd)

    def adjoint_map(self, other:object) -> PauliString|None:
        r"""
        Compute the adjoint map :math:`\mathrm{ad}_{A}(B) = [A, B]`.

        Args:
            other: Pauli string to adjoint map with self.
        Returns:
            PauliString|None: None if the commutator is zero (i.e., if A and B commute).

            Otherwise, returns a PauliString proportional to the commutator.

        Raises:
            ValueError: If the length of `self` is not equal to the length of `other`.
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

    def expand(self, n: int) -> PauliString:
        """
        Expand the Pauli string by tensoring it with identity operators at the end.

        Args:
            n (int): New Pauli string length.
        Returns:
            PauliString: Pauli string of expanded size.
        """
        return self + PauliString(n = n - len(self))

    def gen_all_pauli_strings(self) -> Generator[PauliString, None, None]:
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

    def get_commutants(self, generators: list[PauliString] | None = None) -> list[PauliString]:
        """
        Get a list of Pauli strings that commute with this string.

        Args:
            generators (list[PauliString], optional): Collection of Pauli strings on which commutant
                is searched. If None, then the search area is all Pauli strings of the same size.
        Returns:
            list[PauliString]: List of Pauli strings that commute with this string.
        """
        if generators is None:
            generators = list(self.gen_all_pauli_strings())

        return [g for g in generators if self | g]

    def get_anti_commutants(self, generators: list[PauliString] | None = None) -> list[PauliString]:
        """
        Get a list of Pauli strings that anticommute with this string.

        Args:
            generators (list[PauliString], optional): Collection of Pauli strings on which
                anticommutant is searched. If None, then the search area is all Pauli strings of the
                same size.
        Returns:
            list[PauliString]: List of Pauli strings that anticommute with this string.
        """
        if generators is None:
            generators = list(self.gen_all_pauli_strings())

        return [g for g in generators if not self | g]

    def get_nested(self, generators: list[PauliString] | None = None
                   ) -> list[tuple[PauliString, PauliString]]:
        """
        Get pairs of Pauli strings whose commutator is proportional to `self`.

        Args:
            generators (list[PauliString], optional): Collection of Pauli strings on which commutant
                is to be searched. If not specified, then the search area is all Pauli strings of
                the same size.
        Returns:
            list[tuple[PauliString, PauliString]]:
            List of pairs of Pauli strings whose commutator is proportional to `self`.
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

    def _match_matrix(self, v:str) -> np.ndarray:
        """
        Return the matrix representation of the Pauli operator.

        Args:
            v (str): String representation of the Pauli operator.
        Returns:
            numpy.ndarray: Matrix representation of the Pauli operator.
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
            case _:
                raise ValueError(f"Invalid Pauli: {v}")

    def get_matrix(self) -> np.ndarray:
        """
        Get the matrix representation of the Pauli string.

        Returns:
            numpy.ndarray: Matrix representation for the Pauli string.
        """
        pauli_str = str(self)
        matrix = self._match_matrix(pauli_str[0])
        for v in pauli_str[1:]:
            matrix = np.kron(matrix, self._match_matrix(v))
        return matrix

    def get_count_non_trivially(self) -> int:
        """
        Get the count of non-identity operators in the Pauli string.

        Returns:
            int: Count of non-identity operators.
        """
        return count_or(self.bits_even, self.bits_odd)

    def get_support(self) -> list[int]:
        """
        Get the indices of non-identity operators in the Pauli string.

        Returns:
            list[int]: List of indices where the Pauli operator is not the identity.
        """
        support_bits = self.bits_even | self.bits_odd
        return [i for i, bit in enumerate(support_bits) if bit]
