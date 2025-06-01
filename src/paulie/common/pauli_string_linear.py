"""Representation of a Pauli string as a bitarray."""

from typing import Self, Generator

from six.moves import reduce
import numpy as np
from paulie.common.pauli_string_bitarray import PauliString

class PauliStringLinearException(Exception):
    """
    Exception for the linear combination of Pauli strings class
    """


class PauliStringLinear(PauliString):
    """Representation of a linear combination of Pauli string."""

    def __init__(self, combinations: list[tuple[complex, str|PauliString]]) -> None:
        """Initialize a linear combination of Pauli strings.
        
        Args:
            combination: list of tuple (weight, Pauli string),
            weight - weight of Pauli string in linear combination,
            Pauli string - Pauli string like PauliString or string
        """
        self.nextpos = 0
        self.combinations = [(c[0], PauliString(pauli_str=str(c[1]))) for c in combinations]


    def _gtzero(self, z: complex) -> bool:
        if z.real > 0:
            return True
        if z.real == 0 and z.imag > 0:
            return True
        return False

    def _print_complex(self, z: complex):
        if z.imag == 0:
            return "" if abs(z.real) == 1 else str(abs(z.real)) + "*"
        if z.real == 0:
            return "i*" if abs(z.imag) == 1 else str(abs(z.imag)) + "*i*"
        if z.real > 0:
            return "(" + z.real + "-" if z.imag < 0 else "+" + abs(z.imag) + ")*"

        return "(" + abs(z.real) + "-" if z.imag > 0 else "+" + abs(z.imag) + ")*"

    def __str__(self) -> str:
        """Convert PauliStringLinear to readable string (e.g., 7*"XYZI" + 5*"ZZYX")."""
        str_value = ''
        for i, c in enumerate(self.combinations):
            if i == 0:
                str_value = self._print_complex(c[0]) + str(c[1])
                continue
            if self._gtzero(c[0]):
                str_value += ' + ' + self._print_complex(c[0]) + str(c[1])
            else:
                str_value += ' - ' + self._print_complex(c[0]) + str(c[1])

        return str_value


    def __eq__(self, other: Self) -> bool:
        """Overloading the equality operator relating two linear combination of Pauli strings.
        Args:
             other: The linear combination of Pauli strings to compare with
        Returns the result of the comparison
        """
        if len(self) != len(other):
            return False

        for c in self.combinations:
            is_eq = None
            for o in other:
                if o[1] == c[1]:
                    is_eq = c[0] == o[0]
                    break
            if not is_eq:
                return False
        return True

    def __lt__(self, other:Self) -> bool:
        """
        Overloading < operator for two linear combination of Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        raise PauliStringLinearException("Not implemented")

    def __le__(self, other:Self) -> bool:
        """
        Overloading <= operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        raise PauliStringLinearException("Not implemented")

    def __gt__(self, other:Self) -> bool:
        """
        Overloading > operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        raise PauliStringLinearException("Not implemented")

    def __ge__(self, other:Self) -> bool:
        """
        Overloading >= operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        raise PauliStringLinearException("Not implemented")

    def __ne__(self, other:Self) -> bool:
        """
        Overloading != operator of two Pauli strings
        Args:
             other: The Pauli string to compare with
        Returns the result of the comparison
        """
        raise PauliStringLinearException("Not implemented")

    def __hash__(self) -> int:
        """Make PauliStringLinear hashable so it can be used in sets"""
        return hash("".join([str(c[0]) + str(c[1]) for c in self.combinations]))

    def __len__(self) -> int:
        """
        Returns the lenght of the Pauli string
        """
        return len(self.combinations)

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
        value = self.combinations[self.nextpos]
        self.nextpos += 1
        return value

    def __setitem__(self, position: int, combination: tuple[complex,PauliString]):
        """
        Sets a specified Pauli at a given position in the Paulistring
        """
        self.combinations[position] = combination

    def __getitem__(self, position: int) -> Self:
        """
        Gets the Pauli at specified position
        """
        return self.combinations[position]

    def __copy__(self) -> Self:
        """
        Pauli string linear combination copy operator
        """
        return PauliStringLinear(self.combinations)

    def copy(self) -> Self:
        """ Copy Linear combination of Pauli strings """
        return PauliStringLinear(self.combinations)

    def __add__(self, other:Self):
        """
        Linear combination of Pauli string addition operator
        """
        combinations = self.combinations.copy()
        for o in other:
            is_found = False
            for c in combinations:
                if o[1] == c[1]:
                    c[0] += o[0]
                    break
            if not is_found:
                combinations.append(o)
        return PauliStringLinear(combinations)


    def __or__(self, other:Self)->bool:
        """
        Overloading | operator of two Pauli strings like commutes_with
        """
        return self.commutes_with(other)

    def __xor__(self, other:str|Self):
        """
        Overloading ^ operator of two linear combination of Pauli strings like adjoint_map
        """
        return self.adjoint_map(other)

    def __matmul__(self, other:PauliString|Self):
        """
        Overloading @ operator of two Pauli strings like multiply
        """
        return self.multiply(other)

    def __rmatmul__(self, other:PauliString):
        """
        Overloading @ operator of two Pauli strings like multiply
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0]*other.sign(c[1]), c[1]@other))
        return PauliStringLinear(new_combinations)

    def multiply(self, other:PauliString|Self) -> Self:
        """
        Multiplication operator of two linear combination
        of Pauli strings
        Returns a PauliString proportional to the multiplication 
        """
        new_combinations = []
        if isinstance(other, PauliString):
            for c in self.combinations:
                new_combinations.append((c[0]*c[1].sign(other), c[1]@other))
            return PauliStringLinear(new_combinations)

        for c in self.combinations:
            for o in other:
                new_combinations.append((c[0]*o[0]*c[1].sign(o[1]), c[1]@o[1]))
        return PauliStringLinear(new_combinations)

    def commutes_with(self, other:str|Self) -> bool:
        """
        Check if this Pauli string commutes with another
        Returns True if they commute, False if they anticommute
        """
        # Compute symplectic product mod 2
        # Paulis commute iff the symplectic product is 0
        for o in other:
            for c in self.combinations:
                if c[1]|o[1]:
                    return False
        return True

    def get_substring(self, start: int, length: int = 1) -> Self:
        """
        Get a substring of Paulis inside the Pauli string

        Args:
            start: Index to begin extracting the string.
            length: Length of each substring.

        Returns:
            substring of the Pauli string.
        """
        raise PauliStringLinearException("Not implemented")

    def set_substring(self, start: int, pauli_string:str|Self) -> None:
        """
        Set substring starting at position `start`
        """
        raise PauliStringLinearException("Not implemented")


    def is_identity(self) -> bool:
        """Check if this Pauli string is the identity"""
        raise PauliStringLinearException("Not implemented")

    def tensor(self, other: Self) -> Self:
        """Tensor product of this Pauli string with another"""
        raise PauliStringLinearException("Not implemented")


    def kron(self, other:PauliString):
        """
        Kroniker multiplication pauli string on linear combination
        of Pauli strings
        Returns a linera comination of PauliString 
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0], c[1] + other))
        return PauliStringLinear(new_combinations)

    def rkron(self, other:PauliString):
        """
        Right Kroniker multiplication pauli string on linear combination
        of Pauli strings
        Returns a linera comination of PauliString 
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0], other + c[1]))
        return PauliStringLinear(new_combinations)


    def quadratic(self, basis:PauliString):
        """
        Quadratic form
        Returns a linera comination of PauliString 
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0]*basis.sign(c[1]), c[1] + c[1]@basis))
        return PauliStringLinear(new_combinations)


    def adjoint_map(self, other:str|Self) -> Self:
        """
        Compute the adjoint map ad_A(B) = [A,B]
        Returns None if the commutator is zero (i.e., if A and B commute)
        Otherwise returns a PauliString proportional to the commutator
        """
        raise PauliStringLinearException("Not implemented")

    def inc(self) -> None:
        """
        Pauli string increment operator
        """
        raise PauliStringLinearException("Not implemented")

    def expand(self, n: int) -> Self:
        """
        Increasing the size of the Pauli string by taking the tensor product
        with identities in the end
        Args:
            n (int): New Pauli string length
        Returns the Pauli string of extend length
        """
        raise PauliStringLinearException("Not implemented")

    def gen_all_pauli_strings(self) -> Generator[list[Self], None, None]:
        """
        Generate a list of Pauli strings that commute with this string
        Yields the commutant of the Pauli string
        """
        raise PauliStringLinearException("Not implemented")

    def get_commutants(self, generators:list[Self] = None) -> list[Self]:
        """
        Get a list of Pauli strings that commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        raise PauliStringLinearException("Not implemented")

    def get_anti_commutants(self, generators:list[Self] = None) -> list[Self]:
        """
        Get a list of Pauli strings that no-commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        raise PauliStringLinearException("Not implemented")


    def get_nested(self, generators:list[Self] = None) ->list[tuple[Self, Self]]:
        """
        Get nested of Pauli string
        Args:
            generators: Collection of Pauli strings on which nested is searched
                        If not specified, then the search area is all Pauli strings of the same size
        """

        # Retrieve the Pauli strings that anticommute with self.
        raise PauliStringLinearException("Not implemented")



    def get_matrix(self) -> np.array:
        """
        Get matrix representation for Pauli string
        Returns: Matrix representation for the Pauli string
        """

        return reduce(lambda matrix, c: matrix + c[0] * c[1].get_matrix()
                      if matrix is not None else c[0] * c[1].get_matrix(), self, None)
