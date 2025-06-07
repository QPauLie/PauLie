"""Representation of a Pauli string as a bitarray."""

from copy import deepcopy
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
        self.combinations = [[c[0], PauliString(pauli_str=str(c[1]))] for c in combinations]


    def _gtzero(self, z: complex) -> bool:
        if z.real > 0:
            return True
        if np.isclose(z.real, 0) and z.imag > 0:
            return True
        return False

    def _print_complex(self, z: complex):
        if np.isclose(z.imag, 0):
            return "" if np.isclose(abs(z.real), 1) else str(abs(z.real)) + "*"

        if np.isclose(z.real, 0):
            return "i*" if np.isclose(abs(z.imag), 1) else str(abs(z.imag)) + "i*"

        if z.real > 0:
            return "(" + str(z.real) + ("-" if z.imag < 0 else "+") + str(abs(z.imag)) + "i)*"

        return "(" + str(-z.real) + ("-" if z.imag > 0 else "+") + str(abs(z.imag)) + "i)*"

    def __str__(self) -> str:
        """Convert PauliStringLinear to readable string (e.g., 7*"XYZI" + 5*"ZZYX")."""
        str_value = ''
        for i, c in enumerate(self.combinations):
            if i == 0:
                if not self._gtzero(c[0]):
                    str_value += '- '
                str_value += self._print_complex(c[0]) + str(c[1])
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

        for c in self.combinations:
            is_eq = None
            for o in other:
                if o[1] == c[1]:
                    is_eq = np.isclose(np.real(c[0]), np.real(o[0])) \
                            and np.isclose(np.imag(c[0]), np.imag(o[0]))
                    break
            if is_eq is None and not np.isclose(np.abs(c[0]), 0):
                return False
            if is_eq is False:
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
        combinations = deepcopy(self.combinations)
        otherdict = {c[1]: c[0] for c in other.combinations}
        for c in combinations:
            if c[1] in otherdict:
                c[0] += otherdict[c[1]]
                del otherdict[c[1]]
        for c1, c0 in otherdict.items():
            combinations.append([c0, c1])
        return PauliStringLinear(combinations)

    def __iadd__(self, other:Self):
        """
        Linear combination of Pauli string addition operator
        """
        otherdict = {c[1]: c[0] for c in other.combinations}
        for c in self.combinations:
            if c[1] in otherdict:
                c[0] += otherdict[c[1]]
                del otherdict[c[1]]
        for c1, c0 in otherdict.items():
            self.combinations.append([c0, c1])
        return self

    def __mul__(self, num:complex):
        """
        Multiply by a constant
        """
        combinations = deepcopy(self.combinations)
        for c in combinations:
            c[0] *= num
        return PauliStringLinear(combinations)

    def __rmul__(self, num:complex):
        """
        Multiply by a constant
        """
        combinations = deepcopy(self.combinations)
        for c in combinations:
            c[0] *= num
        return PauliStringLinear(combinations)

    def __imul__(self, num:complex):
        """
        Multiply by a constant
        """
        for c in self.combinations:
            c[0] *= num
        return self

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
            new_combinations.append((c[0]*other.sign(c[1]), other@c[1]))
        return PauliStringLinear(new_combinations)

    def multiply(self, other:PauliString|Self) -> Self:
        """
        Multiplication operator of two linear combination
        of Pauli strings
        Returns a PauliString proportional to the multiplication 
        """
        new_combinations = []
        # isinstance also returns True if other is a subclass of PauliString
        # Therefore check instead if other is NOT of type PauliStringLinear
        if not isinstance(other, PauliStringLinear):
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
        Returns a linear comination of PauliString 
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0]*basis.sign(c[1]), c[1] + basis@c[1]))
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

    def trace(self) -> complex:
        """Get the trace of the matrix representation for Pauli string"""
        return sum(c[0] * c[1].trace() for c in self.combinations)
