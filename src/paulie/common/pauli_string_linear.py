"""Representation of a Pauli string as a bitarray."""
from collections import defaultdict
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

    def __init__(self, combinations: list[tuple[complex, str | PauliString]]) -> None:
        """Initialize a linear combination of Pauli strings.
        
        Args:
            combination: list of tuple (weight, Pauli string),
            weight - weight of Pauli string in linear combination,
            Pauli string - Pauli string like PauliString or string
        """
        num_qubits = len(str(combinations[0][1])) if combinations else 0
        super().__init__(n=num_qubits)
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
        """Converts PauliStringLinear to a readable, sorted string."""
        if self.is_zero():
            return f"0.0*I{self.get_size()}"

        simplified_self = self.simplify()
        # Sort by the Pauli string part for a canonical, consistent output
        sorted_combinations = sorted(simplified_self.combinations, key=lambda term: str(term[1]))

        def _format_term(coeff, pauli_str):
            if np.isclose(coeff.imag, 0):
                val_str = f"{coeff.real:.8g}"
            elif np.isclose(coeff.real, 0):
                if np.isclose(coeff.imag, 1):
                    val_str = "i"
                elif np.isclose(coeff.imag, -1):
                    val_str = "-i"
                else: val_str = f"{coeff.imag:.8g}*i"
            else:
                val_str = f"({coeff.real:.8g}{coeff.imag:+.8g}j)"

            if val_str == "1" and pauli_str != "I"*len(pauli_str):
                return pauli_str
            if val_str == "-1" and pauli_str != "I"*len(pauli_str):
                return f"-{pauli_str}"
            return f"{val_str}*{pauli_str}"

        terms = []
        for i, (coeff, pauli) in enumerate(sorted_combinations):
            term_str = _format_term(coeff, str(pauli))
            if i == 0:
                terms.append(term_str)
            else:
                # Handle leading sign for subsequent terms
                if term_str.startswith('-'):
                    terms.append(f" - {term_str[1:]}")
                else:
                    terms.append(f" + {term_str}")
        return "".join(terms)


    def __eq__(self, other):
        """
        Performs a robust equality check between two PauliStringLinear objects.
        It simplifies both objects and compares their terms using a tolerance
        for floating-point coefficients.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented

        # Simplify both objects to get them into a canonical form
        self_simplified = self.simplify()
        other_simplified = other.simplify()

        # Create dictionaries for easy lookup: {pauli_str: coefficient}
        self_dict = {str(p): c for c, p in self_simplified.combinations}
        other_dict = {str(p): c for c, p in other_simplified.combinations}

        # Check if they have the same Pauli string terms
        if self_dict.keys() != other_dict.keys():
            return False

        # Check if the coefficients for each term are close enough
        for pauli_str, self_coeff in self_dict.items():
            other_coeff = other_dict[pauli_str]
            if not np.isclose(self_coeff, other_coeff):
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

    def __add__(self, other):
        """Performs a robust addition of two linear combinations."""
        # Use a dictionary to correctly sum coefficients
        summed_coeffs = defaultdict(complex)
        for coeff, pauli in self.combinations:
            summed_coeffs[str(pauli)] += coeff
        for coeff, pauli in other.combinations:
            summed_coeffs[str(pauli)] += coeff

        new_combinations = [(c, p) for p, c in summed_coeffs.items() if abs(c) > 1e-12]
        if not new_combinations:
            return PauliStringLinear([])
        return PauliStringLinear(new_combinations)

    def __iadd__(self, other):
        """Performs in-place addition."""
        # This calls our robust __add__ method and reassigns self
        new_self = self + other
        self.combinations = new_self.combinations
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

    def __matmul__(self, other: 'PauliStringLinear') -> 'PauliStringLinear':
        """
        Overloading @ operator of two Pauli strings like multiply
        Performs the distributive multiplication of two linear combinations of Pauli strings.
        This version assumes the PauliStringLinear object is directly iterable.
        """
        # pylint: disable=import-outside-toplevel
        from paulie.common.pauli_string_factory import get_pauli_string as p

        # Ensure the other object is also a PauliStringLinear for uniform processing.
        if not isinstance(other, self.__class__):
            raise TypeError(f"Unsupported operand type(s) for @: "
                            f"'{type(self).__name__}' and '{type(other).__name__}'")

        new_pauli_terms = []

        # Distributive multiplication: Loop through all pairs of terms.
        # We now loop over `self` and `other` directly.
        for self_coeff, self_pauli in self:
            for other_coeff, other_pauli in other:

                # Step 1: Calculate the phase of the Pauli product (P1 * P2)
                phase = self_pauli.sign(other_pauli)

                # Step 2: Calculate the resulting Pauli string (the phase-less product)
                product_pauli_string = self_pauli.multiply(other_pauli)

                # Step 3: Calculate the new term's overall coefficient: a * b * phase
                new_coeff = self_coeff * other_coeff * phase

                # Step 4: Append the new (coefficient, PauliString) tuple
                new_pauli_terms.append((new_coeff, product_pauli_string))

        # If there are no resulting terms, return a zero operator
        if not new_pauli_terms:
            size = self.get_size() if self.get_size() > 0 else other.get_size()
            return p([(0.0, 'I' * size)])

        # Create a new PauliStringLinear and simplify it to collect common terms
        return p(new_pauli_terms).simplify()

    def __rmatmul__(self, other:PauliString):
        """
        Overloading @ operator of two Pauli strings like multiply
        """
        new_combinations = []
        for c in self.combinations:
            new_combinations.append((c[0]*other.sign(c[1]), c[1]@other))
        return PauliStringLinear(new_combinations)

    def __rmul__(self, scalar: complex) -> 'PauliStringLinear':
        """
        Performs reverse scalar multiplication: scalar * self.
        """
        # Scalar multiplication is commutative, so we can just call our existing __mul__
        return self.__mul__(scalar)

    def __mul__(self, scalar: complex) -> 'PauliStringLinear':
        """
        Performs scalar multiplication: self * scalar.
        """
        # pylint: disable=import-outside-toplevel
        from paulie.common.pauli_string_factory import get_pauli_string as p

        # Check that we are multiplying by a number
        if not isinstance(scalar, (int, float, complex)):
            # Let Python know this operation is not implemented for other types
            return NotImplemented

        # Create a new list of terms where each coefficient is scaled by the number
        new_terms = [(coeff * scalar, pauli) for coeff, pauli in self]
        return p(new_terms)

    @property
    def h(self) -> 'PauliStringLinear':
        """
        Returns the Hermitian conjugate of this linear combination.
        This is found by taking the complex conjugate of all coefficients,
        as the Pauli matrices themselves are Hermitian.
        """
        # Use a list comprehension to create the new list of terms
        conjugated_combinations = [
            (np.conj(coeff), pauli) for coeff, pauli in self.combinations
        ]
        return PauliStringLinear(conjugated_combinations)

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


    def quadratic(self, basis: PauliString):
        """
        Computes the quadratic form Q_j = ∑_{S ∈ self} S ⊗ (L_j*S)
        where `self` represents the component C_k = ∑S, and `basis` is the
        linear symmetry L_j.
        """
        new_combinations = []
        # `self` is a linear combination of terms `(coeff, S_pauli)`
        # In our case, `coeff` is always 1.0 for a component.
        for coeff, s_pauli in self.combinations:
            # 1. Calculate the phase of the product Lj * S
            phase = basis.sign(s_pauli)

            # 2. Calculate the Pauli letters of the product Lj @ S
            product_letters = basis @ s_pauli

            # 3. Form the tensor product S ⊗ (product_letters)
            tensor_prod_str = str(s_pauli) + str(product_letters)

            # 4. The new coefficient is the original coefficient times the phase
            new_coeff = coeff * phase
            new_combinations.append((new_coeff, tensor_prod_str))
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

    def simplify(self) -> 'PauliStringLinear':
        """
        Combines terms with the same Pauli string by summing their coefficients.
        Removes terms with coefficients close to zero.
        This version assumes the PauliStringLinear object is directly iterable.
        """
        # pylint: disable=import-outside-toplevel
        from paulie.common.pauli_string_factory import get_pauli_string as p
        if not self.combinations:
            return self

        summed_coeffs = defaultdict(complex)
        for coeff, pauli in self.combinations:
            summed_coeffs[str(pauli)] += coeff

        simplified_list = [
            (coeff, pauli_str) for pauli_str, coeff in summed_coeffs.items()
            if abs(coeff) > 1e-12
        ]

        if not simplified_list:
            return p([(0.0, 'I' * self.get_size())])

        return p(simplified_list)

    def trace(self) -> complex:
        """
        Computes the trace of the operator represented by this linear combination.

        The trace is non-zero only if the Identity operator is present in the sum.
        Tr(self) = (coefficient of Identity) * 2^n.

        Returns:
            The complex value of the trace.
        """
        identity_coeff = 0.0  # Initialize the coefficient of the Identity operator

        # Loop through the terms to find the coefficient of the Identity string
        for coeff, pauli in self:
            # The PauliString class should have an `is_identity()` method
            if pauli.is_identity():
                identity_coeff = coeff
                break  # Found it, no need to look further

        # If there was no identity term, its coefficient is zero, so trace is zero.
        if identity_coeff == 0:
            return 0.0

        # Get the number of qubits, n
        num_qubits = self.get_size()

        # The trace is coeff * 2^n
        return identity_coeff * (2**num_qubits)

    # It should get the qubit count from its first Pauli string.
    def get_size(self) -> int:
        """
        Get the length of the Pauli Strings in this linear combination.
        """
        # Assuming the internal list of terms is accessible via iteration
        try:
            # Get the first term to determine the size
            _, first_pauli = next(iter(self))
            return len(first_pauli)
        except StopIteration:
            # Handle the case of an empty linear combination
            return 0

    def is_zero(self) -> bool:
        """
        Check if the linear combination is effectively zero.
        This can be done by checking if all coefficients are close to zero.
        
        Returns:
            True if the linear combination is zero, False otherwise.
        """
        return all(abs(coeff) < 1e-12 for coeff, _ in self)

    def norm(self) -> float:
        """
        Calculates the Frobenius norm of the coefficient vector.
        The norm is sqrt(∑ |cᵢ|²), where cᵢ are the coefficients of the
        Pauli strings in the linear combination.
        """
        # Sum the squared magnitudes of all coefficients
        sum_of_squares = sum(abs(coeff)**2 for coeff, _ in self.combinations)
        return np.sqrt(sum_of_squares)
