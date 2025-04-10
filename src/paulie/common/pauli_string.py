"""Implementation of a Pauli string on bitarray."""
from bitarray import bitarray
from bitarray.util import count_and
from paulie.common.pauli_string_parser import pauli_string_parser
from itertools import combinations

CODEC = {
    "I": bitarray([0, 0]), 
    "X": bitarray([1, 0]), 
    "Y": bitarray([1, 1]), 
    "Z": bitarray([0, 1]),
}

class PauliString:
    """Implementation of a Pauli string on bitarray."""

    def __init__(self, n: int = None, pauli_str: str = None, bits: bitarray = None):
        """Initialize a Pauli string.
        
        Args:
            n: Pauli string length
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

    def create_instance(self, n: int = None, pauli_str: str = None):
        """Create a Pauli string instance.
           Args:
                n: Pauli string length
                pauli_str: String representation of a Pauli string
           Returns the intensity of a Pauli string
        """
        return PauliString(n=n, pauli_str=pauli_str)
   
    def __str__(self) -> str:
        """Convert PauliString to readable string (e.g., "XYZI")."""
        return "".join(self.bits.decode(CODEC))
    
    def __eq__(self, other) -> bool:
        """Overloading the equality operator of two Pauli strings.
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits == other.bits

    def __lt__(self, other):
        """
        Overloading < operator of two Pauli strings
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits < other.bits

    def __le__(self, other):
        """
        Overloading <= operator of two Pauli strings
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits <= other.bits

    def __gt__(self, other):
        """
        Overloading > operator of two Pauli strings
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits > other.bits

    def __ge__(self, other):
        """
        Overloading >= operator of two Pauli strings
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits >= other.bits

    def __ne__(self, other):
        """
        Overloading != operator of two Pauli strings
        Args:
             other: Comparable Pauli string
        Returns the result of the comparison
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        if not isinstance(other, PauliString):
            return False
        return self.bits != other.bits
   
    def __hash__(self) -> int:
        """Make PauliString hashable so it can be used in sets"""
        return hash(str(self.bits))
    
    def __len__(self) -> int:
        """
        Pauli string length
        """
        return len(self.bits) // 2

    def __iter__(self):
        """
        Pauli String Iterator
        """
        self.nextpos = 0
        return self

    def __next__(self):
        """
        The value of the next position of the Pauli string
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = PauliString(bits=self.bits[2*self.nextpos:2*self.nextpos+2])
        self.nextpos += 1
        return value

    def __setitem__(self, position: int, pauli_string):
        """
        Set value to position
        """
        self.set_subsystem(position, pauli_string)

    def __getitem__(self, position: int):
        """
        Get values at position
        """
        return self.get_subsystem(position)

    def __copy__(self):
        """
        Pauli string copy operator
        """
        return PauliString(bits=self.bits)

    def copy(self):
         """Copy Pauli string"""
         return PauliString(bits=self.bits)
 
    def __add__(self, other): 
        """
        Pauli string addition operator
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)
        return self.tensor(other)
    
    def __or__(self, other)->bool:
        """
        Overloading | operator of two Pauli strings like commutes_with
        """
        return self.commutes_with(other)

    def __xor__(self, other):
        """
        Overloading ^ operator of two Pauli strings like adjoint_map
        """
        return self.adjoint_map(other)

    def __matmul__(self, other):
        """
        Overloading @ operator of two Pauli strings like multiply
        """
        return self.multiply(other)

    def commutes_with(self, other) -> bool:
        """
        Check if this Pauli string commutes with another
        Returns True if they commute, False if they anticommute
        """
        # Compute symplectic product mod 2
        # Paulis commute iff the symplectic product is 0
        if isinstance(other, str):
            other = PauliString(pauli_str=other)

        if len(self) != len(other):
            raise ValueError("Pauli arrays must be of equal length")
        return count_and(self.bits_even, other.bits_odd) % 2 == count_and(other.bits_even, self.bits_odd) % 2
   
        #return sum(self.bits[i] & other.bits[i + 1] for i in range(0, len(self.bits), 2)) % 2 == sum(self.bits[i + 1] & other.bits[i] for i in range(0, len(self.bits), 2)) % 2
      
    def get_subsystem(self, start: int, length: int = 1) -> "PauliString":
        """Get a subsystem of this Pauli string"""
        return PauliString(bits=self.bits[2*start:2*start+2*length])

    def set_subsystem(self, position: int, pauli_string):
        """
        Set subsystem value
        """
        if isinstance(pauli_string, str):
            pauli_string = PauliString(pauli_str=pauli_string)

        for i in range(0, len(pauli_string)):
            self.bits[2*position + 2*i] = pauli_string.bits[2*i]
            self.bits[2*position + 2*i + 1] = pauli_string.bits[2*i + 1]
            self.bits_even[position  + i] = pauli_string.bits_even[i]
            self.bits_odd[position + i] = pauli_string.bits_odd[i]

    def is_identity(self) -> bool:
        """Check if this Pauli string is the identity"""
        return bitarray(len(self.bits)) == self.bits
    
    def tensor(self, other: "PauliString") -> "PauliString":
        """Tensor product of this Pauli string with another"""
        new_bits = bitarray(len(self.bits) + len(other.bits))
        for i in range(len(new_bits)):
            s = self.bits if i < len(self.bits) else other.bits
            j = i if i < len(self.bits) else i - len(self.bits)
            new_bits[i] = s[j]

        return PauliString(bits=new_bits)

    def multiply(self, other) -> "PauliString":
        """
        Proportional multiplication operator of two Pauli strings
        Returns a PauliString proportional to the multiplication 
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)

        if len(self.bits) != len(other.bits):
            raise ValueError("Pauli arrays must have the same length")
        return PauliString(bits = self.bits ^ other.bits)  # Bitwise XOR is equivalent to mod-2 addition
    
    def adjoint_map(self, other) -> "PauliString":
        """
        Compute the adjoint map ad_A(B) = [A,B]
        Returns None if the commutator is zero (i.e., if A and B commute)
        Otherwise returns a PauliString proportional to the commutator
        """
        if isinstance(other, str):
            other = PauliString(pauli_str=other)

        if self.commutes_with(other):
            return None
        
        # For Pauli strings, if they anticommute, [A,B] = 2AB
        # In the context of generating a Lie algebra, we only care about
        # the result up to a constant factor
        
        # For anticommuting Paulis, the product gives a new Pauli
        # XOR of the bit vectors gives the non-phase part of the product
        if len(self.bits) != len(other.bits):
            raise ValueError("Pauli arrays must have the same length")
        return PauliString(bits = self.bits ^ other.bits)  # Bitwise XOR is equivalent to mod-2 addition

    def inc(self):
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

    def expand(self, n: int):
        """
        Increasing the size of the Pauli string
        Args:
            n (int): New Pauli string length
        Returns the Pauli string extension
        """
        return self + PauliString(n = n - len(self))

    def gen_all_pauli_strings(self):
        """
        Generate a list of Pauli strings that commute with this string
        Args:
          Yields the generated Pauli string

        """
        n = len(self)
        pauli_string = PauliString(n=n)

        last = PauliString(bits = bitarray([1] * (2 * n)))

        while pauli_string !=last:
            yield pauli_string.copy()
            pauli_string.inc()
        yield pauli_string.copy()

    def get_commutants(self, generators = None):
        """
        Get a list of Pauli strings that commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if self|g]

    def get_anti_commutants(self, generators = None):
        """
        Get a list of Pauli strings that no-commute with this string
        Args:
            generators: Collection of Pauli strings on which commutant is searched
                        If not specified, then the search area is all Pauli strings of the same size

        """
        if generators is None:
            generators = self.gen_all_pauli_strings()

        return [g for g in generators if not self|g]

    def get_nested(self, generators = None):
        """
        Get nested of Pauli string
        Args:
            generators: Collection of Pauli strings on which nested is searched
                        If not specified, then the search area is all Pauli strings of the same size
        """

        generators = self.get_anti_commutants(generators = generators)
        return [
            (a, b)
            for a, b in combinations(generators, 2)
            if a^b == self
        ]
