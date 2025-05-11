"""
Canonical graph classification
"""
import enum

class ClassificatonException(Exception):
    """
    Classification exclusion
    """

class TypeGraph(enum.Enum):
    """
    Type of canonical graph
    """
    A = 0
    B1 = 1
    B2 = 2
    B3 = 3
    NONE = 4

class TypeAlgebra(enum.Enum):
    """
    Lie algebra type
    """
    U = 0
    SU = 1
    SP = 2
    SO = 3


class Morph:
    """
    Class of canonical form of a graph
    """
    def __init__(self, legs, dependents):
        """
        Constructor
        """
        self.legs = legs # center is zero leg
        self.dependents = dependents

    def is_empty(self):
        """
        Checking for emptiness of a graph
        """
        return len(self.legs) == 0

    def is_empty_legs(self):
        """
        Checking for missing legs in a graph
        """
        return len(self.legs) == 1

    def get_vertices(self):
        """
        Get a list of graph vertices
        """
        return [v for leg in self.legs for v in leg ]

    def get_dependents(self):
        """
        Get a list of dependent Pauli strings
        """
        return self.dependents

    def counts(self):
        """
        Get number of leg types
        """
        one_legs = 0
        two_legs = 0
        long_vertices = 0
        for i, leg in enumerate(self.legs):
            if i == 0:
                continue
            if len(leg) == 1:
                one_legs += 1
            if len(leg) == 2:
                two_legs += 1
            if len(leg) > 2:
                if long_vertices > 0:
                    raise ClassificatonException("Graph of non-canonical type")
                long_vertices += len(leg)
        if long_vertices == 0 and two_legs == 1:
            two_legs = 0
            long_vertices = 2
        if long_vertices > 0 and two_legs == 0:
            long_vertices += 1
        if long_vertices == 0 and two_legs == 0 and one_legs == 1:
            long_vertices = 1
        return one_legs, two_legs, long_vertices

    def get_properties(self):
        """
        Get graph properties
        """
        if self.is_empty():
            raise ClassificatonException("Graph of non-canonical type")
        if self.is_empty_legs():
            return  TypeGraph.NONE, 0, 0, 0
        one_legs, two_legs, long_vertices = self.counts()
        if two_legs == 0:
            return  TypeGraph.A, one_legs, two_legs, long_vertices
        if long_vertices == 0:
            return  TypeGraph.B1, one_legs, two_legs, long_vertices
        if long_vertices == 3:
            return  TypeGraph.B3, one_legs, two_legs, long_vertices
        if long_vertices == 4:
            return  TypeGraph.B2, one_legs, two_legs, long_vertices
        raise ClassificatonException("Graph of non-canonical type")

    def get_type(self):
        """
        Get graph type
        """
        type_graph = self.get_properties()[0]
        return type_graph

    def get_algebra_properties(self):
        """
        Get properties of algebra
        """
        type_graph, one_legs, two_legs, long_vertices = self.get_properties()
        if type_graph == TypeGraph.NONE:
            return TypeAlgebra.U, 1, 1
        if type_graph == TypeGraph.A:
            return TypeAlgebra.SO, one_legs, long_vertices + 2
        if type_graph == TypeGraph.B1:
            return TypeAlgebra.SP, one_legs, 2**two_legs
        if type_graph == TypeGraph.B2:
            return TypeAlgebra.SO, one_legs, 2**(two_legs + 3)
        if type_graph == TypeGraph.B3:
            return TypeAlgebra.SU, one_legs, 2**(two_legs + 2)
        return None, None, None

    def check_algebra_properties(self, type_algebra = None, nc = None, size = None):
        """
        Check properties of algebra
        """
        _type_algebra, _nc, _size = self.get_algebra_properties()
        if type_algebra is None and nc is None and size is None:
            return True
        if type_algebra is not None:
            if type_algebra != _type_algebra:
                return False
        if nc is not None and size is not None:
            return size == _size and nc == _nc
        if size is not None:
            return size == _size
        if nc is not None:
            return nc == _nc
        return False

class Classification:
    """
    Algebra classification class
    """
    def __init__(self):
        """
        Constructor
        """
        self.morphs = set()

    def add(self, morph):
        """
        Add canonical form
        """
        self.morphs.add(morph)

    def get_morphs(self):
        """
        Get canonical form
        """
        return self.morphs

    def get_algebra(self):
        """
        Get algebra
        """
        algebras = {}
        for morph in self.morphs:
            type_algebra, nc, size = morph.get_algebra_properties()
            algebra = ""
            if type_algebra == TypeAlgebra.U:
                algebra = f"u({size})"
            if type_algebra == TypeAlgebra.SO:
                algebra = f"so({size})"
            if type_algebra == TypeAlgebra.SP:
                algebra = f"sp({size})"
            if type_algebra == TypeAlgebra.SU:
                algebra = f"su({size})"
            if algebra in algebras:
                algebras[algebra] += nc if nc == 1 else 2**(nc-1)
            else:
                algebras[algebra] = nc if nc == 1 else 2**(nc-1)
        return "+".join([key if v == 1 else str(v) + "*" + key for key, v in algebras.items()])

    def contains_algebra(self, algebra):
        """
        Algebra inclusion check
        """
        _algebra = self.get_algebra()
        algebra.replace(" ", "")
        return _algebra.find(algebra) > -1

    def _parse_algebra(self, algebra):
        """
        Parse algebra
        """
        algebra = algebra.replace(" ", "")
        algebras = algebra.split("+")

        algs = {}
        for alg in algebras:
            name = alg
            q = 1
            if "*" in alg:
                a = name.split("*")
                name = a[1]
                q = int(a[0])
            if name in algs:
                q += algs[name]
            algs[name] = q
        return "+".join([key if v == 1 else str(v) +
               "*" + key for key, v in algs.items()]).split("+")

    def is_algebra(self, algebra):
        """
        Checking for compliance with a given algebra
        """
        _algebra = self.get_algebra()
        algebras = self._parse_algebra(algebra)
        _algebras = _algebra.split("+")
        if len(algebras) != len(_algebras):
            return False
        algebras.sort()
        _algebras.sort()
        for i, a in enumerate(algebras):
            if a != _algebras[i]:
                return False
        return True

    def get_subalgebras(self, algebra=None):
        """
        Get subalgebras
        """
        if algebra is None:
            algebra = self.get_algebra()
        else:
            algebra.replace(" ", "")
        return algebra.split("+")

    def get_vertices(self):
        """
        Get a list of independent strings of Pauli algebra
        """
        return [v for morph in self.morphs for v in morph.get_vertices() ]

    def get_dependents(self):
        """
        Get a list of dependent strings of Pauli algebra
        """
        return [v for morph in self.morphs for v in morph.get_dependents()]
