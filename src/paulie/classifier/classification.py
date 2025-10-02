"""
Canonical graph classification
"""
import enum
from typing import Generator

from paulie.common.pauli_string_bitarray import PauliString


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
    def __init__(self, legs:list[list[PauliString]], dependents:list[PauliString]) -> None:
        """
        Constructor
        """
        self.legs = legs # center is zero leg
        self.dependents = dependents

    def is_empty(self) -> bool:
        """
        Checking for emptiness of a graph
        """
        return len(self.legs) == 0

    def is_empty_legs(self) -> bool:
        """
        Checking for missing legs in a graph
        """
        return len(self.legs) == 1

    def get_vertices(self) -> list[PauliString]:
        """
        Get a list of graph vertices
        """
        return [v for leg in self.legs for v in leg ]

    def get_dependents(self) -> list[PauliString]:
        """
        Get a list of dependent Pauli strings
        """
        return self.dependents

    def get_legs(self) -> list[list[PauliString]]:
        """
        Get legs
        """
        return self.legs


    def counts(self) -> tuple[int,int,int]:
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
        if long_vertices == 0 and two_legs == 0 and one_legs == 2:
            long_vertices = 1
            one_legs = 2
        return one_legs, two_legs, long_vertices

    def get_properties(self) -> tuple[TypeGraph,int,int,int]:
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

    def get_type(self) -> TypeGraph:
        """
        Get graph type
        """
        type_graph = self.get_properties()[0]
        return type_graph

    def get_algebra_properties(self) -> tuple[TypeAlgebra,int,int]:
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

    def check_algebra_properties(self, type_algebra:TypeAlgebra = None,
                                 nc:int = None, size:int = None) -> bool:
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

    def gen_independent_pair(self
        )-> Generator[list[list[PauliString]], None, None]:
        """
        Generate independent pairs
        """
        for i in range(len(self.legs) - 1, 0, -1):
            leg = self.legs[i]
            for j in range(len(leg) - 1, -1, -1):
                v = leg[j]
                for k in range(j - 2, -1, -1):
                    w = leg[k]
                    yield [v, w]
                for m in range(i - 1, -1, -1):
                    for w in leg[m]:
                        yield [v, w]

    def gen_pq(self
        )-> Generator[list[dict], None, None]:
        """
        Generate pq 
        """
        for i in range(len(self.legs) - 1, 0, -1):
            leg = self.legs[i]
            for j in range(len(leg) - 1, -1, -1):
                v = leg[j]
                for k in range(j - 1, -1, -1):
                    w = leg[k]

                    yield { "v": v,
                            "w": w,
                            "vw": v@w,
                            "neighbour": j - k == 1
                          }
                for m in range(i - 1, -1, -1):
                    for w in self.legs[m]:
                        yield { "v": v,
                                "w": w,
                                "vw": v@w,
                                "neighbour": m == 0 and j == 0
                              }

    def _inc_vertices_generator(self, i, init_vertices, vertix_generators, vertices):
        """incriminate the generator"""
        if i == len(vertix_generators):
            return False
        try:
            vw = next(vertix_generators[i])
            if init_vertices[i] != vw["v"]:
                vertices[i] = init_vertices[i]@vw["vw"]
            else:
                if vw["neighbour"]:
                    vertices[i] = init_vertices[i]@vw["w"]
                else:
                    return self._inc_vertices_generator(i, init_vertices,
                           vertix_generators, vertices)
            return True
        except StopIteration:
            for k in range(0, i+1):
                vertix_generators[k] = self.gen_pq()
            return self._inc_vertices_generator(i + 1, init_vertices, vertix_generators, vertices)

    def gen_generators(self
        )-> Generator[list[list[PauliString]], None, None]:
        """
        Generate generators
        """
        vertices = self.get_vertices()
        yield vertices
        init_vertices = vertices.copy()
        vertix_generators = [self.gen_pq() for i in range(0, len(vertices))]

        while True:
            if not self._inc_vertices_generator(0, init_vertices, vertix_generators, vertices):
                break
            if len(list(set(vertices))) == len(vertices):
                yield vertices



class Classification:
    """
    Algebra classification class
    """
    def __init__(self) -> None:
        """
        Constructor
        """
        self.morphs: set[Morph] = set()

    def add(self, morph: Morph) -> None:
        """
        Add canonical form
        """
        self.morphs.add(morph)

    def get_morphs(self) -> Morph:
        """
        Get canonical form
        """
        return self.morphs

    def get_algebra(self) -> str:
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

    def contains_algebra(self, algebra:str) -> bool:
        """
        Algebra inclusion check
        """
        _algebra = self.get_algebra()
        algebra.replace(" ", "")
        return _algebra.find(algebra) > -1

    def _parse_algebra(self, algebra:str) -> bool:
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

    def is_algebra(self, algebra:str) -> bool:
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
            if a != _algebras[i] and a != self.get_isomorphism(_algebras[i]):
                return False
        return True

    def get_subalgebras(self, algebra:str=None) -> list[str]:
        """
        Get subalgebras
        """
        if algebra is None:
            algebra = self.get_algebra()
        else:
            algebra.replace(" ", "")
        return algebra.split("+")

    def get_vertices(self) -> list[PauliString]:
        """
        Get a list of independent strings of Pauli algebra
        """
        return [v for morph in self.morphs for v in morph.get_vertices() ]

    def get_dependents(self) -> list[PauliString]:
        """
        Get a list of dependent strings of Pauli algebra
        """
        return [v for morph in self.morphs for v in morph.get_dependents()]

    def get_isomorphisms(self) -> dict[str, str]:
        """
        Dictionary of isomorphisms
        """
        return {"2*so(2)":"2*su(2)",
                "so(3)":"su(2)",
                "so(4)":"2*su(2)"
               }

    def get_isomorphism(self, algebra:str)->str:
        """
        Get algebra isomorphism
        """
        n = 1
        core_algebra = ""
        if algebra not in self.get_isomorphisms():
            if '*' in algebra:
                algebras = algebra.split("*")
                n = algebras[0]
                core_algebra = algebras[1]
                if core_algebra not in self.get_isomorphisms():
                    return None
                isomorph_core_algebra = self.get_isomorphisms()[core_algebra]
                isomorph_n = 1
                if "*" in isomorph_core_algebra:
                    isomorph_core_algebras = isomorph_core_algebra.split("*")
                    isomorph_n = int(isomorph_core_algebras[0])
                    isomorph_core_algebra = isomorph_core_algebras[1]
                isomorph_n *= n
                return (f"{isomorph_core_algebra}" if isomorph_n == 1
                       else f"{isomorph_n}*{isomorph_core_algebra}")
            return "None"
        return self.get_isomorphisms()[algebra]

    def get_dla_dim(self) -> int:
        """
        Get the dimension of the classified dynamical Lie algebra
        """
        dim_su = lambda n: n**2-1
        dim_so = lambda n: n*(n-1)/2
        dim_sp = lambda n: n*(2*n+1)
        dim = 0
        for morph in self.morphs:
            result = morph.get_algebra_properties()
            type_algebra = result[0]
            n = result[2]
            if type_algebra == TypeAlgebra.SU:
                dim+= dim_su(n)
            if type_algebra == TypeAlgebra.SP:
                dim+= dim_sp(n)
            if type_algebra == TypeAlgebra.SO:
                dim+= dim_so(n)
        return dim

    def _inc_morph_generator(self, ms, morphs, morph_generators, current_morph_generators):
        """incriminate the generator"""
        if ms == len(morph_generators):
            return False
        try:
            current_morph_generators[ms] = next(morph_generators[ms])
            return True
        except StopIteration:
            for k in range(0, ms+1):
                morph_generators[k] = morphs[k].gen_generators()

            #morph_generators[ms] = morphs[ms].gen_generators()
            current_morph_generators[ms] = next(morph_generators[ms])
            return self._inc_morph_generator(ms + 1, morphs,
                   morph_generators, current_morph_generators)

    def gen_generators(self) -> Generator[list[list[PauliString]], None, None]:
        """
        Get a list of independent strings of Pauli algebra
        """
        generators = []
        current_morph_generators = []
        morph_generators = []
        morphs = list(self.morphs)
        for ms in range(0, len(morphs)):
            morph_generators.append(morphs[ms].gen_generators())
            try:
                current_morph_generators.append(next(morph_generators[ms]))
            except StopIteration:
                return
            generators += current_morph_generators[ms]

        yield generators

        while True:
            if not self._inc_morph_generator(0, morphs, morph_generators, current_morph_generators):
                break
            generators = []
            for ms in range(0, len(morphs)):
                generators += current_morph_generators[ms]
            yield generators
