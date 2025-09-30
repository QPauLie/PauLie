"""
Factory for constructing a canonical graph
"""
import traceback
from typing import Generator, Self
from paulie.helpers.printing import Debug
from paulie.classifier.classification import Morph
from paulie.common.pauli_string_bitarray import PauliString

class AppendedException(Exception):
    """
    The vertix is appended
    """

class CheckAppendedException(Exception):
    """
    The vertix is appended
    """

class DebugException(Exception):
    """
    Debug exception
    """

class DependentException(Exception):
    """
    Vertex dependent
    """

class NotConnectedException(Exception):
    """
    No connection with cononic graph
    """

class RaiseException(Exception):
    """
    External exception
    """
class MorphFactoryException(Exception):
    """
    Morph factory exception
    """

class MorphFactory(Debug):
    """
    Factory for constructing a canonical graph
    """
    def __init__(self, debug:bool = False) -> None:
        """
        Constructor
        """
        super().__init__(debug)
        self.legs = [] # center is zero leg
        self.lighting = None
        self.delayed_vertices = []
        self.debug_lighting = None
        self.debug_break = False
        self.dependents = []
        self.is_check = False

    def set_debug(self, debug:bool) -> None:
        """
        Set debug flag
        """
        self.debug = debug

    def set_lighting(self, lighting:PauliString) -> None:
        """
        Set lighting
        """
        self.lighting = lighting

    def get_lighting(self) -> PauliString:
        """
        Get lighting
        """
        return self.lighting

    def get_morph(self) -> Morph:
        """
        Get canonical graph form
        """
        return Morph(self.legs, self.dependents)


    def lit(self, lighting:PauliString, vertix:PauliString) -> PauliString:
        """
        Lit vertix
        """
        lighting = lighting@vertix
        #lighting = lighting^vertix
        if self.is_included(lighting):
            raise DependentException()
        return lighting

    def get_lits(self, lighting:PauliString,
                 vertices:list[PauliString]=None) ->list[PauliString]:
        """
        Return lit vertices (connected to the selected vertex).
        """
        if vertices is None:
            vertices = self.get_vertices()
        return [v for v in vertices if v != lighting and not lighting|v]


    def is_empty(self) -> bool:
        """
        Checking for emptiness
        """
        return len(self.legs) == 0

    def is_empty_legs(self) -> bool:
        """
        Checking for missing legs
        """
        return len(self.legs) < 3

    def _find_in_leg(self, leg:list[PauliString], v:PauliString) -> int:
        """
        Find vertix in leg
        """
        try:
            index = leg.index(v)
        except ValueError:
            index = -1
        return index


    def find(self, v:PauliString) -> tuple[int,int]:
        """
        Find vertix
        """
        for i, leg in enumerate(self.legs):
            index = self._find_in_leg(leg, v)
            if index > -1:
                return i, index
        return -1, -1

    def is_included(self, v:PauliString) -> bool:
        """
        Checking a vertex for inclusion in the graph
        """
        leg_index = self.find(v)[0]
        return leg_index > -1

    def get_vertices(self) -> list[PauliString]:
        """
        Get graph vertices
        """
        return [v for leg in self.legs for v in leg]

    def get_center(self) -> PauliString|None:
        """
        Get center
        """
        if self.is_empty():
            return None
        return self.legs[0][0]

    def set_center(self, v: PauliString) -> None:
        """
        Set center
        """
        if self.is_empty() is False:
            raise MorphFactoryException("Center is setted")
        self.legs.append([v])

    def get_long_leg(self) -> list[PauliString]:
        """
        Get long leg
        """
        if self.is_empty_legs():
            raise MorphFactoryException("No legs")
        return self.legs[len(self.legs) - 1]

    def get_one_vertix(self) -> PauliString:
        """
        Get one vertix in leg
        """
        if self.is_empty_legs():
            raise MorphFactoryException("No legs")
        return self.legs[1][0]

    def _gen_one_legs(self) -> Generator[list[list[PauliString]], None, None]:
        """
        Generate vertices included in single legs
        """
        if self.is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(1, len(self.legs)):
            if len(self.legs[i]) == 1:
                yield self.legs[i]
            else:
                break

    def get_one_vertices(self)->list[PauliString]:
        """
        Get vertices included in single legs
        """
        vertices = []
        for leg in self._gen_one_legs():
            vertices.append(leg[0].copy())
        return vertices

    def get_pq(self, lighting:PauliString) -> tuple[PauliString|None,PauliString|None]:
        """
        Get pq 
        """
        one_verices = self.get_one_vertices()
        lits = self.get_lits(lighting, one_verices)
        p = None
        q = None
        for v in one_verices:
            if v in lits:
                p = v
            else:
                q = v
            if p is not None and q is not None:
                return p@q, p
        return None, None

    def _gen_two_legs(self) -> Generator[list[list[PauliString]], None, None]:
        """
        Generate vertices included in two legs
        """
        if self.is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(1, len(self.legs)):
            if len(self.legs[i]) == 2:
                yield self.legs[i]
            else:
                if len(self.legs[i]) > 2:
                    break

    def get_two_legs(self) -> list[tuple[PauliString, PauliString]]:
        """
        Get vertices included in two legs
        """
        return [(leg[0].copy(), leg[1].copy()) for leg in self._gen_two_legs()]

    def get_count_two_legs(self) -> int:
        """
        Get the number of legs of length two
        """
        return len(self.get_two_legs())

    def is_two_leg(self) -> bool:
        """
        Checking the leg for length two
        """
        count_two_legs = self.get_count_two_legs()
        if count_two_legs == 0:
            return False

        long_leg = self.get_long_leg()
        if len(long_leg) != 2:
            return True
        return count_two_legs > 1

    def _gen_long_legs(self) -> list[tuple[PauliString, PauliString]]:
        """
        Generate long leg vertices
        """
        if self.is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(len(self.legs)-1, 1, -1):
            if len(self.legs) > 2:
                yield self.legs[i]
            else:
                break

    def get_long_legs(self) -> list[list[PauliString]]:
        """
        Get long leg vertices
        """
        return [leg.copy() for leg in self._gen_long_legs()]

    def append(self, v:PauliString, lit:PauliString) -> None:
        """
        Append vertix to graph
        """
        if self.is_check:
            raise CheckAppendedException()
        leg_index, vertix_index = self.find(lit)
        if leg_index == -1:
            raise MorphFactoryException("No vertix")
        if leg_index == 0:
            self.legs.insert(1, [v])
            return
        if vertix_index !=  len(self.legs[leg_index]) - 1:
            raise MorphFactoryException("The vertix is not the last")
        leg = self.legs[leg_index].copy()
        del self.legs[leg_index]
        leg.append(v)
        if len(leg) >= len(self.legs[len(self.legs) - 1]):
            self.legs.append(leg)
            return
        for i in range(len(self.legs) - 1, 0, -1):
            if len(self.legs[i]) <= len(leg):
                self.legs.insert(i+1, leg)
                return
        raise MorphFactoryException("Can't append")

    def append_to_two_center(self, lighting:PauliString) -> None:
        """
        Append vertix to two wertices graph
        """
        center = self.get_center()
        if len(self.legs) == 1:
            self.append(lighting, center)
            return
        vertices = self.get_vertices()
        lits = self.get_lits(lighting, vertices)

        if len(lits) == 1:
            if center in lits:
                self.append(lighting, center)
                return
            else:
                lighting = self.lit(lighting, lits[0])
                lighting = self.lit(lighting, center)
                self.append(lighting, center)
                return
        if len(lits) == 2:
            lighting = self.lit(lighting, center)
            self.append(lighting, center)
            return
        raise NotConnectedException()

    def check_dependency_one_leg(self, lighting:PauliString) -> None:
        """
        Dependency check when attaching a vertex to the center of the graph
        """
        ones = self.get_one_vertices()
        vertices = set(self.get_vertices())
        for one in ones:
            pq = one @ lighting
            for v in vertices:
                if v == one:
                    continue
                n_v = pq @ v
                if n_v in vertices or n_v == lighting:
                    raise DependentException()

    def append_to_center(self, lighting:PauliString) -> None:
        """
        Joining a vertex to the center of the graph
        """
        self.check_dependency_one_leg(lighting)
        center = self.get_center()
        self.append(lighting, center)

    def remove(self, v:PauliString) -> None:
        """
        Removing a graph vertex
        """
        leg_index, vertix_index = self.find(v)
        if leg_index == -1:
            raise MorphFactoryException("No vertix")
        if leg_index == 0:
            raise MorphFactoryException("Can't delete the center")

        leg = [self.legs[leg_index][i] for i in range(0, vertix_index)]
        del self.legs[leg_index]
        if len(leg) == 0:
            return
        if len(leg) == 1:
            self.legs.insert(1, leg)
            return
        if len(leg) >= len(self.legs[len(self.legs) - 1]):
            self.legs.append(leg)
            return
        for i in range(len(self.legs) - 1, 1, -1):
            if len(self.legs[i]) <= len(leg):
                self.legs.insert(i+1, leg)
                return
        raise MorphFactoryException("Can't remove")

    def replace(self, v:PauliString, v_new:PauliString) -> None:
        """
        Replacing a graph vertex with an equivalent one
        """
        leg_index, vertix_index = self.find(v)
        if leg_index == -1:
            raise MorphFactoryException("No vertix")
        self.legs[leg_index][vertix_index] = v_new

    def print_state(self, lighting:PauliString = None) -> None:
        """
        Debug output of graph state
        """
        if self.debug is False:
            return
        self.print_title(f"state graph {len(self.get_vertices())}")
        if lighting is None:
            for i, leg in enumerate(self.legs):
                if i == 0:
                    self.print_title("center")
                    self.print_vertix(self.legs[i][0])
                    continue
                self.print_vertices(leg)
            return
        self.print_vertix(lighting, "print state")
        lits = self.get_lits(lighting)
        for i, leg in enumerate(self.legs):
            if i == 0:
                self.print_title("center")
                clit = ""
                if self.legs[i][0] in lits:
                    clit = "*"
                self.print_vertix(self.legs[i][0], clit)
                continue
            self.print_lit_vertices(leg, lits)

    def get_lit_indexes(self, vertices:list[PauliString], lits:list[PauliString]) -> list[int]:
        """
        Get the indices of the lited vertices in lits 
        """
        indexes = []
        for i, v in enumerate(vertices):
            if v in lits:
                indexes.append(i)
        return indexes

    def _append_three_graph(self) -> Self:
        """
        Step I. Construct a graph of three vertices
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Check and build graph with three vertices")
        if self.is_empty():
            self.set_center(lighting)
            raise AppendedException
        if self.is_included(lighting):
            raise DependentException
        if self.is_empty_legs():
            self.append_to_two_center(lighting)
            raise AppendedException
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _append_one_legs_in_different_state(self) -> Self:
        """
        Step II. Legs of length 1 in different initial lit states.
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Build graph with different one leg state")
        pq, p = self.get_pq(lighting)
        if pq is not None:
            lits = self.get_lits(lighting)
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    if self.is_included(v):
                        raise DependentException()
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    self.replace(lit, v)
            self.append(lighting, p)
            long_leg = self.get_long_leg()
            if len(long_leg) > 4:
                for i in range(4, len(long_leg)):
                    self.append_delayed(long_leg[i])
                self.remove(long_leg[4])
            raise AppendedException
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _append_fast(self) -> Self:
        """
        Quickly obvious connection of lightning to a graph
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Append fast")
        center = self.get_center()
        two_legs = self.get_two_legs()
        long_leg = self.get_long_leg()
        if len(long_leg) == 2:
            del two_legs[len(two_legs) - 1]
        if len(two_legs) == 0:
            lits = self.get_lits(lighting)
            if len(lits) == 1:
                if center in lits:
                    self.append_to_center(lighting)
                    raise AppendedException
                long_leg = self.get_long_leg()
                if long_leg[len(long_leg) - 1] in lits:
                    self.append(lighting, long_leg[len(long_leg) - 1])
                    raise AppendedException
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _lit_only_long_leg(self) -> Self:
        """
        Step III. Lit only the long leg
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Lit only long leg")
        omega = self.get_one_vertix()
        center = self.get_center()
        center_lits = self.get_lits(lighting, [center])
        lits = self.get_lits(lighting, [omega])
        if omega in lits:
            if center not in center_lits:
                lighting = self.lit(lighting, omega)
            lighting = self.lit(lighting, center)
        two_legs = self.get_two_legs()
        long_leg = self.get_long_leg()
        if len(long_leg) == 2:
            del two_legs[len(two_legs) - 1]
        elif len(long_leg) == 1:
            self.print_state(lighting)
            self.set_lighting(lighting)
            return self
        if len(two_legs) == 0:
            self.print_state(lighting)
            self.set_lighting(lighting)
            return self
        long_lits = self.get_lits(lighting, long_leg)
        if len(long_lits) == 0:
            # find lited two leg
            center_lits = self.get_lits(lighting, [center])
            if center in center_lits:
                lighting = self.lit(lighting, center)
                lighting = self.lit(lighting, long_leg[0])
                lighting = self.lit(lighting, omega)
                lighting = self.lit(lighting, center)
            else:
                for leg in two_legs:
                    lits = self.get_lits(lighting, leg)
                    v0 = leg[0]
                    v1 = leg[1]
                    if v1 in lits and v0 not in lits:
                        lighting = self.lit(lighting, v1)
                        lits.append(v0)
                    if v0 in lits:
                        lighting = self.lit(lighting, v0)
                        lighting = self.lit(lighting, center)
                        lighting = self.lit(lighting, long_leg[0])
                        lighting = self.lit(lighting, omega)
                        lighting = self.lit(lighting, center)
                        break
        # lit second vertix on long leg
        long_lits = self.get_lits(lighting, long_leg)
        lit_indexes = self.get_lit_indexes(long_leg, long_lits)
        if 1 not in lit_indexes:
            if 0 in lit_indexes:
                lighting = self.lit(lighting, long_leg[0])
            else:
                if len(lit_indexes) == 0:
                    raise NotConnectedException()
                first_lit = lit_indexes[0]
                for i in range(first_lit, 1, -1):
                    lighting = self.lit(lighting, long_leg[i])
        long_v0 = long_leg[0]
        long_v1 = long_leg[1]
        for i,leg in enumerate(two_legs):
            lits = self.get_lits(lighting, leg)
            v0 = leg[0]
            v1 = leg[1]
            if v0 not in lits and v1 not in lits:
                continue
            if v0 in lits and v1 not in lits:
                lighting = self.lit(lighting, v0)
                lits.append(v1)
            elif v0 not in lits and v1 in lits:
                lighting = self.lit(lighting, v1)
                lits.append(v0)
            if v0 in lits and v1 in lits:
                center_lits = self.get_lits(lighting, [center])
                if center in center_lits:
                    lighting = self.lit(lighting, center)
                    #omega is lited
                    lighting = self.lit(lighting, v1)
                    lighting = self.lit(lighting, v0)
                    lighting = self.lit(lighting, omega)
                    lighting = self.lit(lighting, center)
                else:
                    long_lits = self.get_lits(lighting, [long_leg[0]])
                    if len(long_lits) == 0:
                        lighting = self.lit(lighting, long_v1)
                    lighting = self.lit(lighting, long_v0)
                    lighting = self.lit(lighting, center)
                    lighting = self.lit(lighting, omega)
                    lighting = self.lit(lighting, v1)
                    lighting = self.lit(lighting, v0)
                    lighting = self.lit(lighting, center)
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _lit_center(self) -> Self:
        """
        Lit center
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Liting center")
        center = self.get_center()
        center_lits = self.get_lits(lighting, [center])
        if center not in center_lits:
            long_leg = self.get_long_leg()
            long_lits = self.get_lits(lighting, long_leg)
            lit_indexes = self.get_lit_indexes(long_leg, long_lits)
            first_lit = lit_indexes[0]
            for i in range(first_lit, -1, -1):
                lighting = self.lit(lighting, long_leg[i])
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self


    def _reduce_long_leg_more_than_one_lits(self) -> Self:
        """
        Step IV. Reducing the long leg lits to standard configurations
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Reduce long leg lits")
        long_leg = self.get_long_leg()
        #Offset the first vertex to the end of the linear leg. Since the leg is finite,
        #we will always reach the highlighting of one vertex,
        #it will be either the last one or the penultimate one
        while True:
            lits = self.get_lits(lighting, long_leg)
            if len(lits) == 0:
                self.append_to_center(lighting)
                raise AppendedException()

            if len(lits) == 2:
                lit_indexes = self.get_lit_indexes(long_leg, lits)
                if lit_indexes[0] == 0 and lit_indexes[1] == len(long_leg) - 1:
                    break
            if len(lits) == 1:
                if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                    break
                if long_leg[0] != lits[0]:
                    lit_indexes = self.get_lit_indexes(long_leg, lits)
                    if lit_indexes[0] < len(long_leg) - 1:
                        for i in range(lit_indexes[0] + 1, len(long_leg)):
                            self.append_delayed(long_leg[i])
                        self.remove(long_leg[lit_indexes[0] + 1])
                    break
            lit_indexes = self.get_lit_indexes(long_leg, lits)
            first = lit_indexes[0]
            second = lit_indexes[1]
            if first > 0 and first + 1 != second:
                for i in range(second, first, -1): ## maybe + 1
                    lighting = self.lit(lighting, long_leg[i])
            else:
                lighting = self.lit(lighting, long_leg[second])
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _append_long_leg_first_and_center_lit(self) -> Self:
        """
        Step V. Append long leg with first lit and center
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Append long leg with first lit and center")
        omega = self.get_one_vertix()
        center = self.get_center()
        lits = self.get_lits(lighting, [center, omega])
        is_center_lit = center in lits
        # lit only long leg
        long_leg = self.get_long_leg()
        lits = self.get_lits(lighting, long_leg)
        if is_center_lit and len(lits) == 0:
            #lit only center
            self.print_state(lighting)
            self.append_to_center(lighting)
            raise AppendedException
        lit_indexes = self.get_lit_indexes(long_leg, lits)
        # only long leg and center are lited
        if len(lit_indexes) == 1 and 0 in lit_indexes:
            # if leg less than 4 or not legs long 2, we can connect to the end of long leg
            is_can_connect_to_end = True
            if self.is_two_leg() and len(long_leg) > 3:
                is_can_connect_to_end = False
            if is_can_connect_to_end:
                for v in long_leg:
                    lighting = self.lit(lighting, v)
                self.append(lighting, long_leg[len(long_leg)-1])
                raise AppendedException
            two_legs = self.get_two_legs()
            two_leg = two_legs[0]
            v0 = two_leg[0]
            v1 = two_leg[1]
            lighting = self.lit(lighting, center)
            lighting = self.lit(lighting, v0)
            lighting = self.lit(lighting, omega)
            lighting = self.lit(lighting, center)
            lighting = self.lit(lighting, long_leg[0])
            lighting = self.lit(lighting, v1)
            lighting = self.lit(lighting, v0)
            lighting = self.lit(lighting, center)
            lighting = self.lit(lighting, long_leg[1])
            lighting = self.lit(lighting, long_leg[0])
            lighting = self.lit(lighting, long_leg[2])
            lighting = self.lit(lighting, long_leg[1])
            lighting = self.lit(lighting, long_leg[3])
            lighting = self.lit(lighting, long_leg[2])
            lighting = self.lit(lighting, omega)
            lighting = self.lit(lighting, center)
            lighting = self.lit(lighting, long_leg[0])
            lighting = self.lit(lighting, long_leg[1])
            lighting = self.lit(lighting, v0)
            lighting = self.lit(lighting, v1)
            lighting = self.lit(lighting, center)
            lighting = self.lit(lighting, long_leg[0])
            lighting = self.lit(lighting, v0)
            lighting = self.lit(lighting, center)
            self.append_to_center(lighting)
            raise AppendedException
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _append_long_leg_only_last_lit(self) -> Self:
        """
        Step VI. Append if long leg last and center are lited
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Append if long leg last and center are lited")
        center = self.get_center()
        long_leg = self.get_long_leg()
        lits = self.get_lits(lighting, long_leg)
        if len(lits) == 1:
            self.check_dependency_one_leg(lighting)
            last_v = long_leg[len(long_leg) - 1]
            if len(long_leg) == 1:
                lighting = self.lit(lighting, last_v)
                self.append(lighting, last_v)
                raise AppendedException
            g = long_leg[len(long_leg) - 2]
            omega = self.get_one_vertix()
            pq = omega@lighting
            new_g = pq@g
            if self.is_included(new_g):
                raise DependentException()
            self.remove(last_v)
            self.append(lighting, center)
            self.replace(g, new_g)
            self.append(last_v, lighting)
            long_leg = self.get_long_leg()
            if len(long_leg) > 4:
                for i in range(4, len(long_leg)):
                    self.append_delayed(long_leg[i])
                self.remove(long_leg[4])
            raise AppendedException
        self.print_state(lighting)
        self.set_lighting(lighting)
        return self

    def _append_long_leg_last_and_first_lit(self) -> None:
        """
        Step VII. Append if long leg last, first and center are lited
        """
        lighting = self.get_lighting()
        self.print_vertix(lighting, "Append if long leg last, first and center are lited")
        omega = self.get_one_vertix()
        center = self.get_center()
        long_leg = self.get_long_leg()
        first_v = long_leg[0]
        for i in range(len(long_leg)-1, 0, -1):
            lighting = self.lit(lighting, long_leg[i])
        lighting = self.lit(lighting, center)
        lighting = self.lit(lighting, omega)
        lighting = self.lit(lighting, first_v)
        lighting = self.lit(lighting, center)
        self.append_to_center(lighting)
        raise AppendedException


    def _pipeline(self, lighting: PauliString) -> None:
        """
        Pipeline
        """
        self.print_vertix(lighting, "Appending vertix to graph")
        # pipeline building
        self.set_lighting(lighting)
        self._append_three_graph()
        self._append_one_legs_in_different_state()
        self._append_fast()
        self._lit_only_long_leg()
        self._lit_center()
        self._reduce_long_leg_more_than_one_lits()
        self._append_long_leg_first_and_center_lit()
        self._append_long_leg_only_last_lit()
        self._append_long_leg_last_and_first_lit()

    def append_delayed(self, v:PauliString) -> None:
        """
        Append to delayed
        """
        self.delayed_vertices.append(v)

    def restore_delayed(self, vertices:list[PauliString]) -> list[PauliString]:
        """
        Restore to delayed
        """
        for i in range(len(self.delayed_vertices) - 1, -1, -1):
            vertices.insert(0, self.delayed_vertices[i])
        self.delayed_vertices = []
        return vertices

    def set_debug_vertix(self, lighting:PauliString) -> None:
        """
        Set debug vertix
        """
        self.debug_lighting = lighting

    def set_debug_break(self, lighting:PauliString) -> None:
        """
        Set debug break
        """
        if lighting == self.debug_lighting:
            self.debug_break = True

    def debugbreak(self, number:int=None, lighting:PauliString=None,
                   append:bool=True) -> None:
        """
        Debug break
        """
        if self.is_break():
            if append:
                self.append(lighting, self.get_center())
            self.print_state(lighting)
            raise DebugException()
        if number is not None:
            if number <= len(self.get_vertices()):
                self.print_state(lighting)
                raise DebugException()
        if lighting is not None:
            self.set_debug_break(lighting)

    def is_break(self) -> bool:
        """
        Check debug
        """
        return self.debug_break

    def _get_anti_commutates(self, pauli_string:PauliString,
                             generators) -> list[PauliString]:
        """
            Get a collection of non-commuting Pauli strings
            Args:
                Pauli string to which commutators are defined
            generators: The area of Pauli strings over which to build a graph.
            If not specified, then collection
        """
        return [g for g in generators
               if g != pauli_string and not pauli_string|g]


    def _get_max_connected(self, generators:list[PauliString]
    ) -> tuple[PauliString|None, list[PauliString]|None]:
        """Get the Pauli string that has the maximum number of non-commutable"""
        if len(generators) == 0:
            return None, None
        pauli_string = generators[0]
        anti_commutates = self._get_anti_commutates(pauli_string, generators)
        for p in generators:
            _anti_commutates = self._get_anti_commutates(p, generators)
            if len(_anti_commutates) > len(anti_commutates):
                pauli_string = p
                anti_commutates = _anti_commutates
        return pauli_string, anti_commutates



    def _append_to_queue(self, queue_pauli_strings:list[PauliString],
                         pauli_strings:list[PauliString]) -> None:
        """Append the next related Pauli string to the queue"""
        for p in pauli_strings:
            if p in queue_pauli_strings:
                pauli_strings.remove(p)
                continue
            anti_commutates = self._get_anti_commutates(p, queue_pauli_strings)
            if len(anti_commutates) == 0:
                continue
            if len(anti_commutates) > 1:
                min_index = len(queue_pauli_strings)
                for anti_commutate in anti_commutates:
                    index = queue_pauli_strings.index(anti_commutate)
                    if index < min_index:
                        min_index = index
                        queue_pauli_strings.insert(min_index + 1, p)
            else:
                queue_pauli_strings.append(p)
            pauli_strings.remove(p)
            return

    def _get_queue(self, generators:list[PauliString])->list[PauliString]:
        """Get associated sequence of Pauli strings"""
        new_generators = generators.copy()
        new_generators.sort()
        queue_pauli_strings = []
        pauli_string, anti_commutates = self._get_max_connected(new_generators)

        new_generators.remove(pauli_string)
        queue_pauli_strings.append(pauli_string)
        for anti_commutate in anti_commutates:
            new_generators.remove(anti_commutate)
            if anti_commutate not in queue_pauli_strings:
                queue_pauli_strings.append(anti_commutate)

        while len(new_generators) > 0:
            self._append_to_queue(queue_pauli_strings, new_generators)
        return queue_pauli_strings

    def build(self, generators:list[PauliString]) -> Self:
        """
        Transform a connected graph to a cononic type.
        """
        #self.debuging()
        if len(generators) == 0:
            return self
        vertices = self._get_queue(generators)
        queue = vertices.copy()
        #self.set_debug(generators.get_debug())
        self.print_vertices(vertices, "init")
        unappended = []
        self.dependents = []
        #self.set_debug_vertix("ZYIIII")
        while len(vertices) > 0:
            lighting = vertices[0]
            vertices.remove(lighting)
            try:
                #self.debugbreak(number=37, lighting = None)
                self._pipeline(lighting)
            except AppendedException:
                vertices = self.restore_delayed(vertices)
                if lighting in unappended:
                    unappended.remove(lighting)
            except DependentException:
                self.dependents.append(lighting)
                self.print_vertix(lighting, "Dependency")
                vertices = self.restore_delayed(vertices)
            except NotConnectedException:
                vertices = self.restore_delayed(vertices)
                #exc_type, exc_obj, exc_tb = sys.exc_info()
                self.print_vertix(lighting, f"Exception {traceback.format_exc()}")
                if lighting not in unappended:
                    unappended.append(lighting)
                    vertices.append(lighting)
            except DebugException:
                #exc_type, exc_obj, exc_tb = sys.exc_info()
                self.print_vertix(lighting, f"Debug exception {traceback.format_exc()}")
                break
            except RaiseException:
                self.debuging()
                self.print_vertices(queue, "init")
                self.restore()
                break
            except Exception as e:
                vertices = self.restore_delayed(vertices)
                if self.debug:
                    #exc_type, exc_obj, exc_tb = sys.exc_info()
                    self.print_vertix(lighting, f"Exception {traceback.format_exc()}")
                else:
                    self.debuging()
                    self.print_vertix(lighting, f"Exception {e}")
                    self.restore()
                unappended.append(lighting)
        self.print_state()
        self.print_vertices(unappended, "unappended")
        self.print_vertices(self.dependents, "dependents")
        return self

    def is_eq(self, legs:list[list[PauliString]], generators:list[PauliString]) -> bool:
        """
        Testing for equivalence of two algebras. 
        All Pauli strings of one algebra are dependent on another.
        return bool
        """
        self.legs = legs.copy()
        for g in generators:
            try:
                self._pipeline(g)
            except AppendedException:
                self.legs = legs.copy()
                return False
            except Exception:
                continue

        self.legs = legs.copy()
        return True

    def select_dependents(self, legs:list[list[PauliString]], generators:list[PauliString]
        ) -> list[PauliString]:
        """
        Select dependent strings
        """
        self.legs = legs.copy()
        self.is_check = True
        dependents = []

        for g in generators:
            self.legs = legs.copy()
            try:
                self._pipeline(g)
            except CheckAppendedException:
                continue
            except DependentException:
                dependents.append(g)
            except Exception:
                continue

        return dependents
