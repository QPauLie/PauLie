"""
Factory for constructing a canonical graph
"""
from collections.abc import Generator
from typing import Self, Callable
from paulie.classifier.classification import Morph
from paulie.common.pauli_string_bitarray import PauliString

class AppendedException(Exception):
    """
    The vertex is appended
    """

class CheckAppendedException(Exception):
    """
    The vertex is appended
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

class MorphFactory:
    """
    Factory for constructing a canonical graph
    """
    def __init__(self) -> None:
        """
        Initialize a MorphFactory
        """
        self.legs = [] # center is zero leg
        self.lighting = None
        self.delayed_vertices = []
        self.dependents = []
        self.is_check = False
        self.all_vertices = set()
        self.qubit_to_vertices = {}


    def _set_lighting(self, lighting:PauliString) -> None:
        """
        Set lighting.

        Args:
            lighting: Paulie string, which is lightning.
        Returns:
            None
        """
        self.lighting = lighting

    def _get_lighting(self) -> PauliString:
        """
        Get lighting.

        Returns:
            Lighting.
        """
        return self.lighting

    def _lit(self, lighting:PauliString, vertex:PauliString) -> PauliString:
        """
        Lit vertex.

        Args:
            lighting: Canonical graph join candidate.
            vertex: Vertex that will be lited by lightning.
        Returns:
            New lighting.
        """
        lighting = lighting@vertex
        if self._is_included(lighting):
            raise DependentException()
        return lighting

    def _get_lits(self, lighting:PauliString,
                 vertices:list[PauliString]=None) ->list[PauliString]:
        """
        Get lited vertices (connected to the selected vertex).

        Args:
            lighting: Canonical graph join candidate.
            vertices: List of vertices.
        Returns:
            list of lited vertices (connected to the selected vertex).
        """
        if vertices is not None:
            return [v for v in vertices if v != lighting and not lighting|v]

        # Use overlap index for efficiency in sparse graphs
        candidates = set()
        for qubit in lighting.get_support():
            if qubit in self.qubit_to_vertices:
                for v in self.qubit_to_vertices[qubit]:
                    if v != lighting:
                        candidates.add(v)

        return [v for v in candidates if not lighting|v]

    def _add_to_index(self, v: PauliString) -> None:
        self.all_vertices.add(v)
        for qubit in v.get_support():
            if qubit not in self.qubit_to_vertices:
                self.qubit_to_vertices[qubit] = []
            self.qubit_to_vertices[qubit].append(v)

    def _remove_from_index(self, v: PauliString) -> None:
        if v in self.all_vertices:
            self.all_vertices.remove(v)
        for qubit in v.get_support():
            if qubit in self.qubit_to_vertices:
                if v in self.qubit_to_vertices[qubit]:
                    self.qubit_to_vertices[qubit].remove(v)

    def _is_empty(self) -> bool:
        """
        Checking for emptiness.

        Returns:
            True if graph is empty.
        """
        return len(self.legs) == 0

    def _is_empty_legs(self) -> bool:
        """
        Checking for missing legs.

        Returns:
            True if graph no core.
        """
        return len(self.legs) < 3

    def _find_in_leg(self, leg:list[PauliString], v:PauliString) -> int:
        """
        Find vertex in leg.

        Args:
            leg: List of vertices.
            v: Required vertex.
        Returns:
            Index of vertex in the leg.
        """
        try:
            index = leg.index(v)
        except ValueError:
            index = -1
        return index


    def _find(self, v:PauliString) -> tuple[int,int]:
        """
        Find vertex.

        Args:
            v: Required vertex.
        Returns:
            Tuple index of leg and index vertex in the leg.
        """
        if v not in self.all_vertices:
            return -1, -1
        for i, leg in enumerate(self.legs):
            index = self._find_in_leg(leg, v)
            if index > -1:
                return i, index
        return -1, -1

    def _is_included(self, v:PauliString) -> bool:
        """
        Checking a vertex for inclusion in the graph.

        Args:
            v: Required vertex.
        Returns:
            True if vertex in the graph.
        """
        return v in self.all_vertices

    def _get_vertices(self) -> list[PauliString]:
        """
        Get graph vertices.

        Returns:
           List of vertices in the graph.
        """
        return [v for leg in self.legs for v in leg]

    def _get_center(self) -> PauliString|None:
        """
        Get center.

        Returns:
           Canter of the graph.
        """
        if self._is_empty():
            return None
        return self.legs[0][0]

    def _set_center(self, v: PauliString) -> None:
        """
        Set center.

        Args:
            v: Vertex of the graph center.
        Returns:
            None
        """
        if self._is_empty() is False:
            raise MorphFactoryException("Center is setted")
        self.legs.append([v])
        self._add_to_index(v)

    def _get_long_leg(self) -> list[PauliString]:
        """
        Get long leg.

        Returns:
            List of vertices in long leg.
        """
        if self._is_empty_legs():
            raise MorphFactoryException("No legs")
        return self.legs[len(self.legs) - 1]

    def _get_one_vertex(self) -> PauliString:
        """
        Get one vertex in leg.

        Returns:
            Control vertex in the graph.
        """
        if self._is_empty_legs():
            raise MorphFactoryException("No legs")
        return self.legs[1][0]

    def _gen_one_legs(self) -> Generator[list[list[PauliString]], None, None]:
        """
        Generate vertices included in single legs.

        Yields:
            Vertices included in legs of length 1.
        """
        if self._is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(1, len(self.legs)):
            if len(self.legs[i]) == 1:
                yield self.legs[i]
            else:
                break

    def _get_one_vertices(self)->list[PauliString]:
        """
        Get vertices included in single legs.

        Returns:
            List of vertices included in legs of length 1.
        """
        vertices = []
        for leg in self._gen_one_legs():
            vertices.append(leg[0].copy())
        return vertices

    def _get_pq(self, lighting:PauliString) -> tuple[PauliString|None,PauliString|None]:
        """
        Get pq.

        Args:
            lighting: Canonical graph join candidate.
        Returns:
            tuple:
                pq: p@q where p and q are vertices in one leg,
                and p is lited, q is unlited.
                p: Lited vertex.
        """
        one_verices = self._get_one_vertices()
        lits = self._get_lits(lighting, one_verices)
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
        Generate vertices included in two legs.

        Yields:
            Legs long 2.
        """
        if self._is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(1, len(self.legs)):
            if len(self.legs[i]) == 2:
                yield self.legs[i]
            else:
                if len(self.legs[i]) > 2:
                    break

    def _get_two_legs(self) -> list[tuple[PauliString, PauliString]]:
        """
        Get vertices included in two legs.

        Returns:
            List of vertices included in two legs.
        """
        return [(leg[0].copy(), leg[1].copy()) for leg in self._gen_two_legs()]

    def _get_count_two_legs(self) -> int:
        """
        Get the number of legs of length two.

        Returns:
            Number of legs long 2.
        """
        return len(self._get_two_legs())

    def _is_two_leg(self) -> bool:
        """
        Checking the leg for length two.

        Returns:
            True if there is a leg of length 2.
        """
        count_two_legs = self._get_count_two_legs()
        if count_two_legs == 0:
            return False

        long_leg = self._get_long_leg()
        if len(long_leg) != 2:
            return True
        return count_two_legs > 1

    def _gen_long_legs(self) -> list[tuple[PauliString, PauliString]]:
        """
        Generate long leg vertices.

        Yields:
           Vertices from long leg.
        """
        if self._is_empty_legs():
            raise MorphFactoryException("No legs")
        for i in range(len(self.legs)-1, 1, -1):
            if len(self.legs) > 2:
                yield self.legs[i]
            else:
                break

    def _append(self, v:PauliString, lit:PauliString) -> None:
        """
        Append vertex to graph.

        Args:
            v: Added vertex.
            lit: Vertex to which is added.
        Returns:
            None
        Raises:
            CheckAppendedException:
                If check mode.
            MorphFactoryException:
                If No vertex.
                or lit is not last in leg
                or can't append.
        """
        if self.is_check:
            raise CheckAppendedException()
        leg_index, vertex_index = self._find(lit)
        if leg_index == -1:
            raise MorphFactoryException("No vertex")

        self._add_to_index(v)

        if leg_index == 0:
            self.legs.insert(1, [v])
            return
        if vertex_index !=  len(self.legs[leg_index]) - 1:
            raise MorphFactoryException("The vertex is not the last")
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

    def _check_dependency_one_leg(self, lighting:PauliString) -> None:
        """
        Dependency check when attaching a vertex to the center of the graph.

        Args:
            lighting: Canonical graph join candidate.
        Returns:
            None
        Raises:
            NotConnectedException:
                If lighting is not connected.
        """
        ones = self._get_one_vertices()
        vertices = self.all_vertices
        for one in ones:
            pq = one @ lighting
            for v in ones:
                if v == one:
                    continue
                n_v = pq @ v
                if n_v in vertices or n_v == lighting:
                    raise DependentException()

    def _remove(self, v:PauliString) -> None:
        """
        Removing a graph vertex.

        Args:
            v: Removed vertex.
        Returns:
            None
        """
        leg_index, vertex_index = self._find(v)
        if leg_index == -1:
            raise MorphFactoryException("No vertex")
        if leg_index == 0:
            raise MorphFactoryException("Can't delete the center")

        # Before removing from self.legs, we should find which vertices are actually removed.
        # Actually, all vertices from vertex_index to the end of the leg are removed.
        # Wait, the original code:
        # leg = [self.legs[leg_index][i] for i in range(0, vertex_index)]
        # del self.legs[leg_index]
        # This removes everything FROM vertex_index to the end.

        for i in range(vertex_index, len(self.legs[leg_index])):
            self._remove_from_index(self.legs[leg_index][i])

        leg = [self.legs[leg_index][i] for i in range(0, vertex_index)]
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

    def _replace(self, v:PauliString, v_new:PauliString) -> None:
        """
        Replacing a graph vertex with an equivalent one.

        Args:
            v: Removed vertices.
            v_new: Added vertices.
        Returns:
            None
        """
        leg_index, vertex_index = self._find(v)
        if leg_index == -1:
            raise MorphFactoryException("No vertex")
        self._remove_from_index(v)
        self._add_to_index(v_new)
        self.legs[leg_index][vertex_index] = v_new


    def _get_lit_indexes(self, vertices:list[PauliString], lits:list[PauliString]) -> list[int]:
        """
        Get the indices of the lited vertices in lits.

        Args:
            vertices: List of vertices.
            lits: List of lited vertices.
        Returns:
            List of indexes of lited vertices in vertices.
        """
        indexes = []
        for i, v in enumerate(vertices):
            if v in lits:
                indexes.append(i)
        return indexes

    def _append_delayed(self, v:PauliString) -> None:
        """
        Append to delayed.

        Args:
            v: Vertex to append in delayed.
        Return:
            None
        """
        self.delayed_vertices.append(v)

    def _restore_delayed(self, vertices:list[PauliString]) -> list[PauliString]:
        """
        Restore to delayed.

        Args:
            vertices: List of vertices.
        Returns:
            List of vertices.
        """
        for i in range(len(self.delayed_vertices) - 1, -1, -1):
            vertices.insert(0, self.delayed_vertices[i])
        self.delayed_vertices = []
        return vertices

    def _append_to_center(self, lighting:PauliString) -> None:
        """
        Joining a vertex to the center of the graph.

        Args:
            lighting: Canonical graph join candidate.
        Returns:
            None
        """
        self._check_dependency_one_leg(lighting)
        center = self._get_center()
        self._append(lighting, center)

    def _append_to_two_center(self, lighting:PauliString) -> None:
        """
        Append vertex to two vertices graph.

        Args:
            lighting: Canonical graph join candidate.
        Returns:
            None
        Raises:
            NotConnectedException:
                If lighting is not connected.
        """
        center = self._get_center()
        if len(self.legs) == 1:
            self._append(lighting, center)
            return
        vertices = self._get_vertices()
        lits = self._get_lits(lighting, vertices)

        if len(lits) == 1:
            if center in lits:
                self._append(lighting, center)
                return
            else:
                lighting = self._lit(lighting, lits[0])
                lighting = self._lit(lighting, center)
                self._append(lighting, center)
                return
        if len(lits) == 2:
            lighting = self._lit(lighting, center)
            self._append(lighting, center)
            return
        raise NotConnectedException()

    def _append_three_graph(self) -> Self:
        """
        Step I. Construct a graph of three vertices.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        if self._is_empty():
            self._set_center(lighting)
            raise AppendedException
        if self._is_included(lighting):
            raise DependentException
        if self._is_empty_legs():
            self._append_to_two_center(lighting)
            raise AppendedException
        self._set_lighting(lighting)
        return self

    def _append_one_legs_in_different_state(self) -> Self:
        """
        Step II. Legs of length 1 in different initial lit states.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        pq, p = self._get_pq(lighting)
        if pq is not None:
            lits = self._get_lits(lighting)
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    if self._is_included(v):
                        raise DependentException()
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    self._replace(lit, v)
            self._append(lighting, p)
            long_leg = self._get_long_leg()
            if len(long_leg) > 4:
                for i in range(4, len(long_leg)):
                    self._append_delayed(long_leg[i])
                self._remove(long_leg[4])
            raise AppendedException
        self._set_lighting(lighting)
        return self

    def _append_fast(self) -> Self:
        """
        Quickly obvious connection of lightning to a graph.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        center = self._get_center()
        two_legs = self._get_two_legs()
        long_leg = self._get_long_leg()
        if len(long_leg) == 2:
            del two_legs[len(two_legs) - 1]
        if len(two_legs) == 0:
            lits = self._get_lits(lighting)
            if len(lits) == 1:
                if center in lits:
                    self._append_to_center(lighting)
                    raise AppendedException
                long_leg = self._get_long_leg()
                if long_leg[len(long_leg) - 1] in lits:
                    self._append(lighting, long_leg[len(long_leg) - 1])
                    raise AppendedException
        self._set_lighting(lighting)
        return self

    def _lit_only_long_leg(self) -> Self:
        """
        Step III. Lit only the long leg.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        omega = self._get_one_vertex()
        center = self._get_center()
        center_lits = self._get_lits(lighting, [center])
        lits = self._get_lits(lighting, [omega])
        if omega in lits:
            if center not in center_lits:
                lighting = self._lit(lighting, omega)
            lighting = self._lit(lighting, center)
        two_legs = self._get_two_legs()
        long_leg = self._get_long_leg()
        if len(long_leg) == 2:
            del two_legs[len(two_legs) - 1]
        elif len(long_leg) == 1:
            self._set_lighting(lighting)
            return self
        if len(two_legs) == 0:
            self._set_lighting(lighting)
            return self
        long_lits = self._get_lits(lighting, long_leg)
        if len(long_lits) == 0:
            # find lited two leg
            center_lits = self._get_lits(lighting, [center])
            if center in center_lits:
                lighting = self._lit(lighting, center)
                lighting = self._lit(lighting, long_leg[0])
                lighting = self._lit(lighting, omega)
                lighting = self._lit(lighting, center)
            else:
                for leg in two_legs:
                    lits = self._get_lits(lighting, leg)
                    v0 = leg[0]
                    v1 = leg[1]
                    if v1 in lits and v0 not in lits:
                        lighting = self._lit(lighting, v1)
                        lits.append(v0)
                    if v0 in lits:
                        lighting = self._lit(lighting, v0)
                        lighting = self._lit(lighting, center)
                        lighting = self._lit(lighting, long_leg[0])
                        lighting = self._lit(lighting, omega)
                        lighting = self._lit(lighting, center)
                        break
        # lit second vertex on long leg
        long_lits = self._get_lits(lighting, long_leg)
        lit_indexes = self._get_lit_indexes(long_leg, long_lits)
        if 1 not in lit_indexes:
            if 0 in lit_indexes:
                lighting = self._lit(lighting, long_leg[0])
            else:
                if len(lit_indexes) == 0:
                    raise NotConnectedException()
                first_lit = lit_indexes[0]
                for i in range(first_lit, 1, -1):
                    lighting = self._lit(lighting, long_leg[i])
        long_v0 = long_leg[0]
        long_v1 = long_leg[1]
        for i,leg in enumerate(two_legs):
            lits = self._get_lits(lighting, leg)
            v0 = leg[0]
            v1 = leg[1]
            if v0 not in lits and v1 not in lits:
                continue
            if v0 in lits and v1 not in lits:
                lighting = self._lit(lighting, v0)
                lits.append(v1)
            elif v0 not in lits and v1 in lits:
                lighting = self._lit(lighting, v1)
                lits.append(v0)
            if v0 in lits and v1 in lits:
                center_lits = self._get_lits(lighting, [center])
                if center in center_lits:
                    lighting = self._lit(lighting, center)
                    #omega is lited
                    lighting = self._lit(lighting, v1)
                    lighting = self._lit(lighting, v0)
                    lighting = self._lit(lighting, omega)
                    lighting = self._lit(lighting, center)
                else:
                    long_lits = self._get_lits(lighting, [long_leg[0]])
                    if len(long_lits) == 0:
                        lighting = self._lit(lighting, long_v1)
                    lighting = self._lit(lighting, long_v0)
                    lighting = self._lit(lighting, center)
                    lighting = self._lit(lighting, omega)
                    lighting = self._lit(lighting, v1)
                    lighting = self._lit(lighting, v0)
                    lighting = self._lit(lighting, center)
        self._set_lighting(lighting)
        return self

    def _lit_center(self) -> Self:
        """
        Lit center.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        center = self._get_center()
        center_lits = self._get_lits(lighting, [center])
        if center not in center_lits:
            long_leg = self._get_long_leg()
            long_lits = self._get_lits(lighting, long_leg)
            lit_indexes = self._get_lit_indexes(long_leg, long_lits)
            first_lit = lit_indexes[0]
            for i in range(first_lit, -1, -1):
                lighting = self._lit(lighting, long_leg[i])
        self._set_lighting(lighting)
        return self


    def _reduce_long_leg_more_than_one_lits(self) -> Self:
        """
        Step IV. Reducing the long leg lits to standard configurations.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        long_leg = self._get_long_leg()
        #Offset the first vertex to the end of the linear leg. Since the leg is finite,
        #we will always reach the highlighting of one vertex,
        #it will be either the last one or the penultimate one
        while True:
            lits = self._get_lits(lighting, long_leg)
            if len(lits) == 0:
                self._append_to_center(lighting)
                raise AppendedException()
            if len(lits) == 3 and len(long_leg) > 10:
                lit_indexes = self._get_lit_indexes(long_leg, lits)
                if lit_indexes[0] == 0 and lit_indexes[1] == lit_indexes[2] - 1:
                    raise DependentException()
            if len(lits) == 2:
                lit_indexes = self._get_lit_indexes(long_leg, lits)
                if lit_indexes[0] == 0 and lit_indexes[1] == len(long_leg) - 1:
                    break
            if len(lits) == 1:
                if long_leg[0] == lits[0] or long_leg[len(long_leg) - 1] == lits[0]:
                    break
                if long_leg[0] != lits[0]:
                    lit_indexes = self._get_lit_indexes(long_leg, lits)
                    if lit_indexes[0] < len(long_leg) - 1:
                        for i in range(lit_indexes[0] + 1, len(long_leg)):
                            self._append_delayed(long_leg[i])
                        self._remove(long_leg[lit_indexes[0] + 1])
                    break
            lit_indexes = self._get_lit_indexes(long_leg, lits)
            first = lit_indexes[0]
            second = lit_indexes[1]
            if first > 0 and first + 1 != second:
                for i in range(second, first, -1): ## maybe + 1
                    lighting = self._lit(lighting, long_leg[i])
            else:
                lighting = self._lit(lighting, long_leg[second])
        self._set_lighting(lighting)
        lits = self._get_lits(lighting, long_leg)

        return self

    def _append_long_leg_first_and_center_lit(self) -> Self:
        """
        Step V. Append long leg with first lit and center.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        omega = self._get_one_vertex()
        center = self._get_center()
        lits = self._get_lits(lighting, [center, omega])
        is_center_lit = center in lits
        # lit only long leg
        long_leg = self._get_long_leg()
        lits = self._get_lits(lighting, long_leg)
        if is_center_lit and len(lits) == 0:
            #lit only center
            self._append_to_center(lighting)
            raise AppendedException
        lit_indexes = self._get_lit_indexes(long_leg, lits)
        # only long leg and center are lited
        if len(lit_indexes) == 1 and 0 in lit_indexes:
            # if leg less than 4 or not legs long 2, we can connect to the end of long leg
            is_can_connect_to_end = True
            if self._is_two_leg() and len(long_leg) > 3:
                is_can_connect_to_end = False
            if is_can_connect_to_end:
                for v in long_leg:
                    lighting = self._lit(lighting, v)
                self._append(lighting, long_leg[len(long_leg)-1])
                raise AppendedException
            two_legs = self._get_two_legs()
            two_leg = two_legs[0]
            v0 = two_leg[0]
            v1 = two_leg[1]
            lighting = self._lit(lighting, center)
            lighting = self._lit(lighting, v0)
            lighting = self._lit(lighting, omega)
            lighting = self._lit(lighting, center)
            lighting = self._lit(lighting, long_leg[0])
            lighting = self._lit(lighting, v1)
            lighting = self._lit(lighting, v0)
            lighting = self._lit(lighting, center)
            lighting = self._lit(lighting, long_leg[1])
            lighting = self._lit(lighting, long_leg[0])
            lighting = self._lit(lighting, long_leg[2])
            lighting = self._lit(lighting, long_leg[1])
            lighting = self._lit(lighting, long_leg[3])
            lighting = self._lit(lighting, long_leg[2])
            lighting = self._lit(lighting, omega)
            lighting = self._lit(lighting, center)
            lighting = self._lit(lighting, long_leg[0])
            lighting = self._lit(lighting, long_leg[1])
            lighting = self._lit(lighting, v0)
            lighting = self._lit(lighting, v1)
            lighting = self._lit(lighting, center)
            lighting = self._lit(lighting, long_leg[0])
            lighting = self._lit(lighting, v0)
            lighting = self._lit(lighting, center)
            self._append_to_center(lighting)
            raise AppendedException
        self._set_lighting(lighting)
        return self

    def _append_long_leg_only_last_lit(self) -> Self:
        """
        Step VI. Append if long leg last and center are lited.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        center = self._get_center()
        long_leg = self._get_long_leg()
        lits = self._get_lits(lighting, long_leg)
        if len(lits) == 1:
            self._check_dependency_one_leg(lighting)
            last_v = long_leg[len(long_leg) - 1]
            if len(long_leg) == 1:
                lighting = self._lit(lighting, last_v)
                self._append(lighting, last_v)
                raise AppendedException
            g = long_leg[len(long_leg) - 2]
            omega = self._get_one_vertex()
            pq = omega@lighting
            new_g = pq@g
            if self._is_included(new_g):
                raise DependentException()
            self._remove(last_v)
            self._append(lighting, center)
            self._replace(g, new_g)
            self._append(last_v, lighting)
            long_leg = self._get_long_leg()
            if len(long_leg) > 4:
                for i in range(4, len(long_leg)):
                    self._append_delayed(long_leg[i])
                self._remove(long_leg[4])
            raise AppendedException
        self._set_lighting(lighting)
        return self

    def _append_long_leg_last_and_first_lit(self) -> None:
        """
        Step VII. Append if long leg last, first and center are lited.

        Returns:
            Self
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self._get_lighting()
        omega = self._get_one_vertex()
        center = self._get_center()
        long_leg = self._get_long_leg()
        first_v = long_leg[0]
        for i in range(len(long_leg)-1, 0, -1):
            lighting = self._lit(lighting, long_leg[i])
        lighting = self._lit(lighting, center)
        lighting = self._lit(lighting, omega)
        lighting = self._lit(lighting, first_v)
        lighting = self._lit(lighting, center)
        self._append_to_center(lighting)
        raise AppendedException


    def _pipeline(self, lighting: PauliString) -> None:
        """
        Pipeline.

        Args:
            lighting: canonical graph join candidate.
        Returns:
            None
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        # pipeline building
        PauliString.set_performance('s1')
        self._set_lighting(lighting)
        self._append_three_graph()
        PauliString.set_performance('s2')
        self._append_one_legs_in_different_state()
        PauliString.set_performance('fast')
        self._append_fast()
        PauliString.set_performance('s3')
        self._lit_only_long_leg()
        self._lit_center()
        PauliString.set_performance('s4')
        self._reduce_long_leg_more_than_one_lits()
        PauliString.set_performance('s4.1')
        self._append_long_leg_first_and_center_lit()
        PauliString.set_performance('s4.2')
        self._append_long_leg_only_last_lit()
        PauliString.set_performance('s4.3')
        self._append_long_leg_last_and_first_lit()
        PauliString.set_performance('all')


    def _get_anti_commutates(self, pauli_string:PauliString,
                             generators) -> list[PauliString]:
        """
        Get a collection of non-commuting Pauli strings.

        Args:
            Pauli string to which commutators are defined.
            generators: The area of Pauli strings over which to build a graph.
            If not specified, then collection.
        Returns:
            a list of non-commuting Pauli strings
        """
        return [g for g in generators
               if g != pauli_string and not pauli_string|g]


    def _get_max_connected(self, generators:list[PauliString]
    ) -> tuple[PauliString|None, list[PauliString]|None]:
        """
        Get the Pauli string that has the maximum number of non-commutable.

        Args:
            generators: List of Pauli strings.
        Returns:
            tuple:
                PauliString
                List of Pauli strings connected with first.
        """
        if len(generators) == 0:
            return None, None
        pauli_string = generators[0]
        anti_commutates = self._get_anti_commutates(pauli_string, generators)
        #for p in generators:
        #    _anti_commutates = self._get_anti_commutates(p, generators)
        #    if len(_anti_commutates) > len(anti_commutates):
        #        pauli_string = p
        #        anti_commutates = _anti_commutates
        return pauli_string, anti_commutates



    def _append_to_queue(self, queue_pauli_strings:list[PauliString],
                         pauli_strings:list[PauliString]) -> None:
        """
        Append the next related Pauli string to the queue.

        Args:
            queue_pauli_strings: Queue of Pauli strings.
            pauli_strings: List of Pauli strings.
        Returns:
            None
        """
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

    def _get_queue(self, generators: list[PauliString]) -> list[PauliString]:
        """
        Get associated sequence of Pauli strings efficiently.

        Args:
            generators: List of Pauli strings.
        Returns:
            Associated sequence of Pauli strings.
        """
        if not generators:
            return []

        # O(M log M)
        new_generators = sorted(generators)
        m = len(new_generators)

        # Precompute all anticommutation pairs for this subset using overlap index
        # O(M * k) for sparse graphs
        adj = [[] for _ in range(m)]
        qubit_to_indices = {}
        for i, g in enumerate(new_generators):
            for qubit in g.get_support():
                if qubit not in qubit_to_indices:
                    qubit_to_indices[qubit] = []
                qubit_to_indices[qubit].append(i)

        for i, a in enumerate(new_generators):
            candidates = set()
            for qubit in a.get_support():
                for j in qubit_to_indices[qubit]:
                    if j > i:
                        candidates.add(j)
            for j in candidates:
                b = new_generators[j]
                if not a | b:
                    adj[i].append(j)
                    adj[j].append(i)

        queue_indices = []
        processed = [False] * m

        # Start with index 0 (consistent with sorted order)
        p0_idx = 0
        processed[p0_idx] = True
        queue_indices.append(p0_idx)

        # Add all neighbors of p0 to start the queue
        for neighbor_idx in adj[p0_idx]:
            if not processed[neighbor_idx]:
                processed[neighbor_idx] = True
                queue_indices.append(neighbor_idx)

        # To maintain index mapping efficiently
        vertex_to_queue_pos = {p0_idx: 0}
        for i, idx in enumerate(queue_indices[1:], 1):
            vertex_to_queue_pos[idx] = i

        # Set of remaining vertices to add
        remaining_indices = [i for i in range(m) if not processed[i]]

        # While there are vertices left, find one connected to the current queue
        while remaining_indices:
            found_idx = -1
            best_p_idx = -1
            min_pos_in_queue = m + 1

            # Find a vertex that connects to the queue at the earliest possible position
            for i, p_idx in enumerate(remaining_indices):
                neighbors_in_queue = [n_idx for n_idx in adj[p_idx] if processed[n_idx]]
                if neighbors_in_queue:
                    # Find min position in queue for its neighbors
                    current_min_pos = min(
                        vertex_to_queue_pos[idx]
                        for idx in neighbors_in_queue
                    )

                    found_idx = i
                    best_p_idx = p_idx
                    best_neighbors_in_queue = neighbors_in_queue
                    min_pos_in_queue = current_min_pos
                    break

            if best_p_idx != -1:
                p_idx = remaining_indices.pop(found_idx)
                processed[p_idx] = True

                if len(best_neighbors_in_queue) > 1:
                    # Insert at min_pos + 1 to maintain connectivity structure
                    queue_indices.insert(min_pos_in_queue + 1, p_idx)
                    # Shift all positions in vertex_to_queue_pos for vertices after insertion
                    for k in range(min_pos_in_queue + 1, len(queue_indices)):
                        vertex_to_queue_pos[queue_indices[k]] = k
                else:
                    queue_indices.append(p_idx)
                    vertex_to_queue_pos[p_idx] = len(queue_indices) - 1
            else:
                # Should not happen for connected component
                break

        return [new_generators[i] for i in queue_indices]

    def _build(self, vertices:list[PauliString],
        call_lighting: Callable[[PauliString], None] = None) -> Self:
        """
        Transform a connected graph to a canonic type.

        Args:
            vertices: List of Pauli strings.
            call_lighting: callback function when adding a new lighting
        Returns:
            Self
        """
        unappended = []
        self.dependents = []
        while len(vertices) > 0:
            lighting = vertices[0]
            if call_lighting:
                call_lighting(lighting)
            vertices.remove(lighting)
            try:
                self._pipeline(lighting)
            except AppendedException:
                vertices = self._restore_delayed(vertices)
                if lighting in unappended:
                    unappended.remove(lighting)
            except DependentException:
                self.dependents.append(lighting)
                vertices = self._restore_delayed(vertices)
            except NotConnectedException:
                vertices = self._restore_delayed(vertices)
                if lighting not in unappended:
                    unappended.append(lighting)
                    vertices.append(lighting)
            except RaiseException:
                self.restore()
                break
            except Exception:
                vertices = self._restore_delayed(vertices)
                unappended.append(lighting)
        return self

    def build(self, generators:list[PauliString]) -> Self:
        """
        Transform a connected graph to a canonic type.

        Args:
            generators (list[PauliString]): List of Pauli strings.
        Returns:
            Self: State of MorphFactory after transformation
        """
        if len(generators) == 0:
            return self
        PauliString.set_performance('queue')
        vertices = self._get_queue(generators)
        PauliString.set_performance('all')
        return self._build(vertices)

    def is_eq(self, legs:list[list[PauliString]], generators:list[PauliString]) -> bool:
        """
        Testing for equivalence of two algebras.
        All Pauli strings of one algebra are dependent on another.

        Args:
            legs (list[list[PauliString]]): List of legs of
                the canonical graph to compare with.
            generators (list[PauliString]): List of Pauli strings to check
                for inclusion in the canonical graph
        Returns:
            bool: The result of checking for the equivalence
                of generators to a given canonical graph.
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
        Selecting generators that are dependent on the canonical graph

        Args:
            legs (list[list[PauliString]]): List of legs of
                the canonical graph to check with.
            generators (list[PauliString]): List of Pauli strings to check
                for inclusion in the canonical graph
        Returns:
            list[PauliString]: List of dependent generators.
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

    def get_morph(self) -> Morph:
        """
        Get the canonical graph built by the factory.

        Returns:
            Morph: Canonical graph.
        """
        return Morph(self.legs, self.dependents)
