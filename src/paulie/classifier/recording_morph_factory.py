"""
Recording factory for constructing a canonical graph
"""
from collections.abc import Generator
from typing import Self
from paulie.helpers._recording import recording_graph, RecordGraph
from paulie.classifier.classification import Morph, Classification
from paulie.classifier.morph_factory import MorphFactory, AppendedException,DependentException,
     NotConnectedException, MorphFactoryException
from paulie.common.pauli_string_bitarray import PauliString


class RecordingMorphFactory(MorphFactory):
    """
    Factory for constructing a canonical graph.
    """
    def __init__(self, record:RecordGraph=None) -> None:
        """
        Constructor.

        Args:
            record: Record.
        Returns:
            None

        """
        super().__init__()
        self.record = record

    def lit(self, lighting:PauliString, vertex:PauliString) -> PauliString:
        """
        Lit vertex.

        Args:
            lighting: Canonical graph join candidate
            vertex: Vertex that will be lited by lightning.
        Returns:
            New lighting.
        """
        lighting = lighting@vertex
        if self.is_included(lighting):
            recording_graph(self.record, lighting=lighting,
                            dependent=lighting, title=f"Dependent: {lighting}")
            raise DependentException()
        return lighting

    def check_dependency_one_leg(self, lighting:PauliString) -> None:
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
        ones = self.get_one_vertices()
        vertices = set(self.get_vertices())
        for one in ones:
            pq = one @ lighting
            for v in vertices:
                if v == one:
                    continue
                n_v = pq @ v
                if n_v in vertices or n_v == lighting:
                    recording_graph(self.record, lighting=lighting,
                    dependent=n_v, title=f"Dependent: {lighting}")
                    raise DependentException()


    def append_to_two_center(self, lighting:PauliString) -> None:
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
        center = self.get_center()
        if len(self.legs) == 1:
            recording_graph(self.record, lighting=lighting, lits=[center],
            title=f"Step I: {lighting}")
            recording_graph(self.record, lighting=lighting, lits=[center],
            appending=center, title=f"Step I: {lighting}")
            self.append(lighting, center)
            return
        vertices = self.get_vertices()
        lits = self.get_lits(lighting, vertices)
        recording_graph(self.record, lighting=lighting, lits=lits, title=f"Step I: {lighting}")
        if len(lits) == 1:
            if center in lits:
                self.append(lighting, center)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                appending=center, title=f"Step I: {lighting}")
                return
            else:
                lighting = self.lit(lighting, lits[0])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=lits[0], title=f"Step I: {lighting}")
                lighting = self.lit(lighting, center)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=center, title=f"Step I: {lighting}")
                self.append(lighting, center)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                appending=center, title=f"Step I: {lighting}")
                return
        if len(lits) == 2:
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step I: {lighting}")
            self.append(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            appending=center, title=f"Step I: {lighting}")
            return
        raise NotConnectedException()

    def _append_three_graph(self) -> Self:
        """
        Step I. Construct a graph of three vertices.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step I: {lighting}")
        if self.is_empty():
            self.set_center(lighting)
            raise AppendedException
        if self.is_included(lighting):
            recording_graph(self.record, lighting=lighting, dependent=lighting,
            title=f"Dependent: {lighting}")
            raise DependentException
        if self.is_empty_legs():
            self.append_to_two_center(lighting)
            raise AppendedException
        self.set_lighting(lighting)
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
        lighting = self.get_lighting()
        pq, p, q = self.get_pq(lighting)
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step II: {lighting}")

        if pq is not None:
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            p=p, q=q, title=f"Step II: {lighting}")
            lits = self.get_lits(lighting)
            replacing = []
            for lit in lits:
                if lit != p:
                    replacing.append(lit)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            replacing_vertices=replacing, title=f"Step II: {lighting}")
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    if self.is_included(v):
                        recording_graph(self.record, lighting=lighting, dependent=lit,
                        title=f"Dependent: {lighting}")
                        raise DependentException()
            for lit in lits:
                if lit != p:
                    v = pq@lit
                    self.replace(lit, v)
                    replacing.append(lit)
            recording_graph(self.record, collection=self.get_vertices(),lighting=lighting,
            lits=self.get_lits(lighting), title=f"Step II: {lighting}")
            self.append(lighting, p)
            long_leg = self.get_long_leg()
            if len(long_leg) > 4:
                removing = []
                for i in range(4, len(long_leg)):
                    self.append_delayed(long_leg[i])
                    removing.append(long_leg[i])
                if len(removing) > 0:
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    removing_vertices=removing, title=f"Step II: {lighting}")
                self.remove(long_leg[4])
            recording_graph(self.record, collection=self.get_vertices(), lighting=lighting,
            lits=self.get_lits(lighting), appending=p, title=f"Step II: {lighting}")
            raise AppendedException
        self.set_lighting(lighting)
        return self

    def _lit_only_long_leg(self) -> Self:
        """
        Step III. Lit only the long leg.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step III: {lighting}")
        omega = self.get_one_vertex()
        center = self.get_center()
        center_lits = self.get_lits(lighting, [center])
        lits = self.get_lits(lighting, [omega])
        if omega in lits:
            if center not in center_lits:
                lighting = self.lit(lighting, omega)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=omega, title=f"Step III: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step III: {lighting}")
        two_legs = self.get_two_legs()
        long_leg = self.get_long_leg()
        if len(long_leg) == 2:
            del two_legs[len(two_legs) - 1]
        elif len(long_leg) == 1:
            self.set_lighting(lighting)
            return self
        if len(two_legs) == 0:
            self.set_lighting(lighting)
            return self
        long_lits = self.get_lits(lighting, long_leg)
        if len(long_lits) == 0:
            # find lited two leg
            center_lits = self.get_lits(lighting, [center])
            if center in center_lits:
                lighting = self.lit(lighting, center)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=center, title=f"Step III: {lighting}")
                lighting = self.lit(lighting, long_leg[0])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=long_leg[0], title=f"Step III: {lighting}")
                lighting = self.lit(lighting, omega)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=omega, title=f"Step III: {lighting}")
                lighting = self.lit(lighting, center)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=center, title=f"Step III: {lighting}")
            else:
                for leg in two_legs:
                    lits = self.get_lits(lighting, leg)
                    v0 = leg[0]
                    v1 = leg[1]
                    if v1 in lits and v0 not in lits:
                        lighting = self.lit(lighting, v1)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=v1, title=f"Step III: {lighting}")
                        lits.append(v0)
                    if v0 in lits:
                        lighting = self.lit(lighting, v0)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=v0, title=f"Step III: {lighting}")
                        lighting = self.lit(lighting, center)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=center, title=f"Step III: {lighting}")
                        lighting = self.lit(lighting, long_leg[0])
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=long_leg[0], title=f"Step III: {lighting}")
                        lighting = self.lit(lighting, omega)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=omega, title=f"Step III: {lighting}")
                        lighting = self.lit(lighting, center)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=center, title=f"Step III: {lighting}")
                        break
        long_lits = self.get_lits(lighting, long_leg)
        lit_indexes = self.get_lit_indexes(long_leg, long_lits)
        if 1 not in lit_indexes:
            if 0 in lit_indexes:
                lighting = self.lit(lighting, long_leg[0])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=long_leg[0], title=f"Step III: {lighting}")
            else:
                if len(lit_indexes) == 0:
                    raise NotConnectedException()
                first_lit = lit_indexes[0]
                for i in range(first_lit, 1, -1):
                    lighting = self.lit(lighting, long_leg[i])
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=long_leg[i], title=f"Step III: {lighting}")
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
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=v0, title=f"Step III: {lighting}")
                lits.append(v1)
            elif v0 not in lits and v1 in lits:
                lighting = self.lit(lighting, v1)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=v1, title=f"Step III: {lighting}")
                lits.append(v0)
            if v0 in lits and v1 in lits:
                center_lits = self.get_lits(lighting, [center])
                if center in center_lits:
                    lighting = self.lit(lighting, center)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=center, title=f"Step III: {lighting}")
                    #omega is lited
                    lighting = self.lit(lighting, v1)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=v1, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, v0)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=v0, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, omega)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=omega, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, center)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=center, title=f"Step III: {lighting}")
                else:
                    long_lits = self.get_lits(lighting, [long_leg[0]])
                    if len(long_lits) == 0:
                        lighting = self.lit(lighting, long_v1)
                        recording_graph(self.record, lighting=lighting,
                        lits=self.get_lits(lighting),
                        contracting=long_v1, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, long_v0)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=long_v0, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, center)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=center, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, omega)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=omega, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, v1)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=v1, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, v0)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=v0, title=f"Step III: {lighting}")
                    lighting = self.lit(lighting, center)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=center, title=f"Step III: {lighting}")
        self.set_lighting(lighting)
        return self

    def _lit_center(self) -> Self:
        """
        Lit center.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step IV: {lighting}")
        center = self.get_center()
        center_lits = self.get_lits(lighting, [center])
        if center not in center_lits:
            long_leg = self.get_long_leg()
            long_lits = self.get_lits(lighting, long_leg)
            lit_indexes = self.get_lit_indexes(long_leg, long_lits)
            first_lit = lit_indexes[0]
            for i in range(first_lit, -1, -1):
                lighting = self.lit(lighting, long_leg[i])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=long_leg[i], title=f"Step IV: {lighting}")
        self.set_lighting(lighting)
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
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step IV: {lighting}")
        long_leg = self.get_long_leg()
        #Offset the first vertex to the end of the linear leg. Since the leg is finite,
        #we will always reach the highlighting of one vertex,
        #it will be either the last one or the penultimate one
        while True:
            lits = self.get_lits(lighting, long_leg)
            if len(lits) == 0:
                self.append_to_center(lighting)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                appending=center, title=f"Step IV: {lighting}")
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
                        removing = []
                        for i in range(lit_indexes[0] + 1, len(long_leg)):
                            self.append_delayed(long_leg[i])
                            removing.append(long_leg[i])
                        if len(removing) > 0:
                            recording_graph(self.record, lighting=lighting,
                            lits=self.get_lits(lighting),
                            removing_vertices=removing, title=f"Step IV: {lighting}")
                        self.remove(long_leg[lit_indexes[0] + 1])
                        recording_graph(self.record, collection=self.get_vertices(),
                        lighting=lighting,
                        lits=self.get_lits(lighting), title=f"Step IV: {lighting}")
                    break
            lit_indexes = self.get_lit_indexes(long_leg, lits)
            first = lit_indexes[0]
            second = lit_indexes[1]
            if first > 0 and first + 1 != second:
                for i in range(second, first, -1): ## maybe + 1
                    lighting = self.lit(lighting, long_leg[i])
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=long_leg[i], title=f"Step IV: {lighting}")
            else:
                lighting = self.lit(lighting, long_leg[second])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=long_leg[second], title=f"Step IV: {lighting}")
        self.set_lighting(lighting)
        return self


    def _append_long_leg_first_and_center_lit(self) -> Self:
        """
        Step V. Append long leg with first lit and center.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step V: {lighting}")
        omega = self.get_one_vertex()
        center = self.get_center()
        lits = self.get_lits(lighting, [center, omega])
        is_center_lit = center in lits
        # lit only long leg
        long_leg = self.get_long_leg()
        lits = self.get_lits(lighting, long_leg)
        if is_center_lit and len(lits) == 0:
            #lit only center
            self.append(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            appending=center, title=f"Step V: {lighting}")
            raise AppendedException
        lit_indexes = self.get_lit_indexes(long_leg, lits)
        if len(lit_indexes) == 1 and 0 in lit_indexes:
            # if leg less than 4 or not legs long 2, we can connect to the end of long leg
            is_can_connect_to_end = True
            if self.is_two_leg() and len(long_leg) > 3:
                is_can_connect_to_end = False
            if is_can_connect_to_end:
                for v in long_leg:
                    lighting = self.lit(lighting, v)
                    recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                    contracting=v, title=f"Step V: {lighting}")
                self.append(lighting, long_leg[len(long_leg)-1])
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                appending=long_leg[len(long_leg)-1], title=f"Step V: {lighting}")
                raise AppendedException
            two_legs = self.get_two_legs()
            two_leg = two_legs[0]
            v0 = two_leg[0]
            v1 = two_leg[1]
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v0)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v0, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, omega)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=omega, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[0])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[0], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v1)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v1, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v0)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v0, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[1])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[1], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[0])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[0], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[2])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[2], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[1])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[1], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[3])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[3], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[2])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[2], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, omega)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=omega, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[0])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[0], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[1])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[1], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v0)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v0, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v1)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v1, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, long_leg[0])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[0], title=f"Step V: {lighting}")
            lighting = self.lit(lighting, v0)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=v0, title=f"Step V: {lighting}")
            lighting = self.lit(lighting, center)
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=center, title=f"Step V: {lighting}")
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            appending=center, title=f"Step V: {lighting}")
            self.append(lighting, center)
            raise AppendedException
        self.set_lighting(lighting)
        return self

    def _append_long_leg_only_last_lit(self) -> Self:
        """
        Step VI. Append if long leg last and center are lited.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step VI: {lighting}")
        center = self.get_center()
        long_leg = self.get_long_leg()
        lits = self.get_lits(lighting, long_leg)
        if len(lits) == 1:
            self.check_dependency_one_leg(lighting)
            last_v = long_leg[len(long_leg) - 1]
            if len(long_leg) == 1:
                lighting = self.lit(lighting, last_v)
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                contracting=last_v, title=f"Step VI: {lighting}")
                recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
                appending=last_v, title=f"Step VI: {lighting}")
                self.append(lighting, last_v)
                raise AppendedException
            g = long_leg[len(long_leg) - 2]
            omega = self.get_one_vertex()
            pq = omega@lighting
            new_g = pq@g
            if self.is_included(new_g):
                recording_graph(self.record, lighting=lighting, dependent=g,
                title=f"Dependent: {lighting}")
                raise DependentException()
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            removing_vertices=[last_v], title=f"Step VI: {lighting}")
            self.remove(last_v)
            recording_graph(self.record, collection=self.get_vertices(), lighting=lighting,
            lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            appending=center, title=f"Step VI: {lighting}")
            self.append(lighting, center)
            recording_graph(self.record, collection=self.get_vertices(), lighting=lighting,
            lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
            recording_graph(self.record, lighting=last_v, lits=self.get_lits(lighting),
            replacing_vertices=[g], title=f"Step VI: {lighting}")
            self.replace(g, new_g)
            recording_graph(self.record, collection=self.get_vertices(), lighting=last_v,
            lits=self.get_lits(lighting), title=f"Step VI: {lighting}")
            long_leg = self.get_long_leg()
            if len(long_leg) > 4:
                removing = []
                for i in range(4, len(long_leg)):
                    self.append_delayed(long_leg[i])
                    removing.append(long_leg[i])
                if len(removing) > 0:
                    recording_graph(self.record, lighting=last_v, lits=self.get_lits(lighting),
                    removing_vertices=removing, title=f"Step VI: {lighting}")
                self.remove(long_leg[4])
            recording_graph(self.record, collection=self.get_vertices(), lighting=last_v,
            lits=self.get_lits(lighting), appending=lighting, title=f"Step VI: {lighting}")
            self.append(last_v, lighting)
            raise AppendedException
        self.set_lighting(lighting)
        return self

    def _append_long_leg_last_and_first_lit(self) -> None:
        """
        Step VII. Append if long leg last, first and center are lited.

        Returns:
            Self.
        Raises:
            AppendedException:
                If success.
            DependentException:
                If added vertex is dependent.
        """
        lighting = self.get_lighting()
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        title=f"Step VII: {lighting}")
        omega = self.get_one_vertex()
        center = self.get_center()
        long_leg = self.get_long_leg()
        first_v = long_leg[0]
        for i in range(len(long_leg)-1, 0, -1):
            lighting = self.lit(lighting, long_leg[i])
            recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
            contracting=long_leg[i], title=f"Step VII: {lighting}")
        lighting = self.lit(lighting, center)
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        contracting=center, title=f"Step VII: {lighting}")
        lighting = self.lit(lighting, omega)
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        contracting=omega, title=f"Step VII: {lighting}")
        lighting = self.lit(lighting, first_v)
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        contracting=first_v, title=f"Step VII: {lighting}")
        lighting = self.lit(lighting, center)
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        contracting=center, title=f"Step VII: {lighting}")
        recording_graph(self.record, lighting=lighting, lits=self.get_lits(lighting),
        appending=center, title=f"Step VII: {lighting}")
        self.append(lighting, center)
        raise AppendedException

    def build(self, generators:list[PauliString]) -> Self:
        """
        Transform a connected graph to a canonic type.

        Args:
            generators: List of Pauli strings.
        Returns:
            Self
        """
        if len(generators) == 0:
            return self
        vertices = self._get_queue(generators)
        recording_graph(self.record, collection=vertices, title="Original graph", init=True)
        def record_lighiting(lighting: PauliString) -> None:
            recording_graph(self.record, collection=self.get_vertices(), lighting=lighting,
            title=f"Adding: {lighting}")

        self._build(vertices, record_lighiting)
        classification = Classification()
        classification.add(self.get_morph())
        recording_graph(self.record, collection=self.get_vertices(),
        title=f"Algebra: {classification.get_algebra()}")
        return self
