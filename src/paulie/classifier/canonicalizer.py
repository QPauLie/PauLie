"""
Canonicalizer of generators
"""
from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.classification import Morph
from paulie.classifier.observer import EventManager

class Canonicalizer:
    """
    Class of canonicalizer of generators.

    The canonicalizer is the *publisher* of the observer pattern: it owns an
    :class:`~paulie.classifier.observer.EventManager` and emits events at every relevant step of
    the algorithm through :meth:`_notify`. When no observer is subscribed, the notifications are
    cheap no-ops, so the behaviour of the algorithm is unchanged.
    """
    def __init__(self) -> None:
        """
        Initialize a Canonicalizer
        """
        self.type = 'A' # B when at least two legs of length at least 2
        self.central_vertex = None
        self.legs = []
        self.vertex_stack = []
        self.events = EventManager()

    def _current_collection(self) -> list[PauliString]:
        """
        Snapshot the vertices of the canonical graph being built.

        The central vertex is placed first, followed by the vertices of each leg. Observers use
        this to reconstruct the anticommutation graph of the current state.

        Returns:
            list[PauliString]: Current vertices, central vertex first.
        """
        collection: list[PauliString] = []
        if self.central_vertex is not None:
            collection.append(self.central_vertex)
        for leg in self.legs:
            collection.extend(leg)
        return collection

    def _notify(self, event: str, **data) -> None:
        """
        Emit an event to subscribed observers.

        Args:
            event (str): Name of the event (used as the frame title).
            **data: Contextual data describing the event (node roles, explicit collection, ...).
        Returns:
            None
        """
        if self.events.has_subscribers():
            self.events.notify(event, self, data)

    def _set_vertex_stack(self, vertex_stack: list[PauliString]) -> None:
        """
        Set initial list of generators
        Args:
            vertex_stack (list[PauliString]): Initial list of generators
        """
        self.vertex_stack = vertex_stack.copy()

    def _representative(self, v: PauliString)-> PauliString:
        """
        Pauli string representative
        Args:
            v (PauliString): Pauli string
        Returns:
            PauliString: representative
        """
        return v

    def _is_lit(self, v: PauliString, w: PauliString) -> bool:
        """
        Testing for non-commutativity of two Pauli strings
        Args:
            v (PauliString): Pauli string
            w (PauliString): Pauli string
        Returns:
            bool: True if the two strings are not commutative
        """
        return not self._representative(v) | self._representative(w)

    def _tracked_multiply(self, v: PauliString, w: PauliString) -> PauliString:
        """
        Tracked multiply of two Pauli strings
        Args:
            v (PauliString): Pauli string
            w (PauliString): Pauli string
        Returns:
            PauliString: multiplication result
        """
        return v @ w

    def _build_core(self, v: PauliString) -> PauliString:
        """
        Build the core of a canonical graph of 3 vertices
        Args:
            v (PauliString): Append Pauli string
        Returns:
            PauliString: Append Pauli string after core installation
        """
        if len(self.legs) == 1:
            if self._is_lit(v, self.legs[0][0]):
                if not self._is_lit(v, self.central_vertex):
                    v = self._tracked_multiply(v, self._representative(self.legs[0][0]))
                v = self._tracked_multiply(v, self._representative(self.central_vertex))
        self.legs.append([v])
        self._notify("Attach to centre", appending=v)
        return v

    def _convert_to_single_lit_state(self, p_index: int, q_index: int,
        vertex_stack: list[PauliString], v: PauliString) -> None:
        """
        Convert to single lit state
        Args:
            p_index (int): Lited single leg index
            q_index (int): Unlited single leg index
            vertex_stack (list[PauliString]): Generator stack
            v (PauliString): Append Pauli string
        """
        self._notify("Single legs in different states (p lit, q unlit)",
            lighting=v, p=self.legs[p_index][0], q=self.legs[q_index][0])
        pq = (self._representative(self.legs[p_index][0]) @
            self._representative(self.legs[q_index][0]))
        replaced: list[PauliString] = []
        if self._is_lit(v, self.central_vertex):
            self.central_vertex = self._tracked_multiply(self.central_vertex, pq)
            replaced.append(self.central_vertex)
        for i in range(len(self.legs)):
            for j in range(len(self.legs[i])):
                if i != p_index and self._is_lit(v, self.legs[i][j]):
                    self.legs[i][j] = self._tracked_multiply(self.legs[i][j], pq)
                    replaced.append(self.legs[i][j])
        self.legs[p_index].append(v)
        self._notify("Append to lit leg", appending=v, replacing=replaced)
        # Truncate longest leg if necessary, this happens at most once
        big_leg_cnt = sum(1 for leg in self.legs if len(leg) >= 2)
        if self.type == 'A' and big_leg_cnt >= 2:
            self.type = 'B'
            removed: list[PauliString] = []
            while len(self.legs[-1]) > 4:
                removed.append(self.legs[-1].pop())
                vertex_stack.append(removed[-1])
            if removed:
                self._notify("Type becomes B, trim longest leg", removing=removed)

    def _transfer_lightning(self, lit_2_leg_index: int, v: PauliString) -> PauliString:
        """
        Transferring lighting to the long leg from leg length 2
        Args:
            lit_2_leg_index (int): Lited leg index long 2
            v (PauliString): Append Pauli string
        Returns:
            PauliString: Pauli string after transferring
        """
        # We need to make self.legs[-1][1] lit
        if not self._is_lit(v, self.legs[-1][1]):
            m = None
            for i in range(len(self.legs[-1])):
                if self._is_lit(v, self.legs[-1][i]):
                    m = i
                    break
            if m is None:
                if not self._is_lit(v, self.central_vertex):
                    if self._is_lit(v, self.legs[lit_2_leg_index][0]):
                        v = self._tracked_multiply(v,
                            self._representative(self.legs[lit_2_leg_index][0]))
                    else:
                        v = self._tracked_multiply(v,
                            self._representative(self.legs[lit_2_leg_index][1]
                            ) @ self._representative(self.legs[lit_2_leg_index][0]))
                v = self._tracked_multiply(v,
                    self._representative(self.legs[-1][0]) @ self._representative(
                    self.legs[0][0]))
            else:
                if m == 0:
                    v = self._tracked_multiply(v, self._representative(self.legs[-1][0]))
                else:
                    for i in range(m, -1, -1):
                        v = self._tracked_multiply(v, self._representative(self.legs[-1][i]))
        # Now handle all legs of length 2
        l1b_is_lit = self._is_lit(v, self.legs[-1][0])
        for i in range(len(self.legs) - 1):
            if len(self.legs[i]) < 2:
                continue
            if len(self.legs[i]) > 2:
                break
            if not self._is_lit(v, self.legs[i][0]) and not self._is_lit(v, self.legs[i][1]):
                continue
            if self._is_lit(v, self.legs[i][0]) and not self._is_lit(v, self.legs[i][1]):
                v = self._tracked_multiply(v, self._representative(self.legs[i][0]))
            elif not self._is_lit(v, self.legs[i][0]):
                v = self._tracked_multiply(v, self._representative(self.legs[i][1]))
            if not self._is_lit(v, self.central_vertex):
                if not l1b_is_lit:
                    v = self._tracked_multiply(v, self._representative(self.legs[-1][1]))
                    l1b_is_lit = True
                v = self._tracked_multiply(v, self._representative(self.legs[-1][0]))
            v = self._tracked_multiply(v,
                self._representative(self.legs[i][1]) @ self._representative(
                self.legs[i][0]) @ self._representative(self.legs[0][0]))
        self._notify("Transfer lightning to long leg", lighting=v)
        return v

    def _reduce_lightning(self, vertex_stack: list[PauliString],
        v: PauliString)->PauliString:
        """
        Reduce lighting on long leg
        Args:
            vertex_stack (list[PauliString]): Generator stack
            v (PauliString): Append Pauli string
        Returns:
            PauliString: Pauli string after reduce
        """
        if not self._is_lit(v, self.central_vertex):
            if len(self.legs[-2]) == 1 and len(self.legs[-1]) > 4:
                #self._debug(v)
                if self._is_lit(v, self.legs[-1][-1]):
                    is_append_to_end = True
                    for i in range(len(self.legs[-1]) - 2, -1, -1):
                        if self._is_lit(v, self.legs[-1][i]):
                            is_append_to_end = False
                            break
                    if is_append_to_end:
                        self.legs[-1].append(v)
                        self._notify("Append to leg end", appending=v)
                        return v
            m = None
            for i in range(len(self.legs[-1])):
                if self._is_lit(v, self.legs[-1][i]):
                    m = i
                    break
            for i in range(m, -1, -1):
                v = self._tracked_multiply(v, self._representative(self.legs[-1][i]))
        # Now we need to reduce the lit vertices on the long leg to one position and
        #a list of contractions
        f, s = None, None
        for i in range(len(self.legs[-1])):
            if self._is_lit(v, self.legs[-1][i]):
                if f is None:
                    f = i
                elif s is None:
                    s = i
                    break
        # Exit if no element of longest leg is lit
        if f is None and s is None:
            self.legs.append([v])
            self._notify("Attach as new leg", appending=v)
            return v
        # Otherwise naively reduce until one element is left
        if s is not None:
            self._notify("Reduce lightning on long leg", lighting=v,
                contracting=self.legs[-1][f])
            # Compute prefix products on the leg to perform operations in O(1) and O(n) overall
            pref = self._representative(self.legs[-1][f])
            for i in range(f, s):
                pref = pref @ self._representative(self.legs[-1][i])
            while s < len(self.legs[-1]):
                pref = pref @ self._representative(self.legs[-1][s])
                v = self._tracked_multiply(v, pref)
                f += 1
                s += 1
                pref = pref @ self._representative(self.legs[-1][f])
                while s < len(self.legs[-1]) and not self._is_lit(v, self.legs[-1][s]):
                    pref = pref @ self._representative(self.legs[-1][s])
                    s += 1
        if f == 0:
            # Case 1: First vertex of long leg is lit
            if self.type == 'B' and len(self.legs[-1]) == 4:
                # Graph is of type B2
                v = self._tracked_multiply(v,
                    self._representative(self.legs[-1][1]) @ self._representative(
                    self.legs[-1][3]))
                self.legs.append([v])
                self._notify("First vertex lit, type B2 split", appending=v)
            else:
                for w in self.legs[-1]:
                    v = self._tracked_multiply(v, self._representative(w))
                self.legs[-1].append(v)
                self._notify("First vertex lit, extend long leg", appending=v)
        else:
            # Now we have to do careful case handling based on the type of the graph
            # Here f is either the middle or last vertex
            # If it is an A type graph, it may or may not become B type after this
            # First we break the legs
            self.legs[-1][f - 1] = self._tracked_multiply(
                self.legs[-1][f - 1], self._representative(
                v) @ self._representative(self.legs[0][0]))
            subleg = self.legs[-1][f:]
            self.legs[-1] = self.legs[-1][:f]
            self.legs.append([v] + subleg)
            self._notify("Break long leg at lit vertex", appending=v,
                replacing=[self.legs[-2][f - 1]])
            # Now if the graph was type B, then nothing else has to be done
            # Let's order the legs in increasing size
            if len(self.legs[-1]) < len(self.legs[-2]):
                self.legs[-1], self.legs[-2] = self.legs[-2], self.legs[-1]
            # If the graph was type A and we violate the conditions,
            # we need to remove vertices from legs
            # Remove more from the larger leg and less from the smaller leg
            # This happens at most once
            if self.type == 'A' and len(self.legs[-1]) >= 2 and len(self.legs[-2]) >= 2:
                self.type = 'B'
                removed: list[PauliString] = []
                while len(self.legs[-1]) > 4:
                    removed.append(self.legs[-1].pop())
                    vertex_stack.append(removed[-1])
                while len(self.legs[-2]) > 2:
                    removed.append(self.legs[-2].pop())
                    vertex_stack.append(removed[-1])
                if removed:
                    self._notify("Type becomes B, trim legs", removing=removed)
        return v

    def _dependency_check(self, length_1_legs: list[list[PauliString]]) -> None:
        """
        Check depending on legs long 1
        Args:
            length_1_legs (list[list[PauliString]]): List of legs length 1
        Returns:
            list[list[PauliString]]: list of independent legs
        """
        # We need to do Gaussian elimination on the legs of length 1
        independent_legs = []
        basis: dict[int, PauliString] = {}
        for leg in length_1_legs:
            p = leg[0].copy()
            while True:
                pos = p.bits.find(1)
                if pos == -1:
                    break
                if pos not in basis:
                    basis[pos] = p
                    independent_legs.append(leg)
                    break
                p = p @ basis[pos]
        return independent_legs

    def _connected_canonical_graph(self, vertex_stack: list[PauliString]) -> None:
        """
        Constructing a canonical graph from a stack
        Args:
            vertex_stack (list[PauliString]): Generator stack
        """
        while vertex_stack:
            confirmed_legs = [leg for leg in self.legs if len(leg) != 1]
            length_1_legs = [leg for leg in self.legs if len(leg) == 1]
            independent_legs = self._dependency_check(length_1_legs)
            if self.events.has_subscribers():
                independent_set = {leg[0] for leg in independent_legs}
                dropped = [leg[0] for leg in length_1_legs if leg[0] not in independent_set]
                if dropped:
                    self._notify("Remove dependent vertices",
                        dependent=dropped[0], removing=dropped)
            confirmed_legs.extend(independent_legs)
            self.legs = confirmed_legs
            self.legs.sort(key=len)
            v = vertex_stack.pop()
            if self.central_vertex is None:
                self.central_vertex = v
                self._notify("Central vertex")
                continue
            self._notify("Lighting vertex", lighting=v)
            # Build the core
            if len(self.legs) < 2:
                v = self._build_core(v)
                continue
            # Check if there are legs of length 1 with different lit states
            lit_index, unlit_index = None, None
            one_indexes = [i for i, leg in enumerate(self.legs) if len(leg) == 1]

            for i in reversed(one_indexes):
                if lit_index is not None and unlit_index is not None:
                    break
                if len(self.legs[i]) > 1:
                    break
                if self._is_lit(v, self.legs[i][0]):
                    lit_index = i
                else:
                    unlit_index = i
            if lit_index is not None and unlit_index is not None:
                self._convert_to_single_lit_state(lit_index, unlit_index, vertex_stack, v)
                continue
            # From here on we work with self.legs[0] as the representative of length 1 legs WLOG
            # We need to handle a special case, if v is only connected to the central vertex then
            # we just connect it and exit
            if self._is_lit(v, self.legs[0][0]):
                if not self._is_lit(v, self.central_vertex):
                    v = self._tracked_multiply(v, self._representative(self.legs[0][0]))
                v = self._tracked_multiply(v, self._representative(self.central_vertex))
            any_lit_leg = False

            for leg in reversed(self.legs):
                if len(leg) == 1:
                    continue
                if any_lit_leg:
                    break
                for w in reversed(leg):
                    if self._is_lit(v, w):
                        any_lit_leg = True
                        break
            if not any_lit_leg:
                self.legs.append([v])
                self._notify("Connect to centre only", appending=v)
                continue

            if self.type == 'B':
                # Check if there is a lit vertex in a leg of length 2
                lit_2_leg_index = None
                for i in range(len(self.legs) - 1):
                    if len(self.legs[i]) < 2:
                        continue
                    if len(self.legs[i]) > 2:
                        break
                    if self._is_lit(v, self.legs[i][0]) or self._is_lit(v, self.legs[i][1]):
                        lit_2_leg_index = i
                        break
                if lit_2_leg_index is not None:
                    v = self._transfer_lightning(lit_2_leg_index, v)
            v = self._reduce_lightning(vertex_stack, v)

        confirmed_legs = [leg for leg in self.legs if len(leg) != 1]
        length_1_legs = [leg for leg in self.legs if len(leg) == 1]
        independent_legs = self._dependency_check(length_1_legs)
        if self.events.has_subscribers():
            independent_set = {leg[0] for leg in independent_legs}
            dropped = [leg[0] for leg in length_1_legs if leg[0] not in independent_set]
            if dropped:
                self._notify("Remove dependent vertices",
                    dependent=dropped[0], removing=dropped)
        confirmed_legs.extend(independent_legs)
        self.legs = confirmed_legs
        self.legs.sort(key=len)

    def _get_morph(self) -> Morph:
        """
        Get the canonical graph built by the canonicalizer.

        Returns:
            Morph: Canonical graph.
        """
        legs = self.legs.copy()
        legs.insert(0, [self.central_vertex])

        return Morph(legs, [], self.vertex_stack)

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> Morph:
        """
        Build a canonical graph from a stack of connected generators

        Args:
            vertex_stack (list[PauliString]): Generator stack
        """
        self._set_vertex_stack(vertex_stack)
        self._notify("Anticommutation graph", collection=vertex_stack.copy(), init=True)
        self._connected_canonical_graph(vertex_stack)
        morph = self._get_morph()
        if self.events.has_subscribers():
            try:
                self._notify(f"Canonical graph of type {morph.get_type().name}")
            except Exception:  # pylint: disable=broad-except
                self._notify("Canonical graph")
        return morph
