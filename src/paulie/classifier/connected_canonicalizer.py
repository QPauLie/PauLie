from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.classification import Morph

class ConnectedCanonicalizer:
    def __init__(self):
        self.type = 'A' # B when at least two legs of length at least 2
        self.central_vertex = None
        self.legs = []
        self.vertex_stack = []

    def _set_vertex_stack(self, vertex_stack):
        self.vertex_stack = vertex_stack.copy()

    def _representative(self, v: PauliString):
        return v

    def _is_lit(self, v: PauliString, w: PauliString):
        return self._representative(v) ^ self._representative(w) is not None

    def _tracked_multiply(self, v: PauliString, w: PauliString):
        return v @ w

    def _build_core(self, v: PauliString) -> PauliString:
        if len(self.legs) == 1:
            if self._is_lit(v, self.legs[0][0]):
                if not self._is_lit(v, self.central_vertex):
                    v = self._tracked_multiply(v, self._representative(self.legs[0][0]))
                v = self._tracked_multiply(v, self._representative(self.central_vertex))
        self.legs.append([v])
        return v

    def _convert_to_single_lit_state(self, p_index: int, q_index: int, vertex_stack: list[PauliString], v: PauliString):
        pq = self._representative(self.legs[p_index][0]) @ self._representative(self.legs[q_index][0])
        if self._is_lit(v, self.central_vertex):
            self.central_vertex = self._tracked_multiply(self.central_vertex, pq)
        for i in range(len(self.legs)):
            for j in range(len(self.legs[i])):
                if i != p_index and self._is_lit(v, self.legs[i][j]):
                    self.legs[i][j] = self._tracked_multiply(self.legs[i][j], pq)
        self.legs[p_index].append(v)
        # Truncate longest leg if necessary, this happens at most once
        big_leg_cnt = sum(1 for leg in self.legs if len(leg) >= 2)
        if self.type == 'A' and big_leg_cnt >= 2:
            self.type = 'B'
            while len(self.legs[-1]) > 4:
                vertex_stack.append(self.legs[-1].pop())

    def _transfer_lightning(self, lit_2_leg_index: int, v: PauliString):
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
                        v = self._tracked_multiply(v, self._representative(self.legs[lit_2_leg_index][0]))
                    else:
                        v = self._tracked_multiply(v, self._representative(self.legs[lit_2_leg_index][1]) @ self._representative(self.legs[lit_2_leg_index][0]))
                v = self._tracked_multiply(v, self._representative(self.legs[-1][0]) @ self._representative(self.legs[0][0]))
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
            v = self._tracked_multiply(v, self._representative(self.legs[i][1]) @ self._representative(self.legs[i][0]) @ self._representative(self.legs[0][0]))
        return v

    def _reduce_lightning(self, vertex_stack: list[PauliString], v: PauliString):
        if not self._is_lit(v, self.central_vertex):
            m = None
            for i in range(len(self.legs[-1])):
                if self._is_lit(v, self.legs[-1][i]):
                    m = i
                    break
            for i in range(m, -1, -1):
                v = self._tracked_multiply(v, self._representative(self.legs[-1][i]))
        # Now we need to reduce the lit vertices on the long leg to one position and a list of contractions
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
            return v
        # Otherwise naively reduce until one element is left
        if s is not None:
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
                v = self._tracked_multiply(v, self._representative(self.legs[-1][1]) @ self._representative(self.legs[-1][3]))
                self.legs.append([v])
            else:
                for w in self.legs[-1]:
                    v = self._tracked_multiply(v, self._representative(w))
                self.legs[-1].append(v)
        else:
            # Now we have to do careful case handling based on the type of the graph
            # Here f is either the middle or last vertex
            # If it is an A type graph, it may or may not become B type after this
            # First we break the legs
            self.legs[-1][f - 1] = self._tracked_multiply(self.legs[-1][f - 1], self._representative(v) @ self._representative(self.legs[0][0]))
            subleg = self.legs[-1][f:]
            self.legs[-1] = self.legs[-1][:f]
            self.legs.append([v] + subleg)
            # Now if the graph was type B, then nothing else has to be done
            # Let's order the legs in increasing size
            if len(self.legs[-1]) < len(self.legs[-2]):
                self.legs[-1], self.legs[-2] = self.legs[-2], self.legs[-1]
            # If the graph was type A and we violate the conditions, we need to remove vertices from legs
            # Remove more from the larger leg and less from the smaller leg
            # This happens at most once
            if self.type == 'A' and len(self.legs[-1]) >= 2 and len(self.legs[-2]) >= 2:
                self.type = 'B'
                while len(self.legs[-1]) > 4:
                    vertex_stack.append(self.legs[-1].pop())
                while len(self.legs[-2]) > 2:
                    vertex_stack.append(self.legs[-2].pop())
        return v
    
    def _dependency_check(self, length_1_legs: list[list[PauliString]]):
        # We need to do Gaussian elimination on the legs of length 1
        independent_legs = []
        basis: dict[int, PauliString] = {}
        for leg in length_1_legs:
            p = leg[0].copy()
            while True:
                if p.bits.find(1) == -1:
                    break
                x = basis.get(p.bits.find(1))
                if x is None:
                    basis[p.bits.find(1)] = p
                    independent_legs.append(leg)
                    break
                p = p @ x
        return independent_legs

    def _connected_canonical_graph(self, vertex_stack: list[PauliString]) -> None:
        while vertex_stack:
            v = vertex_stack.pop()
            # Don't forget to sort self.legs by length before accessing them!
            self.legs.sort(key=len)
            if self.central_vertex is None:
                self.central_vertex = v
                continue
            # Build the core
            if len(self.legs) < 2:
                v = self._build_core(v)
                continue
            # Check if there are legs of length 1 with different lit states
            lit_index, unlit_index = None, None
            for i in range(len(self.legs)):
                if len(self.legs[i]) > 1:
                    break
                if lit_index is not None and unlit_index is not None:
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
            for leg in self.legs:
                if len(leg) == 1:
                    continue
                if any_lit_leg:
                    break
                for w in leg:
                    if self._is_lit(v, w):
                        any_lit_leg = True
                        break
            if not any_lit_leg:
                self.legs.append([v])
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
        confirmed_legs.extend(self._dependency_check(length_1_legs))
        self.legs = confirmed_legs
        self.legs.sort(key=len)

    def build_canonical_graph(self, vertex_stack: list[PauliString]) -> None:
        self.__init__()
        self._set_vertex_stack(vertex_stack)
        self._connected_canonical_graph(vertex_stack)

    def get_morph(self) -> Morph:
        """
        Get the canonical graph built by the canonicalizer.

        Returns:
            Morph: Canonical graph.
        """
        legs = self.legs.copy()
        legs.insert(0, [self.central_vertex])

        return Morph(legs, [], self.vertex_stack)
