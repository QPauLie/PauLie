from paulie import get_pauli_string
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_bitarray import PauliString
from paulie.classifier.classification import Classification, Morph
import networkx as nx

class ConnectedClassifier:
    def __init__(self):
        self.type = 'A' # B when at least two legs of length at least 2
        self.central_vertex: PauliString = None
        self.legs: list[list[PauliString]] = []

    def is_lit(self, v: PauliString, w: PauliString):
        return v ^ w is not None

    def build_core(self, v: PauliString) -> PauliString:
        if len(self.legs) == 1:
            if self.is_lit(v, self.legs[0][0]):
                if not self.is_lit(v, self.central_vertex):
                    v = v @ self.legs[0][0]
                v = v @ self.central_vertex
        self.legs.append([v])
        return v

    def convert_to_single_lit_state(self, p_index: int, q_index: int, vertex_stack: list[PauliString], v: PauliString):
        pq = self.legs[p_index][0] @ self.legs[q_index][0]
        if self.is_lit(v, self.central_vertex):
            self.central_vertex = pq @ self.central_vertex
        for i in range(len(self.legs)):
            for j in range(len(self.legs[i])):
                if self.is_lit(v, self.legs[i][j]) and i != p_index:
                    self.legs[i][j] = pq @ self.legs[i][j]
        self.legs[p_index].append(v)
        # Truncate longest leg if necessary, this happens at most once
        big_leg_cnt = sum(1 for leg in self.legs if len(leg) >= 2)
        if self.type == 'A' and big_leg_cnt >= 2:
            self.type = 'B'
            while len(self.legs[-1]) > 4:
                vertex_stack.append(self.legs[-1].pop())

    def transfer_lightning(self, lit_2_leg_index: int, v: PauliString):
        # We need to make self.legs[-1][1] lit
        if not self.is_lit(v, self.legs[-1][1]):
            m = None
            for i in range(len(self.legs[-1])):
                if self.is_lit(v, self.legs[-1][i]):
                    m = i
                    break
            if m is None:
                if not self.is_lit(v, self.central_vertex):
                    if self.is_lit(v, self.legs[lit_2_leg_index][0]):
                        v = v @ self.legs[lit_2_leg_index][0]
                    else:
                        v = v @ self.legs[lit_2_leg_index][1] @ self.legs[lit_2_leg_index][0]
                v = v @ self.legs[-1][0] @ self.legs[0][0]
            else:
                if m == 0:
                    v = v @ self.legs[-1][0]
                else:
                    for i in range(m, 1, -1):
                        v = v @ self.legs[-1][i]
        # Now handle all legs of length 2
        for i in range(len(self.legs) - 1):
            if len(self.legs[i]) < 2:
                continue
            if len(self.legs[i]) > 2:
                break
            if not self.is_lit(v, self.legs[i][0]) and not self.is_lit(v, self.legs[i][1]):
                continue
            if self.is_lit(v, self.legs[i][0]) and not self.is_lit(v, self.legs[i][1]):
                v = v @ self.legs[i][0]
            elif not self.is_lit(v, self.legs[i][0]) and self.is_lit(v, self.legs[i][1]):
                v = v @ self.legs[i][1]
            if not self.is_lit(v, self.central_vertex):
                if not self.is_lit(v, self.legs[-1][0]):
                    v = v @ self.legs[-1][1]
                v = v @ self.legs[-1][0]
            v = v @ self.legs[i][1] @ self.legs[i][0] @ self.legs[0][0]
        return v

    def reduce_lightning(self, vertex_stack: list[PauliString], v: PauliString):
        if not self.is_lit(v, self.central_vertex):
            m = None
            for i in range(len(self.legs[-1])):
                if self.is_lit(v, self.legs[-1][i]):
                    m = i
                    break
            for i in range(m, -1, -1):
                v = v @ self.legs[-1][i]
        # Now we need to reduce the lit vertices on the long leg to one position and a list of contractions
        f, s = None, None
        for i in range(len(self.legs[-1])):
            if self.is_lit(v, self.legs[-1][i]):
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
            pref = self.legs[-1][f]
            for i in range(f, s):
                pref = pref @ self.legs[-1][i]
            while s < len(self.legs[-1]):
                pref = pref @ self.legs[-1][s]
                v = v @ pref
                f += 1
                s += 1
                pref = pref @ self.legs[-1][f]
                while s < len(self.legs[-1]) and not self.is_lit(v, self.legs[-1][s]):
                    pref = pref @ self.legs[-1][s]
                    s += 1
        if f == 0:
            # Case 1: First vertex of long leg is lit
            if self.type == 'B' and len(self.legs[-1]) == 4:
                # Graph is of type B2
                v = v @ self.legs[-1][1] @ self.legs[-1][3]
                self.legs.append([v])
            else:
                for w in self.legs[-1]:
                    v = v @ w
                self.legs[-1].append(v)
        else:
            # Now we have to do careful case handling based on the type of the graph
            # Here f is either the middle or last vertex
            # If it is an A type graph, it may or may not become B type after this
            # First we break the legs
            self.legs[-1][f - 1] = self.legs[-1][f - 1] @ v @ self.legs[0][0]
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
    
    def dependency_check(self):
        # We need to do Gaussian elimination on the legs of length 1
        confirmed_legs = [leg for leg in self.legs if len(leg) != 1]
        length_1_legs = [leg for leg in self.legs if len(leg) == 1]
        basis: dict[int, PauliString] = {}
        for leg in length_1_legs:
            p = leg[0].copy()
            while True:
                if p.bits.find(1) == -1:
                    break
                x = basis.get(p.bits.find(1))
                if x is None:
                    basis[p.bits.find(1)] = p
                    confirmed_legs.append(leg)
                    break
                p = p @ x
        self.legs = confirmed_legs

    def connected_canonical_graph(self, vertex_stack: list[PauliString]):
        self.__init__()
        while vertex_stack:
            v = vertex_stack.pop()
            # Don't forget to sort self.legs by length before accessing them!
            self.legs.sort(key=len)
            if self.central_vertex is None:
                self.central_vertex = v
                continue
            # Build the core
            if len(self.legs) < 2:
                v = self.build_core(v)
                continue
            # Check if there are legs of length 1 with different lit states
            lit_index, unlit_index = None, None
            for i in range(len(self.legs)):
                if len(self.legs[i]) > 1:
                    break
                if lit_index is not None and unlit_index is not None:
                    break
                if self.is_lit(v, self.legs[i][0]):
                    lit_index = i
                else:
                    unlit_index = i
            if lit_index is not None and unlit_index is not None:
                self.convert_to_single_lit_state(lit_index, unlit_index, vertex_stack, v)
                continue
            # From here on we work with self.legs[0] as the representative of length 1 legs WLOG
            # We need to handle a special case, if v is only connected to the central vertex then
            # we just connect it and exit
            if self.is_lit(v, self.legs[0][0]):
                if not self.is_lit(v, self.central_vertex):
                    v = v @ self.legs[0][0]
                v = v @ self.central_vertex
            any_lit_leg = False
            for leg in self.legs:
                for w in leg:
                    if self.is_lit(v, w):
                        any_lit_leg = True
                        break
                if any_lit_leg:
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
                    if self.is_lit(v, self.legs[i][0]) or self.is_lit(v, self.legs[i][1]):
                        lit_2_leg_index = i
                        break
                if lit_2_leg_index is not None:
                    v = self.transfer_lightning(lit_2_leg_index, v)
            v = self.reduce_lightning(vertex_stack, v)
        self.dependency_check()
        self.legs.sort(key=len)
        return self.central_vertex, self.legs

def canonical_graph(gens: PauliStringCollection):
    verts, edges, _ = gens.get_graph()
    g = nx.Graph()
    g.add_nodes_from(verts)
    g.add_edges_from(edges)
    ccs = nx.connected_components(g)
    classif = ConnectedClassifier()
    classification = Classification()
    for cc in ccs:
        vertex_stack = [get_pauli_string(s) for s in nx.dfs_preorder_nodes(g.subgraph(cc))]
        vertex_stack.reverse()
        central_vertex, legs = classif.connected_canonical_graph(vertex_stack)
        legs.insert(0, [central_vertex])
        classification.add(Morph(legs, []))
    return classification
