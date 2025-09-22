"""
Optimal Pauli Compiler (Algorithm from arXiv:2408.03294) implemented using PauLie types.

API overview:
- OptimalPauliCompiler(k_left: int, n_total: int).compile(V_left: PauliString, W_right: PauliString) -> list[PauliString]
  Returns a sequence G in PauLie orientation: nested_adjoint(G[:-1], G[-1]) = V ⊗ W.
- compile_target(target: PauliString, k_left: int) -> list[PauliString]

This module does not replace the existing simple pauli_compiler; it provides an
algorithmically faithful implementation alongside it.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from collections import deque
from typing import Iterable

from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import get_identity, get_pauli_string


def _single(n: int, i: int, label: str) -> PauliString:
    p = get_identity(n)
    p[i] = label
    return p

def _tensor(left: PauliString, right: PauliString) -> PauliString:
    return left.tensor(right)

def _multiply(a: PauliString, b: PauliString) -> PauliString:
    return a @ b

def _commutes(a: PauliString, b: PauliString) -> bool:
    return a | b

def _left_part(p: PauliString, k: int) -> PauliString:
    return p.get_substring(0, k)

def _right_part(p: PauliString, k: int) -> PauliString:
    return p.get_substring(k, len(p) - k)

def _ad_apply(A: PauliString, B: PauliString | None) -> PauliString | None:
    if B is None:
        return None
    return None if _commutes(A, B) else _multiply(A, B)

def _nested_commutator_result(G: list[PauliString]) -> PauliString | None:
    if not G:
        return None
    cur: PauliString | None = G[0]
    for A in G[1:]:
        cur = _ad_apply(A, cur)
        if cur is None:
            return None
    return cur

def _sequence_to_paulie_orientation(G: list[PauliString]) -> list[PauliString]:
    # PauLie nested_adjoint iterates operators in reversed order.
    # Our G is [base, A1, A2, ..., Am] where main.py applies A1 then A2 ... then Am.
    # To make nested_adjoint(seq[:-1], seq[-1]) apply in the same order, we must
    # return operators reversed: [Am, ..., A2, A1, base].
    if not G:
        return []
    base = G[0]
    ops = G[1:]
    ops_rev = list(reversed(ops))
    return ops_rev + [base]


def left_A_minimal(k: int) -> list[PauliString]:
    # {Xi, Zi} for i=0..k-1 plus Z...Z on all k
    A: list[PauliString] = []
    for i in range(k):
        A.append(_single(k, i, "X"))
        A.append(_single(k, i, "Z"))
    zall = get_pauli_string("Z" * k)
    A.append(zall)
    return A

def choose_U_for_B(k: int) -> PauliString:
    return _single(k, 0, "X")

def _all_left_paulis(k: int) -> list[PauliString]:
    # Enumerate all non-identity Pauli strings of length k
    out: list[PauliString] = []
    labels = ["I", "X", "Y", "Z"]
    def rec(i: int, cur: list[str]):
        if i == k:
            s = "".join(cur)
            if set(s) != {"I"}:
                out.append(get_pauli_string(s))
            return
        for L in labels:
            cur.append(L)
            rec(i + 1, cur)
            cur.pop()
    rec(0, [])
    return out


@dataclass
class SubsystemCompilerConfig:
    k_left: int
    n_total: int

class SubsystemCompiler:
    def __init__(self, cfg: SubsystemCompilerConfig):
        if cfg.k_left < 2:
            raise ValueError("k_left must be >= 2 for the Pauli Compiler algorithm")
        self.k = cfg.k_left
        self.N = cfg.n_total
        self.n_right = self.N - self.k
        self.U_tag = choose_U_for_B(self.k)
        self.left_pool = _all_left_paulis(self.k)

    def extend_left(self, A_left: PauliString) -> PauliString:
        return _tensor(A_left, get_identity(self.n_right))

    def extend_pair(self, U_left: PauliString, B_right: PauliString) -> PauliString:
        return _tensor(U_left, B_right)

    def factor_W_orders(self, W_right: PauliString) -> list[list[tuple[PauliString, PauliString]]]:
        assert len(W_right) == self.n_right
        per_site_opts: list[list[list[PauliString]]] = []
        wstr = str(W_right)
        for j, ch in enumerate(wstr):
            if ch == "Y":
                opt1 = [_single(self.n_right, j, "X"), _single(self.n_right, j, "Z")]
                opt2 = [_single(self.n_right, j, "Z"), _single(self.n_right, j, "X")]
                per_site_opts.append([opt1, opt2])
            elif ch == "X":
                per_site_opts.append([[ _single(self.n_right, j, "X") ]])
            elif ch == "Z":
                per_site_opts.append([[ _single(self.n_right, j, "Z") ]])
            else:
                per_site_opts.append([[]])
        sequences: list[list[PauliString]] = []
        def rec(i: int, acc: list[PauliString]):
            if i == len(per_site_opts):
                sequences.append(list(acc))
                return
            for seg in per_site_opts[i]:
                acc.extend(seg)
                rec(i + 1, acc)
                del acc[-len(seg):]
        rec(0, [])
        Ui = self.U_tag
        return [[(Ui, b) for b in flat] for flat in sequences]

    def _choose_A1_A2(self, U: PauliString) -> tuple[PauliString, PauliString]:
        for a1 in self.left_pool:
            if _commutes(a1, U):
                continue
            for a2 in self.left_pool:
                if a2 is a1:
                    continue
                if _commutes(a2, U):
                    continue
                if _commutes(a1, a2):
                    return a1, a2
        raise RuntimeError("Failed to find A1,A2 in iP_k^*.")

    def _choose_Aprime(self, Ui: PauliString, P_left: PauliString) -> PauliString:
        for a in self.left_pool:
            if (not _commutes(a, Ui)) and _commutes(a, P_left):
                return a
        raise RuntimeError("Failed to find A' in iP_k^*.")

    def _rest_full_after(self, Ui_Bi: list[tuple[PauliString, PauliString]], i: int, H: list[PauliString]) -> tuple[PauliString, PauliString]:
        P_left = get_identity(self.k)
        for j in range(i+1, len(Ui_Bi)):
            P_left = _multiply(P_left, Ui_Bi[j][0])
        for A in H:
            P_left = _multiply(P_left, A)
        P_right = get_identity(self.n_right)
        for j in range(i+1, len(Ui_Bi)):
            P_right = _multiply(P_right, Ui_Bi[j][1])
        return P_left, P_right

    def subsystem_compiler(self, W_right: PauliString) -> list[PauliString]:
        assert len(W_right) == self.n_right
        for Ui_Bi in self.factor_W_orders(W_right):
            if not Ui_Bi:
                return []
            r = len(Ui_Bi)
            i = r - 1
            G: list[PauliString] = [self.extend_pair(Ui_Bi[-1][0], Ui_Bi[-1][1])]
            H: list[PauliString] = []
            helper_used_for_i: dict[int, int] = {}
            while i >= 1:
                Ui, Bi = Ui_Bi[i]
                P_left, P_right = self._rest_full_after(Ui_Bi, i, H)
                if _multiply(P_left, Ui).is_identity():
                    cnt = helper_used_for_i.get(i, 0)
                    if cnt >= 1:
                        G.append(self.extend_pair(Ui, Bi))
                        i -= 1
                        continue
                    A1, A2 = self._choose_A1_A2(Ui)
                    H = [A1, A2]
                    helper_used_for_i[i] = cnt + 1
                    G.append(self.extend_left(A1))
                    G.append(self.extend_left(A2))
                    continue
                current = self.extend_pair(Ui, Bi)
                rest_full = _tensor(P_left, P_right)
                if _commutes(current, rest_full):
                    cnt = helper_used_for_i.get(i, 0)
                    if cnt >= 1:
                        G.append(self.extend_pair(Ui, Bi))
                        i -= 1
                        continue
                    A_prime = self._choose_Aprime(Ui, P_left)
                    H = [A_prime]
                    helper_used_for_i[i] = cnt + 1
                    G.append(self.extend_left(A_prime))
                    continue
                G.append(current)
                i -= 1
            return G
        return []


def _key(p: PauliString) -> str:
    return str(p)

def left_map_over_A(V_from: PauliString, V_to: PauliString, A: list[PauliString]) -> list[PauliString]:
    if _key(V_from) == _key(V_to):
        return []
    start_k = _key(V_from)
    goal_k = _key(V_to)
    q: deque[PauliString] = deque([V_from])
    parent: dict[str, tuple[str, str, PauliString]] = {}
    seen = {start_k}
    while q:
        cur = q.popleft()
        cur_k = _key(cur)
        if cur_k == goal_k:
            seq: list[PauliString] = []
            while cur_k != start_k:
                prev_k, _, a_used = parent[cur_k]
                seq.append(a_used)
                cur_k = prev_k
            seq.reverse()
            return seq
        for a in A:
            if _commutes(a, cur):
                continue
            nxt = _multiply(a, cur)
            nk = _key(nxt)
            if nk in seen:
                continue
            seen.add(nk)
            parent[nk] = (cur_k, nk, a)
            q.append(nxt)
    raise RuntimeError("Left map BFS failed.")


@dataclass
class PauliCompilerConfig:
    k_left: int
    n_total: int
    fallback_depth: int = 8
    fallback_nodes: int = 200000

class OptimalPauliCompiler:
    def __init__(self, cfg: PauliCompilerConfig):
        if cfg.k_left < 2:
            raise ValueError("k_left must be >= 2 for the Pauli Compiler algorithm")
        self.k = cfg.k_left
        self.N = cfg.n_total
        self.n_right = self.N - self.k
        self.A_left = left_A_minimal(self.k)
        self.U_tag = choose_U_for_B(self.k)
        self.sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=self.k, n_total=self.N))
        self.fallback_depth = cfg.fallback_depth
        self.fallback_nodes = cfg.fallback_nodes

    def extend_left(self, A_left: PauliString) -> PauliString:
        return _tensor(A_left, get_identity(self.n_right))

    def _left_factor_from_sequence(self, G: list[PauliString]) -> PauliString:
        res = _nested_commutator_result(G)
        if res is None:
            return _single(self.k, 0, "X")
        return _left_part(res, self.k)

    def _candidate_decompositions(self, W: PauliString) -> list[tuple[PauliString, PauliString]]:
        cand: list[tuple[PauliString, PauliString]] = []
        wstr = str(W)
        n_right = len(wstr)
        for j, ch in enumerate(wstr):
            if ch == "I":
                continue
            labels = ["X", "Z"] if ch == "Y" else (["Z"] if ch == "X" else ["X"])
            for lab in labels:
                W1 = _single(n_right, j, lab)
                W2 = _multiply(W1, W)
                if not _commutes(W1, W2):
                    cand.append((W1, W2))
        uniq: list[tuple[PauliString, PauliString]] = []
        seen = set()
        for a, b in cand:
            key = (str(a), str(b))
            if key not in seen:
                seen.add(key)
                uniq.append((a, b))
        return uniq

    def _all_interleavings_preserving(self, A: list[PauliString], B: list[PauliString], C: list[PauliString], cap: int = 60000) -> Iterable[list[PauliString]]:
        count = 0
        NA, NB, NC = len(A), len(B), len(C)
        def rec(i: int, j: int, k: int, prefix: list[PauliString]):
            nonlocal count
            if count >= cap:
                return
            if i == NA and j == NB and k == NC:
                count += 1
                yield list(prefix)
                return
            if i < NA:
                prefix.append(A[i])
                yield from rec(i + 1, j, k, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if j < NB:
                prefix.append(B[j])
                yield from rec(i, j + 1, k, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if k < NC:
                prefix.append(C[k])
                yield from rec(i, j, k + 1, prefix)
                prefix.pop()
        return rec(0, 0, 0, [])

    def _all_interleavings_preserving4(self, A: list[PauliString], B: list[PauliString], C: list[PauliString], D: list[PauliString], cap: int = 120000) -> Iterable[list[PauliString]]:
        count = 0
        NA, NB, NC, ND = len(A), len(B), len(C), len(D)
        def rec(i: int, j: int, k: int, l_idx: int, prefix: list[PauliString]):
            nonlocal count
            if count >= cap:
                return
            if i == NA and j == NB and k == NC and l_idx == ND:
                count += 1
                yield list(prefix)
                return
            if i < NA:
                prefix.append(A[i])
                yield from rec(i + 1, j, k, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if j < NB:
                prefix.append(B[j])
                yield from rec(i, j + 1, k, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if k < NC:
                prefix.append(C[k])
                yield from rec(i, j, k + 1, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if l_idx < ND:
                prefix.append(D[l_idx])
                yield from rec(i, j, k, l_idx + 1, prefix)
                prefix.pop()
        return rec(0, 0, 0, 0, [])

    def _case3_best_reordering(self, G1: list[PauliString], G2: list[PauliString], Aext: list[PauliString], W: PauliString) -> list[PauliString] | None:
        k = self.k
        R = list
        G1r, G2r = list(reversed(G1)), list(reversed(G2))
        candidates = [
            R(G1) + R(Aext) + R(G2r) + R(Aext),
            R(G1r) + R(Aext) + R(G2)  + R(Aext),
            R(Aext) + R(G1) + R(Aext) + R(G2r),
            R(Aext) + R(G1r) + R(Aext) + R(G2),
            R(G2)  + R(Aext) + R(G1r) + R(Aext),
            R(G2r) + R(Aext) + R(G1)  + R(Aext),
        ]
        for seq in candidates:
            res = _nested_commutator_result(seq)
            if res is not None and str(_left_part(res, k)) == "I"*k and str(_right_part(res, k)) == str(W):
                return seq
        blocks = [G1, G2, Aext]
        for perm in permutations(range(3)):
            B0, B1, B2 = [blocks[i] for i in perm]
            for r0 in (False, True):
                for r1 in (False, True):
                    for r2 in (False, True):
                        seq = (list(reversed(B0)) if r0 else list(B0)) \
                            + (list(reversed(B1)) if r1 else list(B1)) \
                            + (list(reversed(B2)) if r2 else list(B2))
                        res = _nested_commutator_result(seq)
                        if res is not None and str(_left_part(res, k)) == "I"*k and str(_right_part(res, k)) == str(W):
                            return seq
        for g1 in (G1, list(reversed(G1))):
            for g2 in (G2, list(reversed(G2))):
                for a in (Aext, list(reversed(Aext))):
                    for seq in self._all_interleavings_preserving(g1, g2, a):
                        res = _nested_commutator_result(seq)
                        if res is not None and str(_left_part(res, k)) == "I"*k and str(_right_part(res, k)) == str(W):
                            return seq
        Aopts = (Aext, list(reversed(Aext)))
        for g1 in (G1, G1r):
            for g2 in (G2, G2r):
                for a1 in Aopts:
                    for a2 in Aopts:
                        for seq in (
                            list(g1)+list(a1)+list(g2)+list(a2),
                            list(g2)+list(a1)+list(g1)+list(a2),
                            list(a1)+list(g1)+list(a2)+list(g2),
                            list(a1)+list(g2)+list(a2)+list(g1),
                        ):
                            res = _nested_commutator_result(seq)
                            if res is not None and str(_left_part(res, k)) == "I"*k and str(_right_part(res, k)) == str(W):
                                return seq
                        for seq in self._all_interleavings_preserving4(g1, g2, a1, a2):
                            res = _nested_commutator_result(seq)
                            if res is not None and str(_left_part(res, k)) == "I"*k and str(_right_part(res, k)) == str(W):
                                return seq
        return None

    def _bfs_case3(self, W: PauliString, depth_cap: int, node_cap: int) -> list[PauliString] | None:
        S = construct_universal_set(self.N, self.k)
        target_left = get_identity(self.k)
        target_right = W
        def key_of(p: PauliString | None) -> str | None:
            return None if p is None else str(p)
        frontier: list[tuple[PauliString | None, list[int]]] = [(None, [])]
        visited: dict[tuple[int, str], bool] = {}
        nodes = 0
        for depth in range(1, depth_cap+1):
            new_front: list[tuple[PauliString | None, list[int]]] = []
            for res, seq_idx in frontier:
                for j, op in enumerate(S):
                    nodes += 1
                    if nodes > node_cap:
                        return None
                    new_res = op if res is None else _ad_apply(op, res)
                    kres = key_of(new_res)
                    if new_res is None:
                        continue
                    if visited.get((depth, kres), False):
                        continue
                    visited[(depth, kres)] = True
                    new_seq = seq_idx + [j]
                    if depth >= 2:
                        L = _left_part(new_res, self.k)
                        R = _right_part(new_res, self.k)
                        if str(L) == str(target_left) and str(R) == str(target_right):
                            return [S[idx] for idx in new_seq]
                    new_front.append((new_res, new_seq))
            frontier = new_front
        return None

    def compile(self, V_left: PauliString, W_right: PauliString) -> list[PauliString]:
        assert len(V_left) == self.k and len(W_right) == self.n_right
        V, W = V_left, W_right
        if W.is_identity():
            Aset = left_A_minimal(self.k)
            for As in Aset:
                try:
                    seqA = left_map_over_A(As, V, Aset)
                except RuntimeError:
                    continue
                G = [self.extend_left(As)] + [self.extend_left(a) for a in seqA]
                res = _nested_commutator_result(G)
                if res is not None and str(_left_part(res, self.k)) == str(V) and str(_right_part(res, self.k)) == str(W):
                    return _sequence_to_paulie_orientation(G)
            raise RuntimeError("Left-only mapping failed.")
        if not V.is_identity():
            Gp = self.sub.subsystem_compiler(W)
            Vp = self._left_factor_from_sequence(Gp)
            Aset = left_A_minimal(self.k)
            seq = left_map_over_A(Vp, V, Aset)
            candidates = [
                list(Gp) + [self.extend_left(a) for a in seq],
                [self.extend_left(a) for a in seq] + list(Gp),
                list(reversed(Gp)) + [self.extend_left(a) for a in seq],
            ]
            for G in candidates:
                res = _nested_commutator_result(G)
                if res is not None and str(_left_part(res, self.k)) == str(V) and str(_right_part(res, self.k)) == str(W):
                    return _sequence_to_paulie_orientation(G)
            return _sequence_to_paulie_orientation(list(Gp) + [self.extend_left(a) for a in seq])
        Aset = left_A_minimal(self.k)
        for W1, W2 in self._candidate_decompositions(W):
            G1 = self.sub.subsystem_compiler(W1)
            G2 = self.sub.subsystem_compiler(W2)
            V1p = self._left_factor_from_sequence(G1)
            V2p = self._left_factor_from_sequence(G2)
            Aseq = left_map_over_A(V2p, V1p, Aset)
            Aext = [self.extend_left(a) for a in Aseq]
            seq = self._case3_best_reordering(G1, G2, Aext, W)
            if seq is not None:
                return _sequence_to_paulie_orientation(seq)
        seq_fb = self._bfs_case3(W, self.fallback_depth, self.fallback_nodes)
        if seq_fb is not None:
            return _sequence_to_paulie_orientation(seq_fb)
        wstr = str(W)
        j = next(i for i, ch in enumerate(wstr) if ch != "I")
        lab = "X" if wstr[j] == "Z" else ("Z" if wstr[j] == "X" else "X")
        W1 = _single(self.n_right, j, lab)
        W2 = _multiply(W1, W)
        G1 = self.sub.subsystem_compiler(W1)
        G2 = self.sub.subsystem_compiler(W2)
        V1p = self._left_factor_from_sequence(G1)
        V2p = self._left_factor_from_sequence(G2)
        Aseq = left_map_over_A(V2p, V1p, left_A_minimal(self.k))
        Aext = [self.extend_left(a) for a in Aseq]
        return _sequence_to_paulie_orientation(list(reversed(G1)) + Aext + list(reversed(G2)))


def construct_universal_set(N: int, k: int) -> list[PauliString]:
    if not (1 <= k < N):
        raise ValueError("Require 1 <= k < N")
    A_k = left_A_minimal(k)
    n_right = N - k
    U = choose_U_for_B(k)
    right_B = [_single(n_right, j, "X") for j in range(n_right)] + [_single(n_right, j, "Z") for j in range(n_right)]
    A_prime = [_tensor(A, get_identity(n_right)) for A in A_k]
    B_prime = [_tensor(U, b) for b in right_B]
    return A_prime + B_prime


def compile_target(target: PauliString, k_left: int) -> list[PauliString]:
    n = len(target)
    if not (1 <= k_left < n):
        raise ValueError("Require 1 <= k_left < len(target)")
    V = target.get_substring(0, k_left)
    W = target.get_substring(k_left, n - k_left)
    opc = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n))
    return opc.compile(V, W)
