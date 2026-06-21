"""Pauli compiler based on arXiv:2408.03294.
``compiler`` takes a generator set and a target Pauli string and outputs a
:math:`\\mathcal{O}(N)` length sequence of Pauli strings that generates the target
Pauli string via nested commutators.
"""
from collections import deque
from collections.abc import Generator
from dataclasses import dataclass
from typing import Iterable

from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_factory import get_identity, get_pauli_string, get_single

def _evaluate_sequence(sequence: list[PauliString]) -> PauliString | None:
    """Evaluate a sequence stored as ``[base, A1, ..., Am]``."""
    return PauliStringCollection(sequence).evaluate_commutator_sequence()

def _evaluate_paulie_orientation(sequence: list[PauliString]) -> PauliString | None:
    """Evaluate a sequence stored as ``[Am, ..., A1, base]``."""
    if not sequence:
        return None
    return PauliStringCollection(sequence[:-1]).nested_adjoint(sequence[-1])

def _sequence_to_paulie_orientation(sequence: list[PauliString]) -> list[PauliString]:
    """Convert ``[base, A1, ..., Am]`` to PauLie's ``nested_adjoint`` orientation."""
    if not sequence:
        return []
    return list(reversed(sequence[1:])) + [sequence[0]]

def left_a_minimal(k: int) -> list[PauliString]:
    r"""Return the minimal left universal set ``{X_i, Z_i}_i \cup {Z^{\otimes k}}``."""

    if k % 3 != 0:
        a_ops = kn_case(k)
    else:
        a_ops = k3_case(k)

    return a_ops

def kn_case(k: int) -> list[PauliString]:
    """Return the minimal left universal set for k equal 3."""

    if k == 1:
        return [get_single(1, 0, "X"), get_single(1, 0, "Z")]

    a_ops: list[PauliString] = [
        get_single(k, 0, "X"),
        get_single(k, 0, "Z"),
        get_single(k, 1, "X"),
        get_single(k, 1, "Z"),
        get_pauli_string("ZZ" + "I" * (k - 2)),
    ]

    for i in range(1, k - 1):
        xz = ["I"] * k
        xz[i] = "X"
        xz[i + 1] = "Z"
        zx = ["I"] * k
        zx[i] = "Z"
        zx[i + 1] = "X"
        a_ops.append(get_pauli_string("".join(xz)))
        a_ops.append(get_pauli_string("".join(zx)))

    return a_ops

def k3_case(N: int) -> list[PauliString]:
    """Return the minimal left universal set for k not equal 3."""

    if N % 3 != 0:
        raise ValueError(
            f"Example 3 requires N divisible by 3; got N={N}. "
            "Use Example 1 or 2 for general N."
        )
    k = N // 3
    M = 4 * k - 1

    HS_HEIS = {"I": "I", "X": "Z", "Y": "X", "Z": "Y"}
    HS_SCHR = {"I": "I", "X": "Y", "Y": "Z", "Z": "X"}

    PMUL = {(a, b): c for a, b, c in [
        ("I","I","I"),("I","X","X"),("I","Y","Y"),("I","Z","Z"),
        ("X","I","X"),("X","X","I"),("X","Y","Z"),("X","Z","Y"),
        ("Y","I","Y"),("Y","X","Z"),("Y","Y","I"),("Y","Z","X"),
        ("Z","I","Z"),("Z","X","Y"),("Z","Y","X"),("Z","Z","I"),
    ]}

    def cz_pair(c1: str, c2: str) -> tuple[str, str]:
        z1 = "Z" if c1 in ("X", "Y") else "I"
        z2 = "Z" if c2 in ("X", "Y") else "I"
        return PMUL[(c1, z2)], PMUL[(z1, c2)]

    def step_heis(p: list[str]) -> list[str]:
        q = [HS_HEIS[c] for c in p]
        for i in range(len(q) - 1):
            q[i], q[i + 1] = cz_pair(q[i], q[i + 1])
        return q

    def step_schr(p: list[str]) -> list[str]:
        q = p[:]
        for i in range(len(q) - 1):
            q[i], q[i + 1] = cz_pair(q[i], q[i + 1])
        return [HS_SCHR[c] for c in q]

    I_M: list[int] = [-1]
    for j in range(k):
        I_M.extend([4*j, 4*j + 1, 4*j + 2, -4*j - 2, -4*j - 3, -4*j - 4])
    assert len(I_M) == 2 * N + 1

    orbit: dict[int, list[str]] = {0: ["Z"] + ["I"] * (M - 1)}
    p = orbit[0][:]
    for ell in range(1, max(I_M) + 1):
        p = step_heis(p)
        orbit[ell] = p[:]
    p = orbit[0][:]
    for ell in range(-1, min(I_M) - 1, -1):
        p = step_schr(p)
        orbit[ell] = p[:]

    drop = {4 * j - 1 for j in range(1, k)}
    def tilde(s: list[str]) -> str:
        return "".join(c for i, c in enumerate(s) if i not in drop)

    return [get_pauli_string(tilde(orbit[ell])) for ell in I_M]

def choose_u_for_b(k: int) -> PauliString:
    """Choose the fixed left tag used when coupling to right-side generators."""
    return get_single(k, 0, "X")

def _all_left_paulis(k: int) -> list[PauliString]:
    """Enumerate all non-identity Pauli strings on ``k`` qubits."""
    return [p for p in get_identity(k).gen_all_pauli_strings() if not p.is_identity()]

@dataclass
class SubsystemCompilerConfig:
    """Configuration of the subsystem compiler."""

    k_left: int
    n_total: int

class SubsystemCompiler:
    """Compiler for the right subsystem contribution."""

    def __init__(self, cfg: SubsystemCompilerConfig):
        if cfg.k_left < 2:
            raise ValueError("k_left must be >= 2 for the Pauli Compiler algorithm")
        self.k = cfg.k_left
        self.n_total = cfg.n_total
        self.n_right = self.n_total - self.k
        self.u_tag = choose_u_for_b(self.k)
        self.left_pool = _all_left_paulis(self.k)

    def extend_left(self, a_left: PauliString) -> PauliString:
        """Extend a left Pauli string by identities on the right."""
        return a_left + get_identity(self.n_right)

    def extend_pair(self, u_left: PauliString, b_right: PauliString) -> PauliString:
        """Combine left and right parts into a full-length Pauli string."""
        return u_left + b_right

    def factor_w_orders(self, w_right: PauliString) -> list[list[tuple[PauliString, PauliString]]]:
        """Factor ``w_right`` into ordered right-local pieces.

        Each element of the output is a list ``[(U_i, B_i)]`` describing one valid
        order of factors. Only the Y case has two possible local orderings.
        """
        assert len(w_right) == self.n_right
        per_site_options: list[list[list[PauliString]]] = []

        for site, label in enumerate(str(w_right)):
            if label == "Y":
                per_site_options.append(
                    [
                        [get_single(self.n_right, site, "X"), get_single(self.n_right, site, "Z")],
                        [get_single(self.n_right, site, "Z"), get_single(self.n_right, site, "X")],
                    ]
                )
            elif label == "X":
                per_site_options.append([[get_single(self.n_right, site, "X")]])
            elif label == "Z":
                per_site_options.append([[get_single(self.n_right, site, "Z")]])
            else:
                per_site_options.append([[]])

        sequences: list[list[PauliString]] = []

        def rec(index: int, acc: list[PauliString]) -> None:
            if index == len(per_site_options):
                sequences.append(list(acc))
                return
            for segment in per_site_options[index]:
                acc.extend(segment)
                rec(index + 1, acc)
                if segment:
                    del acc[-len(segment) :]

        rec(0, [])
        return [[(self.u_tag, b) for b in flat] for flat in sequences]

    def _choose_a1_a2(self, u_i: PauliString) -> tuple[PauliString, PauliString]:
        """Choose helpers ``A_1, A_2`` with the required commutation pattern."""
        anti_with_u = u_i.get_anti_commutants(self.left_pool)
        for a1 in anti_with_u:
            for a2 in anti_with_u:
                if a2 != a1:
                    if a1 | a2:
                        return a1, a2
        raise RuntimeError("Failed to find A1,A2 in iP_k^*.")

    def _choose_aprime(
        self, u_i: PauliString, prod_uj_a: PauliString, prod_uj_gt_a: PauliString
        ) -> PauliString:
        """Choose a helper ``A'`` that anticommutes
        with ``u_i`` and commutes with ``prod_uj_gt_a``."""
        for aprime in u_i.get_anti_commutants(self.left_pool):
            if not aprime | prod_uj_gt_a:
                if aprime != prod_uj_a:
                    return aprime
        raise RuntimeError("Failed to find A' in iP_k^*.")

    def _rest_full_after(
        self,
        ui_bi: list[tuple[PauliString, PauliString]],
        index: int,
        helpers: list[PauliString],
    ) -> tuple[PauliString, PauliString]:
        """Return the accumulated left and right factors after ``index``."""
        p_left = get_identity(self.k)
        for j in range(index + 1, len(ui_bi)):
            p_left = p_left @ ui_bi[j][0]
        for helper in helpers:
            p_left = p_left @ helper

        p_right = get_identity(self.n_right)
        for j in range(index + 1, len(ui_bi)):
            p_right = p_right @ ui_bi[j][1]
        return p_left, p_right

    def _product_uj_a(
        self,
        index: int,
        ui_bi: list[tuple[PauliString, PauliString]],
        helpers: list[PauliString],
    ) -> PauliString:
        """
        Return product of Uj right multiply with product of A for index j to r-1.
        Args:
            index (int): Index of Uj.
            ui_bi (list[tuple(PauliString, PauliString)]): List of (Uj, Bj) pairs.
            helpers (list[PauliString]): List of A
        Returns:
            prod_uj_a (PauliString): Product of Uj right multiply with product of A
        """
        prod_uj = get_identity(self.k)
        for j in range(index, len(ui_bi)):
            prod_uj = prod_uj @ ui_bi[j][0]

        prod_a = get_identity(self.k)
        for a in helpers:
            prod_a = prod_a @ a

        prod_uj_a = prod_uj @ prod_a

        return prod_uj_a

    def _product_uj_bj(
        self,
        index: int,
        ui_bi: list[tuple[PauliString, PauliString]],
    ) -> PauliString:
        """
        Return product of Uj ⊗ Bj for index j to r-1.
        Args:
            index (int): Index of Uj.
            ui_bi (list[tuple(PauliString, PauliString)]): List of (Uj, Bj) pairs.
        Returns:
            prod_uj_bj (PauliString): Product of Uj ⊗ Bj
        """
        prod_uj_bj = get_identity(self.n_total)
        for j in range(index, len(ui_bi)):
            u_j, b_j = ui_bi[j]
            prod_uj_bj = prod_uj_bj @ self.extend_pair(u_j, b_j)

        return prod_uj_bj

    def _product_a_iden(
        self,
        helpers: list[PauliString],
    ) -> PauliString:
        """
        Return product of A ⊗ I in helpers.
        Args:
            helpers (list[PauliString]): List of A
        Returns:
            prod_a_i (PauliString): Product of A ⊗ I
        """
        prod_a_i = get_identity(self.n_total)
        for a in helpers:
            prod_a_i = prod_a_i @ self.extend_left(a)

        return prod_a_i

    def subsystem_compiler(self, w_right: PauliString) -> list[PauliString]:
        """Compile a target supported only on the right subsystem."""
        assert len(w_right) == self.n_right

        for ui_bi in self.factor_w_orders(w_right):
            if not ui_bi:
                return []

            index = len(ui_bi) - 2
            sequence: list[PauliString] = [self.extend_pair(ui_bi[-1][0], ui_bi[-1][1])]
            helpers: list[PauliString] = []

            while index >= 0:
                u_i, b_i = ui_bi[index]
                ub_i = self.extend_pair(u_i, b_i)
                prod_uj_a = self._product_uj_a(index, ui_bi, helpers)
                prod_uj_gt_a = self._product_uj_a(index + 1, ui_bi, helpers)
                prod_ub_j_gt = self._product_uj_bj(index + 1, ui_bi)
                prod_a_iden = self._product_a_iden(helpers)
                prod_ub_j_gt_prod_a_iden = prod_ub_j_gt @ prod_a_iden

                if prod_uj_a.is_identity():
                    a1, a2 = self._choose_a1_a2(u_i)
                    helpers.append(a1)
                    helpers.append(a2)
                    sequence.append(self.extend_left(a1))
                    sequence.append(self.extend_left(a2))
                    continue

                if ub_i | prod_ub_j_gt_prod_a_iden:
                    a_prime = self._choose_aprime(u_i, prod_uj_a, prod_uj_gt_a)
                    helpers.append(a_prime)
                    sequence.append(self.extend_left(a_prime))
                    continue

                sequence.append(ub_i)
                index -= 1

            return sequence[::-1]

        return []

def left_map_over_a(
    v_from: PauliString,
    v_to: PauliString,
    generators: list[PauliString],
) -> list[PauliString]:
    """
    Find a left-only adjoint path from ``v_from`` to ``v_to`` using BFS.
    Args:
        v_from (PauliString): Starting Pauli string.
        v_to (PauliString): Target Pauli string.
        generators (list[PauliString]): Universal set.
    Returns:
        sequence (list[PauliString]): List of [A1,...,As-1] mapping from v_from to v_to.
    """
    if v_from == v_to:
        return []

    queue: deque[PauliString] = deque([v_from])
    parent: dict[PauliString, tuple[PauliString, PauliString]] = {}
    seen: set[PauliString] = {v_from}

    while queue:
        current = queue.popleft()
        if current == v_to:
            sequence: list[PauliString] = []
            cursor = current
            while cursor != v_from:
                previous, used = parent[cursor]
                sequence.append(used)
                cursor = previous
            return sequence

        for helper in current.get_anti_commutants(generators):
            nxt = helper @ current
            if nxt in seen:
                continue
            seen.add(nxt)
            parent[nxt] = (current, helper)
            queue.append(nxt)

    raise RuntimeError("Left map BFS failed.")

@dataclass
class PauliCompilerConfig:
    """Configuration of the optimal compiler."""

    k_left: int
    n_total: int

class OptimalPauliCompiler:
    """Compiler implementing the construction from arXiv:2408.03294.

    Compiles a target Pauli string into an O(N) length sequence of generators
    that produces the target via nested commutators.

    Args:
        cfg: Compiler configuration specifying the left-right partition
            and fallback search limits.
    """

    def __init__(self, cfg: PauliCompilerConfig):
        if cfg.k_left < 2:
            raise ValueError("k_left must be >= 2 for the Pauli Compiler algorithm")
        self.k = cfg.k_left
        self.n_total = cfg.n_total
        self.n_right = self.n_total - self.k
        self.a_left = left_a_minimal(self.k)
        self.sub = SubsystemCompiler(SubsystemCompilerConfig(k_left=self.k, n_total=self.n_total))

    def extend_left(self, a_left: PauliString) -> PauliString:
        """Extend a left Pauli string to the full system."""
        return a_left + get_identity(self.n_right)

    def _left_factor_from_sequence(self, ops: list[PauliString]) -> PauliString:
        """Extract the left factor of a compiled right-subsystem sequence."""
        result = _evaluate_paulie_orientation(ops)
        assert result is not None
        return result.get_substring(0, self.k)

    def _candidate_decompositions(
            self, w_right: PauliString
    ) -> list[tuple[PauliString, PauliString]]:
        """Return candidate decompositions ``W = W1 @ W2``
        with ``W1`` anti-commuting with ``W2``."""
        candidates: list[tuple[PauliString, PauliString]] = []
        seen: set[tuple[PauliString, PauliString]] = set()

        for w1 in w_right.get_anti_commutants():
            w2 = w1 @ w_right

            w1, w2 = sorted((w1, w2))
            pair = (w1, w2)

            if pair not in seen:
                seen.add(pair)
                candidates.append(pair)

        return candidates

    def _all_interleavings_preserving(
        self,
        a_block: list[PauliString],
        b_block: list[PauliString],
        c_block: list[PauliString],
        cap: int = 60000,
    ) -> Iterable[list[PauliString]]:
        """Yield capped interleavings preserving the order inside each block."""
        count = 0
        len_a, len_b, len_c = len(a_block), len(b_block), len(c_block)

        def rec(
            i: int, j: int, k: int, prefix: list[PauliString]
        ) -> Generator[list[PauliString], None, None]:
            nonlocal count
            if count >= cap:
                return
            if i == len_a and j == len_b and k == len_c:
                count += 1
                yield list(prefix)
                return
            if i < len_a:
                prefix.append(a_block[i])
                yield from rec(i + 1, j, k, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if j < len_b:
                prefix.append(b_block[j])
                yield from rec(i, j + 1, k, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if k < len_c:
                prefix.append(c_block[k])
                yield from rec(i, j, k + 1, prefix)
                prefix.pop()

        return rec(0, 0, 0, [])

    def _all_interleavings_preserving4(
        self,
        a_block: list[PauliString],
        b_block: list[PauliString],
        c_block: list[PauliString],
        d_block: list[PauliString],
        cap: int = 120_000,
    ) -> Iterable[list[PauliString]]:
        """Yield capped interleavings preserving the order inside four blocks."""
        count = 0
        len_a, len_b, len_c, len_d = len(a_block), len(b_block), len(c_block), len(d_block)

        def rec(
            i: int, j: int, k: int, l_idx: int, prefix: list[PauliString]
        ) -> Generator[list[PauliString], None, None]:
            nonlocal count
            if count >= cap:
                return
            if i == len_a and j == len_b and k == len_c and l_idx == len_d:
                count += 1
                yield list(prefix)
                return
            if i < len_a:
                prefix.append(a_block[i])
                yield from rec(i + 1, j, k, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if j < len_b:
                prefix.append(b_block[j])
                yield from rec(i, j + 1, k, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if k < len_c:
                prefix.append(c_block[k])
                yield from rec(i, j, k + 1, l_idx, prefix)
                prefix.pop()
                if count >= cap:
                    return
            if l_idx < len_d:
                prefix.append(d_block[l_idx])
                yield from rec(i, j, k, l_idx + 1, prefix)
                prefix.pop()

        return rec(0, 0, 0, 0, [])

    @staticmethod
    def nested_eval(seq: list[PauliString]) -> PauliString | None:
        """ad_{seq[0]} … ad_{seq[-2]}(seq[-1]), or None if any step gives 0."""
        result: PauliString | None = seq[-1]
        for P in reversed(seq[:-1]):
            result = P.adjoint_map(result)
            if result is None:
                return None
        return result

    def reorder_to_nested(self, L_seq: list[PauliString],
                          R_seq: list[PauliString]) -> list[PauliString]:
        """Given two sequences L, R interpreted as nested commutators, return a
        permutation S of L+R such that nested_eval(S) ∝ [nested_eval(L), nested_eval(R)].
        Applies Lemma G.2 of arXiv:2408.03294 iteratively."""
        if len(R_seq) == 1:
            return [R_seq[0]] + L_seq
        if len(L_seq) == 1:
            return L_seq + R_seq
        R_0 = R_seq[0]
        Left = self.nested_eval(L_seq)
        assert Left is not None
        if Left.commutes_with(R_0):
            return [R_0] + self.reorder_to_nested(L_seq, R_seq[1:])
        else:
            return self.reorder_to_nested([R_0] + L_seq, R_seq[1:])

    def compile(self, v_left: PauliString, w_right: PauliString) -> list[PauliString] | None:
        """Compile a target specified by its left and right factors.

        Args:
            v_left: Left factor of the target (length ``k_left``).
            w_right: Right factor of the target (length ``n_total - k_left``).

        Returns:
            Sequence in PauLie's ``nested_adjoint`` orientation: ``[Am, ..., A1, base]``.

        Raises:
            ValueError: If the lengths of ``v_left`` and ``w_right``
                do not match the configured partition.
            RuntimeError: If no valid sequence is found.
        """

        if len(v_left) != self.k or len(w_right) != self.n_right:
            raise ValueError(
                f"Expected v_left of length {self.k} and w_right of length {self.n_right}, "
                f"got {len(v_left)} and {len(w_right)}."
            )

        sequence: list[PauliString] | None

        if w_right.is_identity() and v_left.is_identity():
            sequence = [v_left+w_right]
            return sequence

        if w_right.is_identity():
            for a_s in self.a_left:
                try:
                    seq_a = left_map_over_a(a_s, v_left, self.a_left)
                    seq_a.append(a_s)
                except RuntimeError:
                    continue
                sequence = [self.extend_left(a) for a in seq_a]
                return sequence
            raise RuntimeError("Left-only mapping failed.")

        if v_left.is_identity():
            for w1, w2 in self._candidate_decompositions(w_right):
                g1 = self.sub.subsystem_compiler(w1)
                v1_prime = self._left_factor_from_sequence(g1)
                g2 = self.sub.subsystem_compiler(w2)
                v2_prime = self._left_factor_from_sequence(g2)
                seq_a = left_map_over_a(v2_prime, v1_prime, self.a_left)
                ext_a = [self.extend_left(a) for a in seq_a]
                sequence = self.reorder_to_nested(g1, [*ext_a, *g2])
                if sequence is not None:
                    return sequence

        else:
            g_right = self.sub.subsystem_compiler(w_right)
            v_prime = self._left_factor_from_sequence(g_right)
            seq_a = left_map_over_a(v_prime, v_left, self.a_left)
            ext_a = [self.extend_left(a) for a in seq_a]
            sequence = [*ext_a,*g_right]
            return sequence

        return None

def construct_universal_set(n_total: int, k: int) -> list[PauliString]:
    """Construct the universal generator set used by the compiler.

    The set consists of left-local generators extended with identities
    and right-local generators tagged with a fixed left Pauli string.

    Args:
        n_total: Total number of qubits.
        k: Number of qubits in the left partition.

    Returns:
        List of ``2 * n_total + 1`` Pauli strings forming the universal set.

    Raises:
        ValueError: If ``k`` is out of range.
    """
    if not 1 <= k < n_total:
        raise ValueError("Require 1 <= k < N")

    a_k = left_a_minimal(k)
    n_right = n_total - k
    u_tag = choose_u_for_b(k)
    right_b = [get_single(n_right, index, "X") for index in range(n_right)] + [
        get_single(n_right, index, "Z") for index in range(n_right)
    ]
    a_prime = [a + get_identity(n_right) for a in a_k]
    b_prime = [u_tag + b for b in right_b]
    return a_prime + b_prime


def compile_target(target: PauliString, k_left: int) -> list[PauliString]:
    """Compile a full target Pauli string into a generator sequence.

    Args:
        target: The target Pauli string to compile.
        k_left: Number of qubits in the left partition (must be >= 2).

    Returns:
        Sequence in PauLie's ``nested_adjoint`` orientation such that
        ``PauliStringCollection(seq[:-1]).nested_adjoint(seq[-1]) == target``.

    Raises:
        ValueError: If ``k_left`` is out of range.
    """
    n_total = len(target)
    if not 2 <= k_left < n_total:
        raise ValueError("Require 2 <= k_left < len(target)")

    v_left = target.get_substring(0, k_left)
    w_right = target.get_substring(k_left, n_total - k_left)
    compiler = OptimalPauliCompiler(PauliCompilerConfig(k_left=k_left, n_total=n_total))
    result = compiler.compile(v_left, w_right)
    if result is None:
        raise RuntimeError("No adjoint map sequence found")

    return result
