"""Pauli compiler based on arXiv:2408.03294.
``compiler`` takes a generator set and a target Pauli string and outputs a
:math:`\\mathcal{O}(N)` length sequence of Pauli strings that generates the target
Pauli string via nested commutators.
"""

from collections import deque
from dataclasses import dataclass
from itertools import permutations
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
    a_ops: list[PauliString] = []
    for index in range(k):
        a_ops.append(get_single(k, index, "X"))
        a_ops.append(get_single(k, index, "Z"))
    a_ops.append(get_pauli_string("Z" * k))
    return a_ops


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

    def _choose_a1_a2(self, u_op: PauliString) -> tuple[PauliString, PauliString]:
        """Choose helpers ``A_1, A_2`` with the required commutation pattern."""
        anti_with_u = u_op.get_anti_commutants(self.left_pool)
        for a1 in anti_with_u:
            for a2 in anti_with_u:
                if a2 == a1:
                    continue
                if a1 | a2:
                    return a1, a2
        raise RuntimeError("Failed to find A1,A2 in iP_k^*.")

    def _choose_aprime(self, u_i: PauliString, p_left: PauliString) -> PauliString:
        """Choose a helper ``A'`` that anticommutes with ``u_i`` and commutes with ``p_left``."""
        for helper in u_i.get_anti_commutants(self.left_pool):
            if helper | p_left:
                return helper
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

    def subsystem_compiler(self, w_right: PauliString) -> list[PauliString]:
        """Compile a target supported only on the right subsystem."""
        assert len(w_right) == self.n_right

        for ui_bi in self.factor_w_orders(w_right):
            if not ui_bi:
                return []

            index = len(ui_bi) - 1
            sequence: list[PauliString] = [self.extend_pair(ui_bi[-1][0], ui_bi[-1][1])]
            helpers: list[PauliString] = []
            helper_uses: dict[int, int] = {}

            while index >= 1:
                u_i, b_i = ui_bi[index]
                p_left, p_right = self._rest_full_after(ui_bi, index, helpers)

                if (p_left @ u_i).is_identity():
                    count = helper_uses.get(index, 0)
                    if count >= 1:
                        sequence.append(self.extend_pair(u_i, b_i))
                        index -= 1
                        continue

                    a1, a2 = self._choose_a1_a2(u_i)
                    helpers = [a1, a2]
                    helper_uses[index] = count + 1
                    sequence.append(self.extend_left(a1))
                    sequence.append(self.extend_left(a2))
                    continue

                current = self.extend_pair(u_i, b_i)
                rest_full = p_left + p_right
                if current | rest_full:
                    count = helper_uses.get(index, 0)
                    if count >= 1:
                        sequence.append(current)
                        index -= 1
                        continue

                    a_prime = self._choose_aprime(u_i, p_left)
                    helpers = [a_prime]
                    helper_uses[index] = count + 1
                    sequence.append(self.extend_left(a_prime))
                    continue

                sequence.append(current)
                index -= 1

            return sequence

        return []


def left_map_over_a(
    v_from: PauliString,
    v_to: PauliString,
    generators: list[PauliString],
) -> list[PauliString]:
    """Find a left-only adjoint path from ``v_from`` to ``v_to`` using BFS."""
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
            sequence.reverse()
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
    fallback_depth: int = 8
    fallback_nodes: int = 200000


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
        self.fallback_depth = cfg.fallback_depth
        self.fallback_nodes = cfg.fallback_nodes

    def extend_left(self, a_left: PauliString) -> PauliString:
        """Extend a left Pauli string to the full system."""
        return a_left + get_identity(self.n_right)

    def _left_factor_from_sequence(self, ops: list[PauliString]) -> PauliString:
        """Extract the left factor of a compiled right-subsystem sequence."""
        result = _evaluate_sequence(ops)
        if result is None:
            return get_single(self.k, 0, "X")
        return result.get_substring(0, self.k)

    def _candidate_decompositions(
            self, w_right: PauliString
    ) -> list[tuple[PauliString, PauliString]]:
        """Return candidate decompositions ``W = W1 @ W2``
        with ``W1`` anti-commuting with ``W2``."""
        candidates: list[tuple[PauliString, PauliString]] = []
        n_right = len(w_right)

        for site, label in enumerate(str(w_right)):
            if label == "I":
                continue
            labels = ["X", "Z"] if label == "Y" else (["Z"] if label == "X" else ["X"])
            for local_label in labels:
                w1 = get_single(n_right, site, local_label)
                w2 = w1 @ w_right
                if not w1 | w2:
                    candidates.append((w1, w2))

        unique: list[tuple[PauliString, PauliString]] = []
        seen: set[tuple[PauliString, PauliString]] = set()
        for pair in candidates:
            if pair in seen:
                continue
            seen.add(pair)
            unique.append(pair)
        return unique

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

        def rec(i: int, j: int, k: int, prefix: list[PauliString]):
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

        def rec(i: int, j: int, k: int, l_idx: int, prefix: list[PauliString]):
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

    def _case3_best_reordering(
        self,
        g1: list[PauliString],
        g2: list[PauliString],
        a_ext: list[PauliString],
        w_right: PauliString,
    ) -> list[PauliString] | None:
        """Search for a valid reordering in the ``V = I`` and ``W != I`` case."""
        g1_rev, g2_rev = list(reversed(g1)), list(reversed(g2))
        candidates = [
            list(g1) + list(a_ext) + list(g2_rev) + list(a_ext),
            list(g1_rev) + list(a_ext) + list(g2) + list(a_ext),
            list(a_ext) + list(g1) + list(a_ext) + list(g2_rev),
            list(a_ext) + list(g1_rev) + list(a_ext) + list(g2),
            list(g2) + list(a_ext) + list(g1_rev) + list(a_ext),
            list(g2_rev) + list(a_ext) + list(g1) + list(a_ext),
        ]
        target_left = get_identity(self.k)
        for sequence in candidates:
            result = _evaluate_sequence(sequence)
            if (
                result is not None
                and result.get_substring(0, self.k) == target_left
                and result.get_substring(self.k, self.n_right) == w_right
            ):
                return sequence

        blocks = [g1, g2, a_ext]
        for perm in permutations(range(3)):
            b0, b1, b2 = [blocks[i] for i in perm]
            for r0 in (False, True):
                for r1 in (False, True):
                    for r2 in (False, True):
                        sequence = (
                            (list(reversed(b0)) if r0 else list(b0))
                            + (list(reversed(b1)) if r1 else list(b1))
                            + (list(reversed(b2)) if r2 else list(b2))
                        )
                        result = _evaluate_sequence(sequence)
                        if (
                            result is not None
                            and result.get_substring(0, self.k) == target_left
                            and result.get_substring(self.k, self.n_right) == w_right
                        ):
                            return sequence

        for g1_block in (g1, list(reversed(g1))):
            for g2_block in (g2, list(reversed(g2))):
                for a_block in (a_ext, list(reversed(a_ext))):
                    for sequence in self._all_interleavings_preserving(g1_block, g2_block, a_block):
                        result = _evaluate_sequence(sequence)
                        if (
                            result is not None
                            and result.get_substring(0, self.k) == target_left
                            and result.get_substring(self.k, self.n_right) == w_right
                        ):
                            return sequence

        a_options = (a_ext, list(reversed(a_ext)))
        for g1_block in (g1, g1_rev):
            for g2_block in (g2, g2_rev):
                for a1 in a_options:
                    for a2 in a_options:
                        for sequence in (
                            list(g1_block) + list(a1) + list(g2_block) + list(a2),
                            list(g2_block) + list(a1) + list(g1_block) + list(a2),
                            list(a1) + list(g1_block) + list(a2) + list(g2_block),
                            list(a1) + list(g2_block) + list(a2) + list(g1_block),
                        ):
                            result = _evaluate_sequence(sequence)
                            if (
                                result is not None
                                and result.get_substring(0, self.k) == target_left
                                and result.get_substring(self.k, self.n_right) == w_right
                            ):
                                return sequence
                        for sequence in self._all_interleavings_preserving4(
                                g1_block, g2_block, a1, a2
                        ):
                            result = _evaluate_sequence(sequence)
                            if (
                                result is not None
                                and result.get_substring(0, self.k) == target_left
                                and result.get_substring(self.k, self.n_right) == w_right
                            ):
                                return sequence
        return None

    def _bfs_case3(
        self, w_right: PauliString, depth_cap: int, node_cap: int
    ) -> list[PauliString] | None:
        """Fallback bounded BFS for the case ``V = I`` and ``W != I``."""
        universal_set = construct_universal_set(self.n_total, self.k)
        target_left = get_identity(self.k)
        nodes = 0
        frontier: list[tuple[PauliString | None, list[int]]] = [(None, [])]
        visited: set[tuple[int, PauliString]] = set()

        for depth in range(1, depth_cap + 1):
            new_frontier: list[tuple[PauliString | None, list[int]]] = []
            for result, seq_idx in frontier:
                for op_index, operator in enumerate(universal_set):
                    nodes += 1
                    if nodes > node_cap:
                        return None
                    new_result = operator if result is None else (operator ^ result)
                    if new_result is None:
                        continue
                    state_key = (depth, new_result)
                    if state_key in visited:
                        continue
                    visited.add(state_key)
                    new_sequence = seq_idx + [op_index]
                    if depth >= 2:
                        if (
                            new_result.get_substring(0, self.k) == target_left
                            and new_result.get_substring(self.k, self.n_right) == w_right
                        ):
                            return [universal_set[idx] for idx in new_sequence]
                    new_frontier.append((new_result, new_sequence))
            frontier = new_frontier
        return None

    def compile(self, v_left: PauliString, w_right: PauliString) -> list[PauliString]:
        """Compile a target specified by its left and right factors.

        Args:
            v_left: Left factor of the target (length ``k_left``).
            w_right: Right factor of the target (length ``n_total - k_left``).

        Returns:
            Sequence in PauLie's ``nested_adjoint`` orientation: ``[Am, ..., A1, base]``.

        Raises:
            ValueError: If the lengths of ``v_left`` and ``w_right`` do not match the configured partition.
            RuntimeError: If no valid sequence is found.
        """
        if len(v_left) != self.k or len(w_right) != self.n_right:
            raise ValueError(
                f"Expected v_left of length {self.k} and w_right of length {self.n_right}, "
                f"got {len(v_left)} and {len(w_right)}."
            )
        if w_right.is_identity():
            for seed in self.a_left:
                try:
                    seq_a = left_map_over_a(seed, v_left, self.a_left)
                except RuntimeError:
                    continue
                sequence = [self.extend_left(seed)] + [self.extend_left(a) for a in seq_a]
                result = _evaluate_sequence(sequence)
                if (
                    result is not None
                    and result.get_substring(0, self.k) == v_left
                    and result.get_substring(self.k, self.n_right) == w_right
                ):
                    return _sequence_to_paulie_orientation(sequence)
            raise RuntimeError("Left-only mapping failed.")

        if not v_left.is_identity():
            g_right = self.sub.subsystem_compiler(w_right)
            v_prime = self._left_factor_from_sequence(g_right)
            seq = left_map_over_a(v_prime, v_left, self.a_left)
            candidates = [
                list(g_right) + [self.extend_left(a) for a in seq],
                [self.extend_left(a) for a in seq] + list(g_right),
                list(reversed(g_right)) + [self.extend_left(a) for a in seq],
            ]
            for sequence in candidates:
                result = _evaluate_sequence(sequence)
                if (
                    result is not None
                    and result.get_substring(0, self.k) == v_left
                    and result.get_substring(self.k, self.n_right) == w_right
                ):
                    return _sequence_to_paulie_orientation(sequence)
            return _sequence_to_paulie_orientation(
                list(g_right) + [self.extend_left(a) for a in seq]
            )

        for w1, w2 in self._candidate_decompositions(w_right):
            g1 = self.sub.subsystem_compiler(w1)
            g2 = self.sub.subsystem_compiler(w2)
            v1_prime = self._left_factor_from_sequence(g1)
            v2_prime = self._left_factor_from_sequence(g2)
            a_seq = left_map_over_a(v2_prime, v1_prime, self.a_left)
            a_ext = [self.extend_left(a) for a in a_seq]
            sequence = self._case3_best_reordering(g1, g2, a_ext, w_right)
            if sequence is not None:
                return _sequence_to_paulie_orientation(sequence)

        seq_fb = self._bfs_case3(w_right, self.fallback_depth, self.fallback_nodes)
        if seq_fb is not None:
            return _sequence_to_paulie_orientation(seq_fb)

        w_str = str(w_right)
        site = next(index for index, label in enumerate(w_str) if label != "I")
        label = "X" if w_str[site] == "Z" else ("Z" if w_str[site] == "X" else "X")
        w1 = get_single(self.n_right, site, label)
        w2 = w1 @ w_right
        g1 = self.sub.subsystem_compiler(w1)
        g2 = self.sub.subsystem_compiler(w2)
        v1_prime = self._left_factor_from_sequence(g1)
        v2_prime = self._left_factor_from_sequence(g2)
        a_seq = left_map_over_a(v2_prime, v1_prime, self.a_left)
        a_ext = [self.extend_left(a) for a in a_seq]
        return _sequence_to_paulie_orientation(list(reversed(g1)) + a_ext + list(reversed(g2)))


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
    return compiler.compile(v_left, w_right)
