"""
    Module to get the generator set with optimal generation rate for su(2^n).

    ``get_optimal_universal_generators(n)`` returns 2n+1 Pauli strings that
    generate su(2^n) and anticommute in a fraction f ~ 0.706 of all pairs, the
    fraction that maximises the generation rate (arXiv:2408.03294).

    Method: start from the minimal universal set of Example 1 of the paper,
    then repeatedly replace a generator by its product with an anticommuting
    neighbour, P -> P*Q, until the target fraction is reached.  Replacing keeps
    the element count fixed (a generator is replaced, not added) and keeps the
    algebra fixed ([P, Q] = 2 P*Q puts P*Q in it, [P*Q, Q] ~ P puts P back);
    only the anticommutation pattern moves.

    Speed comes from working on the bit vectors of
    :mod:`paulie.common.symplectic`.  With v = (x | z) over GF(2), P_i and P_j
    anticommute iff x_i.z_j + x_j.z_i = 1, and P_i P_j is v_i XOR v_j.  So the
    graph is one matrix product, a replacement is one row XOR, and every
    candidate move is scored in closed form (``_candidate_deltas``).
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_collection import PauliStringCollection
from paulie.common.pauli_string_factory import get_pauli_string, get_single
from paulie.common.symplectic import from_symplectic, symplectic_gram, to_symplectic

OPTIMAL_FRACTION = 0.706
"""Anticommutation fraction that maximises the generation rate."""

def _contract(v: np.ndarray, a: np.ndarray, src: int, dst: int) -> int:
    """Replace P_src by P_src * P_dst in place; return the new edge count.

    The product is a row XOR, and bilinearity makes the new graph row the XOR
    of the two old rows.
    """
    v[src] ^= v[dst]
    row = a[src] ^ a[dst]
    row[src] = 0
    a[src] = row
    a[:, src] = row
    return int(a.sum()) // 2


def _minimal_universal_seed(n: int) -> list[PauliString]:
    r"""Example 1 of arXiv:2408.03294: 2n+1 Pauli strings generating su(2^n).

    :math:`\{X_1, Z_1, X_2, Z_2, Z_1 Z_2\} \cup \{X_i Z_{i+1}, Z_i X_{i+1}\}_{i=2}^{n-1}`,
    or :math:`\{X, Z\}` when ``n = 1``.  It is minimal, universal and defined
    for every n -- all three needed, since replacement cannot shrink an
    oversized set nor repair a non-universal one.
    """
    if n < 1:
        raise ValueError(f"Number of qubits must be >= 1; got n={n}")
    if n == 1:
        return [get_single(1, 0, "X"), get_single(1, 0, "Z")]

    ops = [
        get_single(n, 0, "X"),
        get_single(n, 0, "Z"),
        get_single(n, 1, "X"),
        get_single(n, 1, "Z"),
        get_pauli_string("ZZ" + "I" * (n - 2)),
    ]
    for i in range(1, n - 1):
        for first, second in (("X", "Z"), ("Z", "X")):
            label = ["I"] * n
            label[i], label[i + 1] = first, second
            ops.append(get_pauli_string("".join(label)))
    return ops


def _candidate_deltas(a: np.ndarray, delta: int) -> tuple[np.ndarray, np.ndarray]:
    """Leftover gap after every possible replacement, from one matrix product.

    Replacing ``i <- i*j`` XORs row j into row i, giving edge count
    ``E + d_j - 2*c_ij - 1`` for degrees ``d`` and shared-neighbour counts
    ``c = a @ a.T``.  Returns ``(dx, dy)``: the gap after ``i <- i*j`` and
    after ``j <- j*i``.
    """
    deg = a.sum(axis=1)
    c = a @ a.T
    return delta - deg[None, :] + 2 * c + 1, delta - deg[:, None] + 2 * c + 1


def _pick_move(
    a: np.ndarray,
    iu: tuple[np.ndarray, np.ndarray],
    delta: int,
    rng: np.random.Generator,
) -> tuple[int, int] | None:
    """Pick the next replacement: the move that gets strictly closer to the
    target (ties broken at random, overshoot allowed if it lands closer), or a
    random edge if none helps.  Returns ``(src, dst)``, or ``None``
    if the graph has no edges.
    """
    big = np.iinfo(np.int32).max
    dx, dy = _candidate_deltas(a, delta)
    is_edge = a == 1
    sx = np.where(is_edge, np.abs(dx), big)[iu]
    sy = np.where(is_edge, np.abs(dy), big)[iu]
    best = min(int(sx.min(initial=big)), int(sy.min(initial=big)))
    improved = best < abs(delta)

    if improved:
        pool = [(int(k), 0) for k in np.flatnonzero(sx == best)]
        pool += [(int(k), 1) for k in np.flatnonzero(sy == best)]
        k, which = pool[int(rng.integers(len(pool)))]
    else:
        edges = np.flatnonzero(is_edge[iu])
        if edges.size == 0:
            return None
        k = int(edges[int(rng.integers(edges.size))])
        which = int(rng.integers(2))

    i, j = int(iu[0][k]), int(iu[1][k])
    return (i, j) if which == 0 else (j, i)


def _search(v: np.ndarray, target: int, seed: int | None = 0) -> np.ndarray:
    """Replace generators until the graph has ``target`` edges, returning the
    closest matrix found.

    Not every target is reachable (the graph stays connected, so
    ``edges >= m - 1``, and no universal set has every pair anticommuting), so
    the search is capped and keeps the best state seen.
    """
    rng = np.random.default_rng(seed)
    v = v.copy()
    a = symplectic_gram(v)
    m = a.shape[0]
    iu = np.triu_indices(m, k=1)
    max_sweeps = max(1000, 40 * m)

    edges = int(a.sum()) // 2
    best_v, best_delta = v.copy(), abs(target - edges)
    sweeps = 0

    while sweeps < max_sweeps and edges != target:
        sweeps += 1
        move = _pick_move(a, iu, target - edges, rng)
        if move is None:
            break
        src, dst = move
        edges = _contract(v, a, src, dst)
        if abs(target - edges) < best_delta:
            best_delta = abs(target - edges)
            best_v = v.copy()

    return best_v


@lru_cache(maxsize=None)
def _cached(n: int, fraction: float, seed: int | None) -> np.ndarray:
    """Build the set for ``n`` once (seed -> search); returned read-only."""
    base = _minimal_universal_seed(n)
    v = to_symplectic(base)
    if len(base) >= 3:  # a single anticommuting pair is already rigid
        v = _search(v, get_optimal_edges_su_2_n(len(base), fraction), seed)
    v.flags.writeable = False
    return v


def get_optimal_edges_su_2_n(ng: int, fraction: float = OPTIMAL_FRACTION) -> int:
    r"""
    Get the optimal number of edges in the anticommutation graph for :math:`\mathfrak{su}(2^{n})`.

    Args:
        ng (int): Number of generators.
        fraction (float): Desired fraction of anticommuting pairs.
    Returns:
        int: Target number of anticommuting pairs, or ``-1`` if ``ng < 2``.
    """
    if ng < 2:
        return -1
    return math.floor(fraction * ng * (ng - 1) / 2)


def get_optimal_universal_generators(
    n: int, fraction: float = OPTIMAL_FRACTION, seed: int | None = 0
) -> PauliStringCollection:
    r"""
    Get an optimal universal generator set for :math:`\mathfrak{su}(2^{n})`.

    Returns ``2n + 1`` Pauli strings (two when ``n = 1``) generating
    :math:`\mathfrak{su}(2^{n})`, cached per ``(n, fraction, seed)``.

    Args:
        n (int): Number of qubits, i.e. the exponent in
            :math:`\mathfrak{su}(2^{n})`.
        fraction (float): Wanted fraction of anticommuting pairs, in ``[0, 1]``.
            Values below about ``2 / (2n + 1)`` or near 1 are unreachable; the
            closest set found is returned.
        seed (int | None): Tie-breaking seed.  The default is deterministic;
            ``None`` gives a different valid set each call.
    Returns:
        PauliStringCollection: An optimal universal generator set for
        :math:`\mathfrak{su}(2^{n})`.

    Raises:
        ValueError: If ``n < 1`` or ``fraction`` is outside ``[0, 1]``.
    """
    if n < 1:
        raise ValueError(f"Number of qubits must be >= 1; got n={n}")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"fraction must be between 0 and 1; got {fraction}")
    return PauliStringCollection(from_symplectic(_cached(n, fraction, seed)))
