"""
Fixed-unitary OTOC and related quantities (e.g. Pauli instability).
"""

from __future__ import annotations

from typing import Literal, cast

import numpy as np
from numpy.typing import NDArray

from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_factory import pauli_string_from_index

# ``method="exact"`` sums over ``16**n`` Pauli pairs; cap ``n`` to avoid hangs.
_MAX_QUBITS_EXACT_OTOC_UNIFORM = 4


def _abs_trace_times_factor(chained: NDArray[np.complex128], factor: float) -> float:
    """``| factor * Tr(chained) |`` as a plain ``float``."""
    z = cast(complex, np.trace(chained).item())
    return abs(z * factor)


def _pauli_dense_matrices(n: int) -> list[np.ndarray]:
    """All :math:`4^n` Pauli matrices as dense ``(2^n, 2^n)`` arrays."""
    return [
        pauli_string_from_index(i, n).get_matrix().astype(np.complex128, copy=False)
        for i in range(4**n)
    ]


def otoc_fixed_unitary(
    u: np.ndarray,
    p1: PauliString,
    p2: PauliString,
    *,
    check_unitary: bool = True,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> np.complex128:
    r"""
    Out-of-time-order correlator for a fixed unitary :math:`U` and Pauli strings :math:`P_1`, :math:`P_2`:

    .. math::

        \mathrm{OTOC}(U, P_1, P_2) = \frac{1}{2^n}\,\mathrm{Tr}\!\left(
        U^\dagger P_1 U P_2 U^\dagger P_1 U P_2\right).

    Args:
        u: Shape ``(2^n, 2^n)`` unitary.
        p1: :math:`P_1` as an :math:`n`-qubit Pauli string.
        p2: :math:`P_2` as an :math:`n`-qubit Pauli string.
        check_unitary: If True, verify ``u`` is unitary up to ``rtol`` / ``atol``.
        rtol: Relative tolerance for the unitary check.
        atol: Absolute tolerance for the unitary check.

    Returns:
        The OTOC value (complex scalar).

    Raises:
        ValueError: If dimensions mismatch or ``u`` is not unitary (when checked).
    """
    if u.ndim != 2 or u.shape[0] != u.shape[1]:
        raise ValueError("u must be a square 2D array")
    d = u.shape[0]
    if d <= 0 or (d & (d - 1)) != 0:
        raise ValueError("u must have size 2^n for some n >= 0")
    n = d.bit_length() - 1
    if len(p1) != n or len(p2) != n:
        raise ValueError(
            f"Pauli strings must have length n={n} matching log2(dim(u))={n}"
        )
    if d == 1:
        if check_unitary and not np.allclose(
            np.abs(u[0, 0]), 1.0, rtol=rtol, atol=atol
        ):
            raise ValueError("u is not unitary (within rtol/atol)")
        return np.complex128(1.0)
    if check_unitary:
        identity: NDArray[np.complex128] = np.eye(d, dtype=np.complex128)
        if not np.allclose(u.conj().T @ u, identity, rtol=rtol, atol=atol):
            raise ValueError("u is not unitary (within rtol/atol)")
    p1_m = p1.get_matrix().astype(np.complex128, copy=False)
    p2_m = p2.get_matrix().astype(np.complex128, copy=False)
    u_c: NDArray[np.complex128] = u.astype(np.complex128, copy=False)
    ud = u_c.conj().T
    chained = ud @ p1_m @ u_c @ p2_m @ ud @ p1_m @ u_c @ p2_m
    return np.complex128(np.trace(chained) / d)


def mean_abs_otoc_uniform(
    u: np.ndarray,
    *,
    method: Literal["exact", "monte_carlo"] = "exact",
    num_samples: int = 10_000,
    seed: int | None = None,
    check_unitary: bool = True,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> float:
    r"""
    Uniform expectation of :math:`\lvert\mathrm{OTOC}(U,P_1,P_2)\rvert` over independent
    :math:`P_1,P_2 \in \{I,X,Y,Z\}^{\otimes n}` (each uniform on :math:`4^n` strings).

    .. math::

        \mathbb{E}_{P_1,P_2}\left[\left\lvert\mathrm{OTOC}(U, P_1, P_2)\right\rvert\right],

    using the same OTOC convention as :func:`otoc_fixed_unitary`.

    Args:
        u: Unitary of shape ``(2^n, 2^n)``.
        method: ``"exact"`` — sum over all :math:`16^n` pairs, allowed only for
            :math:`n \le 4` qubits; ``"monte_carlo"`` — sample :math:`P_1,P_2` i.i.d. uniform
            (any :math:`n` supported by dense arithmetic).
        num_samples: Number of i.i.d. pairs for ``method="monte_carlo"``.
        seed: If given, ``method="monte_carlo"`` uses ``numpy.random.default_rng(seed)``
            (reproducible); if ``None``, a fresh default generator is used.
        check_unitary: Forwarded to :func:`otoc_fixed_unitary`.
        rtol: Forwarded to :func:`otoc_fixed_unitary`.
        atol: Forwarded to :func:`otoc_fixed_unitary`.

    Returns:
        Estimated or exact expectation in :math:`[0,1]` (typically).

    Raises:
        ValueError: Invalid ``u`` or ``method``, or ``method="exact"`` with more than
            four qubits.

    """
    if u.ndim != 2 or u.shape[0] != u.shape[1]:
        raise ValueError("u must be a square 2D array")
    d = u.shape[0]
    if d <= 0 or (d & (d - 1)) != 0:
        raise ValueError("u must have size 2^n for some n >= 0")
    n = d.bit_length() - 1

    inner_check = check_unitary
    if check_unitary:
        if d == 1:
            if not np.allclose(np.abs(u[0, 0]), 1.0, rtol=rtol, atol=atol):
                raise ValueError("u is not unitary (within rtol/atol)")
        else:
            identity_u: NDArray[np.complex128] = np.eye(d, dtype=np.complex128)
            if not np.allclose(u.conj().T @ u, identity_u, rtol=rtol, atol=atol):
                raise ValueError("u is not unitary (within rtol/atol)")
        inner_check = False

    if method == "exact":
        if n > _MAX_QUBITS_EXACT_OTOC_UNIFORM:
            raise ValueError(
                f"method='exact' supports at most {_MAX_QUBITS_EXACT_OTOC_UNIFORM} qubits "
                f"(got n={n}); use method='monte_carlo' for larger systems."
            )
        total_pairs = 16**n
        if n == 0:
            p0 = PauliString(n=0)
            o = otoc_fixed_unitary(
                u, p0, p0, check_unitary=inner_check, rtol=rtol, atol=atol
            )
            return float(abs(complex(o)))
        u_ce: NDArray[np.complex128] = u.astype(np.complex128, copy=False)
        ud = u_ce.conj().T
        pauli_mats = _pauli_dense_matrices(n)
        sandwiches = [ud @ pauli_mats[i] @ u_ce for i in range(4**n)]
        acc = 0.0
        inv_d = 1.0 / d
        for i in range(4**n):
            mi = sandwiches[i]
            for j in range(4**n):
                pj = pauli_mats[j]
                chained = mi @ pj @ mi @ pj
                acc += _abs_trace_times_factor(chained, inv_d)
        return float(acc / total_pairs)

    if method == "monte_carlo":
        gen = np.random.default_rng(seed)
        u_cm: NDArray[np.complex128] = u.astype(np.complex128, copy=False)
        ud = u_cm.conj().T
        inv_d = 1.0 / d
        acc = 0.0
        if n == 0:
            p0 = PauliString(n=0)
            for _ in range(num_samples):
                o = otoc_fixed_unitary(
                    u, p0, p0, check_unitary=inner_check, rtol=rtol, atol=atol
                )
                acc += abs(complex(o))
            return float(acc / num_samples)
        if n <= _MAX_QUBITS_EXACT_OTOC_UNIFORM:
            pauli_mats = _pauli_dense_matrices(n)
            sandwiches = [ud @ pauli_mats[i] @ u_cm for i in range(4**n)]
            for _ in range(num_samples):
                i = int(gen.integers(0, 4**n))
                j = int(gen.integers(0, 4**n))
                mi = sandwiches[i]
                pj = pauli_mats[j]
                chained = mi @ pj @ mi @ pj
                acc += _abs_trace_times_factor(chained, inv_d)
            return float(acc / num_samples)
        for _ in range(num_samples):
            i = int(gen.integers(0, 4**n))
            j = int(gen.integers(0, 4**n))
            p1 = pauli_string_from_index(i, n)
            p2 = pauli_string_from_index(j, n)
            o = otoc_fixed_unitary(
                u, p1, p2, check_unitary=inner_check, rtol=rtol, atol=atol
            )
            acc += abs(complex(o))
        return float(acc / num_samples)

    raise ValueError(f"Invalid method: {method}")


def pauli_instability(
    u: np.ndarray,
    *,
    method: Literal["exact", "monte_carlo"] = "exact",
    num_samples: int = 10_000,
    seed: int | None = None,
    base: float | None = None,
    check_unitary: bool = True,
    rtol: float = 1e-8,
    atol: float = 1e-8,
) -> float:
    r"""
    Pauli instability (Definition 1):

    .. math::

        I(U) = -\log \mathbb{E}_{P_1,P_2}\left[
        \left\lvert\mathrm{OTOC}(U, P_1, P_2)\right\rvert\right],

    expectation over independent uniform :math:`P_1,P_2` on :math:`\{I,X,Y,Z\}^{\otimes n}`.

    Args:
        u: Unitary of shape ``(2^n, 2^n)``.
        method: Passed to :func:`mean_abs_otoc_uniform`.
        num_samples: Passed to :func:`mean_abs_otoc_uniform` when ``method="monte_carlo"``.
        seed: Passed to :func:`mean_abs_otoc_uniform` when ``method="monte_carlo"``.
        base: Logarithm base; ``None`` means natural log (``numpy.log``).
        check_unitary: Passed to :func:`mean_abs_otoc_uniform` / :func:`otoc_fixed_unitary`.
        rtol: Passed through.
        atol: Passed through.

    Returns:
        Finite :math:`I(U)` when the expectation is strictly positive.

    Raises:
        ValueError: If the expectation is numerically zero (``log`` undefined), ``u`` is invalid,
            or ``method="exact"`` is used with more than four qubits.

    """
    m = mean_abs_otoc_uniform(
        u,
        method=method,
        num_samples=num_samples,
        seed=seed,
        check_unitary=check_unitary,
        rtol=rtol,
        atol=atol,
    )
    if m == 0:
        raise ValueError("Expectation numerically zero")
    if base is None:
        return float(-np.log(m))
    return float(-np.log(m) / np.log(base))
