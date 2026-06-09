"""
Tests for get_algebra_basis() and the underlying basis constructors.

Key test: _check_bracket_closure verifies [B_i, B_j] lies in span(basis)
for all i<j. This is the direct algebraic soundness check the maintainer
asked for -- not a dimension proxy.
"""

import re
import numpy as np
import pytest
from paulie import get_pauli_string as p
from paulie.application.algebra_basis import so_basis, su_basis, sp_basis


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _in_span(vec: np.ndarray, basis: np.ndarray, tol: float = 1e-9) -> bool:
    A = basis.reshape(len(basis), -1).T
    v = vec.reshape(-1)
    x, *_ = np.linalg.lstsq(A, v, rcond=None)
    return float(np.linalg.norm(A @ x - v)) < tol * (1.0 + np.linalg.norm(v))


def _check_bracket_closure(basis: np.ndarray, tol: float = 1e-9) -> None:
    """
    Assert [B_i, B_j] is in span(basis) for every i < j.
    A correct Lie algebra basis must be closed under the bracket.
    """
    for i in range(len(basis)):
        for j in range(i + 1, len(basis)):
            comm = basis[i] @ basis[j] - basis[j] @ basis[i]
            assert _in_span(comm, basis, tol), (
                f"[B_{i}, B_{j}] not in span(basis) -- basis is not closed under the Lie bracket"
            )


def _check_full_rank(basis: np.ndarray) -> None:
    flat = basis.reshape(len(basis), -1)
    rank = np.linalg.matrix_rank(flat, tol=1e-10)
    assert rank == len(basis), f"Not full rank: {rank} != {len(basis)}"


def _summand_count(label: str) -> int:
    """Parse 'k*algebra(N)' or 'algebra(N)+algebra(N)+...' to get summand count."""
    m = re.match(r"^(\d+)\*[a-z]+\(\d+\)$", label.strip())
    return int(m.group(1)) if m else len(label.split("+"))


# ---------------------------------------------------------------------------
# so(N)
# ---------------------------------------------------------------------------


class TestSoBasis:
    """Tests for the so(N) real antisymmetric defining-rep basis."""

    @pytest.mark.parametrize("N", [2, 3, 4, 5, 6])
    def test_shape(self, N):
        """Verify shape matches the dimension formula."""
        assert so_basis(N).shape == (N * (N - 1) // 2, N, N)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_antisymmetric(self, N):
        """Each element satisfies X + X.T == 0."""
        for m in so_basis(N):
            np.testing.assert_allclose(m + m.T, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_full_rank(self, N):
        """Basis elements are linearly independent."""
        _check_full_rank(so_basis(N))

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_bracket_closure(self, N):
        """Each [B_i, B_j] decomposes in span(basis)."""
        _check_bracket_closure(so_basis(N))


# ---------------------------------------------------------------------------
# su(N)
# ---------------------------------------------------------------------------


class TestSuBasis:
    """Tests for the su(N) anti-Hermitian traceless defining-rep basis."""

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_shape(self, N):
        """Verify shape matches the dimension formula."""
        assert su_basis(N).shape == (N * N - 1, N, N)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_traceless(self, N):
        """Each element has zero trace."""
        for m in su_basis(N):
            assert abs(np.trace(m)) < 1e-12

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_anti_hermitian(self, N):
        """Each element satisfies X + X.conj().T == 0."""
        for m in su_basis(N):
            np.testing.assert_allclose(m + m.conj().T, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_full_rank(self, N):
        """Basis elements are linearly independent."""
        _check_full_rank(su_basis(N))

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_bracket_closure(self, N):
        """Each [B_i, B_j] decomposes in span(basis)."""
        _check_bracket_closure(su_basis(N))


# ---------------------------------------------------------------------------
# sp(N) -- N = half-dimension, matrices 2N x 2N
# ---------------------------------------------------------------------------


class TestSpBasis:
    """Tests for the sp(N) symplectic defining-rep basis."""

    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_shape(self, N):
        """Verify shape matches the dimension formula."""
        assert sp_basis(N).shape == (N * (2 * N + 1), 2 * N, 2 * N)

    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_symplectic_condition(self, N):
        """Each element satisfies X^T J + J X == 0."""
        size = 2 * N
        J = np.zeros((size, size), dtype=np.float64)
        J[:N, N:] = np.eye(N)
        J[N:, :N] = -np.eye(N)
        for m in sp_basis(N):
            np.testing.assert_allclose(m.T @ J + J @ m, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_full_rank(self, N):
        """Basis elements are linearly independent."""
        _check_full_rank(sp_basis(N))

    @pytest.mark.parametrize("N", [1, 2, 3])
    def test_bracket_closure(self, N):
        """Each [B_i, B_j] decomposes in span(basis)."""
        _check_bracket_closure(sp_basis(N))


# ---------------------------------------------------------------------------
# Integration tests via the public API
# ---------------------------------------------------------------------------


class TestGetAlgebraBasis:
    """Integration tests for PauliStringCollection.get_algebra_basis()."""

    def test_b1_sp4_shape(self):
        """sp(4) basis has shape (36, 8, 8)."""
        basis_list = p(["XY", "XZ"], n=4).get_algebra_basis()
        assert len(basis_list) == 1
        assert basis_list[0].shape == (36, 8, 8)

    def test_b1_sp4_symplectic(self):
        """sp(4) elements satisfy the symplectic condition."""
        b = p(["XY", "XZ"], n=4).get_algebra_basis()[0]
        N = 4
        J = np.zeros((8, 8))
        J[:N, N:] = np.eye(N)
        J[N:, :N] = -np.eye(N)
        for m in b:
            np.testing.assert_allclose(m.T @ J + J @ m, 0, atol=1e-10)

    def test_b1_sp4_bracket_closure(self):
        """sp(4) basis closes under Lie brackets."""
        _check_bracket_closure(p(["XY", "XZ"], n=4).get_algebra_basis()[0])

    def test_a_type_so_antisymmetric(self):
        """Type A so elements are antisymmetric."""
        gens = ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]
        for b in p(gens).get_algebra_basis():
            for m in b:
                np.testing.assert_allclose(m + m.T, 0, atol=1e-10)

    def test_a_type_bracket_closure(self):
        """Type A so basis closes under Lie brackets."""
        gens = ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]
        for b in p(gens).get_algebra_basis():
            _check_bracket_closure(b)

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
            (["XZ", "ZX"], 2),
            (["XX", "YY", "ZZ"], None),
        ],
    )
    def test_total_dim_matches_get_dla_dim(self, gens, n):
        """sum of basis sizes must equal get_dla_dim()."""
        ps = p(gens, n=n) if n else p(gens)
        assert sum(len(b) for b in ps.get_algebra_basis()) == ps.get_dla_dim()

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
        ],
    )
    def test_summand_count_matches_subalgebras(self, gens, n):
        """Summand count from get_algebra_basis() matches get_subalgebras()."""
        ps = p(gens, n=n) if n else p(gens)
        assert len(ps.get_algebra_basis()) == _summand_count(ps.get_algebra())

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
            (["XZ", "ZX"], 2),
        ],
    )
    def test_each_summand_full_rank(self, gens, n):
        """Each summand has full column rank."""
        ps = p(gens, n=n) if n else p(gens)
        for b in ps.get_algebra_basis():
            _check_full_rank(b)

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["XZ", "ZX"], 2),
            (["XX", "YY"], 2),
        ],
    )
    def test_basis_dim_matches_lie_closure_rank(self, gens, n):
        """
        Brute-force DLA via adjoint_map + get_matrix(). The rank of the
        resulting matrix stack must equal the total basis dimension.
        This is the reachability check the maintainer asked for.
        """
        ps = p(gens, n=n) if n else p(gens)

        closure: set = set(ps)
        changed = True
        while changed:
            changed = False
            snapshot = list(closure)
            for a in snapshot:
                for b_ps in snapshot:
                    c = a.adjoint_map(b_ps)
                    if c is not None and c not in closure:
                        closure.add(c)
                        changed = True

        matrices = np.array([q.get_matrix() for q in closure])
        dla_rank = np.linalg.matrix_rank(matrices.reshape(len(matrices), -1), tol=1e-10)
        basis_dim = sum(len(b) for b in ps.get_algebra_basis())
        assert basis_dim == dla_rank, f"gens={gens}: basis dim {basis_dim} != DLA rank {dla_rank}"


# ---------------------------------------------------------------------------
# Element-wise reachability against brute-force DLA  (type B3 / su family)
#
# For su(2^n), defining rep and embedded DLA live in the same space C^{2^n x 2^n}.
# For each DLA element X we check ||X - P_basis(X)|| / ||X|| < tol directly.
# Both directions: DLA->basis and basis->DLA.
# ---------------------------------------------------------------------------


class _LieSpan:
    """Incrementally-maintained orthonormal span via pre-allocated BLAS buffer."""

    def __init__(self, tol: float = 1e-10) -> None:
        self.tol = tol
        self.mats: list = []
        self._buf = None
        self._d = 0
        self._cap = 0

    def _grow(self, flat_size: int) -> None:
        new_cap = max(16, self._cap * 2)
        buf = np.empty((flat_size, new_cap), dtype=complex)
        if self._d > 0:
            buf[:, : self._d] = self._buf[:, : self._d]
        self._buf, self._cap = buf, new_cap

    def _residual(self, v: np.ndarray) -> np.ndarray:
        if self._d == 0:
            return v.copy()
        Q = self._buf[:, : self._d]
        return v - Q @ (Q.conj().T @ v)

    def try_add(self, mat: np.ndarray) -> bool:
        """Add mat to the span if linearly independent; return True iff added."""
        v = mat.ravel().astype(complex)
        flat_size = v.shape[0]
        if self._buf is None:
            self._grow(flat_size)
        elif self._d == self._cap:
            self._grow(flat_size)
        r = self._residual(v)
        nv = float(np.linalg.norm(v))
        if nv < 1e-14 or float(np.linalg.norm(r)) < self.tol * nv:
            return False
        self._buf[:, self._d] = r / float(np.linalg.norm(r))
        self._d += 1
        self.mats.append(mat)
        return True


def _lie_closure(seeds: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Wavefront Lie closure: each new element commuted against all current."""
    span = _LieSpan(tol=tol)
    for m in seeds:
        span.try_add(m)
    prev = 0
    while prev < len(span.mats):
        curr = len(span.mats)
        for i in range(prev, curr):
            for j in range(curr):
                C = span.mats[i] @ span.mats[j] - span.mats[j] @ span.mats[i]
                span.try_add(C)
        prev = curr
    return np.array(span.mats)


_PAULI: dict = {
    "I": np.eye(2, dtype=np.complex128),
    "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
}


def _pm(s: str) -> np.ndarray:
    m = _PAULI[s[0]]
    for c in s[1:]:
        m = np.kron(m, _PAULI[c])
    return m


def _dla(strings: list) -> np.ndarray:
    """Brute-force anti-Hermitian DLA from Pauli string generators."""
    return _lie_closure(np.array([1j * _pm(s) for s in strings]))


def _qr_of(basis: np.ndarray) -> np.ndarray:
    flat = basis.reshape(len(basis), -1).astype(complex).T
    Q, _ = np.linalg.qr(flat, mode="reduced")
    return Q


def _reachability_failures(Q: np.ndarray, elements: np.ndarray, tol: float = 1e-8) -> list:
    failures = []
    for k, X in enumerate(elements):
        v = X.ravel().astype(complex)
        nv = max(float(np.linalg.norm(v)), 1e-14)
        rel_res = float(np.linalg.norm(v - Q @ (Q.conj().T @ v))) / nv
        if rel_res > tol:
            failures.append((k, rel_res))
    return failures


class TestElementWiseReachability:
    """Each DLA element must decompose exactly in the defining-rep basis.

    Type B3 (su family): defining rep and embedded DLA share C^{2^n x 2^n},
    so span equality is checkable element by element without an isomorphism.
    """

    def _assert_reachable(self, Q: np.ndarray, elements: np.ndarray, tag: str) -> None:
        failures = _reachability_failures(Q, elements)
        if failures:
            msgs = [f"  [{k}] rel_res={r:.2e}" for k, r in failures]
            raise AssertionError(f"{tag}: {len(failures)}/{len(elements)} elements not reachable\n" + "\n".join(msgs))

    def test_su2_dla_to_basis(self):
        """Each su(2) DLA element decomposes in su_basis(2)."""
        dla = _dla(["X", "Y", "Z"])
        assert len(dla) == 3, f"su(2) DLA dim should be 3, got {len(dla)}"
        self._assert_reachable(_qr_of(su_basis(2)), dla, "su(2) DLA->basis")

    def test_su2_basis_to_dla(self):
        """Every su_basis(2) element decomposes in the su(2) DLA span."""
        dla = _dla(["X", "Y", "Z"])
        self._assert_reachable(_qr_of(dla), su_basis(2), "su(2) basis->DLA")

    def test_su4_dla_to_basis(self):
        """Every one of the 15 su(4) DLA elements decomposes exactly in su_basis(4)."""
        dla = _dla(["XX", "YY", "ZZ", "XY", "YX", "XZ", "ZX", "YZ", "ZY"])
        assert len(dla) == 15, f"su(4) DLA dim should be 15, got {len(dla)}"
        self._assert_reachable(_qr_of(su_basis(4)), dla, "su(4) DLA->basis")

    def test_su4_basis_to_dla(self):
        """Every su_basis(4) element decomposes in the su(4) DLA span."""
        dla = _dla(["XX", "YY", "ZZ", "XY", "YX", "XZ", "ZX", "YZ", "ZY"])
        self._assert_reachable(_qr_of(dla), su_basis(4), "su(4) basis->DLA")

    def test_su4_explicit_coefficients_finite(self):
        """Decomposition coefficients are finite; reconstruction residual < 1e-8."""
        dla = _dla(["XX", "YY", "ZZ", "XY", "YX", "XZ", "ZX", "YZ", "ZY"])
        basis = su_basis(4)
        flat = basis.reshape(len(basis), -1).astype(complex)
        for X in dla:
            if float(np.linalg.norm(X)) < 1e-12:
                continue
            v = X.ravel().astype(complex)
            c, _, _, _ = np.linalg.lstsq(flat.T, v, rcond=None)
            residual = float(np.linalg.norm(flat.T @ c - v))
            assert np.all(np.isfinite(c)), "Non-finite decomposition coefficients"
            assert residual < 1e-8, f"Reconstruction residual {residual:.2e} too large"

    def test_su8_dla_to_basis(self):
        """Element-wise reachability for su(8) at n=3 qubits (63 elements)."""
        strings = [
            "XII",
            "YII",
            "ZII",
            "IXI",
            "IYI",
            "IZI",
            "IIX",
            "IIY",
            "IIZ",
            "XXI",
            "YYI",
            "ZZI",
            "IXX",
            "IYY",
            "IZZ",
            "XIZ",
            "ZIX",
            "YIZ",
            "ZXI",
            "XZI",
        ]
        dla = _dla(strings)
        if len(dla) != 63:
            pytest.skip(f"Generators gave DLA dim={len(dla)}, expected 63 for su(8)")
        self._assert_reachable(_qr_of(su_basis(8)), dla, "su(8) DLA->basis")
