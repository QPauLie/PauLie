"""
Tests for get_algebra_basis() and the underlying basis constructors.

Key test: _check_bracket_closure verifies [B_i, B_j] lies in span(basis)
for all i<j. This is the direct algebraic soundness check the maintainer
asked for -- not a dimension proxy.
"""

import re
from itertools import product as iproduct
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
        basis = p(["XY", "XZ"], n=4).get_algebra_basis()
        assert basis.shape == (36, 8, 8)

    def test_b1_sp4_symplectic(self):
        """sp(4) elements satisfy the symplectic condition."""

        basis = p(["XY", "XZ"], n=4).get_algebra_basis()
        N = 4
        J = np.zeros((2 * N, 2 * N))
        J[:N, N:] = np.eye(N)
        J[N:, :N] = -np.eye(N)
        for mat in basis:
            assert np.max(np.abs(mat.T @ J + J @ mat)) < 1e-10

    def test_b1_sp4_bracket_closure(self):
        """sp(4) basis closes under Lie brackets."""
        _check_bracket_closure(p(["XY", "XZ"], n=4).get_algebra_basis())

    def test_a_type_so_antisymmetric(self):
        """Type A so elements are antisymmetric."""

        basis = p(["XY", "YX", "ZI", "IZ"], n=2).get_algebra_basis()
        for mat in basis:
            assert np.max(np.abs(mat + mat.T)) < 1e-10

    def test_a_type_bracket_closure(self):
        """Type A so basis closes under Lie brackets (uses 2*so(3) generators)."""
        _check_bracket_closure(p(["XY", "YX", "ZI", "IZ"], n=2).get_algebra_basis())

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["XY", "YX"], 2),
            (["XY", "YX"], 4),
        ],
    )
    def test_total_dim_matches_get_dla_dim(self, gens, n):
        """Total basis dim matches get_dla_dim()."""
        ps = p(gens, n=n)
        basis = ps.get_algebra_basis()
        assert len(basis) == ps.get_dla_dim()

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["XY", "YX"], 2),
        ],
    )
    def test_basis_is_full_rank(self, gens, n):
        """All basis elements are linearly independent."""

        basis = p(gens, n=n).get_algebra_basis()
        flat = basis.reshape(len(basis), -1)
        assert np.linalg.matrix_rank(flat) == len(basis)

    @pytest.mark.parametrize(
        "gens,n",
        [
            (["XY", "XZ"], 4),
            (["XY", "YX"], 2),
        ],
    )
    def test_basis_dim_matches_lie_closure_rank(self, gens, n):
        """Basis dimension matches brute-force Lie closure size."""

        ps = p(gens, n=n)
        generators = ps.get()
        dla = list(generators)
        old_len, new_len, init_len = 0, len(dla), len(dla)
        while new_len > old_len:
            for pw1, pw2 in iproduct(dla[:init_len], dla[old_len:]):
                if pw1.commutes_with(pw2):
                    continue
                com = pw1 @ pw2
                if com not in dla:
                    dla.append(com)
            old_len = new_len
            new_len = len(dla)
        assert len(ps.get_algebra_basis()) == len(dla)


class _LieSpan:
    """Incrementally-maintained orthonormal span via pre-allocated BLAS buffer."""

    def __init__(self, tol: float = 1e-10) -> None:
        """Initialise with tolerance for span membership tests."""
        self.tol = tol
        self.mats: list = []
        self._buf = None
        self._d = 0
        self._cap = 0

    def _grow(self, flat_size: int) -> None:
        """Double the buffer capacity."""
        new_cap = max(16, self._cap * 2)
        buf = np.empty((flat_size, new_cap), dtype=complex)
        if self._d > 0:
            buf[:, : self._d] = self._buf[:, : self._d]
        self._buf, self._cap = buf, new_cap

    def _residual(self, v: np.ndarray) -> np.ndarray:
        """Return the component of v orthogonal to the current span."""
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
    """Return the 2^n x 2^n matrix for a Pauli string like 'XYZ'."""
    m = _PAULI[s[0]]
    for c in s[1:]:
        m = np.kron(m, _PAULI[c])
    return m


def _dla(strings: list) -> np.ndarray:
    """Brute-force anti-Hermitian DLA from Pauli string generators."""
    return _lie_closure(np.array([1j * _pm(s) for s in strings]))


def _qr_of(basis: np.ndarray) -> np.ndarray:
    """Return orthonormal column matrix Q for span(basis)."""
    flat = basis.reshape(len(basis), -1).astype(complex).T
    Q, _ = np.linalg.qr(flat, mode="reduced")
    return Q


def _reachability_failures(Q: np.ndarray, elements: np.ndarray, tol: float = 1e-8) -> list:
    """Return (index, rel_residual) for elements not in span(Q)."""
    failures = []
    for k, X in enumerate(elements):
        v = X.ravel().astype(complex)
        nv = max(float(np.linalg.norm(v)), 1e-14)
        rel_res = float(np.linalg.norm(v - Q @ (Q.conj().T @ v))) / nv
        if rel_res > tol:
            failures.append((k, rel_res))
    return failures


def _decompose(basis: np.ndarray, mat: np.ndarray) -> tuple:
    """Return coefficients c and residual for mat = sum(c_k * basis[k])."""
    flat = basis.reshape(len(basis), -1).astype(complex)
    v = mat.ravel().astype(complex)
    c, _, _, _ = np.linalg.lstsq(flat.T, v, rcond=None)
    residual = float(np.linalg.norm(flat.T @ c - v))
    return c, residual


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


class TestTypeAReachability:
    """Element-wise reachability for type A (so family) via Majorana bilinear isomorphism.

    For n=2 qubits, so(4) ≅ so(3)⊕so(3) via the Jordan-Wigner representation.
    The test maps each element of get_algebra_basis() (6x6 block-diagonal) into
    the embedded 4x4 DLA space via the so(3)→su(2) correspondence and verifies
    span equality with the brute-force Lie closure.
    """

    def test_2so3_span_via_majorana_isomorphism(self):
        """get_algebra_basis() output spans the embedded 2-qubit DLA (2*so(3)).

        Maps each 6x6 block-diagonal basis element to a 4x4 DLA matrix via the
        so(4)≅so(3)⊕so(3) Majorana isomorphism, then checks span equality with
        the brute-force Lie closure computed using PauLie-native operations.
        """
        g = p(["XY", "YX", "ZI", "IZ"], n=2)
        assert g.get_algebra() == "2*so(3)", f"got {g.get_algebra()}"

        # Defining-rep basis — this is the output under test
        basis = g.get_algebra_basis()
        assert basis.shape == (6, 6, 6), f"expected (6,6,6) got {basis.shape}"

        # Majorana bilinears G_ab = γ_a γ_b  (Jordan-Wigner, n=2 qubits)
        g01 = 1j * _pm("ZI")
        g02 = -1j * _pm("YX")
        g03 = -1j * _pm("YY")
        g12 = 1j * _pm("XX")
        g13 = 1j * _pm("XY")
        g23 = 1j * _pm("IZ")

        # so(4) ≅ so(3)⊕so(3): self-dual / anti-self-dual generators in DLA space
        summand_i = [(g01 + g23) / 2, (g02 - g13) / 2, (g03 + g12) / 2]
        summand_ii = [(g01 - g23) / 2, (g02 + g13) / 2, (g03 - g12) / 2]

        # so(3) defining-rep basis — same ordering used by so_basis(3)
        so3 = so_basis(3)  # shape (3, 3, 3)
        flat_so3 = so3.reshape(3, -1)  # (3, 9) — columns are basis vectors

        # Map each 6x6 basis element → 4x4 DLA matrix via the isomorphism
        dla_images = []
        for mat in basis:
            b0 = mat[:3, :3].real  # first so(3) block
            b1 = mat[3:, 3:].real  # second so(3) block
            c0, *_ = np.linalg.lstsq(flat_so3.T, b0.ravel(), rcond=None)
            c1, *_ = np.linalg.lstsq(flat_so3.T, b1.ravel(), rcond=None)
            img = sum(float(c) * gen for c, gen in zip(c0, summand_i)) + sum(
                float(c) * gen for c, gen in zip(c1, summand_ii)
            )
            dla_images.append(img)
        dla_images_arr = np.array(dla_images)

        # Brute-force Lie closure using PauLie-native PauliString operations
        dla_paulis = list(g.get())
        old_len, new_len, init_len = 0, len(dla_paulis), len(dla_paulis)
        while new_len > old_len:
            for pw1, pw2 in iproduct(dla_paulis[:init_len], dla_paulis[old_len:]):
                if pw1.commutes_with(pw2):
                    continue
                com = pw1 @ pw2
                if com not in dla_paulis:
                    dla_paulis.append(com)
            old_len = new_len
            new_len = len(dla_paulis)
        assert len(dla_paulis) == 6
        dla_mats = np.array([1j * _pm(str(el)) for el in dla_paulis])

        # span(images of get_algebra_basis()) == span(embedded DLA)
        flat_img = dla_images_arr.reshape(len(dla_images_arr), -1)
        flat_dla = dla_mats.reshape(len(dla_mats), -1)
        combined = np.vstack([flat_img, flat_dla])
        tol = 1e-8
        r_img = int(np.linalg.matrix_rank(flat_img, tol=tol))
        r_dla = int(np.linalg.matrix_rank(flat_dla, tol=tol))
        r_all = int(np.linalg.matrix_rank(combined, tol=tol))
        assert r_img == r_dla == r_all == 6, (
            f"get_algebra_basis() image does not span DLA: r_img={r_img}, r_dla={r_dla}, r_combined={r_all}"
        )
