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
from paulie.application.algebra_basis import so_basis, su_basis, sp_basis, u1_basis


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
                f"[B_{i}, B_{j}] not in span(basis) -- "
                "basis is not closed under the Lie bracket"
            )


def _check_full_rank(basis: np.ndarray) -> None:
    flat = basis.reshape(len(basis), -1)
    rank = np.linalg.matrix_rank(flat, tol=1e-10)
    assert rank == len(basis), f"Not full rank: {rank} != {len(basis)}"



def _summand_count(label: str) -> int:
    """Parse 'k*algebra(N)' or 'algebra(N)+algebra(N)+...' to get summand count."""
    m = re.match(r'^(\d+)\*[a-z]+\(\d+\)$', label.strip())
    return int(m.group(1)) if m else len(label.split('+'))

# ---------------------------------------------------------------------------
# so(N)
# ---------------------------------------------------------------------------

class TestSoBasis:
    @pytest.mark.parametrize("N", [2, 3, 4, 5, 6])
    def test_shape(self, N):
        assert so_basis(N).shape == (N * (N - 1) // 2, N, N)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_antisymmetric(self, N):
        for m in so_basis(N):
            np.testing.assert_allclose(m + m.T, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_full_rank(self, N):
        _check_full_rank(so_basis(N))

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_bracket_closure(self, N):
        _check_bracket_closure(so_basis(N))


# ---------------------------------------------------------------------------
# su(N)
# ---------------------------------------------------------------------------

class TestSuBasis:
    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_shape(self, N):
        assert su_basis(N).shape == (N * N - 1, N, N)

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_traceless(self, N):
        for m in su_basis(N):
            assert abs(np.trace(m)) < 1e-12

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_anti_hermitian(self, N):
        for m in su_basis(N):
            np.testing.assert_allclose(m + m.conj().T, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_full_rank(self, N):
        _check_full_rank(su_basis(N))

    @pytest.mark.parametrize("N", [2, 3, 4])
    def test_bracket_closure(self, N):
        _check_bracket_closure(su_basis(N))


# ---------------------------------------------------------------------------
# sp(N) -- N = half-dimension, matrices 2N x 2N
# ---------------------------------------------------------------------------

class TestSpBasis:
    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_shape(self, N):
        assert sp_basis(N).shape == (N * (2 * N + 1), 2 * N, 2 * N)

    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_symplectic_condition(self, N):
        size = 2 * N
        J = np.zeros((size, size), dtype=np.float64)
        J[:N, N:] =  np.eye(N)
        J[N:, :N] = -np.eye(N)
        for m in sp_basis(N):
            np.testing.assert_allclose(m.T @ J + J @ m, 0, atol=1e-12)

    @pytest.mark.parametrize("N", [1, 2, 3, 4])
    def test_full_rank(self, N):
        _check_full_rank(sp_basis(N))

    @pytest.mark.parametrize("N", [1, 2, 3])
    def test_bracket_closure(self, N):
        _check_bracket_closure(sp_basis(N))


# ---------------------------------------------------------------------------
# Integration tests via the public API
# ---------------------------------------------------------------------------

class TestGetAlgebraBasis:

    def test_b1_sp4_shape(self):
        basis_list = p(["XY", "XZ"], n=4).get_algebra_basis()
        assert len(basis_list) == 1
        assert basis_list[0].shape == (36, 8, 8)

    def test_b1_sp4_symplectic(self):
        b = p(["XY", "XZ"], n=4).get_algebra_basis()[0]
        N = 4
        J = np.zeros((8, 8))
        J[:N, N:] = np.eye(N)
        J[N:, :N] = -np.eye(N)
        for m in b:
            np.testing.assert_allclose(m.T @ J + J @ m, 0, atol=1e-10)

    def test_b1_sp4_bracket_closure(self):
        _check_bracket_closure(p(["XY", "XZ"], n=4).get_algebra_basis()[0])

    def test_a_type_so_antisymmetric(self):
        gens = ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]
        for b in p(gens).get_algebra_basis():
            for m in b:
                np.testing.assert_allclose(m + m.T, 0, atol=1e-10)

    def test_a_type_bracket_closure(self):
        gens = ["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]
        for b in p(gens).get_algebra_basis():
            _check_bracket_closure(b)

    @pytest.mark.parametrize("gens,n", [
        (["XY", "XZ"], 4),
        (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
        (["XZ", "ZX"], 2),
        (["XX", "YY", "ZZ"], None),
    ])
    def test_total_dim_matches_get_dla_dim(self, gens, n):
        """sum of basis sizes must equal get_dla_dim()."""
        ps = p(gens, n=n) if n else p(gens)
        assert sum(len(b) for b in ps.get_algebra_basis()) == ps.get_dla_dim()

    @pytest.mark.parametrize("gens,n", [
        (["XY", "XZ"], 4),
        (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
    ])
    def test_summand_count_matches_subalgebras(self, gens, n):
        ps = p(gens, n=n) if n else p(gens)
        assert len(ps.get_algebra_basis()) == _summand_count(ps.get_algebra())

    @pytest.mark.parametrize("gens,n", [
        (["XY", "XZ"], 4),
        (["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"], None),
        (["XZ", "ZX"], 2),
    ])
    def test_each_summand_full_rank(self, gens, n):
        ps = p(gens, n=n) if n else p(gens)
        for b in ps.get_algebra_basis():
            _check_full_rank(b)

    @pytest.mark.parametrize("gens,n", [
        (["XY", "XZ"], 4),
        (["XZ", "ZX"], 2),
        (["XX", "YY"], 2),
    ])
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
        dla_rank = np.linalg.matrix_rank(
            matrices.reshape(len(matrices), -1), tol=1e-10
        )
        basis_dim = sum(len(b) for b in ps.get_algebra_basis())
        assert basis_dim == dla_rank, (
            f"gens={gens}: basis dim {basis_dim} != DLA rank {dla_rank}"
        )
