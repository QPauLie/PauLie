"""
Test of an AI-detected error
"""

import pytest
from paulie import get_pauli_string as p

# ---------------------------------------------------------------------------
# Regression: stars whose canonical graph has only legs of length 1
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "generators, algebra, dim",
    [
        (["ZZ", "XI", "IX"], "2*so(3)", 6),                       # K_{1,2}: so(4)
        (["ZZZ", "XII", "IXI", "IIX"], "4*so(3)", 12),            # K_{1,3}
        (["ZZZ", "XII", "IXI", "IIX", "XXX"], "4*so(3)", 12),     # K_{1,3} + dependent leaf
        (["ZZZZ", "XIII", "IXII", "IIXI", "IIIX"], "8*so(3)", 24),  # K_{1,4}
    ],
)
def test_star_with_only_short_legs(generators: list[str], algebra: str, dim: int) -> None:
    """A star with k >= 3 length-1 legs is 2^(k-1) copies of so(3) = su(2),
    not 2^(k-1) copies of so(2). Previously reported as `4*so(2)` with dim 4."""
    gens = p(generators)
    assert gens.get_algebra() == algebra
    assert gens.get_dla_dim() == dim
