"""
Common types for classifier
"""

import enum

class TypeAlgebra(enum.Enum):
    """
    Lie algebra type.
    """

    U = 0
    SU = 1
    SP = 2
    SO = 3
