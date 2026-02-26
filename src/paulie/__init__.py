"""
Public API
"""
import warnings

# Application
from .application.otoc import average_otoc
from .application.fourpoint import fourpoint
from .application.charges import non_commuting_charges
from .application.second_moment import second_moment
from .application.graph_complexity import average_graph_complexity
from .application.animation import animation_anti_commutation_graph
from .application.average_pauli_weight import (
    quantum_fourier_entropy,
    average_pauli_weight,
    get_pauli_weights
)
from .application.matrix_decomposition import (
    matrix_decomposition,
    matrix_decomposition_diagonal
)
from .application.pauli_compiler import (
    compile_target,
    OptimalPauliCompiler,
    PauliCompilerConfig,
    construct_universal_set,
)

# Classification
from .classifier.classification import (
    Morph,
    Classification
)

# Common
from .common.pauli_string_bitarray import PauliString
from .common.pauli_string_collection import PauliStringCollection
from .common.pauli_string_linear import PauliStringLinear
from .common.pauli_string_factory import (
    get_identity,
    get_pauli_string,
    get_single
)
from .common.two_local_generators import (
    get_lie_algebra,
    G_LIE,
    two_local_algebras
)

# Helpers
from .helpers.recording import RecordGraph
from .helpers.drawing import (
    animation_graph,
    plot_graph,
    plot_graph_by_nodes
)

__all__ = [
    # Application
    "average_otoc",
    "fourpoint",
    "non_commuting_charges",
    "second_moment",
    "average_graph_complexity",
    "animation_anti_commutation_graph",
    "quantum_fourier_entropy",
    "average_pauli_weight",
    "get_pauli_weights",
    "matrix_decomposition",
    "matrix_decomposition_diagonal",
    "compile_target",
    "OptimalPauliCompiler",
    "PauliCompilerConfig",
    "construct_universal_set",

    # Classification
    "Morph",
    "Classification",

    # Common
    "PauliString",
    "PauliStringCollection",
    "PauliStringLinear",
    "get_identity",
    "get_pauli_string",
    "get_single",
    "get_lie_algebra",
    "G_LIE",
    "two_local_algebras",

    # Helpers
    "RecordGraph",
    "animation_graph",
    "plot_graph",
    "plot_graph_by_nodes"
]

deprecated_symbols = {
    "g_lie": "G_LIE" # assume g_lie is old name for depreciation test
    # "old_name": "new_name"
}

def __getattr__(name):

    if name in deprecated_symbols:
        new_name = deprecated_symbols[name]
        warnings.warn(
            f"'{name}' is deprecated, use '{new_name}' instead",
            DeprecationWarning,
            stacklevel=2
        )
        return globals()[new_name]

    raise AttributeError(f"module {__name__} has no attribute {name}")
