"""
Test public API
"""
import warnings
import paulie

allowed_dunders = {
    "__all__",
    "__name__",
    "__doc__",
    "__path__",
    "__file__",
    "__spec__",
    "__package__",
    "__loader__",
    "__cached__",
    "__builtins__",
    "__getattr__",
}
allowed_internals = {
    "deprecated_symbols",
    "warnings",
    "exceptions"
}
subpackages = {
    "application",
    "common",
    "helpers",
    "classifier"
}
public_symbols = [
    # Exceptions
    "PaulieError",
    "ParsingError",
    "ValidationError",
    "PauliStringLinearException",
    "PauliStringCollectionException",

    # Application
    "average_otoc",
    "fourpoint",
    "non_commuting_charges",
    "second_moment",
    "get_optimal_su_2_n_generators",
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

def test_public_api_exports():
    """
    Test public API
    """
    exported = set(dir(paulie))
    expected = set(public_symbols) | allowed_dunders | allowed_internals | subpackages
    assert exported == expected

def test_deprecated_symbols():
    """
    Test depreciation symbols
    """
    deprecated_symbols = ["g_lie"]
    for old_name in deprecated_symbols:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _ = getattr(paulie, old_name)
            assert any(issubclass(warn.category, DeprecationWarning) for warn in w)
