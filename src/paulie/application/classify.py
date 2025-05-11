"""
Get algebra
"""
from paulie.common.pauli_string_collection import PauliStringCollection

def get_algebra(generators: PauliStringCollection):
    """
     Get algebra
    """
    classification = generators.get_class()
    return classification.get_algebra()
