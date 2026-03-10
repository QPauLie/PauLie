"""
Exceptions
"""

class PaulieError(Exception):
    """Base exception for all Paulie errors."""

class ParsingError(PaulieError):
    """Raised for parsing-related errors."""

class ValidationError(PaulieError):
    """Raised for validation errors."""

class PauliStringLinearException(PaulieError):
    """Exception for the linear combination of Pauli strings class."""

class PauliStringCollectionException(PaulieError):
    """Exception for the Pauli string collection class."""
