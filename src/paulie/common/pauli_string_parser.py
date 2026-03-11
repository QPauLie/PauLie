"""
    Module to parse Pauli strings.
"""
LOWCASE = "_"
SIZE = "s"
GATES = {"I", "X", "Y", "Z"}
TOKENS = GATES.copy()
TOKENS.add(LOWCASE)
TOKENS.add(SIZE)


def _is_token(char: str) -> bool:
    """
    Check if a character is a valid token.

    Args:
        char (str): Char for checking in tokens.
    Returns: bool
        True if char is token.
    """
    return char in TOKENS


def _is_number(char: str) -> bool:
    """
    Check if a character is a number, raise exception if invalid token.

    Args:
        char (str): Char to check.
    Returns: bool
        True if char is a number.
    Raises:
        ValueError: If the input char format is invalid.
    """
    try:
        int(char)
        return True
    except ValueError as e:
        if char not in TOKENS:
            raise ValueError("Invalid pauli string: unexpected character") from e
        return False


def _to_int(position: str) -> int:
    """
    Convert string to int, raise exception if invalid.

    Args:
        position (str): String representation of number.
    Returns: int
        Integer representation of number.
    Raises:
        ValueError: If the input position format is invalid.
    """
    try:
        return int(position)
    except ValueError as e:
        raise ValueError("Invalid pauli string: position must be a number") from e


def pauli_string_parser(pauli_string: str) -> str:
    """
    Parse a compact Pauli string representation and return the expanded form.

    For example, `X_4s10` is expanded to `IIIXIIIIII` and `ZYX_4s10` is expanded to `ZYIXIIIIII`.

    Args:
        pauli_string (str): Compact string representation of the Pauli string.
    Returns: str
        Expanded Pauli string with I's in unspecified positions.
    Raises:
        ValueError: If the input string format is invalid.
    """
    new_pauli_string = ""
    i = 0
    size = None
    # Extract size if specified
    try:
        index = pauli_string.find(SIZE)
        if index != -1:
            if index == len(pauli_string) - 1:
                raise ValueError("Invalid pauli string: missing size value after 's'")
            size_string = pauli_string[index+1:]
            size = _to_int(size_string)
            pauli_string = pauli_string[:index]
    except ValueError as e:
        raise e
    # Parse the operations
    while i < len(pauli_string):
        if pauli_string[i] not in GATES:
            raise ValueError("Invalid pauli string: unexpected character")

        token = pauli_string[i]
        if i < len(pauli_string) - 2 and pauli_string[i+1] == LOWCASE:
            # Handle positioned operator (e.g., X_4)
            m = i
            i += 2
            p = ""
            while i < len(pauli_string):
                if _is_token(pauli_string[i]):
                    break
                if _is_number(pauli_string[i]):
                    p += pauli_string[i]
                    i += 1
            position = _to_int(p)
            # Check for valid position
            if position - len(new_pauli_string) - 1 < 0:
                raise ValueError("Invalid pauli string: invalid position order")
            # Create string with I's up to position and add the operator
            token = "".join(["I" for _ in range(position - len(new_pauli_string) - 1)])
            token += pauli_string[m]
        else:
            i += 1
        new_pauli_string += token
    # Add padding I's if size is specified
    if size is not None:
        if size < len(new_pauli_string):
            raise ValueError("Invalid pauli string: size too small for operators")
        padding = size - len(new_pauli_string)
        new_pauli_string += "".join(["I" for _ in range(padding)])
    return new_pauli_string
