import re
import numpy as np
from paulie.helpers._lie_bases import (
    get_so_basis,
    get_su_basis,
    get_sp_basis,
)


def parse_algebra_label(label: str):
    """Parsuje etykiety typu '2^n1 family(N)'."""
    # Wyciągamy liczbe skladnikow prostych (n1)
    n1 = 0
    if "^" in label:
        match_n1 = re.match(r"2\^(\d+)", label)
        if match_n1:
            n1 = int(match_n1.group(1))
        label = label.split()[-1]

    # Parse family and dimension
    match_fam = re.match(r"(so|su|sp)\((\d+)\)", label)
    if not match_fam:
        raise ValueError(f"Unknown algebra format: {label}")

    family = match_fam.group(1)
    N = int(match_fam.group(2))
    return n1, family, N


def get_algebra_basis_impl(label: str) -> list:
    """Implementacja podzialu na skladniki proste."""
    n1, family, N = parse_algebra_label(label)
    num_summands = 2**n1

    # Generujemy baze dla jednego skladnika
    if family == "so":
        single_basis = get_so_basis(N)
    elif family == "su":
        single_basis = get_su_basis(N)
    elif family == "sp":
        single_basis = get_sp_basis(N)
    else:
        raise ValueError(f"Unsupported family: {family}")

    # Zwracamy podział po składnikach prostych
    return [single_basis.copy() for _ in range(num_summands)]
