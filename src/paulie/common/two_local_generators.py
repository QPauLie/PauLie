"""
 Two local algebra
"""
#A dictionary representing 2 local generators
G_LIE: dict[str, list[str]] = {
    "a0": ["XX"],
    "a1": ["XY"],
    "a2": ["XY", "YX"],
    "a3": ["XX", "YZ"],
    "a4": ["XX", "YY"],
    "a5": ["XY", "YZ"],
    "a6": ["XX", "YZ", "ZY"],
    "a7": ["XX", "YY", "ZZ"],
    "a8": ["XX", "XZ"],
    "a9": ["XY", "XZ"],
    "a10": ["XY", "YZ", "ZX"],
    "a11": ["XY", "YX", "YZ"],
    "a12": ["XX", "XY", "YZ"],
    "a13": ["XX", "YY", "YZ"],
    "a14": ["XX", "YY", "XY"],
    "a15": ["XX", "XY", "XZ"],
    "a16": ["XY", "YX", "YZ", "ZY"],
    "a17": ["XX", "XY", "ZX"],
    "a18": ["XX", "XZ", "YY", "ZY"],
    "a19": ["XX", "XY", "ZX", "YZ"],
    "a20": ["XX", "YY", "ZZ", "ZY"],
    "a21": ["XX", "YY", "XY", "ZX"],
    "a22": ["XX", "XY", "XZ", "YX"],
    "b0": ["XI", "IX"],
    "b1": ["XX", "XI", "IX"],
    "b2": ["XY", "XI", "IX"],
    "b3": ["XI", "YI", "IX", "IY"],
    "b4": ["XX", "XY", "XZ", "XI", "IX", "IY", "IZ"],
}

def _a3(n):
    """
    a3 algebra
    """
    if n % 8 == 0:
        return f"4*so({2**(n-2)})"
    if n % 8 == 1 or n % 8 == 7:
        return f"so({2**(n-1)})"
    if n % 8 == 2 or n % 8 == 6:
        return f"2*su({2**(n-2)})"
    if n % 8 == 3 or n % 8 == 5:
        return f"sp({2**(n-2)})"
    if n % 8 == 4:
        return f"4*sp({2**(n-3)})"
    return None

def _a5(n):
    """
    a5 algebra
    """
    if n % 6 == 0:
        return f"4*so({2**(n-2)})"
    if n % 6 == 1 or n % 6 == 5:
        return f"so({2**(n-1)})"
    if n % 6 == 2 or n % 6 == 4:
        return f"2*su({2**(n-2)})"
    if n % 6 == 3:
        return f"sp({2**(n-2)})"
    return None

def _a6(n):
    """
    a6 algebra
    """
    if n % 2 == 1: #odd
        return f"su({2**(n-1)})"
    else: #even
        return f"4*su({2**(n-2)})"
    return None

def _a7(n):
    """
    a7 algebra
    """
    return _a6(n)

def _a10(n):
    """
    a10 algebra
    """
    return _a6(n)

def two_local_algebras(n):
    """
    A dictionary of the DLAs n>=3
    """

    return  {
    "a0": f"{n-1}*u(1)",
    "a1": f"so({n})",
    "a2": f"so({n})+so({n})", #same as a_4
    "a3": _a3(n),
    "a4": f"so({n})+so({n})",
    "a5": _a5(n),
    "a6": _a6(n),
    "a7": _a7(n),
    "a8": f"so({2*n-1})",
    "a9": f"sp({2**(n-2)})",
    "a10": _a10(n),
    "a11": f"so({2**n})", # a_16
    "a12": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "a13": f"su({2**(n-1)}) + su({2**(n-1)})", # a_20, a_15
    "a14": f"so({2*n})",
    "a15": f"su({2**(n-1)}) + su({2**(n-1)})", # a_20, a_15
    "a16": f"so({2**n})", # a_16
    "a17": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "a18": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "a19": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "a20": f"su({2**(n-1)}) + su({2**(n-1)})", # a_20, a_15
    "a21": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "a22": f"su({2**n})", # a_17, a_18, a_19, a_21, a_22
    "b0": f"{n}*u(1)",
    "b1": f"{2*n-1}*u(1)",
    "b2": f"sp({2**(n-2)}) + u(1)",
    "b3": f"{n}*so(3)",
    "b4": f"su({2**(n-1)}) + su({2**(n-1)}) + u(1)" ,
}



def get_lie_algebra(name:str) -> list[str]:
    """Returns generators of Lie algebra by name. 
       Args:
            name (str): name of Lie algebra
    """
    return G_LIE[name]

def get_lie_algebras() -> dict[str, list[str]]:
    """Returns the dictionary of generators of Lie algebras."""
    return G_LIE
