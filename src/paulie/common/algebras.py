G_lie = {
    "a0": ["XX"], #0
    "a1": ["XY"], #1
    "a2": ["XY", "YX"], #2
    "a3": ["XX", "YZ"], #3
    "a4": ["XX", "YY"], #4
    "a5": ["XY", "YZ"], #5
    "a6": ["XX","YZ","ZY"],#6
    "a7": ["XX", "YY","ZZ"],#7
    "a8": ["XX", "XZ"],#8
    "a9": ["XY", "XZ"],#9
    "a10": ["XY", "YZ", "ZX"],#10
    "a11": ["XY", "YX", "YZ"],#11
    "a12": ["XX", "XY", "YZ"],#12
    "a13": ["XX", "YY", "YZ"],#13
    "a14": ["XX", "YY", "XY"],#14
    "a15": ["XX", "XY", "XZ"],#15
    "a16": ["XY", "YX", "YZ", "ZY"],#16
    "a17": ["XX", "XY", "ZX" ], #17
    "a18": ["XX", "XZ", "YY", "ZY"], #18
    "a19": ["XX", "XY", "ZX", "YZ"], #19
    "a20": ["XX", "YY", "ZZ", "ZY"],#20
    "a21": ["XX", "YY", "XY", "ZX"], #21
    "a22": [ "XX", "XY", "XZ", "YX"], #22
    "b0": [ "XI", "IX"], #b0
    "b1": [ "XX", "XI", "IX"], #b1
    "b2": [ "XY", "XI", "IX"],#b2
    "b3": [ "XI", "YI", "IX", "IY"], #b3
    "b4": [ "XX", "XY", "XZ", "XI", "IX", "IY", "IZ"], #b4
}


def get_algebra_generators(name: str) -> list[str]:
    return G_lie[name]


def get_algebra() -> dict[str, list[str]]:
    return G_lie
