from paulie.application.classify import get_algebra


if __name__ == "__main__":
    #1
    print("Example III.1")
    alg = get_algebra(["XY", "XZ"])
    print(f"alg {"XY", "XZ"} = {alg}")
    alg = get_algebra(["IX", "XY"])
    print(f"alg {"IX", "XY"} = {alg}")
    #2
    print("Example III.2")
    alg = get_algebra(["XY", "XI", "IX"])
    print(f"alg {"XY", "XI", "IX"} = {alg}")
    alg = get_algebra(["XY", "XZ"])
    print(f"alg {"XY", "XZ"} = {alg}")
    #3
    print("Example III.3")
    alg = get_algebra(["XX", "YY"])
    print(f"alg {"XX", "YY"} = {alg}")
    alg = get_algebra(["XX", "YX"])
    print(f"alg {"XX", "YX"} = {alg}")

    print("Example III.4")
    alg = get_algebra(["XX", "YZ"])
    print(f"alg {"XX", "YZ"} = {alg}")
    alg = get_algebra(["YY", "ZX"])
    print(f"alg {"YY", "ZX"} = {alg}")
    alg = get_algebra(["XX", "YY"])
    print(f"alg {"XX", "YY"} = {alg}")

    print("Example III.5")
    alg = get_algebra(["ZZ", "YX", "XY"])
    print(f"alg {"ZZ", "YX", "XY"} = {alg}")
    alg = get_algebra(["XX", "YZ", "ZY"])
    print(f"alg {"XX", "YZ", "ZY"} = {alg}")
    alg = get_algebra(["YY", "ZX", "XZ"])
    print(f"alg {"YY", "ZX", "XZ"} = {alg}")

    print("Example III.6")
    alg = get_algebra(["XX", "XZ", "IY"])
    print(f"alg {"XX", "XZ", "IY"} = {alg}")
    alg = get_algebra(["XY", "XZ", "IX"])
    print(f"alg {"XY", "XZ", "IX"} = {alg}")

    print("Example III.7")
    alg = get_algebra(["XY", "YX"])
    print(f"alg {"XY", "YX"} = {alg}")
    alg = get_algebra(["XY", "YZ"])
    print(f"alg {"XY", "YZ"} = {alg}")

    print("Example III.8")
    alg = get_algebra(["XX", "YY", "ZZ", "ZY"])
    print(f"alg {"XX", "YY", "ZZ", "ZY"} = {alg}")
