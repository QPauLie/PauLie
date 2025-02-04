from paulie.application.classify import get_algebra


def classify_size_generators(generators):
    for size in range(3, 4):
        algebra = get_algebra(generators, size=size)
        print(f"size = {size} algebra = {algebra}")


if __name__ == "__main__":
    classify_size_generators(["XYI", "IXY", "YIX"])
    classify_size_generators(["XYI", "IXY"])





