from paulie.common.algebras import get_lie_algebra
from paulie.application.classify import get_algebra


def classify_all_algebras(size: int):
    print(f"Classification of dynamic Lie algebras size = {size}")
    print("--------------------------------------------------")

    for name, generators in get_lie_algebra().items():
        algebra = get_algebra(generators, size=size)
        print(f"name={name} algebra={algebra}")

if __name__ == '__main__':
    classify_all_algebras(8)
