from paulie.common.algebras import get_algebra_generators, get_algebra_dict
from paulie.application.classify import get_algebra


def classify_all_algebras(size):
    print(f"Classification of dynamic Lia algebras size = {size}")
    print("--------------------------------------------------")
    for name in get_algebra_dict():
        generators = get_algebra_generators(name)
        algebra = get_algebra(generators, size=size)
        print(f"name={name} algebra={algebra}")

if __name__ == '__main__':
    classify_all_algebras(8)






