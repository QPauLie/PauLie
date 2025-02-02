from paulie.common.algebras import get_algebra_generators
from paulie.helpers.drawing import get_algebra

def classifyAllAlgebras(size):
    print(f"Classification of dynamic Lia algebras size = {size}")
    print("--------------------------------------------------")
    for name in get_algebra():
        generators = get_algebra_generators(name)
        algebra = get_algebra(generators, size=size)
        print(f"name={name} algebra={algebra}")

if __name__ == '__main__':
    classifyAllAlgebras(8)






