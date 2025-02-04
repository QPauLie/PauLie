from paulie.helpers.transform import app_transform_to_canonics


# Get algebra
# generators - list of generators
# size - Generator extensions to size size
def get_algebra(generators, size=0):
    canonics = app_transform_to_canonics(generators, size=size)
    algebras = []
    for canonic in canonics:
        algebras.append(canonic["shape"].get_algebra())
    return " + ".join(algebras)



