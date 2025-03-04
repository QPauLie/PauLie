from paulie.helpers.classify_generators import classify_generators


# Get algebra
# generators - list of generators
# size - Generator extensions to size size
def get_algebra(generators, size=0):
    classification = classify_generators(generators, size=size)
    return classification.get_algebra()



