from paulie.common.pauli import get_I, get_Y, get_pauli_array, get_pauli_string, inc_pauli_array, is_sub_in_array
from paulie.common.algebras import get_algebra_generators
from paulie.common.nested import get_nested_nodes_in_algebra


def gen_k_local(n, p, converter = None, used = []):
    if n < len(p)//2:
        raise ValueError(f"Size must be greater than {len(p)//2}")

    #if len(p) != 4:
    #    raise ValueError("Pauli string should be equal to 2")

    np = n - len(p)//2
    k = 0
    while k <= np:
       left = get_I(k)
       left.extend(p)
       right = get_I(np-k)
       left.extend(right)
       k = k + 1
       if left in used:
           continue
       if converter is None:
           yield left
       else:
           yield converter(left)
       used.append(left)


def gen_k_local_by_string(n, pauli_string, used = []):
    yield from gen_k_local(n, get_pauli_array(pauli_string), used=used)


def gen_k_local_string(n, pauli_string, used = []):
    yield from gen_k_local(n, get_pauli_array(pauli_string), get_pauli_string, used=used)


def gen_k_local_ext(n, p, converter = None, used = []):

    if n < 2:
        raise ValueError("Size must be greater than 1")

    #if len(p) != 4:
    #    raise ValueError("Pauli string should be equal to 2")

    np = n - len(p)//2
    k = 0
    while(k <= np):
        left = get_I(k)
        left_one = get_Y(k)
        right = get_I(np-k)
        right_one = get_Y(np-k)
        isFinish = False  
        while isFinish is False:
             gen = get_I(0)
             gen.extend(left)
             gen.extend(p)
             gen.extend(right)
             if right == right_one:
                 if left == left_one:
                     isFinish = True
                 else:  
                     left = inc_pauli_array(left)
             else:
                right = inc_pauli_array(right)
             if len(used) > 0:
                 isUsed = False
                 for u in used:
                     if is_sub_in_array(u, gen):
                        isUsed = True
                        break
                 if isUsed:
                     continue
             if converter is None:
                 yield gen
             else:
                 yield converter(gen)
        k = k + 1


def gen_k_local_ext_by_string(n, pauli_string):
    yield from gen_k_local_ext(n, get_pauli_array(pauli_string))


def gen_k_local_ext_string(n, pauli_string):
    yield from gen_k_local_ext(n, get_pauli_array(pauli_string), get_pauli_string)


def gen_k_local_algebra_generators(n, name):
    generators = get_algebra_generators(name)
    used = []
    for g in generators:
        yield from gen_k_local_by_string(n, g, used=used)


def get_k_local_algebra_generators(n, name):
    generators = []
    for g in gen_k_local_algebra_generators(n, name):
        generators.append(g)
    return generators


def gen_k_local_string_algebra_generators(n, name):
    generators = get_algebra_generators(name)
    used = []
    for g in generators:
        yield from gen_k_local_string(n, g, used=used)


def get_k_local_string_algebra_generators(n: int, name: str) -> list:
    generators = []
    for g in gen_k_local_string_algebra_generators(n, name):
        generators.append(g)
    return generators


def gen_k_local_generators(n, generators):
    used = []
    for g in generators:
        yield from gen_k_local_by_string(n, g, used=used)


def get_k_local_generators(n, generators):
    gens = []
    for g in gen_k_local_generators(n, generators):
        gens.append(g)
    return gens


def gen_k_local_string_generators(n, generators):
    used = []
    for g in generators:
        yield from gen_k_local_string(n, g, used=used)


def get_k_local_string_generators(n, generators):
    gens = []
    for g in gen_k_local_string_generators(n, generators):
        gens.append(g)
    return gens


def gen_k_local_nested_nodes_in_agebra_generator(n, name):
    nested = get_nested_nodes_in_algebra(name)
    used = []
    for node in nested:
        yield from gen_k_local_ext(n, node, used = used)
        used.append(node)


def get_k_local_nested_nodes_in_algebra_generator(n, name):
    nodes = []
    for node in gen_k_local_nested_nodes_in_agebra_generator(n, name):
        nodes.append(node)
    return nodes