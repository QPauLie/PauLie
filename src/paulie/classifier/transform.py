import queue
import random

from paulie.classifier.shape import Shape
from paulie.common.pauli import is_commutate, multi_pauli_arrays
from paulie.graphs.subgraphs import get_subgraphs
from paulie.helpers.printing import print_node, print_nodes
from paulie.helpers.recording import recording_graph


def get_lits(lighting, nodes):
    """Return highlighted vertices (connected to the selected vertex)."""
    lits = []
    for node in nodes:
        if node != lighting:
           if is_commutate(lighting, node) is False:
               lits.append(node)
    return lits


def find_center_and_lits(nodes):
    """Find a center with maximum connections. And bring back these connections."""
    center = nodes[0]
    center_lits = []
    for node in nodes:
        lits = get_lits(node, nodes)
        if len(lits) > len(center_lits):
            center = node
            center_lits = lits
    return center, center_lits


def append_next_nodes(prepared_nodes, nodes):
    """Adding the next node according to the principle of having a connection with the previous one.
    
    If there are several connected ones, then we insert them after the first one. (To reduce the risk of graph
    reassembly).
    """
    for node in nodes:
        if node in prepared_nodes:
            nodes.remove(node)
            continue

        lits = get_lits(node, prepared_nodes)
        if len(lits) == 0:
            continue
        if len(lits) > 1:
           min_index = len(prepared_nodes)
           for lit in lits:
               index = prepared_nodes.index(lit)
               if index < min_index:
                   min_index = index
           prepared_nodes.insert(min_index + 1, node)
        else:
            prepared_nodes.append(node)
        nodes.remove(node)
        return


def prepare_nodes(nodes, debug=False):
    """Sorting nodes in order from the center and then by connections."""
    nodes = sorted(nodes)
    prepared_nodes = []
    center, lits =  find_center_and_lits(nodes)
    print_node(debug, center,"center")
    print_nodes(debug, lits, "lits")
    nodes.remove(center)
    prepared_nodes.append(center)

    for lit in lits:
        nodes.remove(lit)
        if lit not in prepared_nodes:
            prepared_nodes.append(lit)

    while len(nodes) > 0:
         append_next_nodes(prepared_nodes, nodes)

    return prepared_nodes


def get_lits_or_append(shape, center, lighting, canonic, debug):
    """Get lits nodes and if possible then append."""
    lits = get_lits(lighting, canonic)
    if len(lits) != 1:
        return canonic, False, lits

    if lits[0] == center:
        return canonic, shape.append_is_canonical(lighting, center), lits
    else:
        if shape.is_end(lits[0]):
            return canonic, shape.append_is_canonical(lighting, lits[0]), lits
    return canonic, False, lits


def reconstruct_if_not_connection(shape, lighting, canonic, debug):
    """Reconstruct cononic graph if there are no connections."""
    lits = get_lits(lighting, canonic)
    is_attachable = len(lits) > 0

    while is_attachable is False:
        v, lit = shape.pop_vertix()
        if v is None:
            shape.clear_prohibiteds_for_pop()
            return canonic, False

        canonic.remove(v)
        nodes, is_appended = append_node(shape, v, canonic, debug)

        if is_appended:
            last = canonic[len(canonic)]
            if is_commutate(lighting, last) is False:
                shape.clear_prohibiteds_for_pop()
                return  canonic, is_attachable
            else:
                shape.append_prohibiteds_for_pop(v)
        else:
            shape.append_prohibiteds_for_pop(v)
        shape.reset_prohibited()
        shape.append_is_canonical(v, lit)

    shape.clear_prohibiteds_for_pop()
    return canonic, is_attachable


def append_to_canonic(shape, lighting, canonic, debug, record=None):
    """Append to cononic graph."""

    if len(canonic) == 0:
        canonic.append(lighting)
        return canonic, True

    center = canonic[0]
    if len(canonic) < 2:
        shape.append_is_canonical(lighting, center)
        canonic.append(lighting)
        return canonic, True

    canonic, is_attachable = reconstruct_if_not_connection(shape, lighting, canonic, debug)
    if is_attachable is False:
        return canonic, False

    used = []
    q = queue.Queue()
    q.put(lighting)

    while q.empty() is False:
        lighting = q.get()
        if record is not None:
            r = canonic.copy()
            r.append(lighting)
            recording_graph(record, r)

        canonic, is_appended, lits = get_lits_or_append(shape, center, lighting, canonic, debug)
        if is_appended:
            canonic.append(lighting)
            return canonic, True
        for lit in lits:
             contractor = multi_pauli_arrays(lighting, lit)
             if contractor in canonic:
                 return canonic, True

             if contractor not in used:
                 used.append(contractor)
                 q.put(contractor)

    return canonic, False


def get_untoggleables(lighting, lits):
    """Find untoggleable nodes in lits."""
    if len(lits) == 1:
        return None, None
    for litA in lits:
        contractor = multi_pauli_arrays(lighting, litA)
        for litB in lits:
             if litB != litA:
                 if is_commutate(contractor, litB) is False:
                     return litA, litB
    return None, None


def deconstruct_canonic_until_connect(shape, canonic, untoggleableA, untoggleableB, debug=False):
    """Deconstruct the canonical graph up to the first connection."""
    route = shape.get_route_to_end(untoggleableA)
    print_nodes(debug, route, "route A")
    if len(route) > 1 or len(route) == 0:
        routeB = shape.get_route_to_end(untoggleableB)
        print_nodes(debug, routeB, "route B")
        if len(routeB) != 0 and len(routeB) < len(route) or len(route) == 0:
            route = routeB
    reverseNodes = route[::-1]
    while len(route) > 0:
        v = route[0]
        route.remove(v)
        shape.remove_vertix(v)
        canonic.remove(v)
    return reverseNodes


def reconstruct_and_append_to_toggleable(shape, lighting, canonic, debug, record):
    """Reconstruct the canonical graph if there are untoggleable nodes and append."""
    print_node(debug, lighting, "reconstruct and append to toggleable")
    lenNodes = len(canonic)
    lits = get_lits(lighting, canonic)

    untoggleableA, untoggleableB = get_untoggleables(lighting, lits)
    if untoggleableA is not None:
        print_node(debug, untoggleableA, "untoggleableA")
        print_node(debug, untoggleableB, "untoggleableB")
        deconstructed_nodes = deconstruct_canonic_until_connect(shape, canonic, untoggleableA, untoggleableB, debug)
        print_nodes(debug, deconstructed_nodes, "deconstucted nodes")
        print_nodes(debug, canonic, "remain")
        if len(deconstructed_nodes) == 0:
            return canonic, False
        canonic, is_appended = append_node(shape, lighting, canonic, debug, record)
        if is_appended is False:
            return canonic, False
        for node in deconstructed_nodes:
            canonic, is_appended = append_node(shape, node, canonic, debug, record)
            if is_appended is False:
                return canonic, False

    if lenNodes < len(canonic):
        return canonic, True
    return canonic, False


def resort_nodes(center, nodes, debug):
    new_nodes = []
    random.shuffle(nodes)
    lits = get_lits(center, nodes)
    for node in nodes:
        if node in lits:
           new_nodes.append(node)
           nodes.remove(node)
    while len(nodes) > 0:
        for node in nodes:
            lits = get_lits(node, new_nodes)
            if len(lits) > 0:
                new_nodes.append(node)
                nodes.remove(node)
                break
    return new_nodes


def reconfing_and_append(shape, lighting, canonic, debug, record):
    """Reconfigure the canonical graph and append a node."""
    print_node(debug, lighting, "reconfig")
    canonic, is_appended = reconstruct_and_append_to_toggleable(shape, lighting, canonic, debug, record)
    if is_appended:
        return canonic, True

    nodes = shape.decrease_len_line()
    print_nodes(debug, nodes, "decreased")
    if len(nodes) > 1:
        center = shape.get_center()
        shape.set_prohibited_line(True)
        canonic = [center]
        return canonic, False
    else:
        for node in nodes:
            canonic.remove(node)

    mix = False
    for node in nodes:
        if mix is False:
            lits = get_lits(lighting, canonic)
            if len(lits) > 0:
                 canonic, is_appended = append_node(shape, lighting, canonic, debug, record)
                 mix = True
        canonic, is_appended = append_node(shape, node, canonic, debug, record)
        if is_appended is False:
            return canonic, False

    if mix is False:
        canonic, is_appended = append_node(shape, lighting, canonic, debug, record)

    return canonic, is_appended


def append_node(shape, lighting, canonic, debug=False, record=None):
    """Appending a node to a cononic graph."""
    print_node(debug, lighting, "append node")
    if lighting in canonic:
        return canonic, True

    canonic, is_appended = append_to_canonic(shape, lighting, canonic, debug, record)
    if is_appended is False:
        canonic, is_appended = reconfing_and_append(shape, lighting, canonic, debug, record)
    
    return canonic, is_appended


def transform_to_canonic(nodes, debug=False, record=None):
    """Transform a connected graph to a cononic type."""

    shape = Shape(debug)
    if len(nodes) == 0:
        return shape, nodes

    canonic = []
    nodes = prepare_nodes(nodes, debug)
    center = nodes[0]
    shape.set_center(nodes[0])
    original = nodes.copy()
    original.remove(center)

    while len(nodes) > 0:
        lighting = nodes[0]
        nodes.remove(lighting)
        is_exception = False
        try:
            canonic, is_appended = append_node(shape, lighting, canonic, debug, record)
            if is_appended is False:
                nodes = resort_nodes(center, original, debug)
                shape.reset()
                canonic = [center]
        except Exception as ex:
            print_node(debug, lighting, f"exception {ex}")
            nodes = resort_nodes(center, original, debug)
            shape.reset()
            canonic = [center]
            is_exception = True
        if len(canonic) < 3:
            recording_graph(record, canonic)
        if is_exception:
            print_nodes(debug, nodes, "nodes after exception")

    return shape, canonic


def transform_to_canonics(nodes, debug=False, record=None):
    """Split the original graph into connected subgraphs and transform them to a cononical type."""
    subgraphs = get_subgraphs(nodes)
    canonics = []

    for subgraph in subgraphs:
        shape, canonic = transform_to_canonic(subgraph, debug, record)
        canonics.append({"shape": shape, "canonic": canonic})
  
    return canonics


def merge_canonics(canonics):
    """Merge cononic graphs."""
    nodes = []
    for canonic in canonics:
        for node in canonic["canonic"]:
            nodes.append(node)
    return nodes