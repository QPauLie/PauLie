"""
Collection of Pauli strings
"""
import re
import networkx as nx
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.get_graph import get_graph
from paulie.classifier.classification import Classification
from paulie.classifier.morph_factory import MorphFactory
from paulie.classifier.recording_morph_factory import RecordingMorphFactory

class PauliStringCollectionException(Exception):
    """
    Pauli string exception
    """

class PauliStringCollection:
    """
    Collection of Pauli strings
    """
    def __init__(self, generators:list[PauliString]=None, debug: bool = False):
        """
        Initializing a collection
        Args:
            generators: List of Pauli strings
        """
        self.nextpos = 0
        self.generators = []
        self.classification = None
        self.debug = debug
        self.record = None
        if not generators:
            return
        self.sv_queue = None
        self.debug_queue = None

        longest = len(max(generators, key=len))

        for g in generators:
            if len(g) < longest:
                g = g.expand(longest)
            self.generators.append(g)

    def get(self)->list[PauliString]:
        """Get an array of Pauli strings"""
        return self.generators

    def set_debug(self, debug):
        """ Set debug flag"""
        self.debug = debug

    def get_debug(self):
        """Get debug flag"""
        return self.debug

    def set_record(self, record):
        """ Set record 
            Args:
                 record - record of the classification algorithm
        """
        self.record = record

    def get_record(self):
        """
        Get record
        """
        return self.record

    def __str__(self) -> str:
        """Convert PauliStringCollection to readable string (e.g., "XYZI, YYYS")"""
        return ",".join([str(g) for g in self.generators])

    def __len__(self) -> int:
        """ Collection length"""
        return len(self.generators)

    def __iter__(self):
        """Collection Iterator"""
        self.nextpos = 0
        return self

    def __next__(self):
        """Next iterator value"""
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = self.generators[self.nextpos] #.copy()
        self.nextpos += 1
        return value

    def __delitem__(self, key):
        """Overloading the remove operator from a collection"""
        self.classification = None
        del self.generators[key]

    def __copy__(self):
        """Overloading the collection copy operator"""
        return PauliStringCollection(self.generators, debug=self.debug)

    def copy(self):
        """ Copy collection"""
        return PauliStringCollection(self.generators, debug=self.debug)

    def __add__(self, p: PauliString):
        """Overloading the addition operator with a collection"""
        self.classification = None
        new_generators = []
        for g in self.generators:
            new_generators.append(g + p)
        return PauliStringCollection(new_generators, debug=self.debug)

    def mul(self, a, b):
        """ multiplication on collection"""
        self.classification = None
        new_generators = []
        for ga in a.generators:
            for gb in b.generators:
                new_generators.append(ga + gb)
        return PauliStringCollection(new_generators, debug=self.debug)

    def __mul__(self, other):
        """Overloading the multiplication operator with a collection"""
        return self.mul(self, other)

    def __rmul__(self, other):
        """Overloading the right multiplication operator with a collection"""
        return self.mul(other, self)

    def expand(self, n: int):
        """Expand collection string to specified size"""
        self.classification = None
        new_generators = []
        for g in enumerate(self.generators):
            g = g.expand(n)
        self.generators = new_generators

    def _processing(self, p: PauliString):
        """Preprocessing adding a row to a collection. Aligns the size of Pauli strings"""
        if len(self.generators) == 0:
            return p
        longest = len(max(self.generators, key=len))
        if len(p) < longest:
            p = p.expand(longest)
        elif len(p) > longest:
            self.expand(len(p))
        return p

    def append(self, p: PauliString):
        """Append Pauli string"""
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.append(p)

    def insert(self, i: int, p: PauliString):
        """Insert Pauli string"""
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.insert(i, p)

    def remove(self, p: PauliString):
        """Remove Pauli string"""
        self.classification = None
        if p in self.generators:
            self.generators.remove(p)

    def index(self, p: PauliString):
        """Index Pauli string"""
        return self.generators.index(p)

    def get_size(self):
        """Get length of Pauli String in collection"""
        return 0 if len(self.generators) == 0 else len(self.generators[0])

    def create_instance(self, n: int = None, pauli_str: str = None):
        """
        Create a new instance of the same type as the rest of the Pauli strings in the collection
        """
        if len(self.generators) == 0:
            raise PauliStringCollectionException("Empty generator")
        return self.generators[0].create_instance(n=n, pauli_str=pauli_str)

    def sort(self):
        """Sort collection"""
        self.generators.sort()

    def get_commutants(self, generators=None):
        """
           Get Pauli strings that commute with the entire collection
           Args:
               generators: Generators, Pauli string search list.
               If empty, then all lines are the same length
        """
        if len(self) == 0:
            return PauliStringCollection([], debug=self.debug)
        for p in self:
            generators = PauliStringCollection(p.get_commutants(generators=generators),
                         debug=self.debug)
        return generators

    def get_anti_commutants(self, generators=None):
        """
           Get Pauli strings that no-commute with the entire collection
           Args:
               generators: Generators, Pauli string search list.
               If empty, then all lines are the same length
        """
        if len(self) == 0:
            return PauliStringCollection([], debug=self.debug)
        for p in self:
            generators = PauliStringCollection(p.get_anti_commutants(generators=generators),
                         debug=self.debug)
        return generators


    def get_graph(self, generators: list["PauliStringCollection"]=None):
        """
        Get graph
        Args:
            generators: The area of Pauli strings over which to build a graph.
            If not specified, then that's it
        Returns the vertices, edges, and labels of edges
        """
        return get_graph(self.generators, commutators=generators)

    def get_commutator_graph(self):
        """
        Get commutator graph
        """
        n = self.get_size()
        i = PauliString(n=n)
        return get_graph(i.get_commutants(), commutators=self, flag_labels=False)

    def frame_potential(self):
        """
        Get frame potential
        """
        vertices, edges = self.get_commutator_graph()
        graph = nx.Graph()
        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)
        n_comp= nx.number_connected_components(graph)
        n_iso =len(list(nx.isolates(graph)))
        return n_comp*n_iso

    def _convert(self, generators: list[str])->"PauliStringCollection":
        """ Convert list of string to List of PauliString"""
        return PauliStringCollection([self.create_instance(pauli_str=g)
               for g in generators],
               debug=self.debug)

    def get_subgraphs(self):
        """Get related subsets of a collection"""
        vertices, edges, _ = self.get_graph()

        g = nx.Graph()
        g.add_nodes_from(vertices)
        g.add_edges_from(edges)
        return [self._convert(subgaph) for subgaph in
               sorted(nx.connected_components(g), key=len, reverse=True)]

    def debug_classify(self):
        """
        debug classify
        """
        self.classification = Classification()
        morph_factory = MorphFactory()
        self.classification.add(morph_factory.build(self.debug_queue).get_morph())
        return self.classification

    def classify(self):
        """Classify collection"""
        if self.debug_queue:
            return self.debug_classify()

        subgraphs = self.get_subgraphs()
        self.classification = Classification()
        for subgraph in subgraphs:
            if not self.record:
                morph_factory = MorphFactory()
            else:
                morph_factory = RecordingMorphFactory(record = self.record)

            self.classification.add(morph_factory.build(subgraph).get_morph())
            sv_queue = subgraph.get_sv_queue()
            if sv_queue:
                self.set_sv_queue(sv_queue)

        return self.classification

    def get_class(self):
        """Get collection class"""
        if self.classification is None:
            self.classification = self.classify()
        return self.classification

    def get_algebra(self):
        """
         Get algebra
        """
        classification = self.get_class()
        return classification.get_algebra()

    def is_algebra(self, algebra):
        """
        Checking for compliance with a given algebra
        """
        classification = self.get_class()
        return classification.is_algebra(algebra)

    def reset_algebra(self):
        """
        Reset classification algebra calculation
        """
        self.classification = None

    def dla_dim(self):
        """
        Get dim of DLA
        """
        dim_su = lambda n: n**2-1
        dim_so = lambda n: n*(n-1)/2
        dim_sp = lambda n: n*(2*n+1)
        subgraphs = self.get_class().get_subalgebras()
        dim = 0
        for s in subgraphs:
            n = [int(x) for x in re.split('(|)', s) if x.isdigit()][0]
            if "su" in s:
                dim+= dim_su(n)
            if "sp" in s:
                dim+= dim_sp(n)
            if "so" in s:
                dim+= dim_so(n)
        return dim

    def get_dependents(self):
        """Get a list of dependent strings in the collection"""
        return PauliStringCollection(self.get_class().get_dependents())

    def get_canonic_graph(self):
        """Get the canonical representation of a graph"""
        classification = self.get_class()
        generators = PauliStringCollection(classification.get_vertices(), debug=self.debug)
        return generators.get_graph()

    def get_anti_commutates(self, pauli_string, generators = None):
        """
            Get a collection of non-commuting Pauli strings
            Args:
                Pauli string to which commutators are defined
            generators: The area of Pauli strings over which to build a graph.
            If not specified, then collection

        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection([g for g in generators
               if g != pauli_string and not pauli_string|g],
               debug=self.debug)

    def get_commutates(self, pauli_string, generators = None):
        """
            Get a collection of non-commuting Pauli strings
            Args:
                Pauli string to which commutators are defined
            generators: The area of Pauli strings over which to build a graph.
            If not specified, then collection

        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection([g for g in generators
               if g != pauli_string and g|pauli_string],
               debug=self.debug)

    def get_max_connected(self):
        """Get the Pauli string that has the maximum number of non-commutable"""
        if len(self) == 0:
            return None, None
        pauli_string = self.generators[0]
        anti_commutates = self.get_anti_commutates(pauli_string)
        for p in self.generators:
            _anti_commutates = self.get_anti_commutates(p)
            if len(_anti_commutates) > len(anti_commutates):
                pauli_string = p
                anti_commutates = _anti_commutates
        return pauli_string, anti_commutates

    def _append_to_queue(self, queue_pauli_strings, pauli_strings):
        """Append the next related Pauli string to the queue"""
        for p in pauli_strings:
            if p in queue_pauli_strings:
                pauli_strings.remove(p)
                continue
            anti_commutates = self.get_anti_commutates(p, generators = queue_pauli_strings)
            if len(anti_commutates) == 0:
                continue
            if len(anti_commutates) > 1:
                min_index = len(queue_pauli_strings)
                for anti_commutate in anti_commutates:
                    index = queue_pauli_strings.index(anti_commutate)
                    if index < min_index:
                        min_index = index
                        queue_pauli_strings.insert(min_index + 1, p)
            else:
                queue_pauli_strings.append(p)
            pauli_strings.remove(p)
            return

    def get_queue(self):
        """Get associated sequence of Pauli strings"""
        if self.debug_queue:
            return self.debug_queue.get()

        new_generators = self.copy()
        new_generators.sort()
        queue_pauli_strings = PauliStringCollection(debug=self.debug)
        pauli_string, anti_commutates = self.get_max_connected()

        new_generators.remove(pauli_string)
        queue_pauli_strings.append(pauli_string)
        for anti_commutate in anti_commutates:
            new_generators.remove(anti_commutate)
            if anti_commutate not in queue_pauli_strings:
                queue_pauli_strings.append(anti_commutate)

        while len(new_generators) > 0:
            self._append_to_queue(queue_pauli_strings, new_generators)
        if len(queue_pauli_strings) > 1:
            self.sv_queue = PauliStringCollection(queue_pauli_strings)
        return queue_pauli_strings

    def get_sv_queue(self):
        """
        Get saved queue
        """
        return self.sv_queue

    def set_sv_queue(self, sv_queue):
        """
        Set saved queue
        """
        self.sv_queue = sv_queue

    def set_debug_queue(self, debug_queue):
        """
         Set debug queue
        """
        self.debug_queue = debug_queue
