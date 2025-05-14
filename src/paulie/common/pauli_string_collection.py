"""
Class for a set/collection of Pauli strings with various features
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
    Exception for the Pauli string collection class
    """

class PauliStringCollection:
    """
    Class for a collection of Pauli strings with various features
    """
    def __init__(self, generators:list[PauliString]=None):
        """
        Initializing a collection of Pauli strings
        Args:
            generators: List of Pauli strings of type PauliString
        """
        self.nextpos = 0
        self.generators = []
        self.classification = None
        self.record = None
        if not generators:
            return

        longest = len(max(generators, key=len))

        for g in generators:
            if len(g) < longest:
                g = g.expand(longest)
            self.generators.append(g)

    def get(self)->list[PauliString]:
        """ Get an array of Pauli strings of type PauliString
         corresponding to the generator elements """
        return self.generators

    def set_record(self, record):
        """ In order to animate the transformations that lead to the canonical graph
        and thus to the classifcation, set record of type RecordGraph
            Args:
                 record - record of type RecordGraph
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
        """ Returns the number of generators in the collection"""
        return len(self.generators)

    def __iter__(self):
        """ Iterator over the generators """
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
        """Overloading the remove operator for a collection"""
        self.classification = None
        del self.generators[key]

    def __copy__(self):
        """Overloading the collection copy operator"""
        return PauliStringCollection(self.generators)

    def copy(self):
        """ Copy collection"""
        return PauliStringCollection(self.generators)

    def __add__(self, p: PauliString):
        """Overloading the addition operator with a collection"""
        self.classification = None
        new_generators = []
        for g in self.generators:
            new_generators.append(g + p)
        return PauliStringCollection(new_generators)

    def mul(self, a, b):
        """ multiplication on collection"""
        self.classification = None
        new_generators = []
        for ga in a.generators:
            for gb in b.generators:
                new_generators.append(ga + gb)
        return PauliStringCollection(new_generators)

    def __mul__(self, other):
        """Overloading the multiplication operator with a collection"""
        return self.mul(self, other)

    def __rmul__(self, other):
        """Overloading the right multiplication operator with a collection"""
        return self.mul(other, self)

    def expand(self, n: int):
        """ Expands each string in the collection to specificed length n by taking
        the tensor product with identities"""
        self.classification = None
        new_generators = []
        for g in enumerate(self.generators):
            g = g.expand(n)
        self.generators = new_generators

    def _processing(self, p: PauliString):
        """ Enforcing that each string in the collection
        is of the same size.  Each string will be expanded
        with identities to have the length of the longest Pauli string."""
        if len(self.generators) == 0:
            return p
        longest = len(max(self.generators, key=len))
        if len(p) < longest:
            p = p.expand(longest)
        elif len(p) > longest:
            self.expand(len(p))
        return p

    def append(self, p: PauliString):
        """Append a specified Pauli string to the collection to the end """
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.append(p)

    def insert(self, i: int, p: PauliString):
        """Insert a specificed Pauli string to the collection at a specificed position"""
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.insert(i, p)

    def remove(self, p: PauliString):
        """Remove a specificed Pauli string from the collection"""
        self.classification = None
        if p in self.generators:
            self.generators.remove(p)

    def index(self, p: PauliString):
        """Returns the index of a given Pauli string inside the collection"""
        return self.generators.index(p)

    def get_size(self):
        """ Get the length of the Pauli Strings in the collection"""
        return 0 if len(self.generators) == 0 else len(self.generators[0])

    def create_instance(self, n: int = None, pauli_str: str = None):
        """
        Create a new instance of the same type as the rest of the Pauli strings in the collection
        """
        if len(self.generators) == 0:
            raise PauliStringCollectionException("Empty generator")
        return self.generators[0].create_instance(n=n, pauli_str=pauli_str)

    def sort(self):
        """ Sort the collection Pauli strings according to their bit value
         given by the bitarray representation """
        self.generators.sort()
        return self

    def get_commutants(self, generators=None):
        """
           Get Pauli strings that commute with the entire collection
           Args:
               generators: Generators, Pauli string search list.
               If empty, then all lines are the same length
           Returns:
                 commutant of the set of Pauli strings
        """
        if len(self) == 0:
            return PauliStringCollection([])
        for p in self:
            generators = PauliStringCollection(p.get_commutants(generators=generators))
        return generators

    def get_anti_commutants(self, generators=None):
        """
           Get Pauli strings that do not commute with the entire collection
           Args:
               generators: Generators, Pauli string search list.
               If empty, then all lines are the same length
           Returns:
               anticommutant of the set of Pauli strings
        """
        if len(self) == 0:
            return PauliStringCollection([])
        for p in self:
            generators = PauliStringCollection(p.get_anti_commutants(generators=generators))
        return generators


    def get_graph(self, generators: list["PauliStringCollection"]=None):
        """
        Get the anticommutation graph whose vertices are the generators and edges are
        determined by the commutator between the vertices
        Args:
            generators: The area of Pauli strings over which to build a graph.
        Returns:
              vertices, edges, and labels of edges of the anticommutation graph
        """
        return get_graph(self.generators, commutators=generators)

    def get_commutator_graph(self):
        """
        Get the commutator graph whose vertices are all Paulistrings of a given dimension
        and an edge between two vertices exist if there is a element in the generator
        to which the one vertex anticommutes with to the other vertex.
        """
        n = self.get_size()
        i = PauliString(n=n)
        return get_graph(i.get_commutants(), commutators=self, flag_labels=False)

    def get_frame_potential(self):
        """
        Returns the frame potential of the system
        generated by the collection.
        The frame potential is a measure of quantum choas.
        """
        vertices, edges = self.get_commutator_graph()
        graph = nx.Graph()
        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)
        n_comp= nx.number_connected_components(graph)
        n_iso =len(list(nx.isolates(graph)))
        return n_comp*n_iso

    def _convert(self, generators: list[str])->"PauliStringCollection":
        """ Convert a list of type string to List of type PauliString """
        return PauliStringCollection([self.create_instance(pauli_str=g)
               for g in generators])

    def get_subgraphs(self):
        """ Get the subgraphs of the anticommutation graph induced by the connected components"""
        vertices, edges, _ = self.get_graph()

        g = nx.Graph()
        g.add_nodes_from(vertices)
        g.add_edges_from(edges)
        return [self._convert(subgaph) for subgaph in
               sorted(nx.connected_components(g), key=len, reverse=True)]

    def classify(self):
        """ Returns a set of canonical graphs corresponding to the direct sum Lie algebra"""
        subgraphs = self.get_subgraphs()
        self.classification = Classification()
        for subgraph in subgraphs:
            if not self.record:
                morph_factory = MorphFactory()
            else:
                morph_factory = RecordingMorphFactory(record = self.record)

            self.classification.add(morph_factory.build(subgraph.get()).get_morph())
        return self.classification

    def get_class(self):
        """ Get the set of canonical graphs used for the classification """
        if self.classification is None:
            self.classification = self.classify()
        return self.classification

    def get_algebra(self):
        """
         Get the dynamical Lie algebra generated by the set of Pauli strings
        """
        classification = self.get_class()
        return classification.get_algebra()

    def is_algebra(self, algebra):
        """
        Checks whether the classifed algebra complies with a given algebra
        """
        classification = self.get_class()
        return classification.is_algebra(algebra)

    def get_dla_dim(self):
        """
        Get the dimension of the classified dynamical Lie algebra
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
        generators = PauliStringCollection(classification.get_vertices())
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
               if g != pauli_string and not pauli_string|g])

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
               if g != pauli_string and g|pauli_string])
