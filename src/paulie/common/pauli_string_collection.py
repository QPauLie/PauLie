"""
Class for a set/collection of Pauli strings with various features
"""
from random import randint
from itertools import combinations
from typing import Self, Generator
import numpy as np
import networkx as nx
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_linear import PauliStringLinear
from paulie.common.get_graph import get_graph
from paulie.classifier.classification import Classification
from paulie.classifier.morph_factory import MorphFactory
from paulie.classifier.recording_morph_factory import RecordingMorphFactory
from paulie.helpers.recording import RecordGraph



class PauliStringCollectionException(Exception):
    """
    Exception for the Pauli string collection class.
    """


# class PauliStringCollection(object):
#    """
#     Dummy PauliStringCollection
#    """


class PauliStringCollection:
    """
    Class for a collection of Pauli strings with various features.
    """

    def __init__(self, generators: list[PauliString] | Self = None) -> None:
        """
        Initializing a collection of Pauli strings.
        Args:
            generators: List of Pauli strings of type PauliString.
        """
        self.nextpos: int = 0
        self.generators: list[PauliString] = []
        self.classification: Classification = None
        self.record: RecordGraph = None
        if not generators:
            return

        longest: int = len(max(generators, key=len))

        for g in generators:
            if len(g) < longest:
                g = g.expand(longest)
            self.generators.append(g)

    def get(self) -> list[PauliString]:
        """
        Get an array of Pauli strings of type PauliString corresponding to the generator elements.
        Returns:
            Array of Pauli strings of type PauliString corresponding to the generator elements.
        """
        return self.generators

    def get_len(self) -> int:
        """
        Get the length of string in collection.
        Returns:
            Length of string in collection.
        """
        if len(self) == 0:
            return 0
        return len(self.generators[0])

    def set_record(self, record:RecordGraph) -> None:
        """
        In order to animate the transformations that lead to the canonical graph and thus to the classification, set record of type RecordGraph.
        Args:
            record: Record of type RecordGraph, recording the construction of canonical graphics.
        Returns: None
        """
        self.record = record

    def get_record(self) -> RecordGraph:
        """
        Get the record of graph construction.
        Returns: Record of graph construction.
        """
        return self.record

    def __str__(self) -> str:
        """
        Convert PauliStringCollection to readable string (e.g., "XYZI, YYYS").
        Returns:
            Readable string (e.g., "XYZI, YYYS").
        """
        return ",".join([str(g) for g in self.generators])

    def __len__(self) -> int:
        """ 
        Get length of collection.
        Returns:
            Number of generators in the collection.
        """
        return len(self.generators)

    def __iter__(self) -> Self:
        """
        Iterator over the generators.
        Returns:
           Iterator over the generators.
        """
        self.nextpos = 0
        return self

    def __next__(self) -> PauliString:
        """
        Next iterator value.
        Returns:
           next Iterator value.
        Raises:
           StopIteration:
               End of collection.
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = self.generators[self.nextpos] #.copy()
        self.nextpos += 1
        return value

    def __delitem__(self, key) -> PauliString:
        """
        Overloading the remove operator for a collection.
        Args: 
            key: Index key in collection.
        Returns:
            Result of the removed operator.
        """
        self.classification = None
        del self.generators[key]

    def __copy__(self) -> Self:
        """
        Overloading the collection copy operator.
        Returns:
            Copy of self.
        """
        return PauliStringCollection(self.generators)

    def copy(self) -> Self:
        """
        Copy collection.
        Returns:
            Copy of self.
        """
        return PauliStringCollection(self.generators)

    def __add__(self, p: PauliString) -> Self:
        """
        Overloading the addition operator with a collection.
        Args: 
            p:
            Pauli string to be added.
        Returns:
            PauliStringCollection with appended PauliString.

        """
        self.classification = None
        new_generators = []
        for g in self.generators:
            new_generators.append(g + p)
        return PauliStringCollection(new_generators)

    def mul(self, a:Self, b:Self) -> Self:
        """
        Multiplication on collection.
        Args:
            a: First Pauli string collection to multiplication.
            b: Second Pauli string collection to multiplication.
        Returns:
            Result of multiplication a * b.
        """
        self.classification = None
        new_generators = []
        for ga in a.generators:
            for gb in b.generators:
                new_generators.append(ga + gb)
        return PauliStringCollection(new_generators)

    def __mul__(self, other:Self) -> Self:
        """
        Overloading the multiplication operator with a collection.
        Args: 
            other: Collection for multiplication on self.
        Returns:
            Result of multiplication self * other.
        """
        return self.mul(self, other)

    def __rmul__(self, other:Self) -> Self:
        """
        Overloading the right multiplication operator with a collection.
        Args: 
            other: Collection for right multiplication on self.
        Returns:
            Result of multiplication other * self.
        """
        return self.mul(other, self)

    def expand(self, n: int) -> None:
        """
        Expands each string in the collection to specified length n by taking the tensor product with identities.
        Args:
            n: New length of Pauli string in collection.
        Returns:
            None
        """
        self.classification = None
        new_generators = []
        for g in enumerate(self.generators):
            g = g.expand(n)
        self.generators = new_generators

    def _processing(self, p: PauliString) -> PauliString:
        """
        Enforcing that each string in the collection is of the same size. Each string will be expanded with identities to have the length of the longest Pauli string.
        Args: 
            p: Pauli string arbitrary length.
        Returns:
            Pauli string with the same length as the collection element
        """
        if len(self.generators) == 0:
            return p
        longest = len(max(self.generators, key=len))
        if len(p) < longest:
            p = p.expand(longest)
        elif len(p) > longest:
            self.expand(len(p))
        return p

    def append(self, p: PauliString) -> None:
        """
        Append a specified Pauli string to the collection to the end.
        Args: 
            p: Pauli string arbitrary for appending to the end.
        Returns:
            None
        """
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.append(p)

    def insert(self, i: int, p: PauliString) -> None:
        """
        Insert a specified Pauli string to the collection at a specified position.
        Args:
            i: Position to insert
            p: Pauli string arbitrary length for inserting in position.
        Returns:
            None
        """
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.insert(i, p)

    def remove(self, p: PauliString) -> None:
        """
        Remove a specified Pauli string from the collection.
        Args:
            p: Pauli string for removing from collection.
        Returns:
            None
        """
        self.classification = None
        if p in self.generators:
            self.generators.remove(p)

    def index(self, p: PauliString) -> int:
        """
        Get the index of a given Pauli string inside the collection.
        Args:
            p: Pauli string for searching in collection.
        Returns:
            Index of a given Pauli string inside the collection.
        """
        return self.generators.index(p)

    def get_size(self) -> int:
        """
        Get the length of string in collection.
        Returns:
            Length of string in collection.
        """
        return 0 if len(self.generators) == 0 else len(self.generators[0])

    def create_instance(self, n: int = None, pauli_str: str = None) -> PauliString:
        """
        Create a new instance of the same type as the rest of the Pauli strings in the collection.
        Args:
            n: Length of new Pauli string.
            pauli_str: String representation of the Pauli string for new instance, if None then identity.
        Returns:
               New created PauliString object.
        Raises:
              PauliStringCollectionException: Empty collection
        """
        if len(self.generators) == 0:
            raise PauliStringCollectionException("Empty generator")
        return self.generators[0].create_instance(n=n, pauli_str=pauli_str)

    def sort(self) -> Self:
        """
        Sort the collection Pauli strings according to their bit value given by the bitarray representation
        Returns:
            Self
        """
        self.generators.sort()
        return self

    def get_anticommutation_fraction(self) -> float:
        """
        Computes the fraction of anticommuting pairs of generators.
        Returns:
            Fraction of anticommuting pairs of generators.
        """
        anti_commute_count = 0
        pair = 0
        for x,y in combinations(self.generators, r=2):
            pair += 1
            if not x|y:
                anti_commute_count += 1
        return anti_commute_count / pair

    def get_pair(self) -> int:
        """
        Get the maximum possible number of pairs of vertices in the graph.
        Returns:
            Maximum possible number of pairs of vertices in the graph.
        """
        return len(self)*(len(self) - 1)/2

    def get_anticommutation_pair(self) -> int:
        """
        Get the number of anticommutation pair.
        Returns:
            Number of anticommutation pair.
        """
        anti_commute_count = 0
        for x,y in combinations(self.generators, r=2):
            if not x|y:
                anti_commute_count += 1
        #n = len(self.generators)
        #n_com = n*(n-1)/2
        return anti_commute_count

    def get_commutants(self) -> Self:
        """
        Get the set of Pauli strings that commute with ALL generators in this collection. This finds the linear symmetries L_j of the system.
        Returns:
            PauliStringCollection of the set of Pauli strings that commute with ALL generators in this collection.

        """
        if not self.generators:
            # If there are no generators, all Paulis are symmetries by definition.
            # However, a group with no generators is trivial, so we can return empty.
            return PauliStringCollection([])

        num_qubits = self.get_size()

        # Start with a list of all possible Pauli strings.
        # Note: gen_all_pauli_strings is a method on the PauliString object.
        identity = self.create_instance(n=num_qubits)
        candidate_symmetries = list(identity.gen_all_pauli_strings())

        # For each generator in our system, filter the candidate list.
        for g in self.generators:
            # Keep only the candidates that commute with the current generator g.
            candidate_symmetries = [
                p for p in candidate_symmetries if g.commutes_with(p)
            ]

        # Return the final filtered list as a new collection.
        return PauliStringCollection(candidate_symmetries)

    def get_anti_commutants(self, generators: list[PauliString] | Self = None) -> Self:
        """
        Get Pauli strings that do not commute with the entire collection.
        Args:
            generators: Generators, Pauli string search list. If empty, then all lines are the same length
        Returns:
            Anticommutant of the set of Pauli strings.
        """
        if len(self) == 0:
            return PauliStringCollection([])
        for p in self:
            generators = PauliStringCollection(
                p.get_anti_commutants(generators=generators)
            )
        return generators

    def get_graph(
        self, generators: list[PauliString] | Self = None
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the anticommutation graph whose vertices are the generators and edges are determined by the commutator between the vertices.
        Args:
            generators: Area of Pauli strings over which to build a graph.
        Returns:
            Vertices, edges, and labels of edges of the anticommutation graph.
        """
        return get_graph(self.generators, commutators=generators)

    def get_commutator_graph(self
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the commutator graph whose vertices are all Paulistrings of a given dimension and an edge between two vertices exist if there is a element in the generator to which the one vertex anticommutes with to the other vertex.

        Returns:
            Vertices, edges, and labels of edges of graph.

        """
        n = self.get_size()
        i = PauliString(n=n)
        return get_graph(i.get_commutants(), commutators=self, flag_labels=False)

    def get_frame_potential(self) -> int:
        """
        Get the frame potential of the system. The frame potential is a measure of quantum chaos.
        Returns:
            Frame potential of the system generated by the collection.
        """
        vertices, edges = self.get_commutator_graph()
        graph = nx.Graph()
        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)
        n_comp= nx.number_connected_components(graph)
        n_iso =len(list(nx.isolates(graph)))
        return n_comp*n_iso

    def _convert(self, generators: list[str]) -> Self:
        """
        Convert a list of type string to List of type PauliString
        Args: 
            generators: a list of string representation of PauliString
        Returns: 
            PauliStringCollection of the list of generators
        """
        return PauliStringCollection(
            [self.create_instance(pauli_str=g) for g in generators]
        )

    def get_subgraphs(self) -> list[Self]:
        """
        Get the subgraphs of the anticommutation graph induced by the connected components.
        Returns: 
            List of connected subgraphs.
        """
        vertices, edges, _ = self.get_graph()

        g = nx.Graph()
        g.add_nodes_from(vertices)
        g.add_edges_from(edges)
        return [self._convert(subgaph) for subgaph in
               sorted(nx.connected_components(g), key=len, reverse=True)]

    def classify(self) -> Classification:
        """
        Build canonical graph.
        Returns:
            Result of building canonical graphs corresponding to the direct sum Lie algebra.
        """
        subgraphs = self.get_subgraphs()
        self.classification = Classification()
        for subgraph in subgraphs:
            if not self.record:
                morph_factory = MorphFactory()
            else:
                morph_factory = RecordingMorphFactory(record = self.record)

            self.classification.add(morph_factory.build(subgraph.get()).get_morph())
        return self.classification

    def get_class(self) -> Classification:
        """
        Get the set of canonical graphs used for the classification.
        Returns:
            Set of canonical graphs used for the classification.
        """
        if self.classification is None:
            self.classification = self.classify()
        return self.classification

    def get_algebra(self) -> str:
        """
        Get the dynamical Lie algebra generated by the set of Pauli strings.
        Returns:
            Dynamical Lie algebra generated by the set of Pauli strings.
        """
        classification = self.get_class()
        return classification.get_algebra()

    def is_algebra(self, algebra:str) -> bool:
        """
        Checks whether the classified algebra complies with a given algebra.
        Args: 
            algebra: Name of the algebra.
        Returns:
            True if the collection forms the specified algebra.
        """
        classification = self.get_class()
        return classification.is_algebra(algebra)

    def is_in(self, generators: Self) -> bool:
        """
        Testing generators in algebra. All Pauli strings of one algebra are dependent on another.
        Args: 
            generators:
            Generator collection for checking.
        Returns:
            True when all generators can be obtained from DLA self.
        """
        if len(self) == 0:
            return False

        classification = self.get_class()
        subgraphs = generators.get_subgraphs()
        for subgraph in subgraphs:
            is_eq = False
            morphs = classification.get_morphs()
            for morph in morphs:
                legs = morph.get_legs()
                morph_factory = MorphFactory()
                is_eq = morph_factory.is_eq(legs, subgraph.get())
                if is_eq:
                    break
            if not is_eq:
                return False

        return True

    def is_eq(self, generators: Self) -> bool:
        """
        Testing for equivalence of two algebras. All Pauli strings of one algebra are dependent on another.
        Args: 
            generators: Generator collection for checking.
        Returns:
            True when two collections of generators are equivalent.
        """
        return self.is_in(generators) and generators.is_in(self)

    def select_dependents(self, generators: Self) -> Self:
        """
        Select dependents from self.
        Args: 
            generators: Generator collection for selecting.
        Returns:
           Collection of dependent Pauli strings.
        """
        if len(self) == 0:
            return False
        dependents = []
        classification = self.get_class()
        subgraphs = generators.get_subgraphs()
        for subgraph in subgraphs:
            morphs = classification.get_morphs()
            for morph in morphs:
                legs = morph.get_legs()
                morph_factory = MorphFactory()
                morph_dependents = morph_factory.select_dependents(legs, subgraph.get())
                #print(f"morph = {PauliStringCollection(morph_dependents)}")
                dependents += morph_dependents

        return PauliStringCollection(dependents)

    def get_space(self) -> Self:
        """
        Get all space.
        Returns:
          Collection of all Pauli strings formed by a collection of generators.
        """
        all_space = PauliStringCollection([g
                    for g in PauliString(n=self.get_len()).gen_all_pauli_strings()
                    if g != PauliString(n=self.get_len())])
        return self.select_dependents(all_space)

    def gen_generators(self) -> Generator[list[list[PauliString]], None, None]:
        """
        Get generators.
        Returns:
            Generator of collections of generators of the same algebra as self.
        """
        classification = self.get_class()
        gen = classification.gen_generators()
        while True:
            try:
                g = next(gen)
                cg = PauliStringCollection(g)
                if cg.get_algebra() == self.get_algebra():
                    yield cg
            except StopIteration:
                break

    def get_dla_dim(self) -> int:
        """
        Get the dimension of the classified dynamical Lie algebra.
        Returns:
            Dimension of the classified dynamical Lie algebra.
        """
        return self.get_class().get_dla_dim()

    def get_dependents(self) -> Self:
        """
        Get a list of dependent strings in the collection
        Returns:
            List of dependent strings in the collection.
        """
        return PauliStringCollection(self.get_class().get_dependents())

    def get_independents(self) -> Self:
        """
        Get independents
        Returns:
            List of independent strings in the collection.
        """
        dependents = self.get_class().get_dependents()
        return PauliStringCollection([v for v in self if v not in dependents])

    def get_canonic_vertices(self) -> Self:
        """
        Get a list of canonic strings in the collection.
        Returns:
            List of canonic strings in the collection.
        """
        return PauliStringCollection(self.get_class().get_vertices())

    def get_canonic_graph(
        self,
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the canonical representation of a graph.
        Returns:
            Vertices, edges, and labels of edges.
        """
        classification = self.get_class()
        generators = PauliStringCollection(classification.get_vertices())
        return generators.get_graph()

    def get_anti_commutates(
        self, pauli_string: PauliString, generators: list[PauliString] | Self = None
    ) -> Self:
        """
        Get a collection of non-commuting Pauli strings.
        Args:
           pauli_string: Pauli string to which commutators are defined.
           generators: Area of Pauli strings over which to build a graph. If not specified, then collection.
        Returns:
            Collection of non-commuting Pauli strings.
        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection([g for g in generators
               if g != pauli_string and not pauli_string|g])

    def get_commutates(self, pauli_string:PauliString,
                       generators:list[PauliString]|Self = None) -> Self:
        """
        Get a collection of non-commuting Pauli strings.
        Args:
           pauli_string: Pauli string to which commutators are defined.
           generators: Area of Pauli strings over which to build a graph. If not specified, then collection.
        Returns:
            Collection of non-commuting Pauli strings.
        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection(
            [g for g in generators if g != pauli_string and g | pauli_string]
        )

    def get_graph_components(self,
                             graph_type: str = 'anticommutator'
                             ) -> list['PauliStringCollection']:
        """
        Computes the connected components of the specified graph (anticommutator or commutator) constructed from the Pauli strings in the collection.
        Args:
            graph_type (str): Type of graph to use for finding components. Must be either 'anticommutator' (default) or 'commutator'.
            list[PauliStringCollection]: List of PauliStringCollection objects, each representing a connected component of the selected graph.
        Raises:
            ValueError: If graph_type is not 'anticommutator' or 'commutator'.
        Returns: Vertices, edges, and labels of edges.
        """
        if graph_type == 'anticommutator':
            nodes, edges = self.get_graph() # The anti-commutation graph
        elif graph_type == 'commutator':
            nodes, edges = self.get_commutator_graph() # The commutator graph
        else:
            raise ValueError("graph_type must be 'anticommutator' or 'commutator'")

        return self._get_connected_components(nodes, edges)

    def _get_connected_components(self, nodes: list[str],
                                  edges: list[tuple[str, str]]
                                  ) -> list[Self]:
        """
        Helper method to compute connected components from nodes and edges.
        Args:
            nodes (list[str]): List of node strings.
            edges (list[tuple[str, str]]): List of edges as tuples of node strings.
        Returns:
            list[PauliStringCollection]: List of PauliStringCollection objects, each representing a connected component.
        """
        # Create a graph from the nodes and edges
        graph = nx.Graph()
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

        # Find connected components
        connected_components = list(nx.connected_components(graph))

        # Convert sets of node strings back into PauliStringCollection objects
        return [self._convert(subgraph) for subgraph in connected_components]

    def get_symmetries_for_component(
        self, linear_symmetries: Self
    ) -> list[PauliStringLinear]:
        """
        Private helper: Computes the quadratic symmetries for THIS collection, assuming THIS collection is a single connected component (C_k) and is provided with a pre-computed list of linear symmetries (L_j). This method performs the inner loop of the full calculation.
        Args: 
            linear_symmetries: Collection of linear symmetries.
        Returns:
            List of linear combinations.
        """
        # Convert this component collection into a single PauliStringLinear object
        component_as_linear = PauliStringLinear([(1.0, str(s)) for s in self])

        quadratic_symmetries_for_this_component = []

        # For each linear symmetry Lj...
        for lj_pauli in linear_symmetries:
            # ...compute Q_kj using the low-level .quadratic() method.
            q_kj = component_as_linear.quadratic(lj_pauli)
            quadratic_symmetries_for_this_component.append(q_kj)

        return quadratic_symmetries_for_this_component


    def get_full_quadratic_basis(self, normalized: bool = False) -> list[PauliStringLinear]:
        """
        Public Method: Calculates the full basis of quadratic symmetries for this system.

        This method performs the complete, high-level operation:
        1. Finds the linear symmetries (L_j) for this system.
        2. Finds the connected components (C_k) of this system's commutator graph.
        3. Loops through each component and calls the internal helper to get the Q_kj.
        
        Args:
            normalized (bool): If True, returns an ORTHONORMAL basis. If False (default), returns an ORTHOGONAL basis.
        Returns:
            List of PauliStringLinear objects representing the basis.
        """
        # Step 1: Find the linear symmetries for this system
        linear_symmetries = self.get_commutants()

        # Step 2: Find the connected components of the commutator graph
        connected_components = self.get_graph_components(graph_type='commutator')

        # Step 3: Loop through components and get symmetries for each
        full_basis = []
        for component in connected_components:
            # Call the private helper for each component
            symmetries = component.get_symmetries_for_component(linear_symmetries)
            full_basis.extend(symmetries)

        # Filter out any zero vectors that may have been created
        orthogonal_basis = [q for q in full_basis if not q.is_zero()]

        if not normalized:
            return orthogonal_basis

        normalized_basis = []
        for q_vector in orthogonal_basis:
            # The squared norm is Tr(Q†Q). Use the .H property.
            squared_norm_trace = (q_vector.h @ q_vector).trace()

            # The trace should be real and positive for a non-zero operator
            if squared_norm_trace.real > 1e-12:
                # The Hilbert-Schmidt norm is the sqrt of this trace
                norm = np.sqrt(squared_norm_trace.real)
                normalized_basis.append(q_vector * (1.0 / norm))

        return normalized_basis

    def find(self, pauli_string: PauliString) -> int:
        """
        Find index pauli string in collection.
        Args:
            pauli_string: Pauli string for searching.
        Returns:
            Index pauli string in collection.
        """
        for i, p in enumerate(self.generators):
            if p == pauli_string:
                return i
        return -1

    def replace(self, pauli_string: PauliString, new_pauli_string: PauliString) -> None:
        """
        Replacing one Pauli string with another.
        Args:
            pauli_string: Pauli string that will be replaced.
            new_pauli_string: New Pauli string on place the old Pauli string.
        Returns:
            None
        """
        index = self.find(pauli_string)
        if index != -1:
            self.generators[index] = new_pauli_string.copy()

    def contract(self, pauli_string: PauliString, contracted_pauli_string: PauliString) -> None:
        """
        Contracting pauli string with other pauli string, pauli_string will be replaced on pauli_string @ contracted_pauli_string.
        Args:
            pauli_string: Pauli string that will be replaced with the contract.
            contracted_pauli_string: Contractable Pauli string.
        Returns:
            None
        """
        self.replace(pauli_string, pauli_string@contracted_pauli_string)

    def list_connections(self):
        """ 
        Get list of connections.
        Returns:
            List of pairs of connected vertices in an anti-commutative graph.
        """
        return [(x, y, len(x.get_anti_commutants(self)), len(y.get_anti_commutants(self)))
               for x,y in combinations(self.generators, r=2) if not x|y]

    def _get_delta(self, generators: Self, number_connections:int)->int:
        """
        Get the difference between the required number and the number of connections in the graph.
        Args: 
            generators: Generator collection of graph.
            number_connections: Required number of connections.
        Returns:
            Difference between the required number and the number of connections in the graph.
        """
        return number_connections - generators.get_anticommutation_pair()

    def find_generators_with_connection(self, number_connections:int)->Self:
        """
        Get generators with connections.
        Args:
            number_connections: Required number of connections.
        Returns: 
            Collection of generators with the required number of connections.
        """
        generators = self.copy().get_canonic_vertices()
        max_iter = number_connections//2
        i = 0
        while i < max_iter:
            delta = self._get_delta(generators, number_connections)
            if delta == 0:
                break
            delta_min = delta
            current_generators = generators.copy()
            for connections in generators.list_connections():
                gx = generators.copy()
                gy = generators.copy()
                x = connections[0]
                y = connections[1]

                gx.contract(x,y)
                gy.contract(y,x)
                delta_x = self._get_delta(gx, number_connections)
                delta_y = self._get_delta(gy, number_connections)
                if abs(delta_x) < abs(delta_y):
                    if delta_x < 0:
                        continue

                    if abs(delta_min) > abs(delta_x):
                        current_generators = gx
                        delta_min = abs(delta_x)
                else:
                    if delta_y < 0:
                        continue

                    if abs(delta_min) > abs(delta_y):
                        current_generators = gy
                        delta_min = abs(delta_y)
            if delta == delta_min:
                list_connections = generators.list_connections()
                while True:
                    index = randint(0, i)
                    connections = list_connections[index]
                    gx = generators.copy()
                    gy = generators.copy()
                    x = connections[0]
                    y = connections[1]

                    gx.contract(x,y)
                    gy.contract(y, x)
                    delta_x = self._get_delta(gx, number_connections)
                    delta_y = self._get_delta(gy, number_connections)
                    if delta_x > 0:
                        current_generators = gx
                        break
                    if delta_x > 0:
                        current_generators = gy
                        break

            generators = current_generators

            i += 1
        return generators
