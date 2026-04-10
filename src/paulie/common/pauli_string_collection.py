"""
    Class for a set/collection of Pauli strings with various features
"""
from __future__ import annotations

from random import randint
from itertools import combinations
from collections.abc import Generator
from typing import Self
import numpy as np
import networkx as nx
from paulie.common.pauli_string_bitarray import PauliString
from paulie.common.pauli_string_linear import PauliStringLinear
from paulie.common.get_graph import get_graph
from paulie.classifier.classification import Classification
from paulie.classifier.morph_factory import MorphFactory
from paulie.classifier.recording_morph_factory import RecordingMorphFactory
from paulie.helpers._recording import RecordGraph
from paulie.exceptions import PauliStringCollectionException

class PauliStringCollection:
    """
    Class for a collection of Pauli strings with various features.
    """

    def __init__(self, generators: list[PauliString] | Self | None = None) -> None:
        """
        Initializing a collection of Pauli strings.

        Args:
            generators (list[PauliString] | Self, optional): List of Pauli strings of type
                PauliString.
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
        Get an array of Pauli strings corresponding to the generator elements.

        Returns:
            list[PauliString]: Array of Pauli strings corresponding to the generator elements.
        """
        return self.generators

    def get_len(self) -> int:
        """
        Get the length of Pauli strings in the collection.

        Returns:
            int: Length of Pauli strings in the collection.
        """
        if len(self) == 0:
            return 0
        return len(self.generators[0])

    def set_record(self, record:RecordGraph) -> None:
        """
        In order to animate the transformations that lead to the canonical graph
        and thus to the classification, set record of type RecordGraph.

        Args:
            record (RecordGraph): Record of type RecordGraph, recording the construction of
                the canonical graph.
        Returns:
            None
        """
        self.record = record

    def get_record(self) -> RecordGraph:
        """
        Get the record of graph construction.

        Returns:
            RecordGraph: Record of graph construction.
        """
        return self.record

    def __repr__(self) -> str:
        """
        Convert PauliStringCollection to readable string
        (e.g., PauliStringCollection([XYZI, YYYS])).

        Returns:
            str: String representation.
        """
        return f"PauliStringCollection([{', '.join(str(g) for g in self.generators)}])"

    def __str__(self) -> str:
        """
        Convert PauliStringCollection to readable string
        (e.g., [XYZI, YYYS]).

        Returns:
            str: String representation.
        """
        return f"[{', '.join(str(g) for g in self.generators)}]"

    def __len__(self) -> int:
        """
        Get the number of Pauli strings in the collection.

        Returns:
            int: Number of Pauli strings in the collection.
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
        Next Pauli string.

        Returns:
           PauliString: Next Pauli string.

        Raises:
           StopIteration: End of collection.
        """
        if self.nextpos >= len(self):
            # we are done
            raise StopIteration
        value = self.generators[self.nextpos]
        self.nextpos += 1
        return value

    def __delitem__(self, key:int) -> PauliString:
        """
        Delete a Pauli string from the collection.

        Args:
            key (int): Index key in collection.
        Returns:
            PauliString: The removed Pauli string.
        """
        self.classification = None
        del self.generators[key]

    def __copy__(self) -> PauliStringCollection:
        """
        Overloading the collection copy operator.

        Returns:
            PauliStringCollection: Copy of self.
        """
        return PauliStringCollection(self.generators)

    def copy(self) -> PauliStringCollection:
        """
        Copy collection.

        Returns:
            PauliString: Copy of self.
        """
        return PauliStringCollection(self.generators)

    def __add__(self, p: PauliString) -> PauliStringCollection:
        """
        Tensor product a Pauli string at the end of every element of the collection.

        Args:
            p (PauliString): Pauli string to be tensored.
        Returns:
            PauliStringCollection:
            PauliStringCollection with the Pauli string tensored at the end to every element.
        """
        self.classification = None
        new_generators = []
        for g in self.generators:
            new_generators.append(g + p)
        return PauliStringCollection(new_generators)

    def mul(self, a: Self, b: Self) -> PauliStringCollection:
        """
        Tensor product of two collections.

        Args:
            a: First Pauli string collection to multiply.
            b: Second Pauli string collection to multiply.
        Returns:
            Collection with elements as tensor products of all pairs of Pauli strings from the two
            collections.
        """
        self.classification = None
        new_generators = []
        for ga in a.generators:
            for gb in b.generators:
                new_generators.append(ga + gb)
        return PauliStringCollection(new_generators)

    def __mul__(self, other: Self) -> PauliStringCollection:
        """
        Overloading the multiplication operator with a collection.

        Args:
            other: Collection for multiplication on self.
        Returns:
            PauliStringCollection: Result of multiplication self * other.
        """
        return self.mul(self, other)

    def __rmul__(self, other: Self) -> PauliStringCollection:
        """
        Overloading the right multiplication operator with a collection.

        Args:
            other: Collection for right multiplication on self.
        Returns:
            PauliStringCollection: Result of multiplication other * self.
        """
        return self.mul(other, self)

    def expand(self, n: int) -> None:
        """
        Expands each string in the collection to specified length n
        by taking the tensor product at the end with identities.

        Args:
            n (int): New length of Pauli string in collection.
        Returns:
            None
        """
        self.classification = None
        self.generators = [g.expand(n) for g in self.generators]

    def _processing(self, p: PauliString) -> PauliString:
        """
        Enforcing that each string in the collection is of the same size.
        Each string will be expanded with identities to have the length of the longest Pauli string.

        Args:
            p (PauliString): A Pauli string of arbitrary length.
        Returns:
            PauliString:
            A Pauli string with the same length as the longest element of the collection.
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
        Add a Pauli string to the collection at the end.

        Args:
            p (PauliString): Pauli string to be added to the collection.
        Returns:
            None
        """
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.append(p)

    def insert(self, i: int, p: PauliString) -> None:
        """
        Insert a Pauli string to the collection at a specified position.

        Args:
            i (int): Position to insert the Pauli string at.
            p (PauliString): Pauli string arbitrary length for inserting in position.
        Returns:
            None
        """
        self.classification = None
        p = self._processing(p)
        if p not in self.generators:
            self.generators.insert(i, p)

    def remove(self, p: PauliString) -> None:
        """
        Remove a Pauli string from the collection.

        Args:
            p (PauliString): Pauli string to be removed from the collection.
        Returns:
            None
        """
        self.classification = None
        if p in self.generators:
            self.generators.remove(p)

    def index(self, p: PauliString) -> int:
        """
        Get the index of a Pauli string inside the collection.

        Args:
            p (PauliString): Pauli string to search for in collection.
        Returns:
            Index of the Pauli string inside the collection.
        """
        return self.generators.index(p)

    def get_size(self) -> int:
        """
        Get the length of the Pauli strings in the collection.

        Returns:
            Length of the Pauli strings in the collection.
        """
        return 0 if len(self.generators) == 0 else len(self.generators[0])

    def create_instance(self, n: int = None, pauli_str: str = None) -> PauliString:
        """
        Create a new instance of the same type as the rest of the Pauli strings in the collection.

        Args:
            n (int, optional): Length of the new Pauli string.
            pauli_str (str, optional): String representation of the new Pauli string instance.
                Defaults to None, in which case the identity is returned.
        Returns:
            PauliString: The created Pauli string instance.

        Raises:
              PauliStringCollectionException: If the collection is empty.
        """
        if len(self.generators) == 0:
            raise PauliStringCollectionException("Empty generator")
        return self.generators[0].create_instance(n=n, pauli_str=pauli_str)

    def sort(self) -> Self:
        """
        Sort the Pauli strings in the collection according to their bit value given by the bitarray
        representation.

        Returns:
            Self
        """
        self.generators.sort()
        return self

    def get_anticommutation_fraction(self) -> float:
        """
        Computes the fraction of anticommuting pairs of generators.

        Returns:
            float: Fraction of anticommuting pairs of generators.
        """
        anti_commute_count = 0
        pair = 0
        for x, y in combinations(self.generators, r=2):
            pair += 1
            if not x | y:
                anti_commute_count += 1
        return anti_commute_count / pair

    def get_pair(self) -> int:
        """
        Get the maximum possible number of pairs of vertices in the graph.

        Returns
            int: Maximum possible number of pairs of vertices in the graph.
        """
        return len(self) * (len(self) - 1) // 2

    def get_anticommutation_pair(self) -> int:
        """
        Get the number of anticommuting pairs in the collection.

        Returns:
            int: Number of anticommuting pairs in the collection.
        """
        anti_commute_count = 0
        for x, y in combinations(self.generators, r=2):
            if not x | y:
                anti_commute_count += 1
        return anti_commute_count

    def get_commutants(self) -> PauliStringCollection:
        """
        Get the set of Pauli strings that commute with all generators in the collection.
        This finds the linear symmetries :math:`L_{j}` of the system.

        Returns:
            PauliStringCollection:
            Collection of Pauli strings that commute with all generators in the collection.
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

    def get_anti_commutants(self, generators: PauliStringCollection) -> PauliStringCollection:
        """
        Get the set of Pauli strings that anticommute with the entire collection.

        Args:
            generators (PauliStringCollection): The set of Pauli strings to search over.
        Returns:
            PauliStringCollection:
            Subset of the generators that belong to the anti-commutant of the entire collection.
        """
        if len(self) == 0:
            return PauliStringCollection([])
        for p in self:
            generators = PauliStringCollection(
                p.get_anti_commutants(generators=generators)
            )
        return generators

    def get_graph(
        self, generators: list[PauliString] | PauliStringCollection | None = None
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the anticommutation graph whose vertices are the generators
        and edges are determined by the commutator between the vertices.

        Args:
            generators (list[PauliString] | PauliStringCollection): Area of Pauli strings over which
                to build a graph.
        Returns:
            Vertices, edges, and labels of edges of the anticommutation graph.
        """
        return get_graph(self.generators, commutators=generators)

    def get_commutator_graph(self) -> tuple[list[str], list[tuple[str, str]]]:
        """
        Get the commutator graph whose vertices are all Pauli strings
        of a given dimension and an edge between two vertices exist
        if there is an element in the generators with which
        one vertex anticommutes to give the other vertex.

        Returns:
            Vertices and edges of the commutator graph.
        """
        n = self.get_size()
        i = PauliString(n=n)
        return get_graph(i.get_commutants(), commutators=self, flag_labels=False)

    def get_frame_potential(self) -> int:
        """
        Get the frame potential of the system. The frame potential is a measure of quantum chaos.

        Returns:
            int: Frame potential of the system generated by the collection.
        """
        vertices, edges = self.get_commutator_graph()
        graph: nx.Graph[str] = nx.Graph()
        graph.add_nodes_from(vertices)
        graph.add_edges_from(edges)
        n_comp = nx.number_connected_components(graph)
        n_iso = len(list(nx.isolates(graph)))
        return n_comp * n_iso

    def _convert(
            self,
            generators: set[PauliString] | list[PauliString] | set[str] | list[str],
    ) -> PauliStringCollection:
        """
        Convert a set of Pauli strings or strings to a Pauli string collection.

        Args:
            generators: A set or list of PauliString objects or their string representations.
        Returns:
            A Pauli string collection formed from the Pauli strings.
        """
        processed_generators = []
        for g in generators:
            if isinstance(g, str):
                processed_generators.append(self.create_instance(pauli_str=g))
            else:
                processed_generators.append(g)
        return PauliStringCollection(processed_generators)

    def get_subgraphs(self) -> list[PauliStringCollection]:
        """
        Get the subgraphs of the anticommutation graph induced by the connected components.

        Returns:
            list[PauliStringCollection]: List of connected subgraphs.
        """
        # We don't need labels or edges labels here, just the graph structure.
        # Modified get_graph returns (vertices, edges) if flag_labels=False.
        # But let's just use the objects directly for efficiency.

        # Build adjacency using the same efficient logic as in get_graph
        # or just call get_graph and use the results.

        # Let's rebuild the graph with objects as nodes.
        g: nx.Graph[PauliString] = nx.Graph()
        g.add_nodes_from(self.generators)

        # Build overlap index
        qubit_to_indices: dict[int, list[int]] = {}
        for i, pauli in enumerate(self.generators):
            for qubit in pauli.get_support():
                if qubit not in qubit_to_indices:
                    qubit_to_indices[qubit] = []
                qubit_to_indices[qubit].append(i)

        for i, a in enumerate(self.generators):
            candidates = set()
            for qubit in a.get_support():
                if qubit in qubit_to_indices:
                    for j in qubit_to_indices[qubit]:
                        if j > i:
                            candidates.add(j)
            for j in candidates:
                b = self.generators[j]
                if not a | b: # Anticommuting
                    g.add_edge(a, b)

        return [self._convert(subgraph) for subgraph in
                sorted(nx.connected_components(g), key=len, reverse=True)]

    def classify(self) -> Classification:
        """
        Build the canonical graph of the generators and therefore classify its Lie algebra.

        Returns:
            Classification:
            Result of building canonical graphs corresponding to the direct sum Lie algebra.
        """
        subgraphs = self.get_subgraphs()
        self.classification = Classification()
        for subgraph in subgraphs:
            if not self.record:
                morph_factory = MorphFactory()
            else:
                morph_factory = RecordingMorphFactory(record=self.record)

            self.classification.add(morph_factory.build(subgraph.get()).get_morph())
        return self.classification

    def get_class(self) -> Classification:
        """
        Get the result of the classification of the generator set.

        Returns:
            Classification: Classification object of the generator set.
        """
        if self.classification is None:
            self.classification = self.classify()
        return self.classification

    def get_algebra(self) -> str:
        """
        Get the dynamical Lie algebra generated by the set of Pauli strings.

        Returns:
            str: Dynamical Lie algebra generated by the set of Pauli strings.
        """
        classification = self.get_class()
        return str(classification.get_algebra())

    def is_algebra(self, algebra: str) -> bool:
        """
        Checks whether the classified algebra is equal to the given algebra.

        Args:
            algebra (str): Name of the algebra.
        Returns:
            bool: True if the collection forms the specified algebra.
        """
        classification = self.get_class()
        return bool(classification.is_algebra(algebra))

    def is_in(self, generators: Self) -> bool:
        """
        Testing generators in algebra. All Pauli strings of one algebra are dependent on another.

        Args:
            generators (Self): Generator collection for checking.
        Returns:
            bool: True when all generators can be obtained from DLA self.
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
        Testing for equivalence of two algebras. All Pauli strings of one algebra are dependent on
        another.

        Args:
            generators (Self): Generator collection for checking.
        Returns:
            bool: True when two collections of generators are equivalent.
        """
        return self.is_in(generators) and generators.is_in(self)

    def select_dependents(self, generators: Self) -> PauliStringCollection:
        """
        Select dependents from self.

        Args:
            generators (Self): Generator collection for selecting.
        Returns:
           PauliStringCollection: Collection of dependent Pauli strings.
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
                dependents += morph_dependents

        return PauliStringCollection(dependents)

    def get_space(self) -> PauliStringCollection:
        """
        Get all space.

        Returns:
          PauliStringCollection:
          Collection of all Pauli strings formed by a collection of generators.
        """
        all_space = PauliStringCollection(
            [g for g in PauliString(n=self.get_len()).gen_all_pauli_strings()
            if g != PauliString(n=self.get_len())])
        return self.select_dependents(all_space)

    def gen_generators(self) -> Generator[PauliStringCollection, None, None]:
        """
        Get generators.

        Yields:
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
            int: Dimension of the classified dynamical Lie algebra.
        """
        return int(self.get_class().get_dla_dim())

    def get_dependents(self) -> PauliStringCollection:
        """
        Get a list of dependent strings in the collection.

        Returns:
            PauliStringCollection: List of dependent strings in the collection.
        """
        return PauliStringCollection(self.get_class().get_dependents())

    def get_independents(self) -> PauliStringCollection:
        """
        Get a list of independent strings in the collection.

        Returns:
            PauliStringCollection: List of independent strings in the collection.
        """
        dependents = self.get_class().get_dependents()
        return PauliStringCollection([v for v in self if v not in dependents])

    def get_canonic_vertices(self) -> PauliStringCollection:
        """
        Get a collection of vertices of the canonical graph of the collection.

        Returns:
            PauliStringCollection: Collection of vertices of the canonical graph of the collection.
        """
        return PauliStringCollection(self.get_class().get_vertices())

    def get_canonic_graph(
            self,
    ) -> tuple[list[str], list[tuple[str, str]], dict[tuple[str, str], str]]:
        """
        Get the canonical graph of the collection.

        Returns:
            Vertices, edges, and labels of edges.
        """
        classification = self.get_class()
        generators = PauliStringCollection(classification.get_vertices())
        return generators.get_graph()

    def get_anti_commutates(
            self, pauli_string: PauliString, generators: PauliStringCollection=None
    ) -> PauliStringCollection:
        """
        Get a collection of Pauli strings which anticommute with the given Pauli string.

        Args:
            pauli_string (PauliString): Pauli string against which anticommutation is checked.
            generators (PauliStringCollection, optional): Collection of Pauli strings to check.
                Default is None, in which case the generators of this collection are used.

        Returns:
            PauliStringCollection: Collection of anticommuting Pauli strings.
        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection([g for g in generators
                                      if g != pauli_string and not pauli_string | g])

    def get_commutates(self, pauli_string: PauliString,
                       generators: PauliStringCollection) -> PauliStringCollection:
        """
        Get a collection of Pauli strings which commute with the given Pauli string.

        Args:
            pauli_string (PauliString): Pauli string against which commutation is checked.
            generators (PauliStringCollection, optional): Collection of Pauli strings to check.
                Default is None, in which case the generators of this collection are used.

        Returns:
            PauliStringCollection: Collection of commuting Pauli strings.
        """
        if generators is None:
            generators = self.generators
        return PauliStringCollection(
            [g for g in generators if g != pauli_string and g | pauli_string]
        )

    def get_non_commuting_charges(self) -> PauliStringCollection:
        """
        Finds the non-commuting charges of the DLA generated by this collection.

        Non-commuting charges describe non-Abelian symmetries, i.e. elements of the
        stabilizer that do not commute with each other.

        Returns:
            PauliStringCollection: Set of non-commuting charges.
        """
        non_commuting_charges = PauliStringCollection()
        commutants = self.get_commutants()

        for c, q in combinations(commutants, 2):
            if not c | q:
                non_commuting_charges.append(c)
                non_commuting_charges.append(q)

        return non_commuting_charges

    def get_graph_components(self,
                             graph_type: str = 'anticommutator'
                             ) -> list[PauliStringCollection]:
        """
        Computes the connected components of the specified graph (anticommutator
        or commutator) constructed from the Pauli strings in the collection.

        Args:
            graph_type (str, optional): Type of graph to use for finding components.
                Must be either "anticommutator" (default) or "commutator".
        Returns:
            list[PauliStringCollection]:
            List of Pauli string collections, each representing a connected component of the
            selected graph.

        Raises:
            ValueError: If graph_type is not "anticommutator" or "commutator".
        """
        if graph_type == 'anticommutator':
            nodes, edges, _edge_labels = self.get_graph(self)  # The anti-commutation graph
        elif graph_type == 'commutator':
            nodes, edges = self.get_commutator_graph()  # The commutator graph
        else:
            raise ValueError("graph_type must be 'anticommutator' or 'commutator'")

        return self._get_connected_components(nodes, edges)

    def _get_connected_components(self, nodes: list[str],
                                  edges: list[tuple[str, str]]
                                  ) -> list[PauliStringCollection]:
        """
        Helper method to compute connected components from nodes and edges.

        Args:
            nodes (list[str]): List of node strings.
            edges (list[tuple[str, str]]): List of edges as tuples of node strings.
        Returns:
            list[PauliStringCollection]:
            List of Pauli string collections, each representing a connected component.
        """
        # Create a graph from the nodes and edges
        graph: nx.Graph[str] = nx.Graph()
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
        Private helper: Computes the quadratic symmetries for THIS collection,
        assuming THIS collection is a single connected component :math:`(C_{k})`
        and is provided with a pre-computed list of linear symmetries :math:`(L_{j})`.
        This method performs the inner loop of the full calculation.

        Args:
            linear_symmetries: Collection of linear symmetries.
        Returns:
            list[PauliStringLinear]: List of linear combinations of Pauli strings.
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
            1. Finds the linear symmetries :math:`(L_{j})` for this system.
            2. Finds the connected components :math:`(C_{k})` of this system's commutator graph.
            3. Loops through each component and calls the internal helper to get the :math:`Q_{kj}`.

        Args:
            normalized (bool): If True, returns an ORTHONORMAL basis.
            If False (default), returns an ORTHOGONAL basis.
        Returns:
            list[PauliStringLinear]:
            List of linear combinations of Pauli strings representing the basis.
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
        Find the index of a Pauli string in the collection.

        Args:
            pauli_string (PauliString): Pauli string to search for.
        Returns:
            int: Index of the Pauli string in the collection. Returns ``-1`` if not found.
        """
        for i, p in enumerate(self.generators):
            if p == pauli_string:
                return i
        return -1

    def replace(self, pauli_string: PauliString, new_pauli_string: PauliString) -> None:
        """
        Replace a Pauli string in the collection with another.

        Args:
            pauli_string (PauliString): Pauli string that will be replaced.
            new_pauli_string (PauliString): New Pauli string to replace the old Pauli string.
        Returns:
            None
        """
        index = self.find(pauli_string)
        if index != -1:
            self.generators[index] = new_pauli_string.copy()

    def contract(self, pauli_string: PauliString, contracted_pauli_string: PauliString) -> None:
        """
        Replace a Pauli string in the collection with its contraction with another Pauli string,

        Args:
            pauli_string (PauliString): Pauli string that will be replaced with the contracted Pauli
                string.
            contracted_pauli_string (PauliString): Pauli string to contract with.
        Returns:
            None
        """
        self.replace(pauli_string, pauli_string @ contracted_pauli_string)

    def list_connections(self) -> list[tuple[PauliString, PauliString, int, int]]:
        """
        Get list of edges in the anticommutation graph.

        Returns:
            list[tuple[PauliString, PauliString, int, int]]:
            List of pairs of connected vertices in the anticommutation graph.
        """
        return [(x, y, len(x.get_anti_commutants(self)), len(y.get_anti_commutants(self)))
                for x, y in combinations(self.generators, r=2) if not x | y]

    def _get_delta(self, generators: PauliStringCollection, number_connections: int) -> int:
        """
        Get the delta between the desired and actual number of edges in the anticommutation graph.

        Args:
            generators (PauliString): Generating collection of graph.
            number_connections (int): Required number of connections.
        Returns:
            Difference between the desired and actual number of edges in the anticommutation graph.
        """
        return number_connections - generators.get_anticommutation_pair()

    def find_generators_with_connection(self, number_connections: int) -> PauliStringCollection:
        """
        Get a subset of the generators with the desired number of anticommuting pairs.

        Args:
            number_connections (int): Required number of anticommuting pairs.
        Returns:
            PauliStringCollection:
            Collection of generators with the required number of anticommuting pairs.
        """
        generators = self.copy().get_canonic_vertices()
        max_iter = number_connections // 2
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

                gx.contract(x, y)
                gy.contract(y, x)
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

                    gx.contract(x, y)
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

    def nested_adjoint(self, target: PauliString | None) -> PauliString | None:
        """
        Return the nested commutator between the given Pauli string and all the generators.

        Args:
            target (PauliString|None): Given Pauli string.
        Returns:
            PauliString|None:
            Returns the nested commutator if it is non-zero, or None if it is zero.
        """
        current = target
        for op in reversed(self.generators):
            if current is None:
                return None
            current = op ^ current
        return current

    def evaluate_commutator_sequence(self) -> PauliString | None:
        """
        Evaluate the nested commutator between the first element of the collections and all others.

        Returns:
            PauliString|None:
            Returns the nested commutator if it is non-zero, or None if it is zero.
        """
        if len(self.generators) == 0:
            return None
        return PauliStringCollection(self.generators[1:]).nested_adjoint(self.generators[0])
