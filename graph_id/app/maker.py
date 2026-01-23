from __future__ import annotations

from collections import Counter

import networkx as nx
import numpy as np
from graph_id_cpp import GraphIDGenerator as CppGraphIDGenerator
from graph_id_cpp import MinimumDistanceNN as CppMinimumDistanceNN
from pymatgen.analysis.local_env import MinimumDistanceNN

from graph_id.core.graph_id import GraphIDGenerator as PyGraphIDGenerator

# from pymatgen.analysis.local_env import CrystalNN, MinimumDistanceNN


class GraphIDMaker:

    """A simple, high-level interface for generating Graph IDs.

    GraphIDMaker provides an easy-to-use API for generating unique identifiers
    for crystal and molecular structures. It automatically selects the appropriate
    backend (C++ or Python) and handles common configuration.

    Examples
    --------
    >>> from pymatgen.core import Structure, Lattice
    >>> from graph_id import GraphIDMaker
    >>> structure = Structure.from_spacegroup(
    ...     "Fm-3m", Lattice.cubic(5.692),
    ...     ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    ... )
    >>> maker = GraphIDMaker()
    >>> maker.get_id(structure)
    'NaCl-88c8e156db1b0fd9'
    """

    def __init__(
        self,
        nn=None,
        depth: int | None = None,
        reduce: bool = False,
        engine: str = "c++",
    ) -> None:
        """Initialize the GraphIDMaker.

        Parameters
        ----------
        nn : NearNeighbors, optional
            A neighbor-finding strategy object. When using the C++ engine,
            you must supply a C++ implementation (e.g., `graph_id_cpp.MinimumDistanceNN`).
            If None, uses the default MinimumDistanceNN for the selected engine.
        depth : int, optional
            Fixed traversal depth for the graph. If None (default), the depth
            is dynamically determined based on the graph diameter.
        reduce : bool, default False
            If True, reduces site sequences by their GCD. This is useful for
            comparing primitive and conventional cells of the same structure.
        engine : str, default "c++"
            The computation engine to use. Options are:

            - ``"c++"`` : Fast C++ implementation (recommended)
            - ``"py"`` or ``"python"`` : Pure Python implementation

        Examples
        --------
        >>> maker = GraphIDMaker()  # Default C++ engine
        >>> maker = GraphIDMaker(engine="py")  # Python engine
        >>> maker = GraphIDMaker(depth=5)  # Fixed depth
        >>> maker = GraphIDMaker(reduce=True)  # Enable reduction
        """
        self.reduce = reduce

        if "py" in engine.lower():
            self.engine = "python"
            if nn is None:
                nn = MinimumDistanceNN()
        elif "c" in engine.lower():
            self.engine = "c++"
            if nn is None:
                nn = CppMinimumDistanceNN()

        diameter_factor = 2
        additional_depth = 1
        if depth is not None:
            diameter_factor = 0
            additional_depth = depth

        if engine == "py":
            self.generator = PyGraphIDGenerator(
                nn=nn,
                diameter_factor=diameter_factor,
                additional_depth=additional_depth,
                prepend_composition=False,
                prepend_dimensionality=False,
            )

        elif engine == "c++":
            self.generator = CppGraphIDGenerator(
                nn=nn,
                diameter_factor=diameter_factor,
                additional_depth=additional_depth,
            )

    def get_id(self, structure) -> str:
        """Generate a Graph ID for the given structure.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object representing the crystal or molecule.

        Returns
        -------
        str
            The Graph ID in the format ``"{formula}-{hash}"``, where:

            - ``formula`` is the reduced chemical formula
            - ``hash`` is a 16-character hexadecimal topological fingerprint

        Examples
        --------
        >>> from pymatgen.core import Structure
        >>> structure = Structure.from_file("NaCl.cif")
        >>> maker = GraphIDMaker()
        >>> maker.get_id(structure)
        'NaCl-88c8e156db1b0fd9'
        """
        graph_id = self.get_id_reducing_site_sequences(structure) if self.reduce else self.generator.get_id(structure)

        return f"{structure.composition.reduced_formula}-{graph_id}"

    def get_id_reducing_site_sequences(self, structure):
        """Generate a Graph ID with reduced site sequences.

        This method divides repeated compositional sequences by their GCD,
        which allows primitive and conventional cells of the same structure
        to produce the same Graph ID.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.

        Returns
        -------
        str
            The reduced Graph ID hash (without formula prefix).

        Notes
        -----
        This is an internal method called when ``reduce=True``.
        """
        sg = self.generator.prepare_structure_graph(structure)

        gcd_list = []
        components_counters = []

        for component in sg.cc_cs:
            _counter = Counter(component["cs_list"])
            _gcd = np.gcd.reduce(list(_counter.values()))
            gcd_list.append(_gcd)
            components_counters.append(_counter)

        divider = min(gcd_list)

        labels_list = []
        for counter in components_counters:
            labels = []
            for label, count in counter.items():
                labels += [label] * int(count / divider)

            labels_list.append(self.generator._join_cs_list(labels))
        return self.generator._component_strings_to_whole_id(labels_list)

    def get_site_ids(self, structure):
        """Get unique compositional sequence identifiers for each atomic site.

        Each site in the structure receives a unique identifier based on its
        local chemical environment (compositional sequence). Sites with identical
        environments will have identical IDs.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.

        Returns
        -------
        dict
            A dictionary mapping site indices (int) to their compositional
            sequence identifiers (str). The format is ``{site_index: "Element_hash"}``.

        Examples
        --------
        >>> maker = GraphIDMaker()
        >>> site_ids = maker.get_site_ids(nacl_structure)
        >>> for idx, cs in site_ids.items():
        ...     print(f"Site {idx}: {cs[:20]}...")
        Site 0: Na_Na-8ac4127fe153ff...
        Site 1: Cl_Cl-0d9c88475b88c4...

        Notes
        -----
        For the C++ engine, constructs site IDs from ``cc_nodes`` and ``cc_cs``.
        For the Python engine, uses NetworkX node attributes directly.
        """
        sg = self.generator.prepare_structure_graph(structure)

        if self.engine == "c++":
            # For C++ engine, construct site IDs from cc_nodes and cc_cs
            site_ids = {}
            labels = sg.labels
            cc_nodes = sg.cc_nodes
            # cc_cs_labels is the list of lists of compositional sequences
            cc_cs = sg.cc_cs_labels

            for cc_idx, nodes in enumerate(cc_nodes):
                cs_list = cc_cs[cc_idx]
                for node_idx, cs_string in zip(nodes, cs_list):
                    # Format: label + "_" + compositional_sequence
                    site_ids[node_idx] = f"{labels[node_idx]}_{cs_string}"

            return site_ids

        # For Python engine, use NetworkX node attributes
        return nx.get_node_attributes(sg.graph, "compositional_sequence")
