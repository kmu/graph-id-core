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
    def __init__(
        self,
        nn=None,
        depth: int | None = None,
        reduce: bool = False,
        engine: str = "c++",
    ) -> None:
        """
        A simple interface to make GraphID.

        nn: NearNeighbor object to use for neighbor finding.
            You must supply C++ implementation if you use `c++` engine.
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
        graph_id = self.get_id_reducing_site_sequences(structure) if self.reduce else self.generator.get_id(structure)

        return f"{structure.composition.reduced_formula}-{graph_id}"

    def get_id_reducing_site_sequences(self, structure):
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
        """
        Get site IDs for a structure.
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
