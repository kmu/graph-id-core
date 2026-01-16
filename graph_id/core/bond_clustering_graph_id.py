from hashlib import blake2b

import networkx as nx
import numpy as np
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.core import Element

from graph_id.analysis.graphs import StructureGraph
from graph_id.analysis.local_env import BondClusteringNN
from graph_id.core.graph_id import GraphIDGenerator

__version__ = "0.1.0"


def blake(s):
    return blake2b(s.encode()).hexdigest()


class BondClusteringGraphID(GraphIDGenerator):
    def __init__(
        self,
        nn=None,
        wyckoff=False,
        diameter_factor=2,
        additional_depth=1,
        symmetry_tol=0.1,
        topology_only=False,
        loop=False,
        digest_size=8,
    ) -> None:
        super().__init__(
            nn,
            wyckoff,
            diameter_factor,
            additional_depth,
            symmetry_tol,
            topology_only,
            loop,
            digest_size,
        )

        self.digest_size = digest_size

        if nn is None:
            self.nn = BondClusteringNN()
        else:
            self.nn = nn

    def _join_cs_list(self, cs_list):
        return blake("-".join(sorted(cs_list)))

    def _component_strings_to_whole_id(self, component_strings):
        long_str = ":".join(np.sort(component_strings))
        return blake2b(long_str.encode("ascii"), digest_size=self.digest_size).hexdigest()

    def get_id(self, structure):
        sg = self.prepare_structure_graph(structure)
        n = len(sg.cc_cs)
        array = np.empty(
            [
                n,
            ],
            dtype=object,
        )
        for i, component in enumerate(sg.cc_cs):
            array[i] = self._join_cs_list(component["cs_list"])
        gid = self._component_strings_to_whole_id(array)

        return self.elaborate_comp_dim(sg, gid)

    def elaborate_comp_dim(self, sg, gid):
        if self.prepend_dimensionality:
            dim = get_dimensionality_larsen(sg)
            gid = f"{dim}D-{gid}"

        if self.prepend_composition and not self.topology_only:
            gid = f"{sg.structure.composition.reduced_formula}-{gid}"

        return gid

    @property
    def version(self):
        return __version__

    def prepare_structure_graph(self, structure):
        sg = StructureGraph.with_local_env_strategy(structure, self.nn, weights=True)
        use_previous_cs = False

        compound = sg.structure
        prev_num_uniq = len(compound.composition)

        if self.topology_only:
            for site_i in range(len(sg.structure)):
                sg.structure.replace(site_i, Element("H"))

        if self.wyckoff:
            sg.set_wyckoffs(symmetry_tol=self.symmetry_tol)

            # remove nx?
            prev_num_uniq = len(list(set(nx.get_node_attributes(sg.graph, "compositional_sequence").values())))

        elif self.loop:
            sg.set_loops(
                diameter_factor=self.diameter_factor,
                additional_depth=self.additional_depth,
            )

        else:
            sg.set_elemental_labels()

        while True:
            sg.set_compositional_sequence_node_attr(
                hash_cs=True,
                wyckoff=self.wyckoff,
                additional_depth=self.additional_depth,
                diameter_factor=self.diameter_factor,
                use_previous_cs=use_previous_cs or self.wyckoff,
            )

            num_unique_nodes = len(list(set(nx.get_node_attributes(sg.graph, "compositional_sequence").values())))
            use_previous_cs = True

            if prev_num_uniq == num_unique_nodes:
                return sg

            prev_num_uniq = num_unique_nodes
