from __future__ import annotations

import multiprocessing as multi
from hashlib import blake2b
from multiprocessing import Pool
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.core import Element
from tqdm import tqdm

from graph_id.analysis.graphs import StructureGraph

if TYPE_CHECKING:
    from pymatgen.core.structure import Structure


__version__ = "0.1.0"


def blake(s):
    return blake2b(s.encode()).hexdigest()


class GraphIDGenerator:
    def __init__(
        self,
        nn=None,
        wyckoff=False,
        depth_factor=2,
        additional_depth=1,
        symmetry_tol=0.1,
        topology_only=False,
        loop=False,
        digest_size=8,
    ):
        """
        comp_dim: include composition and dimensionality as the prefix
        """

        if wyckoff and loop:
            msg = "wyckoff and loop cannot be True at the same time"
            raise ValueError(msg)

        if loop and topology_only:
            msg = "loop and topology_only cannot be True at the same time"
            raise ValueError(msg)

        if nn is None:
            self.nn = MinimumDistanceNN()
        else:
            self.nn = nn

        self.wyckoff = wyckoff
        self.additional_depth = additional_depth
        self.depth_factor = depth_factor
        self.symmetry_tol = symmetry_tol
        self.topology_only = topology_only
        self.loop = loop
        self.digest_size = digest_size

    # def get_graph_I#

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
            array[i] = blake("-".join(sorted(component["cs_list"])))
        long_str = ":".join(np.sort(array))
        gid = blake2b(long_str.encode("ascii"), digest_size=self.digest_size).hexdigest()

        return self.elaborate_comp_dim(sg, gid)

    def elaborate_comp_dim(self, sg, gid):
        dim = get_dimensionality_larsen(sg)
        gid = f"{dim}D-{gid}"

        if not self.topology_only:
            gid = f"{sg.structure.composition.reduced_formula}-{gid}"

        return gid

    @property
    def version(self):
        return __version__

    def get_id_catch_error(self, structure):
        try:
            return self.get_id(structure)
        except Exception:  # noqa: BLE001
            return ""

    def get_many_ids(self, structures, parallel=False):
        if parallel:
            n_cores = multi.cpu_count()
            # ctx = multi.get_context("spawn")
            # p = ctx.Pool(n_cores)
            p = Pool(n_cores)
            imap = p.imap(self.get_id_catch_error, structures)
            # ids = p.map(self.get_id, structures)
            return list(tqdm(imap, total=len(structures)))

        return [self.get_id(s) for s in structures]

    def get_component_ids(self, structure):
        sg = self.prepare_structure_graph(structure)
        cc_gid = np.empty(
            [
                len(sg.cc_cs),
            ],
            dtype=object,
        )
        for i, component in enumerate(sg.cc_cs):
            each_long_str = blake("-".join(sorted(component["cs_list"])))
            gid = blake2b(each_long_str.encode("ascii"), digest_size=16).hexdigest()
            # cc_gid[] = gid
            cc_gid[i] = {"site_i": component["site_i"], "graph_id": gid}

        return cc_gid

    def are_same(self, structure1, structure2):
        return self.get_id(structure1) == self.get_id(structure2)

    def prepare_structure_graph(self, structure):
        sg = StructureGraph.with_local_env_strategy(structure, self.nn)
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
                depth_factor=self.depth_factor,
                additional_depth=self.additional_depth,
            )

        else:
            sg.set_elemental_labels()

        while True:
            sg.set_compositional_sequence_node_attr(
                hash_cs=True,
                wyckoff=self.wyckoff,
                additional_depth=self.additional_depth,
                depth_factor=self.depth_factor,
                use_previous_cs=use_previous_cs or self.wyckoff,
            )

            num_unique_nodes = len(list(set(nx.get_node_attributes(sg.graph, "compositional_sequence").values())))
            use_previous_cs = True

            if prev_num_uniq == num_unique_nodes:
                break

            prev_num_uniq = num_unique_nodes

        return sg

    def get_unique_structures(self, structures: list[Structure]) -> list[Structure]:
        unique_structures = []
        graph_ids = set()

        for strct in structures:
            new_graph_id = self.get_id(strct)
            if new_graph_id not in graph_ids:
                graph_ids.add(new_graph_id)
                unique_structures.append(strct)

        return unique_structures
