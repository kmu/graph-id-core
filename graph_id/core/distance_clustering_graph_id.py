from copy import deepcopy
from hashlib import blake2b

import networkx as nx
import numpy as np
from graph_id.analysis.local_env import DistanceClusteringNN
from pymatgen.analysis.local_env import MinimumDistanceNN
from graph_id.analysis.graphs import StructureGraph
from graph_id.core.graph_id import GraphIDGenerator
from pymatgen.core import Element

__version__ = "0.1.0"


def blake(s):
    return blake2b(s.encode()).hexdigest()


class DistanceClusteringGraphID(GraphIDGenerator):
    def __init__(  # noqa: PLR0913
        self,
        nn=None,
        wyckoff=False,
        depth_factor=2,
        additional_depth=1,
        symmetry_tol=0.1,
        topology_only=False,
        loop=False,
        rank_k=3,
        cutoff=6.0,
        digest_size=8,
    ) -> None:
        super().__init__(
            nn,
            wyckoff,
            depth_factor,
            additional_depth,
            symmetry_tol,
            topology_only,
            loop,
            digest_size,
        )

        self.rank_k = rank_k
        self.cutoff = cutoff
        self.digest_size = digest_size

        if nn is None:
            self.nn = DistanceClusteringNN()
        else:
            self.nn = nn

    def get_id(self, structure):
        gid_list = []
        _sg = StructureGraph.with_local_env_strategy(structure, MinimumDistanceNN())
        for cluster_idx in range(self.rank_k):
            long_str_list = []
            # _sg = StructureGraph.with_local_env_strategy(structure, MinimumDistanceNN())
            for idx in range(len(structure)):
                copied_sg = deepcopy(_sg)
                # まず原子idxが含まれる結合を削除する
                for from_index, to_index, dct in _sg.graph.edges(keys=False, data=True):
                    if from_index == idx or to_index == idx:
                        copied_sg.break_edge(from_index, to_index, dct["to_jimage"], allow_reverse=True)
                sg = self.prepare_structure_graph(structure, copied_sg, idx, cluster_idx)
                n = len(sg.cc_cs)
                array = np.empty(
                    [
                        n,
                    ],
                    dtype=object,
                )
                for i, component in enumerate(sg.cc_cs):
                    array[i] = blake("-".join(sorted(component["cs_list"])))
                    # array[i] = blake2b("-".join(sorted(component["cs_list"])).encode("ascii"), digest_size=16).hexdigest()
                long_str_tmp = ":".join(np.sort(array))
                # long_str_tmp = blake2b(":".join(np.sort(array)).encode("ascii"), digest_size=16).hexdigest()
                long_str_list.append(long_str_tmp)
            long_str = ":".join(np.sort(long_str_list))
            gid = blake2b(long_str.encode("ascii"), digest_size=self.digest_size).hexdigest()
            gid_list.append(gid)

        long_gid = "".join(gid_list)
        # return self.elaborate_comp_dim(sg, blake2b(long_gid.encode("ascii"), digest_size=16).hexdigest())
        return blake2b(long_gid.encode("ascii"), digest_size=self.digest_size).hexdigest()

    def prepare_structure_graph(self, structure, _sg, n, rank_k):
        
        sg = StructureGraph.with_indivisual_state_comp_strategy(
            structure=structure,
            strategy=self.nn,
            _sg=_sg,
            n=n,
            rank_k=rank_k,
            cutoff=self.cutoff,
        )

        use_previous_cs = False

        compound = sg.structure
        prev_num_uniq = len(compound.composition)

        if self.topology_only:
            for site_i in range(len(sg.structure)):
                sg.structure.replace(site_i, Element("H"))

        if self.wyckoff:
            sg.set_wyckoffs(symmetry_tol=self.symmetry_tol)
            prev_num_uniq = len(list(set(nx.get_node_attributes(sg.graph, "compositional_sequence").values())))

        elif self.loop:
            sg.set_loops_as_starting_labels(
                depth_factor=self.depth_factor,
                additional_depth=self.additional_depth,
            )

        else:
            sg.set_elemental_labels()

        while True:
            sg.set_indivisual_compositional_sequence_node_attr(
                n=n,
                hash_cs=False,
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
