import functools

import networkx as nx
import numpy as np
from graph_id.analysis.compositional_sequence import CompositionalSequence
from networkx.algorithms.distance_measures import diameter
from pymatgen.analysis.graphs import StructureGraph as PmgStructureGraph
from pymatgen.core import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from itertools import combinations
from hashlib import blake2b


def standardize_loop(lst):
    lst2 = list(reversed(lst))
    starting_point = lst2.pop(-1)
    lst2.insert(0, starting_point)

    return sorted([lst, lst2], key=lambda x: "".join(x))[-1]


class SiteOnlySpeciesString:
    def __init__(self, species_string):
        self.species_string = species_string


class ConnectedSiteLight:
    def __init__(
        self,
        site,
        jimage,
        index,
        weight,
        dist,
    ):
        self.site = SiteOnlySpeciesString(site.species_string)
        self.jimage = jimage
        self.index = index
        self.weight = weight
        self.dist = dist


class StructureGraph(PmgStructureGraph):  # type: ignore
    @staticmethod
    def from_pymatgen_structure_graph(sg: PmgStructureGraph):
        graph_data = sg.as_dict()["graphs"]

        return StructureGraph(sg.structure, graph_data)

    # Copied from original pymatgen with modifications
    @staticmethod
    def with_local_env_strategy(structure, strategy, weights=False):
        """
        Constructor for StructureGraph, using a strategy
        from :Class: `pymatgen.analysis.local_env`.

        :param structure: Structure object
        :param strategy: an instance of a
            :Class: `pymatgen.analysis.local_env.NearNeighbors` object
        :param weights: if True, use weights from local_env class
            (consult relevant class for their meaning)
        :return:
        """

        if not strategy.structures_allowed:
            raise ValueError(
                "Chosen strategy is not designed for use with structures! " "Please choose another strategy."
            )

        sg = StructureGraph.from_empty_graph(structure, name="bonds")

        for n, neighbors in enumerate(strategy.get_all_nn_info(structure)):
            for neighbor in neighbors:
                # local_env will always try to add two edges
                # for any one bond, one from site u to site v
                # and another form site v to site u: this is
                # harmless, so warn_duplicates=False
                sg.add_edge(
                    from_index=n,
                    from_jimage=(0, 0, 0),
                    to_index=neighbor["site_index"],
                    to_jimage=neighbor["image"],
                    weight=neighbor["weight"] if weights else None,
                    warn_duplicates=False,
                )

        return sg
    
    @staticmethod
    def with_indivisual_state_comp_strategy(structure, strategy, _sg, n, weights=False, rank_k=1, cutoff=6.0):
        """
        Constructor for StructureGraph, using a StateCompNN strategy
        from :Class: `chemsys.pymatgen.analysis.local_env`.
        :param structure: Structure object
        :param strategy: an instance of StateCompNN
        :param n: (int) an index of focused site
        :param weights: if True, use weights from local_env class
            (consult relevant class for their meaning)
        :rank_k: (int) cluster_idx
        :cutoff: (float)
        :return:
        """

        if not strategy.structures_allowed:
            raise ValueError(  # noqa: TRY003
                "Chosen strategy is not designed for use with structures!",  # noqa: EM101
            )

        nn_info = strategy.get_nn_info(structure, n, rank_k, cutoff)

        for neighbor in nn_info:
            # local_env will always try to add two edges
            # for any one bond, one from site u to site v
            # and another form site v to site u: this is
            # harmless, so warn_duplicates=False
            _sg.add_edge(
                from_index=n,
                from_jimage=(0, 0, 0),
                to_index=neighbor["site_index"],
                to_jimage=neighbor["image"],
                weight=neighbor["weight"] if weights else None,
                warn_duplicates=False,
                edge_properties=neighbor["edge_properties"],
            )

        return _sg

    def set_elemental_labels(self):
        self.starting_labels = [site.species_string for site in self.structure]

    def get_connected_sites_light(self, n, jimage=(0, 0, 0)):
        """
        A light version of get_connected_sites.
        periodic_site -> SiteOnlySpeciesString
        """

        connected_sites = set()
        connected_site_images = set()

        out_edges = [(u, v, d, "out") for u, v, d in self.graph.out_edges(n, data=True)]
        in_edges = [(u, v, d, "in") for u, v, d in self.graph.in_edges(n, data=True)]

        for u, v, d, dir in out_edges + in_edges:

            to_jimage = d["to_jimage"]

            if dir == "in":
                u, v = v, u
                to_jimage = np.multiply(-1, to_jimage)

            to_jimage = tuple(map(int, np.add(to_jimage, jimage)))

            if (v, to_jimage) not in connected_site_images:
                connected_site = ConnectedSiteLight(
                    site=self.structure[v], jimage=to_jimage, index=v, weight=None, dist=None
                )

                connected_sites.add(connected_site)
                connected_site_images.add((v, to_jimage))

        _connected_sites = list(connected_sites)

        return _connected_sites

    def set_wyckoffs(self, symmetry_tol: float = 0.1) -> None:
        siteless_strc = self.structure.copy()

        for site_i in range(len(self.structure)):
            siteless_strc.replace(site_i, Element("H"))

        sga = SpacegroupAnalyzer(siteless_strc)
        sym_dataset = sga.get_symmetry_dataset()

        if sym_dataset is None:
            self.set_elemental_labels()
            return None

        wyckoffs = sym_dataset["wyckoffs"]
        number = sym_dataset["number"]

        attribute_values = {}

        self.starting_labels = []
        for site_i, w in enumerate(wyckoffs):
            attribute_values[site_i] = f"{self.structure[site_i].species_string}_{w}_{number}"
            self.starting_labels.append(f"{self.structure[site_i].species_string}_{w}_{number}")

    def set_compositional_sequence_node_attr(
        self,
        hash_cs: bool = False,
        wyckoff: bool = False,
        additional_depth: int = 0,
        depth_factor: int = 2,
        use_previous_cs: bool = False,
    ) -> None:

        node_attributes = {}
        self.cc_cs = []
        get_connected_sites_light = functools.lru_cache(maxsize=None)(self.get_connected_sites_light)

        ug = self.graph.to_undirected()

        for cc in nx.connected_components(ug):
            cs_list = []

            d = diameter(ug.subgraph(cc))

            for focused_site_i in cc:

                depth = depth_factor * d + additional_depth

                cs = CompositionalSequence(
                    focused_site_i=focused_site_i,
                    starting_labels=self.starting_labels,
                    hash_cs=hash_cs,
                    use_previous_cs=use_previous_cs or wyckoff,
                )

                for _ in range(depth):
                    for c_site in cs.get_current_starting_sites():
                        nsites = get_connected_sites_light(c_site[0], c_site[1])
                        cs.count_composition_for_neighbors(nsites)

                    cs.finalize_this_depth()

                this_cs = str(cs)

                node_attributes[focused_site_i] = self.starting_labels[focused_site_i] + "_" + this_cs
                cs_list.append(this_cs)

            self.cc_cs.append({"site_i": cc, "cs_list": cs_list})

        nx.set_node_attributes(self.graph, values=node_attributes, name="compositional_sequence")

    def get_loops(self, depth: int, index: int, shortest: bool = True):  # noqa: C901
        """
        各原子を起点としてループを計算し、そのインデックス情報を返す。

        Parameters:
            indices: ループの起点としたいインデックス
            depth: ループの最大の大きさ

        Returns:
            [[(index, image), ...], ...]
        """

        get_connected_sites = functools.lru_cache(maxsize=None)(self.get_connected_sites)

        def find_all_rings(index, ring_list):
            neighbors = get_connected_sites(index, (0, 0, 0))
            for n0, n1 in combinations(neighbors, 2):
                found = False
                for ring in ring_list:
                    term0 = ring[1]
                    term1 = ring[-2]

                    if all(
                        (
                            n0.index == term0[0],
                            n0.jimage == term0[1],
                            n1.index == term1[0],
                            n1.jimage == term1[1],
                        ),
                    ):
                        found = True
                        break

                    if all(
                        (
                            n1.index == term0[0],
                            n1.jimage == term0[1],
                            n0.index == term1[0],
                            n0.jimage == term1[1],
                        ),
                    ):
                        found = True
                        break

                if found is False:
                    return False

            return True

        def get_further_lines_from_lines(lines):
            new_lines = []
            for line in lines:
                ind, image = line[-1]
                neighbors = get_connected_sites(ind, image)

                for n in neighbors:
                    new_line = [*line, (n.index, n.jimage)]

                    # 戻らない場合のみ。
                    if len(new_line[:-1]) == len(set(new_line[:-1])):
                        new_lines.append(new_line)

            return new_lines

        lines = []
        lines.append([(index, (0, 0, 0))])

        ring_list = []

        for depth_i in range(depth):
            next_lines = []
            lines = get_further_lines_from_lines(lines)

            for line in lines:
                # 前と後ろが同じ
                if line[0] == line[-1]:
                    if depth_i > 1 and list(reversed(line)) not in ring_list:
                        ring_list.append(line)
                else:
                    next_lines.append(line)

            lines = next_lines

            # ここで理論上の値に達したら探索を打ち切る
            if shortest and find_all_rings(index, ring_list):
                return ring_list

        return list(ring_list)

    def set_loops(self, depth_factor: int, additional_depth: int) -> None:
        self.starting_labels = []

        undirected_graph = self.graph.to_undirected()

        max_diameter = 0
        for cc in nx.connected_components(undirected_graph):
            d = diameter(undirected_graph.subgraph(cc))
            if d > max_diameter:
                max_diameter = d

        depth = max_diameter * depth_factor + additional_depth

        for site_i in range(len(self.graph.nodes)):
            all_loops = self.get_loops(depth=depth, index=site_i)
            all_loop_strings = []
            # print(all_loops)
            for loop in all_loops:
                loop_elements = []
                for site_i_jimage in loop:
                    loop_species_string = self.structure[site_i_jimage[0]].species_string
                    # print(loop_species_string)
                    loop_elements.append(loop_species_string)

                loop_elements = standardize_loop(loop_elements)

                seed_str = "-".join(loop_elements)
                hashed_loop = blake2b(seed_str.encode(), digest_size=8).hexdigest()

                all_loop_strings.append(hashed_loop)

            seed_str_all_loops = ":".join(sorted(all_loop_strings))
            hashed_all_loops = blake2b(seed_str_all_loops.encode(), digest_size=8).hexdigest()

            self.starting_labels.append(hashed_all_loops)

    def set_indivisual_compositional_sequence_node_attr(
        self,
        n: int,
        hash_cs: bool = False,
        wyckoff: bool = False,
        additional_depth: int = 0,
        depth_factor: int = 2,
        use_previous_cs: bool = False,
    ) -> None:
        node_attributes = {}
        self.cc_cs = []
        get_connected_sites_light = functools.lru_cache(maxsize=None)(self.get_connected_sites_light)

        ug = self.graph.to_undirected()

        for cc in nx.connected_components(ug):
            cs_list = []

            d = diameter(ug.subgraph(cc))

            if n in cc:
                depth = depth_factor * d + additional_depth

                cs = CompositionalSequence(
                    focused_site_i=n,
                    starting_labels=self.starting_labels,
                    hash_cs=hash_cs,
                    use_previous_cs=use_previous_cs or wyckoff,
                )

                for _this_depth in range(depth):
                    for c_site in cs.get_current_starting_sites():
                        nsites = get_connected_sites_light(c_site[0], c_site[1])
                        cs.count_composition_for_neighbors(nsites)

                    cs.finalize_this_depth()

                this_cs = str(cs)

                node_attributes[n] = self.starting_labels[n] + "_" + this_cs
                cs_list.append(this_cs)

                self.cc_cs.append({"site_i": cc, "cs_list": cs_list})

        nx.set_node_attributes(self.graph, values=node_attributes, name="compositional_sequence")
