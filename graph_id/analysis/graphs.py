import functools

import networkx as nx
import numpy as np
from graph_id.analysis.compositional_sequence import CompositionalSequence
from networkx.algorithms.distance_measures import diameter
from pymatgen.analysis.graphs import StructureGraph as PmgStructureGraph
from pymatgen.core import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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

        sg = StructureGraph.with_empty_graph(structure, name="bonds")

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
