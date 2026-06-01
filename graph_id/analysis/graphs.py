import functools
from hashlib import blake2b
from itertools import combinations

import networkx as nx
import numpy as np
from networkx.algorithms.distance_measures import diameter
from pymatgen.analysis.graphs import StructureGraph as PmgStructureGraph
from pymatgen.core import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from graph_id.analysis.compositional_sequence import CompositionalSequence


def standardize_loop(lst):
    """Standardize a loop to its canonical form.

    Compares the loop with its reverse and returns the lexicographically
    larger version to ensure consistent representation.

    Parameters
    ----------
    lst : list
        A list of element symbols representing a loop.

    Returns
    -------
    list
        The standardized loop.

    """
    lst2 = list(reversed(lst))
    starting_point = lst2.pop(-1)
    lst2.insert(0, starting_point)

    return sorted([lst, lst2], key=lambda x: "".join(x))[-1]


class SiteOnlySpeciesString:

    """Lightweight wrapper containing only the species string of a site.

    Used to reduce memory when only element information is needed.

    Parameters
    ----------
    species_string : str
        The species string (e.g., "Na", "Cl").

    """

    def __init__(self, species_string):
        """Initialize with the species string."""
        self.species_string = species_string


class ConnectedSiteLight:

    """Lightweight representation of a connected site.

    A memory-efficient alternative to pymatgen's ConnectedSite that stores
    only essential information for compositional sequence computation.

    Parameters
    ----------
    site : Site
        The pymatgen Site object (only species_string is retained).
    jimage : tuple
        The periodic image indices (i, j, k).
    index : int
        The site index in the structure.
    weight : float or None
        The bond weight (if applicable).
    dist : float or None
        The distance to the connected site.

    Attributes
    ----------
    site : SiteOnlySpeciesString
        Lightweight site wrapper.
    jimage : tuple
        The periodic image.
    index : int
        The site index.
    weight : float or None
        The bond weight.
    dist : float or None
        The distance.

    """

    def __init__(
        self,
        site,
        jimage,
        index,
        weight,
        dist,
    ):
        """Initialize a connected site with lightweight data."""
        self.site = SiteOnlySpeciesString(site.species_string)
        self.jimage = jimage
        self.index = index
        self.weight = weight
        self.dist = dist


class StructureGraph(PmgStructureGraph):  # type: ignore

    """Extended StructureGraph with methods for Graph ID computation.

    This class extends pymatgen's StructureGraph with additional functionality
    for computing compositional sequences and handling loops/rings in the
    structure graph.

    Attributes
    ----------
    starting_labels : list of str
        Labels for each site used as starting points for compositional sequences.
    cc_cs : list of dict
        Compositional sequences for each connected component.
        Each dict contains ``site_i`` (set of site indices) and ``cs_list``
        (list of compositional sequence strings).

    See Also
    --------
    pymatgen.analysis.graphs.StructureGraph : Base class

    """

    @staticmethod
    def from_pymatgen_structure_graph(sg: PmgStructureGraph):
        """Create a StructureGraph from a pymatgen StructureGraph.

        Parameters
        ----------
        sg : PmgStructureGraph
            A pymatgen StructureGraph object.

        Returns
        -------
        StructureGraph
            A new StructureGraph instance.

        """
        graph_data = sg.as_dict()["graphs"]

        return StructureGraph(sg.structure, graph_data)

    @staticmethod
    def with_local_env_strategy(structure, strategy, weights=False):
        """Create a StructureGraph using a neighbor-finding strategy.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.
        strategy : NearNeighbors
            A neighbor-finding strategy from pymatgen.analysis.local_env,
            such as MinimumDistanceNN, CrystalNN, etc.
        weights : bool, default False
            If True, include bond weights from the strategy.

        Returns
        -------
        StructureGraph
            A new StructureGraph with edges representing bonds.

        Raises
        ------
        ValueError
            If the strategy does not support structures.

        Examples
        --------
        >>> from pymatgen.analysis.local_env import MinimumDistanceNN
        >>> sg = StructureGraph.with_local_env_strategy(
        ...     structure, MinimumDistanceNN()
        ... )

        """
        if not strategy.structures_allowed:
            msg = "Chosen strategy is not designed for use with structures! Please choose another strategy."
            raise ValueError(msg)

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
        """Add edges for a specific site using a distance clustering strategy.

        This method is used by DistanceClusteringGraphID to add bonds
        for a specific site and distance cluster.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.
        strategy : DistanceClusteringNN
            A distance clustering neighbor-finding strategy.
        _sg : StructureGraph
            An existing StructureGraph to modify.
        n : int
            The site index to add edges for.
        weights : bool, default False
            If True, include bond weights from the strategy.
        rank_k : int, default 1
            The distance cluster index (0-based).
        cutoff : float, default 6.0
            Maximum distance cutoff in Angstroms.

        Returns
        -------
        StructureGraph
            The modified StructureGraph with new edges.

        Raises
        ------
        ValueError
            If the strategy does not support structures.

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
        """Set element symbols as starting labels for compositional sequences.

        This is the default labeling scheme where each site is labeled
        by its element symbol (e.g., "Na", "Cl").
        """
        self.starting_labels = [site.species_string for site in self.structure]

    def get_connected_sites_light(self, n, jimage=(0, 0, 0)):
        """Get connected sites with minimal memory footprint.

        A lightweight version of get_connected_sites that returns
        ConnectedSiteLight objects instead of full ConnectedSite objects.

        Parameters
        ----------
        n : int
            The site index to get neighbors for.
        jimage : tuple, default (0, 0, 0)
            The periodic image of the site.

        Returns
        -------
        list of ConnectedSiteLight
            List of connected sites with minimal information.

        """
        connected_sites = set()
        connected_site_images = set()

        out_edges = [(u, v, d, "out") for u, v, d in self.graph.out_edges(n, data=True)]
        in_edges = [(u, v, d, "in") for u, v, d in self.graph.in_edges(n, data=True)]

        for u, v, d, direction in out_edges + in_edges:
            to_jimage = d["to_jimage"]

            if direction == "in":
                u, v = v, u  # noqa: PLW2901
                to_jimage = np.multiply(-1, to_jimage)

            to_jimage = tuple(map(int, np.add(to_jimage, jimage)))

            if (v, to_jimage) not in connected_site_images:
                connected_site = ConnectedSiteLight(
                    site=self.structure[v],
                    jimage=to_jimage,
                    index=v,
                    weight=None,
                    dist=None,
                )

                connected_sites.add(connected_site)
                connected_site_images.add((v, to_jimage))

        return list(connected_sites)

    def set_wyckoffs(self, symmetry_tol: float = 0.01) -> None:
        """Set Wyckoff position labels for each site.

        Labels each site with its element, Wyckoff letter, and space group
        number in the format ``"{element}_{wyckoff}_{spacegroup}"``.

        Parameters
        ----------
        symmetry_tol : float, default 0.01
            Tolerance for symmetry detection in Angstroms.

        Notes
        -----
        If symmetry detection fails, falls back to elemental labels.

        """
        siteless_strc = self.structure.copy()

        for site_i in range(len(self.structure)):
            siteless_strc.replace(site_i, Element("H"))

        sga = SpacegroupAnalyzer(siteless_strc, symprec=symmetry_tol)
        sym_dataset = sga.get_symmetry_dataset()

        if sym_dataset is None:
            self.set_elemental_labels()
            return

        wyckoffs = sym_dataset.wyckoffs
        number = sym_dataset.number

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
        diameter_factor: int = 2,
        use_previous_cs: bool = False,
    ) -> None:
        """Compute and set compositional sequences as node attributes.

        This is the core method that computes the local environment
        fingerprint for each site by traversing the graph and counting
        neighbors at each depth.

        Parameters
        ----------
        hash_cs : bool, default False
            If True, hash the compositional sequence incrementally
            during computation for memory efficiency.
        wyckoff : bool, default False
            If True, use Wyckoff labels in the computation.
        additional_depth : int, default 0
            Extra traversal depth to add.
        diameter_factor : int, default 2
            Multiplier for graph diameter to determine traversal depth.
        use_previous_cs : bool, default False
            If True, use previous compositional sequence as starting labels.

        Notes
        -----
        After calling this method:

        - Node attributes are set with key ``"compositional_sequence"``
        - ``self.cc_cs`` contains compositional sequences per component

        """
        node_attributes = {}
        self.cc_cs = []
        get_connected_sites_light = functools.lru_cache(maxsize=None)(self.get_connected_sites_light)

        ug = self.graph.to_undirected()

        for cc in nx.connected_components(ug):
            cs_list = []

            d = diameter(ug.subgraph(cc))

            for focused_site_i in cc:
                depth = diameter_factor * d + additional_depth

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
        """Find all loops/rings starting from a given atom.

        Traverses the graph to find closed loops that start and end at the
        specified atom index.

        Parameters
        ----------
        depth : int
            Maximum loop size to search for.
        index : int
            The starting atom index.
        shortest : bool, default True
            If True, stop searching when all theoretically possible shortest
            loops are found.

        Returns
        -------
        list of list of tuple
            A list of loops, where each loop is a list of ``(index, image)``
            tuples representing the path.

        Notes
        -----
        Loops are found by breadth-first traversal and tracking when paths
        return to their starting point.

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

    def set_loops(self, diameter_factor: int, additional_depth: int) -> None:
        """Set loop-based labels for each site.

        Computes all loops for each site and creates a hashed label
        representing the ring topology around that site.

        Parameters
        ----------
        diameter_factor : int
            Multiplier for graph diameter to determine search depth.
        additional_depth : int
            Extra depth to add to the search.

        Notes
        -----
        Sets ``self.starting_labels`` to hashed loop representations.
        Used when ``loop=True`` in GraphIDGenerator.

        """
        self.starting_labels = []

        undirected_graph = self.graph.to_undirected()

        max_diameter = 0
        for cc in nx.connected_components(undirected_graph):
            d = diameter(undirected_graph.subgraph(cc))
            max_diameter = max(max_diameter, d)

        depth = max_diameter * diameter_factor + additional_depth

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
        diameter_factor: int = 2,
        use_previous_cs: bool = False,
    ) -> None:
        """Compute compositional sequence for a single site.

        Similar to set_compositional_sequence_node_attr, but only computes
        the sequence for the specified site. Used by DistanceClusteringGraphID.

        Parameters
        ----------
        n : int
            The site index to compute the sequence for.
        hash_cs : bool, default False
            If True, hash the sequence incrementally.
        wyckoff : bool, default False
            If True, use Wyckoff labels.
        additional_depth : int, default 0
            Extra traversal depth.
        diameter_factor : int, default 2
            Multiplier for graph diameter.
        use_previous_cs : bool, default False
            If True, use previous sequence as starting labels.

        """
        node_attributes = {}
        self.cc_cs = []
        get_connected_sites_light = functools.lru_cache(maxsize=None)(self.get_connected_sites_light)

        ug = self.graph.to_undirected()

        for cc in nx.connected_components(ug):
            cs_list = []

            d = diameter(ug.subgraph(cc))

            if n in cc:
                depth = diameter_factor * d + additional_depth

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
