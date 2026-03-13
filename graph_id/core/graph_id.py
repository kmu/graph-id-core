from __future__ import annotations

import multiprocessing as multi
from hashlib import blake2b
from multiprocessing import Pool

import networkx as nx
import numpy as np
from ase import Atoms
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.core import Element, Lattice, Molecule, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from graph_id.analysis.graphs import StructureGraph

__version__ = "0.1.0"


def blake(s):
    """Hash a string using BLAKE2b."""
    return blake2b(s.encode()).hexdigest()


class GraphIDGenerator:

    """Core Python implementation of Graph ID generation.

    GraphIDGenerator converts atomic structures into unique, deterministic identifiers
    by analyzing the topological and compositional properties of the structure graph.

    The algorithm:

    1. Constructs a graph where atoms are nodes and bonds are edges
    2. Computes compositional sequences for each atom (local environment fingerprints)
    3. Iteratively refines sequences until convergence
    4. Hashes the sequences to produce the final ID

    Examples
    --------
    >>> from pymatgen.core import Structure
    >>> from graph_id.core.graph_id import GraphIDGenerator
    >>> structure = Structure.from_file("NaCl.cif")
    >>> gen = GraphIDGenerator()
    >>> gen.get_id(structure)
    'NaCl-3D-88c8e156db1b0fd9'

    See Also
    --------
    GraphIDMaker : High-level interface with simpler API
    DistanceClusteringGraphID : Variant using distance clustering

    """

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
        prepend_composition=True,
        prepend_dimensionality=True,
    ):
        """Initialize the GraphIDGenerator.

        Parameters
        ----------
        nn : NearNeighbors, optional
            A neighbor-finding strategy from pymatgen.analysis.local_env.
            If None, defaults to MinimumDistanceNN().
        wyckoff : bool, default False
            If True, include Wyckoff position information in the ID.
            Cannot be used together with ``loop=True``.
        diameter_factor : int, default 2
            Multiplier for the graph diameter to determine traversal depth.
            The total depth is ``diameter_factor * diameter + additional_depth``.
        additional_depth : int, default 1
            Extra depth added to the calculated traversal depth.
        symmetry_tol : float, default 0.1
            Tolerance for symmetry operations when detecting Wyckoff positions.
            Only used when ``wyckoff=True``.
        topology_only : bool, default False
            If True, generate topology-only IDs that ignore element types.
            Useful for finding isostructural materials.
            Cannot be used together with ``loop=True``.
        loop : bool, default False
            If True, use loop/ring-based identification algorithm.
            Cannot be used together with ``wyckoff=True`` or ``topology_only=True``.
        digest_size : int, default 8
            Size of the BLAKE2b hash digest in bytes.
            The output will be ``2 * digest_size`` hexadecimal characters.
        prepend_composition : bool, default True
            If True, prepend the reduced chemical formula to the ID.
        prepend_dimensionality : bool, default True
            If True, prepend the dimensionality (0D, 1D, 2D, 3D) to the ID.

        Raises
        ------
        ValueError
            If incompatible options are specified:

            - ``wyckoff=True`` and ``loop=True``
            - ``loop=True`` and ``topology_only=True``

        Examples
        --------
        >>> gen = GraphIDGenerator()  # Default settings
        >>> gen = GraphIDGenerator(topology_only=True)  # Topology-only
        >>> gen = GraphIDGenerator(wyckoff=True, symmetry_tol=0.01)  # With Wyckoff
        >>> gen = GraphIDGenerator(diameter_factor=3, additional_depth=2)  # Deeper traversal

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
        self.diameter_factor = diameter_factor
        self.symmetry_tol = symmetry_tol
        self.topology_only = topology_only
        self.loop = loop
        self.digest_size = digest_size
        self.prepend_composition = prepend_composition
        self.prepend_dimensionality = prepend_dimensionality

    def _join_cs_list(self, cs_list):
        """Join and hash a list of compositional sequences."""
        return blake("-".join(sorted(cs_list)))

    def _component_strings_to_whole_id(self, component_strings):
        """Combine component hashes into a single ID."""
        long_str = ":".join(np.sort(component_strings))
        return blake2b(long_str.encode("ascii"), digest_size=self.digest_size).hexdigest()

    def get_id(self, structure):
        """Generate a Graph ID for the given structure.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object representing the crystal or molecule.

        Returns
        -------
        str
            The Graph ID. Format depends on configuration:

            - Default: ``"{formula}-{dim}D-{hash}"``
            - ``prepend_composition=False``: ``"{dim}D-{hash}"``
            - ``prepend_dimensionality=False``: ``"{formula}-{hash}"``
            - Both False: ``"{hash}"``

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> gen.get_id(nacl_structure)
        'NaCl-3D-88c8e156db1b0fd9'

        >>> gen = GraphIDGenerator(prepend_composition=False)
        >>> gen.get_id(nacl_structure)
        '3D-88c8e156db1b0fd9'

        """
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
        """Add composition and dimensionality prefixes to a Graph ID.

        Parameters
        ----------
        sg : StructureGraph
            The prepared structure graph.
        gid : str
            The base Graph ID hash.

        Returns
        -------
        str
            The elaborated Graph ID with prefixes.

        """
        if self.prepend_dimensionality:
            dim = get_dimensionality_larsen(sg)
            gid = f"{dim}D-{gid}"

        if self.prepend_composition and not self.topology_only:
            gid = f"{sg.structure.composition.reduced_formula}-{gid}"

        return gid

    @property
    def version(self):
        """Str : The version of the GraphIDGenerator."""
        return __version__

    def get_id_catch_error(self, structure):
        """Generate a Graph ID with error handling.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.

        Returns
        -------
        str
            The Graph ID, or an empty string if an error occurs.

        Notes
        -----
        This method catches all exceptions silently and returns an empty string.
        Useful for batch processing where some structures may fail.

        """
        try:
            return self.get_id(structure)
        except Exception:  # noqa: BLE001
            return ""

    def get_many_ids(self, structures, parallel=False):
        """Generate Graph IDs for multiple structures.

        Parameters
        ----------
        structures : list of Structure
            A list of pymatgen Structure objects.
        parallel : bool, default False
            If True, use parallel processing with all available CPU cores.
            Shows a progress bar via tqdm.

        Returns
        -------
        list of str
            A list of Graph IDs corresponding to each input structure.
            Failed structures will have empty string IDs.

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> structures = [Structure.from_file(f) for f in cif_files]
        >>> ids = gen.get_many_ids(structures, parallel=True)

        """
        if parallel:
            n_cores = multi.cpu_count()

            p = Pool(n_cores)
            imap = p.imap(self.get_id_catch_error, structures)

            return list(tqdm(imap, total=len(structures)))

        return [self.get_id(s) for s in structures]

    def get_component_ids(self, structure):
        """Get Graph IDs for each connected component in the structure.

        For structures with multiple disconnected fragments (e.g., molecular
        crystals), this returns a separate ID for each component.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.

        Returns
        -------
        numpy.ndarray
            Array of dictionaries, each containing:

            - ``site_i``: Set of site indices in this component
            - ``graph_id``: The Graph ID for this component

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> components = gen.get_component_ids(molecular_crystal)
        >>> for comp in components:
        ...     print(f"Sites {comp['site_i']}: {comp['graph_id']}")

        """
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

    def _molecule_to_structure(self, mol: Molecule, vacuum: float = 10.0) -> Structure:
        coords = mol.cart_coords
        species = mol.species

        min_c = coords.min(axis=0)
        max_c = coords.max(axis=0)
        lengths = max_c - min_c + 2 * vacuum

        lattice = Lattice.from_parameters(lengths[0], lengths[1], lengths[2], 90, 90, 90)

        shifted_coords = coords - min_c + vacuum

        return Structure(lattice, species, shifted_coords, coords_are_cartesian=True)

    def get_merged_id(self, materials_list: list[Structure, Molecule, Atoms]):
        """Generate a merged Graph ID for multiple materials.

        This method computes Graph IDs for all connected components found in
        the input materials and merges them into a single identifier. The
        connected components extracted from each material are converted into
        canonical strings, sorted, concatenated with ``:`` separators, and
        hashed using BLAKE2b to produce the final merged Graph ID.

        The input materials may be provided as pymatgen ``Structure`` objects,
        pymatgen ``Molecule`` objects, or ASE ``Atoms`` objects. Molecules and
        ASE atoms are internally converted to ``Structure`` objects before
        graph analysis.

        Parameters
        ----------
        materials_list : list of Structure or Molecule or Atoms
            A list of materials to be included in the merged Graph ID. Each
            item must be one of the following types:

            - ``pymatgen.core.Structure``
            - ``pymatgen.core.Molecule``
            - ``ase.Atoms``

        Returns
        -------
        str
            A hexadecimal string representing the merged Graph ID of all
            connected components found in the input materials.

        Raises
        ------
        TypeError
            If an item in ``materials_list`` is not a supported material type.

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> gid = gen.get_merged_id([structure1, structure2])
        >>> print(gid)

        The method can also accept mixed object types:

        >>> gid = gen.get_merged_id([structure, molecule, ase_atoms])
        >>> print(gid)

        """
        array_list = []
        for material in materials_list:
            if isinstance(material, Structure):
                structure = material
            elif isinstance(material, Molecule):
                structure = self._molecule_to_structure(material)
            elif isinstance(material, Atoms):
                structure = AseAtomsAdaptor.get_structure(material)
            else:
                error_message = f"Item of materials_list must be pymatgen.core.Strucuture or pymatgen.core.Molecule \
                or ase.Atoms, got {type(material.__name__)}"
                raise TypeError(error_message)

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
            array_list.extend(array)

        long_str = ":".join(np.sort(array_list))

        return blake2b(long_str.encode("ascii"), digest_size=self.digest_size).hexdigest()

    def are_same(self, structure1, structure2):
        """Check if two structures have the same Graph ID.

        Parameters
        ----------
        structure1 : Structure
            The first pymatgen Structure object.
        structure2 : Structure
            The second pymatgen Structure object.

        Returns
        -------
        bool
            True if both structures have identical Graph IDs, False otherwise.

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> if gen.are_same(struct1, struct2):
        ...     print("Structures are topologically equivalent")

        """
        return self.get_id(structure1) == self.get_id(structure2)

    def prepare_structure_graph(self, structure):
        """Build and prepare the structure graph with compositional sequences.

        This method constructs a graph representation of the structure,
        computes compositional sequences for each site, and iteratively
        refines them until convergence.

        Parameters
        ----------
        structure : Structure
            A pymatgen Structure object.

        Returns
        -------
        StructureGraph
            The prepared structure graph with compositional sequence node
            attributes. The graph also has a ``cc_cs`` attribute containing
            the compositional sequences for each connected component.

        Notes
        -----
        This is primarily an internal method, but can be useful for
        advanced analysis of the structure graph.

        """
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

    def get_unique_structures(self, structures: list[Structure]) -> list[Structure]:
        """Filter a list of structures to keep only unique ones.

        Removes duplicate structures based on their Graph IDs. When duplicates
        are found, only the first occurrence is kept.

        Parameters
        ----------
        structures : list of Structure
            A list of pymatgen Structure objects, possibly containing duplicates.

        Returns
        -------
        list of Structure
            A list containing only unique structures (first occurrence of each).

        Examples
        --------
        >>> gen = GraphIDGenerator()
        >>> all_structures = load_many_cifs()
        >>> unique = gen.get_unique_structures(all_structures)
        >>> print(f"Reduced {len(all_structures)} to {len(unique)} unique")

        """
        unique_structures = []
        graph_ids = set()

        for strct in structures:
            new_graph_id = self.get_id(strct)
            if new_graph_id not in graph_ids:
                graph_ids.add(new_graph_id)
                unique_structures.append(strct)

        return unique_structures
