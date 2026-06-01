from __future__ import annotations

from collections import Counter
from hashlib import blake2b
from typing import TYPE_CHECKING

from pymatgen.util.string import formula_double_format

if TYPE_CHECKING:
    from pymatgen.core.structure import Neighbor


def blake(s):
    """Hash a string using BLAKE2b."""
    return blake2b(s.encode()).hexdigest()


class CompositionalSequence:

    """Compute the compositional sequence for a site in a structure graph.

    A compositional sequence is a fingerprint of the local chemical environment
    around an atom, computed by traversing the graph in shells and counting
    the elements encountered at each depth.

    For example, for Na in NaCl rock salt structure:

    - Depth 0: Na (the central atom)
    - Depth 1: Cl6 (6 nearest Cl neighbors)
    - Depth 2: Na12 (12 next-nearest Na neighbors)

    The sequence ``"Na-Cl6-Na12-..."`` uniquely identifies the local environment.

    Parameters
    ----------
    focused_site_i : int
        The index of the central atom.
    starting_labels : list of str
        Labels for each site in the structure.
    hash_cs : bool, default False
        If True, hash the sequence incrementally to save memory.
    use_previous_cs : bool, default False
        If True, use previous compositional sequences as labels
        (for iterative refinement).

    Attributes
    ----------
    focused_site_i : int
        The central site index.
    first_element : str
        The label of the central site.
    compositional_seq : list of str
        The composition at each depth (if hash_cs=False).
    cs_for_hashing : str
        The incrementally hashed sequence (if hash_cs=True).

    Examples
    --------
    >>> cs = CompositionalSequence(0, ["Na", "Cl", "Na", "Cl"])
    >>> # ... add neighbors at each depth ...
    >>> print(str(cs))
    'Na-Cl6-Na12...'

    """

    def __init__(self, focused_site_i, starting_labels, hash_cs=False, use_previous_cs=False):
        """Initialize the compositional sequence computation."""
        self.hash_cs = hash_cs
        if hash_cs:
            self.cs_for_hashing = ""
        else:
            self.compositional_seq = []

        self.focused_site_i = focused_site_i
        self.new_sites = [(focused_site_i, (0, 0, 0))]

        self.seen_sites = set(self.new_sites)
        self.use_previous_cs = use_previous_cs
        self.labels = starting_labels
        self.composition_counter: Counter = Counter()
        self.first_element = starting_labels[focused_site_i]

    def __str__(self):
        """Return the string representation of the compositional sequence.

        Returns
        -------
        str
            Format: ``"{first_element}-{depth1}-{depth2}-..."``

        """
        if self.hash_cs:
            return f"{self.first_element}-{self.cs_for_hashing}"  # type: ignore

        return f"{self.first_element}-{'-'.join(self.compositional_seq)}"  # type: ignore

    def get_current_starting_sites(self):
        """Get the sites to expand from for the next depth.

        Returns
        -------
        list of tuple
            List of ``(site_index, jimage)`` tuples for the frontier sites.

        """
        new_sites = self.new_sites
        self.new_sites = []
        return [*new_sites]

    def count_composition_for_neighbors(
        self,
        nsites: list[Neighbor],
    ) -> None:
        """Count the composition of neighboring sites.

        Adds new neighbors to the frontier and counts their labels
        for the current depth.

        Parameters
        ----------
        nsites : list of Neighbor
            The neighboring sites to count.

        """
        for neighbor in nsites:
            neighbor_info = (neighbor.index, neighbor.jimage)

            if neighbor_info not in self.seen_sites:
                self.seen_sites.add(neighbor_info)

                self.new_sites.append(neighbor_info)

                if self.use_previous_cs:
                    cs = self.labels[neighbor.index]
                    self.composition_counter[cs] += 1
                else:
                    self.composition_counter[self.labels[neighbor.index]] += 1

    def finalize_this_depth(self):
        """Finalize counting for the current depth.

        Converts the composition counter to a formula string and
        either appends it to the sequence or hashes it incrementally.
        Resets the counter for the next depth.
        """
        formula = self.get_sorted_composition_list_from(self.composition_counter)

        if self.hash_cs:
            self.cs_for_hashing = blake(f"{self.cs_for_hashing}-{''.join(formula)}")
        else:
            self.compositional_seq.append("".join(formula))

    def get_sorted_composition_list_from(self, composition_counter: Counter) -> list[str]:
        """Convert a composition counter to a sorted formula list.

        Parameters
        ----------
        composition_counter : Counter
            Counts of each element/label.

        Returns
        -------
        list of str
            Sorted list of ``"{element}{count}"`` strings.

        """
        sorted_symbols = sorted(composition_counter.keys())
        return [s + str(formula_double_format(composition_counter[s], ignore_ones=False)) for s in sorted_symbols]
