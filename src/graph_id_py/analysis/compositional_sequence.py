from collections import Counter
from hashlib import blake2b
from typing import List

from pymatgen.core.structure import Neighbor
from pymatgen.util.string import formula_double_format


def blake(s):
    return blake2b(s.encode()).hexdigest()


class CompositionalSequence:
    def __init__(self, focused_site_i, starting_labels, hash_cs=False, use_previous_cs=False):
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
        if self.hash_cs:
            return f"{self.first_element}-{self.cs_for_hashing}"  # type: ignore

        else:
            return f"{self.first_element}-{'-'.join(self.compositional_seq)}"  # type: ignore

    def get_current_starting_sites(self):
        new_sites = self.new_sites
        self.new_sites = []
        return [*new_sites]

    def count_composition_for_neighbors(
        self,
        nsites: List[Neighbor],
        # graph: nx.Graph,
        # labels: List[str],
    ) -> None:

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
        formula = self.get_sorted_composition_list_from(self.composition_counter)

        if self.hash_cs:
            self.cs_for_hashing = blake(f"{self.cs_for_hashing}-{''.join(formula)}")
        else:
            self.compositional_seq.append("".join(formula))

    def get_sorted_composition_list_from(self, composition_counter: Counter) -> List[str]:
        sorted_symbols = sorted(composition_counter.keys())
        formula = [s + str(formula_double_format(composition_counter[s], False)) for s in sorted_symbols]
        return formula
