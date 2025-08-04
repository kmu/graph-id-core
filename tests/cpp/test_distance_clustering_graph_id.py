import glob
import os
import unittest

from graph_id_py.core.distance_clustering_graph_id import DistanceClusteringGraphID
from pymatgen.core import Structure

import graph_id


test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../py/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        # if s.num_sites <= max_sites:
        if s.num_sites <= max_sites and name != "VSbO4":
            res.append((name, s))
    return res


class TestDistanceClusteringGraphIDGenerator(unittest.TestCase):
    # @pytest.mark.limit_leaks("1 MB")
    def test_get_distance_clustering_id(self):
        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)
        for name, s in small_test_structure():
            with self.subTest(name):
                print(name)

                aid = a.get_id(s)
                bid = b.get_distance_clustering_id(s)

                self.assertEqual(aid, bid)

    def test_carbon_allotrope(self):
        """
        Assert that graphite and diamond have different IDs
        """
        graphite = Structure.from_file(f"{test_file_dir}/graphite.cif")
        diamond = Structure.from_file(f"{test_file_dir}/diamond.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        graphite_aid = a.get_id(graphite)
        graphite_bid = b.get_distance_clustering_id(graphite)

        diamond_aid = a.get_id(diamond)
        diamond_bid = b.get_distance_clustering_id(diamond)

        self.assertEqual(graphite_aid, graphite_bid)
        self.assertEqual(diamond_aid, diamond_bid)

        self.assertEqual(graphite_aid, "06f29c886936a611")
        self.assertEqual(diamond_aid, "0ff44c1761e1940e")

    def test_one_site_structure(self):
        """
        Structures with only one site
        If there is no bond, all structures have the same ID
        """
        s1 = Structure.from_file(f"{test_file_dir}/mp-121.cif")
        s2 = Structure.from_file(f"{test_file_dir}/mp-611219.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s1_aid = a.get_id(s1)
        s1_bid = b.get_distance_clustering_id(s1)

        s2_aid = a.get_id(s2)
        s2_bid = b.get_distance_clustering_id(s2)

        self.assertEqual(s1_aid, s1_bid)
        self.assertEqual(s2_aid, s2_bid)

        self.assertEqual(s1_aid, "5b39d303b4c6149a")
        self.assertEqual(s2_aid, "be82ee60bfd614a3")

    def test_empty_compositional_sequence(self):
        """
        Structures with an empty CompositionalSequence
        If there is no bond, all structures have the same ID
        """
        s_298 = Structure.from_file(f"{test_file_dir}/298 K.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_298_aid = a.get_id(s_298)
        s_298_bid = b.get_distance_clustering_id(s_298)

        self.assertEqual(s_298_aid, s_298_bid)
        self.assertEqual(s_298_aid, "3f8e7842377ffa36")

    def test_sort_jstrs(self):
        """
        If the CompositionalSequence of each site is not sorted,
        the ID will change.
        """
        s_1078 = Structure.from_file(f"{test_file_dir}/1078 K.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_1078_aid = a.get_id(s_1078)
        s_1078_bid = b.get_distance_clustering_id(s_1078)

        self.assertEqual(s_1078_aid, s_1078_bid)
        self.assertEqual(s_1078_aid, "39b5c4ae1aa4aa4c")

    def test_break_edge(self):
        """
        If the reverse edge is not deleted when break_edge is called,
        the ID will change.
        """
        s_ca = Structure.from_file(f"{test_file_dir}/mp-1067285.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_ca_aid = a.get_id(s_ca)
        s_ca_bid = b.get_distance_clustering_id(s_ca)

        self.assertEqual(s_ca_aid, s_ca_bid)
        self.assertEqual(s_ca_aid, "35d25fa388dfaed6")

    def test_no_supercell_structure(self):
        """
        If only a supercell is made,
        all bonds up to 6.0 angstroms are not found.
        """

        s_si = Structure.from_file(f"{test_file_dir}/mp-1056579.cif")

        a = DistanceClusteringGraphID(rank_k=3, cutoff=6.0)
        b = graph_id.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_si_aid = a.get_id(s_si)
        s_si_bid = b.get_distance_clustering_id(s_si)

        self.assertEqual(s_si_aid, s_si_bid)
        self.assertEqual(s_si_aid, "fc7064789538834f")


if __name__ == "__main__":
    unittest.main()
