import glob
import os
import unittest

from graph_id.core.long_distance_graph_id import LongDistanceGraphID
from pymatgen.core import Structure

from .imports import graph_id_cpp


test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        # if s.num_sites <= max_sites:
        if s.num_sites <= max_sites and name != "VSbO4":
            res.append((name, s))
    return res


class TestLongDistanceGraphIDGenerator(unittest.TestCase):
    # @pytest.mark.limit_leaks("1 MB")
    def test_get_long_distance_id(self):
        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)
        for name, s in small_test_structure():
            with self.subTest(name):
                print(name)

                aid = a.get_id(s)
                bid = b.get_long_distance_id(s)

                self.assertEqual(aid, bid)
    
    def test_carbon_allotrope(self):
        """
        グラファイトとダイアモンドは次元を考慮することで区別していた
        """
        graphite = Structure.from_file(f"{test_file_dir}/graphite.cif")
        diamond = Structure.from_file(f"{test_file_dir}/diamond.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        graphite_aid = a.get_id(graphite)
        graphite_bid = b.get_long_distance_id(graphite)

        diamond_aid = a.get_id(diamond)
        diamond_bid = b.get_long_distance_id(diamond)

        self.assertEqual(graphite_aid, graphite_bid)
        self.assertEqual(diamond_aid, diamond_bid)

        self.assertEqual(graphite_aid, "eb8e457303fed0a236ae6814ed491acb")
        self.assertEqual(diamond_aid, "628ceb4fa82fce02232d846c0d8960e6")

if __name__ == "__main__":
    unittest.main()
