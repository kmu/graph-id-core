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

    def test_one_site_structure(self):
        """
        サイトが1つだけの構造
        結合ができていないと全て同じIDになる
        """
        s1 = Structure.from_file(f"{test_file_dir}/mp-121.cif")
        s2 = Structure.from_file(f"{test_file_dir}/mp-611219.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s1_aid = a.get_id(s1)
        s1_bid = b.get_long_distance_id(s1)

        s2_aid = a.get_id(s2)
        s2_bid = b.get_long_distance_id(s2)

        self.assertEqual(s1_aid, s1_bid)
        self.assertEqual(s2_aid, s2_bid)

        self.assertEqual(s1_aid, "4c1d6581dacf78dfc8d3da5bdc07b0d5")
        self.assertEqual(s2_aid, "d701ce2ac1a778013c49e0d4d1fe8645")

    def test_empty_compositional_sequence(self):
        """
        空のCompositionalSequenceを持つサイトができる構造
        空のCompositionalSequenceをハイフンでつながない処理にしている
        """
        s_298 = Structure.from_file(f"{test_file_dir}/298 K.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_298_aid = a.get_id(s_298)
        s_298_bid = b.get_long_distance_id(s_298)

        self.assertEqual(s_298_aid, s_298_bid)
        self.assertEqual(s_298_aid, "18a8dce917bce31f14aa96352ec2d24b")

    def test_sort_jstrs(self):
        """
        各サイトのCompositionalSequenceをソートしていないと
        IDが変わってしまう構造
        """
        s_1078 = Structure.from_file(f"{test_file_dir}/1078 K.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_1078_aid = a.get_id(s_1078)
        s_1078_bid = b.get_long_distance_id(s_1078)

        self.assertEqual(s_1078_aid, s_1078_bid)
        self.assertEqual(s_1078_aid, "9c59aafb0dec3cae9a54db39544cf762")

    def test_break_edge(self):
        """
        break_edgeする際に逆向きの辺も削除しないとIDが変わってしまう構造
        """
        s_ca = Structure.from_file(f"{test_file_dir}/mp-1067285.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_ca_aid = a.get_id(s_ca)
        s_ca_bid = b.get_long_distance_id(s_ca)

        self.assertEqual(s_ca_aid, s_ca_bid)
        self.assertEqual(s_ca_aid, "6c022438a0f73d7122716edb769f2959")

    def test_no_supercell_structure(self):
        """
        スーパーセルを作るだけだと6.0angまでの結合を
        全ては見つけられない構造
        """

        s_si = Structure.from_file(f"{test_file_dir}/mp-1056579.cif")

        a = LongDistanceGraphID(rank_k=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)

        s_si_aid = a.get_id(s_si)
        s_si_bid = b.get_long_distance_id(s_si)

        self.assertEqual(s_si_aid, s_si_bid)
        self.assertEqual(s_si_aid, "e140bcc85363351dba38629faacf2257")


if __name__ == "__main__":
    unittest.main()
