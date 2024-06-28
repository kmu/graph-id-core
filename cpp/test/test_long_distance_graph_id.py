import glob
import os
import unittest
import pytest

import graph_id
from graph_id.core.long_distance_graph_id import LongDistanceGraphID
from pymatgen.core import Structure

from .imports import graph_id_cpp

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        # if s.num_sites <= max_sites and "Fe" in name:
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


class TestLongDistanceGraphIDGenerator(unittest.TestCase):
    # @pytest.mark.limit_leaks("1 MB")
    def test_get_long_distance_id(self):
        a = LongDistanceGraphID(max_cluster_num=3, cutoff=6.0)
        b = graph_id_cpp.GraphIDGenerator(rank_k=3, cutoff=6.0)
        for name, s in small_test_structure():
            with self.subTest(name):
                print(name)

                aid = a.get_id(s)
                bid = b.get_long_distance_id(s)

                self.assertEqual(aid, bid)

if __name__ == "__main__":
    unittest.main()
