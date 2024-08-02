import os
from unittest import TestCase

from graph_id.core.long_distance_graph_id import LongDistanceGraphID
from graph_id.analysis.local_env import LongDistanceNN
from pymatgen.core import Structure

TEST_FILES = os.path.dirname(os.path.abspath(__file__)) + "/test_files"


class TestLongDistanceGraphID(TestCase):
    def test_small(self):
        """
        Sc単体構造についてのテスト
        """
        s = Structure.from_file(f"{TEST_FILES}/mp-36.cif")

        ldgid = LongDistanceGraphID(nn=LongDistanceNN(), max_cluster_num=3, cutoff=6.0)

        id_1 = ldgid.get_id(s)

        assert id_1 == "52c959e1543449ed5d47fcbfe0bec2ff"
        
