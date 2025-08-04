import os
from unittest import TestCase

from graph_id_py.core.distance_clustering_graph_id import DistanceClusteringGraphID
from graph_id_py.analysis.local_env import DistanceClusteringNN
from pymatgen.core import Structure

TEST_FILES = os.path.dirname(os.path.abspath(__file__)) + "/test_files"


class TestDistanceClusteringGraphID(TestCase):
    def test_small(self):
        """
        Sc単体構造についてのテスト
        """
        s = Structure.from_file(f"{TEST_FILES}/mp-36.cif")

        ldgid = DistanceClusteringGraphID(nn=DistanceClusteringNN(), rank_k=3, cutoff=6.0)

        id_1 = ldgid.get_id(s)

        assert id_1 == "3b2bcf29df1ce648"
        
