from pathlib import Path
from unittest import TestCase

from pymatgen.core import Structure

from graph_id.analysis.local_env import DistanceClusteringNN
from graph_id.core.distance_clustering_graph_id import DistanceClusteringGraphID

TEST_FILES = str(Path(__file__).resolve().parent / "test_files")


class TestDistanceClusteringGraphID(TestCase):
    def test_small(self):
        """
        Single site structures test
        """
        s = Structure.from_file(f"{TEST_FILES}/mp-36.cif")

        ldgid = DistanceClusteringGraphID(nn=DistanceClusteringNN(), rank_k=3, cutoff=6.0)

        id_1 = ldgid.get_id(s)

        assert id_1 == "3b2bcf29df1ce648"
