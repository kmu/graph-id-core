import os
import unittest

from pymatgen.core import Structure

from graph_id.analysis.graphs import StructureGraph
from imports import graph_id_cpp

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


class TestStructureGraph(unittest.TestCase):
    def test_set_elemental_labels(self):
        s = Structure.from_file(test_file_dir + "/mp-1299593.cif")
        sg_py = StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        sg_cpp = graph_id_cpp.StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        sg_py.set_elemental_labels()
        sg_cpp.set_elemental_labels()
        self.assertListEqual(sg_py.starting_labels, sg_cpp.labels)

    def test_set_wyckoff_labels(self):
        s = Structure.from_file(test_file_dir + "/mp-1299593.cif")
        sg_py = StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        sg_cpp = graph_id_cpp.StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        sg_py.set_wyckoffs()
        sg_cpp.set_wyckoffs()
        self.assertListEqual(sg_py.starting_labels, sg_cpp.labels)


if __name__ == '__main__':
    unittest.main()
