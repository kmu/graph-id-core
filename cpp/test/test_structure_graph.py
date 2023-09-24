import glob
import os
import unittest

from pymatgen.core import Structure, Lattice

from graph_id.analysis.graphs import StructureGraph
from imports import graph_id_cpp
import networkx as nx

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure():
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= 30:
            res.append((name, s))
    return res


class TestStructureGraph(unittest.TestCase):
    def test_graph(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id_cpp.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id_cpp.StructureGraph.with_local_env_strategy(s, nn)
                py = list(sg_py.graph.to_undirected().edges())
                py = sorted(py + [(b, a) for a, b in py])
                cpp = sg_cpp.get_connected_site_index()
                self.assertSetEqual(set(cpp), {(b, a) for a, b in cpp})
                self.assertSetEqual(set(py), set(cpp))
                self.assertListEqual(py, cpp)

    def test_graph_diameter(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id_cpp.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id_cpp.StructureGraph.with_local_env_strategy(s, nn)
                ug = sg_py.graph.to_undirected()
                diameter = max([nx.diameter(ug.subgraph(cc)) for cc in nx.connected_components(ug)])
                self.assertEqual(diameter, sg_cpp.graph_diameter)

    def test_graph_diameter_not_strongly_connected(self):
        s = Structure(Lattice([[10, 0, 0], [0, 10, 0], [0, 0, 10]]), ["H"]*4, [[0, 0, 0], [0, 0.001, 0], [0.5, 0, 0], [0.5, 0.001, 0]])
        sg_py = StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        sg_cpp = graph_id_cpp.StructureGraph.with_local_env_strategy(s, graph_id_cpp.MinimumDistanceNN())
        ug = sg_py.graph.to_undirected()
        diameter = max([nx.diameter(ug.subgraph(cc)) for cc in nx.connected_components(ug)])
        self.assertEqual(diameter, sg_cpp.graph_diameter)


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
