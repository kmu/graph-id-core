import unittest
from pathlib import Path

import networkx as nx
from graph_id_py.analysis.graphs import StructureGraph
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
from pymatgen.core import Lattice, Structure

import graph_id

test_file_dir = Path(__file__).parent.parent.parent / "tests" / "py" / "test_files"


def small_test_structure(max_sites=30):
    res = []
    for p in test_file_dir.glob("*.cif"):
        name = p.stem.replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


class TestStructureGraph(unittest.TestCase):
    def test_graph(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                py = list(sg_py.graph.to_undirected().edges())
                py = sorted(py + [(b, a) for a, b in py])
                cpp = sg_cpp.get_connected_site_index()
                self.assertSetEqual(set(cpp), {(b, a) for a, b in cpp})
                self.assertSetEqual(set(py), set(cpp))
                self.assertListEqual(py, cpp)

    def test_graph_diameter(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                ug = sg_py.graph.to_undirected()
                diameter = [
                    nx.diameter(ug.subgraph(cc)) for cc in nx.connected_components(ug)
                ]
                self.assertListEqual(diameter, sg_cpp.cc_diameter)

    def test_graph_diameter_not_strongly_connected(self):
        for name, s in small_test_structure(1000):
            with self.subTest(name):
                s = Structure(
                    Lattice([[10, 0, 0], [0, 10, 0], [0, 0, 10]]),
                    ["H"] * 4,
                    [[0, 0, 0], [0, 0.001, 0], [0.5, 0, 0], [0.5, 0.001, 0]],
                )
                sg_py = StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                ug = sg_py.graph.to_undirected()
                diameter = [
                    nx.diameter(ug.subgraph(cc)) for cc in nx.connected_components(ug)
                ]
                self.assertEqual(diameter, sg_cpp.cc_diameter)

    def test_set_elemental_labels(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                sg_py = StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                sg_py.set_elemental_labels()
                sg_cpp.set_elemental_labels()
                self.assertListEqual(sg_py.starting_labels, sg_cpp.labels)

    def test_set_wyckoff_labels(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                sg_py = StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(
                    s, graph_id.MinimumDistanceNN()
                )
                sg_py.set_wyckoffs()
                sg_cpp.set_wyckoffs()
                self.assertListEqual(sg_py.starting_labels, sg_cpp.labels)

    def test_set_compositional_sequence_node_attr(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                sg_py.set_elemental_labels()
                sg_cpp.set_elemental_labels()
                sg_py.set_compositional_sequence_node_attr()
                sg_cpp.set_compositional_sequence_node_attr()
                self.assertListEqual(sg_py.cc_cs, sg_cpp.cc_cs)

    def test_set_compositional_sequence_node_attr_hash(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                sg_py.set_elemental_labels()
                sg_cpp.set_elemental_labels()
                sg_py.set_compositional_sequence_node_attr(hash_cs=True)
                sg_cpp.set_compositional_sequence_node_attr(hash_cs=True)
                self.assertListEqual(sg_py.cc_cs, sg_cpp.cc_cs)

    def test_set_compositional_sequence_node_attr_twice(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                sg_py.set_elemental_labels()
                sg_cpp.set_elemental_labels()
                sg_py.set_compositional_sequence_node_attr()
                sg_cpp.set_compositional_sequence_node_attr()
                sg_py.set_compositional_sequence_node_attr(use_previous_cs=True)
                sg_cpp.set_compositional_sequence_node_attr(use_previous_cs=True)
                self.assertListEqual(sg_py.cc_cs, sg_cpp.cc_cs)

    def test_get_dimensionality_larsen(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                self.assertEqual(
                    get_dimensionality_larsen(sg_py), sg_cpp.get_dimensionality_larsen()
                )

    def test_get_dimensionality_larsen_corner_case(self):
        s = Structure.from_file(
            test_file_dir / "structure_for_dimensionality_larsen.cif"
        )
        nn = graph_id.CutOffDictNN({("Si", "O"): 2.0})
        sg_py = StructureGraph.with_local_env_strategy(s, nn)
        sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
        self.assertEqual(
            get_dimensionality_larsen(sg_py), sg_cpp.get_dimensionality_larsen()
        )

    def test_from_py(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                nn = graph_id.MinimumDistanceNN()
                sg_py = StructureGraph.with_local_env_strategy(s, nn)
                sg_cpp = graph_id.StructureGraph.with_local_env_strategy(s, nn)
                self.assertEqual(sg_py, sg_cpp.to_py())
                sg_cpp_from_py = graph_id.StructureGraph.from_py(sg_py)
                n_conn_cpp = len(sg_cpp.get_connected_site_index())
                n_conn_cpp_from_py = len(sg_cpp_from_py.get_connected_site_index())
                self.assertEqual(n_conn_cpp, n_conn_cpp_from_py)

    def test_from_py_zeolite(self):
        s = Structure.from_file(test_file_dir / "ABW.cif")
        sg_py = StructureGraph.with_local_env_strategy(
            s, graph_id.CutOffDictNN({("Si", "O"): 2})
        )
        sg_cpp = graph_id.StructureGraph.with_local_env_strategy(
            s, graph_id.CutOffDictNN({("Si", "O"): 2})
        )
        sg_cpp_from_py = graph_id.StructureGraph.from_py(sg_py)
        n_conn_cpp = len(sg_cpp.get_connected_site_index())
        n_conn_cpp_from_py = len(sg_cpp_from_py.get_connected_site_index())
        self.assertEqual(n_conn_cpp, n_conn_cpp_from_py)


if __name__ == "__main__":
    unittest.main()
