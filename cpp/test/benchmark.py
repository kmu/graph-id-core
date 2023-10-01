import timeit
import unittest

import glob
import os.path

from pymatgen.core import Structure

from graph_id import GraphIDGenerator
from graph_id.analysis.graphs import StructureGraph
from imports import graph_id_cpp
from pymatgen.analysis.local_env import MinimumDistanceNN

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


class TestBenchmark(unittest.TestCase):
    def test_minimum_distance(self):
        print("MinimumDistanceNN:")
        a = MinimumDistanceNN()
        b = graph_id_cpp.MinimumDistanceNN()
        for name, s in small_test_structure():
            try:
                at = timeit.timeit("a.get_all_nn_info(s)", number=10, globals=locals()) * 100
                bt = timeit.timeit("b.get_all_nn_info(s)", number=10, globals=locals()) * 100
                print("{: 3d} site. Python: {: 8.3f}ms, C++: {: 7.3f}ms, {: 4.1f} times faster [{}]".format(s.num_sites, at, bt, at / bt, name))
            except Exception as e:
                print(e)

    def test_structure_graph(self):
        print("StructureGraph.set_compositional_sequence_node_attr:")
        py = StructureGraph
        cpp = graph_id_cpp.StructureGraph
        for name, s in small_test_structure():
            nn = graph_id_cpp.MinimumDistanceNN()
            def f(cls):
                sg = cls.with_local_env_strategy(s, nn)
                sg.set_elemental_labels()
                sg.set_compositional_sequence_node_attr(hash_cs=True)
                sg.set_compositional_sequence_node_attr(use_previous_cs=True, hash_cs=True)
            try:
                at = timeit.timeit("f(py)", number=10, globals=locals()) * 100
                bt = timeit.timeit("f(cpp)", number=10, globals=locals()) * 100
                print("{: 3d} site. Python: {: 8.3f}ms, C++: {: 7.3f}ms, {: 4.1f} times faster [{}]".format(s.num_sites, at, bt, at / bt, name))
            except Exception as e:
                print(e)

    def test_graph_id(self):
        print("GraphIDGenerator.get_id:")
        a = GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure():
            try:
                N = 1
                at = timeit.timeit("a.get_id(s)", number=N, globals=locals()) * 1000 / N
                bt = timeit.timeit("b.get_id(s)", number=N, globals=locals()) * 1000 / N
                print("{: 3d} site. Python: {: 8.3f}ms, C++: {: 7.3f}ms, {: 4.1f} times faster [{}]".format(s.num_sites, at, bt, at / bt, name))
            except Exception as e:
                print(e)

    def test_graph_id_cpp_only(self):
        print("GraphIDGenerator.get_id:")
        import numpy as np
        g = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure(1000):
            try:
                N = 1000
                if s.num_sites > 20: N = 10
                t = timeit.repeat("g.get_id(s)", number=N, repeat=5, globals=locals())
                mean = np.mean(t) * 1000 / N
                std = np.std(t) * 1000 / N
                print("{: 3d} site. C++: {: 8.3f}ms+={:.1f}% [{}]".format(s.num_sites, mean, std / mean * 100, name))
            except Exception as e:
                print(e)


if __name__ == "__main__":
    unittest.main()
