import glob
import os.path
import timeit
import unittest

from pymatgen.analysis.local_env import (
    BrunnerNN_real,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    MinimumDistanceNN,
    MinimumOKeeffeNN,
    VoronoiNN,
)
from pymatgen.core import Structure

from graph_id import GraphIDGenerator
from graph_id.analysis.graphs import StructureGraph

from .imports import graph_id_cpp

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


def run_benchmark(pymatgen, our_nn):
    for name, s in small_test_structure():
        at = timeit.timeit("pymatgen.get_all_nn_info(s)", number=10, globals=locals()) * 100
        bt = timeit.timeit("our_nn.get_all_nn_info(s)", number=10, globals=locals()) * 100
        print(
            f"{s.num_sites: 3d} site. Python: {at: 8.3f}ms, C++: {bt: 7.3f}ms, {at / bt: 4.1f} times faster [{name}]",
        )


class TestBenchmark(unittest.TestCase):
    def test_voronoi(self):
        print("VoronoiNN:")
        run_benchmark(VoronoiNN(), graph_id_cpp.VoronoiNN())

    def test_minimum_distance(self):
        print("MinimumDistanceNN:")
        run_benchmark(MinimumDistanceNN(), graph_id_cpp.MinimumDistanceNN())

    def test_minimum_okeefee(self):
        print("MinimumOKeeffeNN:")
        run_benchmark(MinimumOKeeffeNN(), graph_id_cpp.MinimumOKeeffeNN())

    def test_crystal(self):
        print("CrystalNN:")
        run_benchmark(CrystalNN(), graph_id_cpp.CrystalNN())

    def test_cut_off_dict(self):
        print("CutOffDictNN:")
        run_benchmark(
            CutOffDictNN.from_preset("vesta_2019"),
            graph_id_cpp.CutOffDictNN.from_preset("vesta_2019"),
        )

    def test_brunner(self):
        print("BrunnerNN_real:")
        run_benchmark(BrunnerNN_real(), graph_id_cpp.BrunnerNN_real())

    def test_econ(self):
        print("EconNN:")
        run_benchmark(EconNN(), graph_id_cpp.EconNN())

    def test_structure_graph(self):
        print("StructureGraph.set_compositional_sequence_node_attr:")
        py = StructureGraph
        cpp = graph_id_cpp.StructureGraph
        for name, s in small_test_structure():
            nn = graph_id_cpp.MinimumDistanceNN()

            def f(cls):
                sg = cls.with_local_env_strategy(s, nn)  # noqa: B023
                sg.set_elemental_labels()
                sg.set_compositional_sequence_node_attr(hash_cs=True)
                sg.set_compositional_sequence_node_attr(use_previous_cs=True, hash_cs=True)

            at = timeit.timeit("f(py)", number=10, globals=locals()) * 100
            bt = timeit.timeit("f(cpp)", number=10, globals=locals()) * 100
            print(
                "{: 3d} site. Python: {: 8.3f}ms, C++: {: 7.3f}ms, {: 4.1f} times faster [{}]".format(
                    s.num_sites,
                    at,
                    bt,
                    at / bt,
                    name,
                ),
            )

    def test_graph_id(self):
        print("GraphIDGenerator.get_id:")
        a = GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure():
            n = 1
            at = timeit.timeit("a.get_id(s)", number=n, globals=locals()) * 1000 / n
            bt = timeit.timeit("b.get_id(s)", number=n, globals=locals()) * 1000 / n
            print(
                "{: 3d} site. Python: {: 8.3f}ms, C++: {: 7.3f}ms, {: 4.1f} times faster [{}]".format(
                    s.num_sites,
                    at,
                    bt,
                    at / bt,
                    name,
                ),
            )

    def test_graph_id_cpp_only(self):
        print("GraphIDGenerator.get_id:")
        import numpy as np

        g = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure(1000):
            n = 1000
            if s.num_sites > 20:
                n = 10
            t = timeit.repeat("g.get_id(s)", number=n, repeat=5, globals=locals())
            mean = np.mean(t) * 1000 / n
            std = np.std(t) * 1000 / n
            print(f"{s.num_sites: 3d} site. C++: {mean: 8.3f}ms+={std / mean * 100:.1f}% [{name}]")


if __name__ == "__main__":
    unittest.main()
