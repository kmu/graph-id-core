import glob
import os
import unittest
import timeit
import graph_id

from imports import graph_id_cpp

from pymatgen.core import Structure

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


class TestGraphIDGenerator(unittest.TestCase):
    def test_get_id(self):
        a = graph_id.GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
            with self.subTest(p.split("/")[-1]):
                s = Structure.from_file(p)
                self.assertEqual(a.get_id(s), b.get_id(s))

    def test_benchmark(self):
        a = graph_id.GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
            with self.subTest(p.split("/")[-1]):
                s = Structure.from_file(p)
                at = timeit.timeit("a.get_id(s)", number=10, globals=locals()) * 100
                bt = timeit.timeit("b.get_id(s)", number=10, globals=locals()) * 100
                print("Python: {:.3f}ms, C++: {:.3f}ms, {:.1f} times faster".format(at, bt, at / bt))
