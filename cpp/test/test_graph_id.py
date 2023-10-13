import glob
import os
import unittest
import timeit
import graph_id

from .imports import graph_id_cpp

from pymatgen.core import Structure

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


class TestGraphIDGenerator(unittest.TestCase):
    def test_get_id(self):
        a = graph_id.GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure():
            with self.subTest(name):
                try:
                    aid = a.get_id(s)
                except Exception as e:
                    self.skipTest('pymatgen error')
                self.assertEqual(aid, b.get_id(s))


if __name__ == "__main__":
    unittest.main()
