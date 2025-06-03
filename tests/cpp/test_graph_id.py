import glob
import os
import unittest

import graph_id
from pymatgen.core import Structure

from .imports import graph_id_cpp

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites and name != "VSbO4":
            res.append((name, s))
    return res


class TestGraphIDGenerator(unittest.TestCase):
    def test_get_id(self):
        a = graph_id.GraphIDGenerator()
        b = graph_id_cpp.GraphIDGenerator()
        for name, s in small_test_structure():
            with self.subTest(name):
                # try:
                print(name)

                aid = a.get_id(s)
                # except Exception:
                #     self.skipTest("pymatgen error")
                self.assertEqual(aid, b.get_id(s))

    def test_digest_size(self):
        a_4 = graph_id_cpp.GraphIDGenerator(digest_size=4)
        s = Structure.from_file(os.path.join("tests/py/test_files/Fe.cif"))
        aid_4 = a_4.get_id(s)
        self.assertEqual(len(aid_4), 8)


if __name__ == "__main__":
    unittest.main()
