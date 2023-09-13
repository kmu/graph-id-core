import unittest

import glob
import os.path
from pymatgen.core import Structure

from imports import graph_id_cpp
from pymatgen.analysis.local_env import MinimumDistanceNN


test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


class TestNN(unittest.TestCase):
    def assert_nn_info(self, a, b):
        self.assertEqual(len(a), len(b), 'mismatch array length')
        self.assertListEqual([len(x) for x in a], [len(x) for x in b], 'mismatch bonds count')

        for i in range(len(a)):
            for j in range(len(a[i])):
                self.assert_nn_info_single(a[i][j], b[i][j], f"i={i}, j={j}")

    def assert_nn_info_single(self, a, b, msg):
        self.assertAlmostEqual(a["weight"], b["weight"], msg=msg)
        self.assertEqual(a["site_index"], b["site_index"], msg=msg)
        self.assertListEqual(list(a["image"]), list(b["image"]), msg=msg)
        self.assertEqual(a["site"], b["site"], msg=msg)


class TestMinimumDistanceNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumDistanceNN().structures_allowed)

    def test_structures(self):
        for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
            with self.subTest(p.split("/")[-1]):
                s = Structure.from_file(p)
                cpp_result = graph_id_cpp.MinimumDistanceNN().get_all_nn_info(s)
                pymatgen_result = MinimumDistanceNN().get_all_nn_info(s)
                self.assert_nn_info(cpp_result, pymatgen_result)


