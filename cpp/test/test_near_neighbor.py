import timeit
import unittest

import glob
import os.path

import numpy as np
from pymatgen.core import Structure, Lattice, Molecule

from imports import graph_id_cpp
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.optimization.neighbors import find_points_in_spheres

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files"))


def small_test_structure(max_sites=30):
    res = []
    for p in glob.glob(os.path.join(test_file_dir, "*.cif")):
        name = p.split("/")[-1].replace(".cif", "").replace("-", "_")
        s = Structure.from_file(p)
        if s.num_sites <= max_sites:
            res.append((name, s))
    return res


class TestNN(unittest.TestCase):
    def assert_nn_info(self, a, b):
        self.assertEqual(len(a), len(b), 'mismatch array length')
        self.assertListEqual([len(x) for x in a], [len(x) for x in b], 'mismatch bonds count')

        for i in range(len(a)):
            ai = self.sort(a[i])
            bi = self.sort(b[i])
            for j in range(len(a[i])):
                self.assert_nn_info_single(ai[j], bi[j], f"i={i}, j={j}")

    def assert_nn_info_single(self, a, b, msg):
        self.assertAlmostEqual(a["weight"], b["weight"], msg=msg)
        self.assertEqual(a["site_index"], b["site_index"], msg=msg)
        if a["image"] is None:
            self.assertIsNone(b["image"], msg=msg)
        else:
            self.assertListEqual(list(a["image"]), list(b["image"]), msg=msg)
        # self.assertEqual(a["site"], b["site"], msg=msg) site は未対応

    def sort(self, a):
        return sorted(a, key=lambda x: ((x["site_index"], *x["image"]) if x["image"] else x["site_index"]))


class TestNNHelper(unittest.TestCase):
    def test_find_near_neighbors(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                for r in [0.1, 1.0, 10.0]:
                    indices_a, indices_b, images, distances = find_points_in_spheres(
                        all_coords=s.cart_coords,
                        center_coords=s.cart_coords,
                        r=r,
                        pbc=np.ascontiguousarray(s.lattice.pbc, dtype=int),
                        lattice=s.lattice.matrix,
                        tol=1e-8,
                    )
                    a = self.sort(indices_a, indices_b, images, distances)
                    indices_a2, indices_b2, images2, distances2 = graph_id_cpp.find_near_neighbors(
                        s.cart_coords,
                        s.cart_coords,
                        r,
                        np.ascontiguousarray(s.lattice.pbc, dtype=int),
                        s.lattice.matrix,
                        1e-8,
                        1.0,
                    )
                    b = self.sort(indices_a2, indices_b2, images2, distances2)

                    self.assert_result(a, b)

    def sort2(self, a, b, c, d):
        return sorted(list(zip(a, b, c, d)), key=lambda x: (x[2], x[3][0], x[3][1], x[3][2]))

    def sort(self, a, b, c, d):
        return sorted(list(zip(a, b, c, d)), key=lambda x: (x[0], x[1], x[2][0], x[2][1], x[2][2]))

    def assert_result(self, a, b):
        sa = {(x[0], x[1], tuple(x[2])) for x in a}
        sb = {(x[0], x[1], tuple(x[2])) for x in b}
        self.assertSetEqual(sa, sb, 'mismatch bonds set')
        self.assertEqual(len(a), len(b), 'mismatch array length')
        for i in range(len(a)):
            self.assertEqual(a[i][0], b[i][0], 'mismatch all_index')
            self.assertEqual(a[i][1], b[i][1], 'mismatch center_index')
            self.assertTrue(all(np.equal(a[i][2], b[i][2])), 'mismatch jindex')
            self.assertAlmostEqual(a[i][3], b[i][3], msg='mismatch distance')


class TestMinimumDistanceNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumDistanceNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumDistanceNN().molecules_allowed)

    def test_structures(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                cpp_result = graph_id_cpp.MinimumDistanceNN().get_all_nn_info(s)
                try:
                    pymatgen_result = MinimumDistanceNN().get_all_nn_info(s)
                except Exception as e:
                    print(e)
                    self.skipTest("pymatgen error")
                self.assert_nn_info(cpp_result, pymatgen_result)

    def test_molecules(self):
        m = Molecule(['H', 'H'], [[0, 0, 0], [0, 0, 1]])
        cpp_result = graph_id_cpp.MinimumDistanceNN().get_all_nn_info(m)
        try:
            pymatgen_result = MinimumDistanceNN().get_all_nn_info(m)
        except Exception as e:
            print(e)
            self.skipTest("pymatgen error")
        self.assert_nn_info(cpp_result, pymatgen_result)

    def test_structures_get_all_sites(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                cpp_result = graph_id_cpp.MinimumDistanceNN(get_all_sites=True).get_all_nn_info(s)
                pymatgen_result = MinimumDistanceNN(get_all_sites=True).get_all_nn_info(s)
                self.assert_nn_info(cpp_result, pymatgen_result)


if __name__ == "__main__":
    unittest.main()
