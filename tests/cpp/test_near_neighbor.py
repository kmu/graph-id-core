import glob
import os.path
import unittest

import numpy as np
import pytest
from pymatgen.analysis.local_env import (
    BrunnerNN_real,
    BrunnerNN_reciprocal,
    BrunnerNN_relative,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    MinimumDistanceNN,
    MinimumOKeeffeNN,
    VoronoiNN,
)
from pymatgen.core import Lattice, Molecule, Structure
from pymatgen.core.periodic_table import Specie
from pymatgen.optimization.neighbors import find_points_in_spheres

from .imports import graph_id_cpp

test_file_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../py/test_files"))


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
        self.assertEqual(len(a), len(b), "mismatch array length")
        self.assertListEqual([len(x) for x in a], [len(x) for x in b], "mismatch bonds count")

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
        return sorted(
            a,
            key=lambda x: ((x["site_index"], *x["image"]) if x["image"] else x["site_index"]),
        )

    def run_for_small_structures(self, pymatgen_nn, out_nn, skip_names=[]):  # noqa: B006
        for name, s in small_test_structure():
            print(name)
            if name not in skip_names:
                with self.subTest(name):
                    pymatgen_result = pymatgen_nn.get_all_nn_info(s)

                    cpp_result = out_nn.get_all_nn_info(s)
                    self.assert_nn_info(cpp_result, pymatgen_result)


class TestNNHelper(unittest.TestCase):
    def test_find_near_neighbors(self):
        for name, s in small_test_structure():
            with self.subTest(name):
                for r in [0.1, 1.0, 10.0]:
                    indices_a, indices_b, images, distances = find_points_in_spheres(
                        all_coords=s.cart_coords,
                        center_coords=s.cart_coords,
                        r=r,
                        pbc=np.ascontiguousarray(s.lattice.pbc, dtype=np.int64),
                        lattice=s.lattice.matrix,
                        tol=1e-8,
                    )
                    a = self.sort(indices_a, indices_b, images, distances)
                    (
                        indices_a2,
                        indices_b2,
                        images2,
                        distances2,
                    ) = graph_id_cpp.find_near_neighbors(
                        s.cart_coords,
                        s.cart_coords,
                        r,
                        s.lattice.pbc,
                        s.lattice.matrix,
                        1e-8,
                        1.0,
                    )
                    b = self.sort(indices_a2, indices_b2, images2, distances2)

                    self.assert_result(a, b)

    def sort2(self, a, b, c, d):
        return sorted(list(zip(a, b, c, d)), key=lambda x: (x[2], x[3][0], x[3][1], x[3][2]))  # noqa: C414

    def sort(self, a, b, c, d):
        return sorted(list(zip(a, b, c, d)), key=lambda x: (x[0], x[1], x[2][0], x[2][1], x[2][2]))  # noqa: C414

    def assert_result(self, a, b):
        sa = {(x[0], x[1], tuple(x[2])) for x in a}
        sb = {(x[0], x[1], tuple(x[2])) for x in b}
        self.assertSetEqual(sa, sb, "mismatch bonds set")
        self.assertEqual(len(a), len(b), "mismatch array length")
        for i in range(len(a)):
            self.assertEqual(a[i][0], b[i][0], "mismatch all_index")
            self.assertEqual(a[i][1], b[i][1], "mismatch center_index")
            self.assertTrue(all(np.equal(a[i][2], b[i][2])), "mismatch jindex")
            self.assertAlmostEqual(a[i][3], b[i][3], msg="mismatch distance")


class TestVoronoiNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.VoronoiNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertFalse(graph_id_cpp.VoronoiNN().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(VoronoiNN(), graph_id_cpp.VoronoiNN(), skip_names=["pt_slab", "VSbO4"])

    def test_get_voronoi_polyhedra(self):
        pmg = VoronoiNN()
        cpp = graph_id_cpp.VoronoiNN()
        for name, s in small_test_structure():
            if name != "VSbO4":
                with self.subTest(name):
                    for i in range(s.num_sites):
                        pmg_res = pmg.get_voronoi_polyhedra(s, i)
                        cpp_res = cpp.get_voronoi_polyhedra(s, i)
                        self.assert_voronoi_polyhedra(pmg_res, cpp_res)

    def test_get_all_voronoi_polyhedra(self):
        pmg = VoronoiNN(allow_pathological=True)
        cpp = graph_id_cpp.VoronoiNN(allow_pathological=True)
        for name, s in small_test_structure():
            with self.subTest(name):
                pmg_res = pmg.get_all_voronoi_polyhedra(s)

                self.skipTest("pymatgen error")
                cpp_res = cpp.get_all_voronoi_polyhedra(s)
                self.assertEqual(len(pmg_res), len(cpp_res))
                for i in range(s.num_sites):
                    self.assert_voronoi_polyhedra(pmg_res[i], cpp_res[i])

    def assert_voronoi_polyhedra(self, pmg_res, cpp_res):
        self.assertEqual(len(pmg_res), len(cpp_res))
        # 返り値の順番が違うので、site の座標を基準に順番を揃える
        pmg_sites = np.array([s["site"].frac_coords for _, s in sorted(pmg_res.items())])
        cpp_sites = np.array([s["site"].frac_coords for _, s in sorted(cpp_res.items())])
        pmg_keys = sorted(pmg_res.keys())
        cpp_keys = sorted(cpp_res.keys())
        a = np.argmin(
            np.linalg.norm(pmg_sites.reshape(1, -1, 3) - cpp_sites.reshape(-1, 1, 3), axis=-1),
            axis=-1,
        )
        cpp2pmg = {cpp_keys[i]: pmg_keys[a[i]] for i in range(len(pmg_res))}

        for cpp_key in cpp_keys:
            pmg_key = cpp2pmg[cpp_key]
            self.assertAlmostEqual(
                np.linalg.norm(pmg_res[pmg_key]["site"].frac_coords - cpp_res[cpp_key]["site"].frac_coords),
                0,
                msg="site mismatch",
            )
            self.assertAlmostEqual(
                np.linalg.norm(pmg_res[pmg_key]["normal"] - cpp_res[cpp_key]["normal"]),
                0,
                msg="site mismatch",
            )
            self.assertAlmostEqual(
                pmg_res[pmg_key]["solid_angle"],
                cpp_res[cpp_key]["solid_angle"],
                msg="solid_angle mismatch",
            )
            self.assertAlmostEqual(pmg_res[pmg_key]["area"], cpp_res[cpp_key]["area"], msg="area mismatch")
            self.assertAlmostEqual(
                pmg_res[pmg_key]["face_dist"],
                cpp_res[cpp_key]["face_dist"],
                msg="face_dist mismatch",
            )
            self.assertAlmostEqual(
                pmg_res[pmg_key]["volume"],
                cpp_res[cpp_key]["volume"],
                msg="volume mismatch",
            )
            self.assertEqual(
                pmg_res[pmg_key]["n_verts"],
                cpp_res[cpp_key]["n_verts"],
                msg="n_verts mismatch",
            )
            self.assertEqual(
                set(pmg_res[pmg_key]["adj_neighbors"]),
                {cpp2pmg[v] for v in cpp_res[cpp_key]["adj_neighbors"]},
                msg="adj_neighbors mismatch",
            )


class TestMinimumDistanceNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumDistanceNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumDistanceNN().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(MinimumDistanceNN(), graph_id_cpp.MinimumDistanceNN(), skip_names=["VSbO4"])

    def test_molecules(self):
        m = Molecule(["H", "H"], [[0, 0, 0], [0, 0, 1]])
        cpp_result = graph_id_cpp.MinimumDistanceNN().get_all_nn_info(m)

        pymatgen_result = MinimumDistanceNN().get_all_nn_info(m)

        self.assert_nn_info(cpp_result, pymatgen_result)

    def test_structures_get_all_sites(self):
        self.run_for_small_structures(
            MinimumDistanceNN(get_all_sites=True),
            graph_id_cpp.MinimumDistanceNN(get_all_sites=True),
        )


class TestMinimumOKeeffeNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumOKeeffeNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertTrue(graph_id_cpp.MinimumOKeeffeNN().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(
            MinimumOKeeffeNN(),
            graph_id_cpp.MinimumOKeeffeNN(),
            skip_names=["pt_slab", "VSbO4"],
        )


class TestCrystalNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.CrystalNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertFalse(graph_id_cpp.CrystalNN().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(CrystalNN(), graph_id_cpp.CrystalNN())


class TestPmgCrystalNN(TestNN):
    # class TestCustomCrystalNN(unittest.TestCase):
    def setUp(self):
        self.lattice = Lattice.cubic(10)
        self.nacl_structure = Structure(
            self.lattice,
            [Specie("Na", 1), Specie("Cl", -1)],
            [[0, 0, 0], [0.5, 0, 0]],  # Distance 5, relative to cell this is 0.05
        )
        self.fe_o_structure = Structure(
            self.lattice,
            [Specie("Fe", 2), Specie("O", -2)],  # Ensure species have oxi_state and X
            [[0, 0, 0], [2, 0, 0]],
        )  # Dist 2A

    def test_get_nn_data_cation_anion(self):
        """Test CrystalNN.get_nn_data with cation_anion=True (covers lines 348-355)."""
        nn = CrystalNN(cation_anion=True)
        nn_data = nn.get_nn_data(self.nacl_structure, 0)  # Na+ at site 0
        assert nn_data is not None
        if nn_data.all_nninfo:  # Check if any neighbors found
            for entry in nn_data.all_nninfo:
                assert entry["site"].specie.symbol == "Cl"
                assert entry["site"].specie.oxi_state * self.nacl_structure[0].specie.oxi_state < 0

        # Test ValueError if no valid targets
        s_no_valid_targets = Structure(self.lattice, [Specie("Na", 1), Specie("K", 1)], [[0, 0, 0], [0.5, 0, 0]])
        with pytest.raises(
            ValueError,
            match="No valid targets for site within cation_anion constraint!",
        ):
            nn.get_nn_data(s_no_valid_targets, 0)


class TestCutoffDictNN(TestNN):
    def test_structures(self):
        self.run_for_small_structures(
            CutOffDictNN.from_preset("vesta_2019"),
            graph_id_cpp.CutOffDictNN.from_preset("vesta_2019"),
        )

    def test_molecules(self):
        m = Molecule(["H", "H"], [[0, 0, 0], [0, 0, 1]])
        cpp_result = graph_id_cpp.CutOffDictNN.from_preset("vesta_2019").get_all_nn_info(m)

        pymatgen_result = CutOffDictNN.from_preset("vesta_2019").get_all_nn_info(m)

        self.assert_nn_info(cpp_result, pymatgen_result)

    def test_structures_with_dict(self):
        d = CutOffDictNN.from_preset("vesta_2019").cut_off_dict
        self.run_for_small_structures(CutOffDictNN(cut_off_dict=d), graph_id_cpp.CutOffDictNN(cut_off_dict=d))

    def test_structures_with_empty_dict(self):
        self.run_for_small_structures(CutOffDictNN(), graph_id_cpp.CutOffDictNN())


class TestBrunnerNNReciprocal(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.BrunnerNN_reciprocal().structures_allowed)

    def test_molecule_allowed(self):
        self.assertFalse(graph_id_cpp.BrunnerNN_reciprocal().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(
            BrunnerNN_reciprocal(),
            graph_id_cpp.BrunnerNN_reciprocal(),
            skip_names=["VSbO4"],
        )


class TestBrunnerNNRelative(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.BrunnerNN_relative().structures_allowed)

    def test_molecule_allowed(self):
        self.assertFalse(graph_id_cpp.BrunnerNN_relative().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(BrunnerNN_relative(), graph_id_cpp.BrunnerNN_relative(), skip_names=["VSbO4"])


class TestBrunnerNNReal(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.BrunnerNN_real().structures_allowed)

    def test_molecule_allowed(self):
        self.assertFalse(graph_id_cpp.BrunnerNN_real().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(BrunnerNN_real(), graph_id_cpp.BrunnerNN_real(), skip_names=["VSbO4"])


class TestEconNN(TestNN):
    def test_structure_allowed(self):
        self.assertTrue(graph_id_cpp.EconNN().structures_allowed)

    def test_molecule_allowed(self):
        self.assertTrue(graph_id_cpp.EconNN().molecules_allowed)

    def test_structures(self):
        self.run_for_small_structures(EconNN(), graph_id_cpp.EconNN(), skip_names=["VSbO4"])

    def test_structures_using_fir(self):
        self.run_for_small_structures(
            EconNN(use_fictive_radius=True),
            graph_id_cpp.EconNN(use_fictive_radius=True),
            skip_names=["VSbO4"],
        )


if __name__ == "__main__":
    unittest.main()
