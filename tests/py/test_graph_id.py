from pathlib import Path
from unittest import TestCase

import pytest
from pymatgen.analysis.local_env import CrystalNN, MinimumDistanceNN
from pymatgen.core import Element, Lattice, Structure

from graph_id.core.graph_id import FixedDepthGraphIDGenerator, GraphIDGenerator

TEST_FILES = (Path(__file__).resolve().parent / "test_files").as_posix()


class TestGraphIDGenerator(TestCase):
    def test_diameter_0(self):
        si = Structure.from_file(f"{TEST_FILES}/mp-1056579.cif")
        sr = Structure.from_file(f"{TEST_FILES}/mp-1056418.cif")
        gid = GraphIDGenerator()

        assert not gid.are_same(si, sr)

    def test_version(self):
        gid = GraphIDGenerator()
        self.assertTrue(gid.version > "0.0.0")

    def test_nacl(self):
        nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        cscl = nacl.copy()
        cscl.replace(0, Element("Cs"))

        gid = GraphIDGenerator()
        gid_topo = GraphIDGenerator(
            topology_only=True,
        )
        gid_topo_wyckoff = GraphIDGenerator(
            topology_only=True,
            wyckoff=True,
        )

        self.assertEqual("NaCl-3D-88c8e156db1b0fd9", gid.get_id(nacl))
        self.assertEqual("CsNa3Cl4-3D-1a85e9c5247fb74a", gid.get_id(cscl))

        self.assertEqual("3D-c144aa1ceffbd9af", gid_topo.get_id(nacl))
        self.assertEqual("3D-c144aa1ceffbd9af", gid_topo.get_id(cscl))

        self.assertEqual("3D-e40be9333fa6f8ae", gid_topo_wyckoff.get_id(nacl))
        self.assertEqual("3D-e40be9333fa6f8ae", gid_topo_wyckoff.get_id(cscl))

    def test_limnteo(self):
        s1 = Structure.from_file(f"{TEST_FILES}/mp-1299593.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1307172.cif")

        gid = GraphIDGenerator(diameter_factor=1)

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_vsb_o_4(self):
        """
        MinimumDistanceNN does not work for this.
        """
        s1 = Structure.from_file(f"{TEST_FILES}/VSbO4.cif")

        gid = GraphIDGenerator(nn=CrystalNN(), diameter_factor=1)
        id_1 = gid.get_id(s1)

        self.assertEqual(id_1, "VSbO4-0D-e50201525efe4cd5")

    def test_quartz(self):
        alpha = Structure.from_file(f"{TEST_FILES}/298 K.cif")
        beta = Structure.from_file(f"{TEST_FILES}/1078 K.cif")

        gid = GraphIDGenerator()

        id_a = gid.get_id(alpha)
        id_b = gid.get_id(beta)

        self.assertEqual(id_a, id_b)
        self.assertEqual(id_a, "SiO2-3D-aab30488e82694e7")

        gid = GraphIDGenerator(wyckoff=True)

        id_a_w = gid.get_id(alpha)
        id_b_w = gid.get_id(beta)

        self.assertNotEqual(id_a_w, id_a)
        self.assertNotEqual(id_a_w, id_b_w)
        self.assertEqual(id_a_w, "SiO2-3D-a37c720fadbbd609")

    def test_carbon(self):
        layer = Structure.from_file(f"{TEST_FILES}/mp-48.cif")
        bulk = Structure.from_file(f"{TEST_FILES}/mp-1018088.cif")

        gid = GraphIDGenerator(nn=CrystalNN(), diameter_factor=1)

        id_1 = gid.get_id(layer)
        id_2 = gid.get_id(bulk)

        self.assertNotEqual(id_1, id_2)

        ids = gid.get_component_ids(layer)
        self.assertEqual(len(ids), 2)

        assert not gid.are_same(layer, bulk)

        ids = gid.get_many_ids([layer, bulk], parallel=True)
        self.assertEqual(ids[0], id_1)
        self.assertEqual(ids[1], id_2)

        ids2 = gid.get_many_ids([layer, bulk], parallel=False)
        self.assertListEqual(ids, ids2)

    def test_simple_same_composition(self):
        """
        Graphs with diameter of 0.
        Compositions are identical.
        Structures are different.
        """
        s1 = Structure.from_file(f"{TEST_FILES}/mp-36.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1008681.cif")

        gid = GraphIDGenerator(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_calcium(self):
        s1 = Structure.from_file(f"{TEST_FILES}/mp-1008498.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1067285.cif")

        gid = GraphIDGenerator(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_connected_components(self):
        s1 = Structure.from_file(f"{TEST_FILES}/mp-121.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-611219.cif")

        gid = GraphIDGenerator(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_get_unique_structures(self):
        alpha = Structure.from_file(f"{TEST_FILES}/298 K.cif")
        beta = Structure.from_file(f"{TEST_FILES}/1078 K.cif")

        gid = GraphIDGenerator()

        unique_structures = gid.get_unique_structures([alpha, beta])
        self.assertEqual(len(unique_structures), 1)

    def test_get_id_catch_error(self):
        s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        gid = GraphIDGenerator()
        self.assertEqual(gid.get_id_catch_error(None), "")
        self.assertEqual(gid.get_id_catch_error(s), "NaCl-3D-88c8e156db1b0fd9")

    def test_inappropriate_combinations(self):
        with pytest.raises(ValueError):  # noqa: PT011
            GraphIDGenerator(wyckoff=True, loop=True)
        with pytest.raises(ValueError):  # noqa: PT011
            GraphIDGenerator(topology_only=True, loop=True)

    def test_reduce_symmetry(self):
        nacl_conventional = Structure.from_spacegroup(
            "Fm-3m",
            Lattice.cubic(5.692),
            ["Na", "Cl"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )

        nacl_primitive = nacl_conventional.get_primitive_structure()

        generator = GraphIDGenerator()
        conventional_id = generator.get_id(nacl_conventional)
        primitive_id = generator.get_id(nacl_primitive)
        assert conventional_id != primitive_id

        generator = GraphIDGenerator(diameter_factor=0, additional_depth=3)
        generator_fixed = FixedDepthGraphIDGenerator(depth=3)
        assert generator.get_id(nacl_primitive) == generator_fixed.get_id(nacl_primitive)

        generator_fixed = FixedDepthGraphIDGenerator(depth=6)
        assert generator_fixed.get_id(nacl_conventional) != generator_fixed.get_id(nacl_primitive)

        generator_fixed = FixedDepthGraphIDGenerator(depth=6, reduce_symmetry=True)
        assert generator_fixed.get_id(nacl_conventional) == generator_fixed.get_id(nacl_primitive)

    def test_one_site_reduction(self):
        """
        Graph ID should be compatible for reduce_symmetry =True/False
        """
        one_site_structure = Structure.from_file(f"{TEST_FILES}/mp-36.cif")
        generator = FixedDepthGraphIDGenerator(depth=6, reduce_symmetry=False)
        generator_reduce = FixedDepthGraphIDGenerator(depth=6, reduce_symmetry=True)
        assert generator.get_id(one_site_structure) == generator_reduce.get_id(one_site_structure)
