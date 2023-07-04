import os
from unittest import TestCase

from pymatgen.analysis.local_env import CrystalNN, MinimumDistanceNN
from pymatgen.core import Element, Lattice, Structure

from graph_id.core.graph_id import GraphID

TEST_FILES = os.path.dirname(os.path.abspath(__file__)) + "/test_files"


class TestGraphID(TestCase):
    def test_diameter_0(self):
        si = Structure.from_file(f"{TEST_FILES}/mp-1056579.cif")
        sr = Structure.from_file(f"{TEST_FILES}/mp-1056418.cif")
        gid = GraphID(depth_factor=1, force_supercell=2)

        assert not gid.are_same(si, sr)

    def test_version(self):
        gid = GraphID()
        self.assertTrue(gid.version > "0.0.0")

    def test_NaCl(self):
        nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        cscl = nacl.copy()
        cscl.replace(0, Element("Cs"))

        gid = GraphID(comp_dim=True)
        gid_topo = GraphID(topology_only=True, comp_dim=True)
        gid_topo_wyckoff = GraphID(topology_only=True, wyckoff=True, comp_dim=True)

        self.assertEqual("NaCl-3D-85844dcbb16110e5c27bb042ea794cf9", gid.get_id(nacl))
        self.assertEqual("CsNa3Cl4-3D-6f29b7e87826fe2a433e502eb1e6e670", gid.get_id(cscl))

        self.assertEqual("3D-ddab0ca91146bfcb0b6029c2690dd616", gid_topo.get_id(nacl))
        self.assertEqual("3D-ddab0ca91146bfcb0b6029c2690dd616", gid_topo.get_id(cscl))

        self.assertEqual("3D-866249dfed8626b6c900e533beb58a57", gid_topo_wyckoff.get_id(nacl))
        self.assertEqual("3D-866249dfed8626b6c900e533beb58a57", gid_topo_wyckoff.get_id(cscl))

    def test_cristobalite_tridymite(self):
        """
        実はCristobaliteとTridymiteは、拡張してあげればdで見分けられてしまう
        """
        s1 = Structure.from_file(f"{TEST_FILES}/SiO2_mp-6945_computed.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/SiO2_mp-7087_computed.cif")

        gid = GraphID(depth_factor=1)

        assert gid.are_same(s1, s2)

        gid_super = GraphID(depth_factor=1, force_supercell=2)
        assert not gid_super.are_same(s1, s2)

    def test_LiMnTeO(self):
        s1 = Structure.from_file(f"{TEST_FILES}/mp-1299593.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1307172.cif")

        gid = GraphID(depth_factor=1, comp_dim=True)

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_VSbO4(self):
        """
        MinimumDistanceNN does not work for this.
        """
        s1 = Structure.from_file(f"{TEST_FILES}/VSbO4.cif")

        gid = GraphID(nn=CrystalNN(), depth_factor=1, comp_dim=True)
        id_1 = gid.get_id(s1)

        self.assertEqual(id_1, "VSbO4-0D-1351c328b38b4ad4737264fdc4be6e47")

    def test_quartz(self):
        alpha = Structure.from_file(f"{TEST_FILES}/298 K.cif")
        beta = Structure.from_file(f"{TEST_FILES}/1078 K.cif")

        gid = GraphID(comp_dim=True)

        id_a = gid.get_id(alpha)
        id_b = gid.get_id(beta)

        self.assertEqual(id_a, id_b)
        self.assertEqual(id_a, "SiO2-3D-20a961a2a8c4c132946f8d7329f3960e")

        gid = GraphID(wyckoff=True, comp_dim=True)

        id_a_w = gid.get_id(alpha)
        id_b_w = gid.get_id(beta)

        self.assertNotEqual(id_a_w, id_a)
        self.assertNotEqual(id_a_w, id_b_w)
        self.assertEqual(id_a_w, "SiO2-3D-c0771f6df4439f53243ff2607a78bb4a")

    def test_carbon(self):
        layer = Structure.from_file(f"{TEST_FILES}/mp-48.cif")
        bulk = Structure.from_file(f"{TEST_FILES}/mp-1018088.cif")

        gid = GraphID(nn=CrystalNN(), depth_factor=1, comp_dim=True)

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
        Graphs with diamter of 0.
        Compositions are identical.
        Structures are different.
        """
        s1 = Structure.from_file(f"{TEST_FILES}/mp-36.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1008681.cif")

        gid = GraphID(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_calcium(self):
        """
        CrystalNNでは問題ないが、MinmumDistanceNNでは同一視されてしまう構造。
        セルを拡張すれば大丈夫。
        あるノードから同じノードに複数個エッジが伸びていた場合、
        スーパーセルにするというロジックを使えば良さそう。
        """
        s1 = Structure.from_file(f"{TEST_FILES}/mp-1008498.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-1067285.cif")

        # sg1 = StructureGraph.with_local_env_strategy(s1, MinimumDistanceNN())
        # sg2 = StructureGraph.with_local_env_strategy(s2, MinimumDistanceNN())

        # vw = VestaWriter(sg1, False)
        # vw.write_file(filename="/home/mrok/sandbox/1.vesta")

        # vw = VestaWriter(sg2, False)
        # vw.write_file(filename="/home/mrok/sandbox/2.vesta")

        # print(diameter(sg1.graph.to_undirected()))

        # self.assertNotEqual(sg1.get_graph_id4(), sg2.get_graph_id4())

        gid = GraphID(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

    def test_connected_components(self):
        """
        セルを拡張したときにはじめてconnected componentsにわかれる例。
        そうした時はセルの拡張が必要。
        """
        s1 = Structure.from_file(f"{TEST_FILES}/mp-121.cif")
        s2 = Structure.from_file(f"{TEST_FILES}/mp-611219.cif")

        gid = GraphID(MinimumDistanceNN())

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)

        gid = GraphID(MinimumDistanceNN(), comp_dim=True)

        id_1 = gid.get_id(s1)
        id_2 = gid.get_id(s2)

        self.assertNotEqual(id_1, id_2)
