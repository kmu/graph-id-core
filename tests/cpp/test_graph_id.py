import glob
import os
import unittest

from graph_id_py.core.graph_id import GraphIDGenerator as PyGraphIDGenerator
from pymatgen.core import Structure
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CutOffDictNN, CrystalNN, MinimumDistanceNN
import graph_id

test_file_dir = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files")
)


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
        a = PyGraphIDGenerator()
        b = graph_id.GraphIDGenerator()
        for name, s in small_test_structure():
            with self.subTest(name):
                # try:
                print(name)

                aid = a.get_id(s)
                # except Exception:
                #     self.skipTest("pymatgen error")
                self.assertEqual(aid, b.get_id(s))

    def test_digest_size(self):
        a_4 = graph_id.GraphIDGenerator(digest_size=4)
        s = Structure.from_file(os.path.join("tests/py/test_files/Fe.cif"))
        aid_4 = a_4.get_id(s)
        self.assertEqual(len(aid_4), 8)

    def test_default_digest_size(self):
        a = graph_id.GraphIDGenerator()  # default: digest_size=8
        s = Structure.from_file(os.path.join("tests/py/test_files/Fe.cif"))
        aid = a.get_id(s)
        self.assertEqual(len(aid), 16)

    def test_zeolite(self):
        structure = Structure.from_file(os.path.join("tests/py/test_files/ABW.cif"))
        gid_gen_cpp = graph_id.GraphIDGenerator()  # default: digest_size=8
        gid_gen_py = PyGraphIDGenerator()
        sg_py_cutoff = StructureGraph.from_local_env_strategy(
            structure, CutOffDictNN({("Si", "O"): 2})
        )
        sg_cpp_cutoff = graph_id.StructureGraph.with_local_env_strategy(
            structure, graph_id.CutOffDictNN({("Si", "O"): 2})
        ).to_py()
        sg_py_cn = StructureGraph.from_local_env_strategy(structure, CrystalNN())
        sg_cpp_cn = graph_id.StructureGraph.with_local_env_strategy(
            structure, graph_id.CrystalNN()
        ).to_py()
        sg_py_mdn = StructureGraph.from_local_env_strategy(structure, MinimumDistanceNN())
        sg_cpp_mdn = graph_id.StructureGraph.with_local_env_strategy(
            structure, graph_id.MinimumDistanceNN()
        ).to_py()
        assert gid_gen_cpp.get_id(structure) == gid_gen_py.get_id(structure).split("-")[2]

        assert gid_gen_cpp.get_id_with_structure_graph(sg_py_cutoff) == gid_gen_cpp.get_id_with_structure_graph(sg_cpp_cutoff)
        assert gid_gen_cpp.get_id_with_structure_graph(sg_py_cn) == gid_gen_cpp.get_id_with_structure_graph(sg_cpp_cn)
        assert gid_gen_cpp.get_id_with_structure_graph(sg_py_mdn) == gid_gen_cpp.get_id_with_structure_graph(sg_cpp_mdn)
    


if __name__ == "__main__":
    unittest.main()
