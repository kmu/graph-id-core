from ase.build import bulk
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor

from graph_id.app.maker import GraphIDMaker


def test_maker():
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    maker = GraphIDMaker()
    assert maker.get_id(nacl) == "NaCl-88c8e156db1b0fd9"

    pymaker = GraphIDMaker(engine="py")
    assert pymaker.get_id(nacl) == "NaCl-88c8e156db1b0fd9"

    maker2 = GraphIDMaker(engine="c++", depth=2)
    assert maker2.get_id(nacl) == "NaCl-b54d6c6b1662ccdb"

    site_ids = maker.get_site_ids(nacl)
    assert len(site_ids) == len(nacl)
    na_site_id = "Na_Na-8ac4127fe153fffb8164eef4dcb4b522d48e2b3c4e89ad8e257f3bf9bec2448d5d8831a3c432a741cf1c77ac3b0ce1bcfa33682aadb3ec4e5f7425857fa86ff0"  # noqa: E501
    cl_site_id = "Cl_Cl-0d9c88475b88c4eedf874aa1caab521d72e324459b530292299e66452a4c1732c7039268e190c887d2ad06fd8720eaccfd5ae981b43d05fc329c1866c2160d57"  # noqa: E501
    assert site_ids[0] == na_site_id
    assert site_ids[1] == na_site_id
    assert site_ids[2] == na_site_id
    assert site_ids[3] == na_site_id
    assert site_ids[4] == cl_site_id
    assert site_ids[5] == cl_site_id
    assert site_ids[6] == cl_site_id
    assert site_ids[7] == cl_site_id

    site_ids = maker.get_site_ids(nacl)
    assert len(site_ids) == len(nacl)
    na_site_id = "Na_Na-8ac4127fe153fffb8164eef4dcb4b522d48e2b3c4e89ad8e257f3bf9bec2448d5d8831a3c432a741cf1c77ac3b0ce1bcfa33682aadb3ec4e5f7425857fa86ff0"  # noqa: E501
    cl_site_id = "Cl_Cl-0d9c88475b88c4eedf874aa1caab521d72e324459b530292299e66452a4c1732c7039268e190c887d2ad06fd8720eaccfd5ae981b43d05fc329c1866c2160d57"  # noqa: E501
    assert site_ids[0] == na_site_id
    assert site_ids[1] == na_site_id
    assert site_ids[2] == na_site_id
    assert site_ids[3] == na_site_id
    assert site_ids[4] == cl_site_id
    assert site_ids[5] == cl_site_id
    assert site_ids[6] == cl_site_id
    assert site_ids[7] == cl_site_id


def test_reduce():
    ase_structure_primitive = bulk("NaCl", "rocksalt", a=5.692)
    ase_structure_conventional = bulk("NaCl", "rocksalt", a=5.692, cubic=True)

    structure_primitive = AseAtomsAdaptor.get_structure(ase_structure_primitive)
    structure_conventional = AseAtomsAdaptor.get_structure(ase_structure_conventional)

    maker = GraphIDMaker()
    assert maker.get_id(structure_primitive) != maker.get_id(structure_conventional)

    maker_cpp = GraphIDMaker(engine="c++", depth=4, reduce=True)
    cpp_id_primitive = maker_cpp.get_id(structure_primitive)
    cpp_id_conventional = maker_cpp.get_id(structure_conventional)

    assert cpp_id_primitive == cpp_id_conventional

    maker = GraphIDMaker(engine="py", depth=4, reduce=True)
    py_id_primitive = maker.get_id(structure_primitive)
    py_id_conventional = maker.get_id(structure_conventional)
    assert py_id_primitive == py_id_conventional

    assert py_id_primitive == cpp_id_primitive

    site_ids = maker_cpp.get_site_ids(structure_primitive)
    assert len(site_ids) == len(structure_primitive)
    na_site_id = "Na-ee0d78f1"
    cl_site_id = "Cl-77752e7a"
    assert na_site_id in site_ids[0]
    assert cl_site_id in site_ids[1]
