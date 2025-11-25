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

    maker = GraphIDMaker(engine="c++", depth=2)
    assert maker.get_id(nacl) == "NaCl-b54d6c6b1662ccdb"


def test_reduce_symmetry():
    ase_structure_primitive = bulk("NaCl", "rocksalt", a=5.692)
    ase_structure_conventional = bulk("NaCl", "rocksalt", a=5.692, cubic=True)

    structure_primitive = AseAtomsAdaptor.get_structure(ase_structure_primitive)
    structure_conventional = AseAtomsAdaptor.get_structure(ase_structure_conventional)

    maker = GraphIDMaker()
    assert maker.get_id(structure_primitive) != maker.get_id(structure_conventional)

    maker_cpp = GraphIDMaker(engine="c++", depth=4, reduce_symmetry=True)
    cpp_id_primitive = maker_cpp.get_id(structure_primitive)
    cpp_id_conventional = maker_cpp.get_id(structure_conventional)

    assert cpp_id_primitive == cpp_id_conventional

    maker = GraphIDMaker(engine="py", depth=4, reduce_symmetry=True)
    py_id_primitive = maker.get_id(structure_primitive)
    py_id_conventional = maker.get_id(structure_conventional)
    assert py_id_primitive == py_id_conventional

    assert py_id_primitive == cpp_id_primitive
