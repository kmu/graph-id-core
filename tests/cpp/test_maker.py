from pymatgen.core import Lattice, Structure

from graph_id.app.maker import GraphIDMaker


def test_maker():
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    maker = GraphIDMaker()
    assert maker.get_id(nacl) == "NaCl-88c8e156db1b0fd9"

    pymaker = GraphIDMaker(engine="py")
    assert pymaker.get_id(nacl) == "NaCl-88c8e156db1b0fd9"

    maker = GraphIDMaker(engine="c++", depth=2)
    assert maker.get_id(nacl) == "NaCl-b54d6c6b1662ccdb"
