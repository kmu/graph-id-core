"""Simple test to check if find_near_neighbors causes segfault"""
import graph_id_cpp
import numpy as np
from pymatgen.core import Lattice, Structure


def test_segfault_simple():
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    pbc = np.ascontiguousarray(s.lattice.pbc, dtype=int)

    print("Testing find_near_neighbors...", flush=True)

    result = graph_id_cpp.find_near_neighbors(s.cart_coords, s.cart_coords, 1.0, pbc, s.lattice.matrix, 1e-8, 1.0)
    print("SUCCESS: No segfault")
    print(f"Result type: {type(result)}")
    assert True
