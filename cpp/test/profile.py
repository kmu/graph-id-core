from pymatgen.core import Structure
from imports import graph_id_cpp
from pymatgen.analysis.dimensionality import get_dimensionality_larsen
import yep

g = graph_id_cpp.GraphIDGenerator()
s = Structure.from_file("graph_id/tests/test_files/YFI.cif")  # 360 sites
# s = Structure.from_file("graph_id/tests/test_files/mp-1008498.cif")  # 4 sites
# s = Structure.from_file("graph_id/tests/test_files/mp-36.cif")  # 1 sites
yep.start("profile.prof")
for i in range(10):
    g.get_id(s)
yep.stop()
