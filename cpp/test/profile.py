from pymatgen.core import Structure
from .imports import graph_id_cpp
import yep

g = graph_id_cpp.GraphIDGenerator()
s_large = Structure.from_file("graph_id/tests/test_files/YFI.cif")  # 360 sites
s_medium = Structure.from_file("graph_id/tests/test_files/mp-1307172.cif")  # 30 sites
s_small = Structure.from_file("graph_id/tests/test_files/mp-1008498.cif")  # 4 sites

yep.start("build/profile-large.prof")
for i in range(10):
    g.get_id(s_large)
yep.stop()

yep.start("build/profile-medium.prof")
for i in range(1000):
    g.get_id(s_medium)
yep.stop()

yep.start("build/profile-small.prof")
for i in range(100000):
    g.get_id(s_small)
yep.stop()
