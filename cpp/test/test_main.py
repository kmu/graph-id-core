import os.path
import unittest

from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.core import Structure

from cpp.test.imports import graph_id_cpp


class TestMail(unittest.TestCase):
    def test_main(self):
        p = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../graph_id/tests/test_files/mp-48.cif"))
        s = Structure.from_file(p)
        self.assertEqual(
            graph_id_cpp.MinimumDistanceNN().get_all_nn_info(s),
            MinimumDistanceNN().get_all_nn_info(s),
        )

