import os.path
import unittest

from cpp.test.imports import graph_id_cpp

class TestMail(unittest.TestCase):
    def test_main(self):
        print(graph_id_cpp.test)
