# ruff: noqa: D101, D102
"""Test main module import."""

import unittest

from .imports import graph_id_cpp


class TestMail(unittest.TestCase):
    def test_import(self):
        print(graph_id_cpp.__version__)
