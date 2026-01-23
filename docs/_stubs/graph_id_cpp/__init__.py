# ruff: noqa: ANN101, ANN003, ANN001, ANN201, ANN401, ARG002
"""Stub module for graph_id_cpp to allow documentation generation."""


class GraphIDGenerator:
    """C++ implementation of Graph ID generator."""

    def __init__(self, **kwargs):
        pass

    def get_id(self, structure) -> str:
        return ""

    def prepare_structure_graph(self, structure):
        pass


class MinimumDistanceNN:
    """C++ implementation of MinimumDistanceNN."""


class CrystalNN:
    """C++ implementation of CrystalNN."""


class CutOffDictNN:
    """C++ implementation of CutOffDictNN."""

    def __init__(self, cutoff_dict):
        pass


class StructureGraph:
    """C++ implementation of StructureGraph."""

    @staticmethod
    def with_local_env_strategy(structure, strategy):
        pass

    def to_py(self):
        pass
