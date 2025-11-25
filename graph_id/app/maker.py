from __future__ import annotations

from graph_id_cpp import GraphIDGenerator as CppGraphIDGenerator
from graph_id_cpp import MinimumDistanceNN

from graph_id.core.graph_id import GraphIDGenerator as PyGraphIDGenerator

# from pymatgen.analysis.local_env import CrystalNN, MinimumDistanceNN


class GraphIDMaker:
    def __init__(
        self,
        nn=None,
        depth: int | None = None,
        reduce_symmetry: bool = False,
        engine: str = "c++",
    ) -> None:
        """
        A simple interface to make GraphID.
        """
        self.engine = engine
        self.reduce_symmetry = reduce_symmetry
        if nn is None:
            nn = MinimumDistanceNN()

        diameter_factor = 2
        additional_depth = 1
        if depth is not None:
            diameter_factor = 0
            additional_depth = depth

        if engine == "py":
            self.generator = PyGraphIDGenerator(
                nn=nn,
                diameter_factor=diameter_factor,
                additional_depth=additional_depth,
                prepend_composition=False,
                prepend_dimensionality=False,
            )

        elif engine == "c++":
            self.generator = CppGraphIDGenerator(
                nn=nn,
                diameter_factor=diameter_factor,
                additional_depth=additional_depth,
            )

    def get_id(self, structure) -> str:
        graph_id = self.generator.get_id(structure)

        return f"{structure.composition.reduced_formula}-{graph_id}"
