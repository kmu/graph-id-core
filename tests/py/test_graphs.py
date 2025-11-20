import pytest
from pymatgen.analysis.local_env import CovalentBondNN
from pymatgen.core import Lattice, Structure

from graph_id.analysis.graphs import StructureGraph


def test_disallowed_strategy():
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    with pytest.raises(ValueError):  # noqa: PT011
        StructureGraph.with_local_env_strategy(s, CovalentBondNN())
