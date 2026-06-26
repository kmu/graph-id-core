import pytest
from pymatgen.analysis.local_env import CovalentBondNN
from pymatgen.core import Lattice, Structure

from graph_id.analysis.graphs import StructureGraph


def test_disallowed_strategy():
    """Test that disallowed strategies raise ValueError."""
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    with pytest.raises(ValueError):  # noqa: PT011
        StructureGraph.from_local_env_strategy(s, CovalentBondNN())


def test_indivisual_strategy_disallowed_strategy():
    """with_indivisual_state_comp_strategy also rejects non-structure strategies."""
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    sg = StructureGraph.from_empty_graph(s, name="bonds")
    with pytest.raises(ValueError):  # noqa: PT011
        StructureGraph.with_indivisual_state_comp_strategy(s, CovalentBondNN(), sg, 0)
