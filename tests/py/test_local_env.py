# ruff: noqa: D103, PT011, B017
"""Tests for graph_id.analysis.local_env helpers."""

import pytest
from pymatgen.core import Lattice, Structure

from graph_id.analysis.local_env import DistanceClusteringNN, _get_original_site


def _nacl() -> Structure:
    return Structure.from_spacegroup(
        "Fm-3m",
        Lattice.cubic(5.692),
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )


def test_get_original_site_structure_branch():
    """A periodic-image neighbor resolves back to a valid original site index."""
    s = _nacl()
    neighbors = s.get_neighbors(s[0], 4.0)
    assert neighbors  # sanity: NaCl has neighbors within 4 A

    for nb in neighbors:
        idx = _get_original_site(s, nb)
        assert 0 <= idx < len(s)
        assert nb.is_periodic_image(s[idx])


def test_get_original_site_iterable_branch():
    """A plain (non-Structure) iterable uses direct site equality."""
    s = _nacl()
    sites = list(s)  # plain list -> exercises the non-Structure branch
    assert _get_original_site(sites, s[2]) == 2


def test_get_original_site_not_found_raises():
    """A site that is absent from the container raises."""
    s = _nacl()
    with pytest.raises(Exception):
        _get_original_site([], s[0])


def test_distance_clustering_nn_rank_k_out_of_range():
    """Requesting a cluster beyond the available clusters yields no neighbors."""
    s = _nacl()
    nn = DistanceClusteringNN()
    assert nn.structures_allowed is True
    # A very high rank_k exceeds the number of distance clusters -> empty list.
    assert nn.get_nn_info(s, n=0, rank_k=999, cutoff=6.0) == []
