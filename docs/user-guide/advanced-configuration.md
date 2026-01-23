# Advanced Configuration

This guide covers advanced options for fine-tuning Graph ID generation.

## Neighbor Detection Strategies



### MinimumDistanceNN (Default)

```python
from pymatgen.analysis.local_env import MinimumDistanceNN
from graph_id import GraphIDMaker

maker = GraphIDMaker(nn=MinimumDistanceNN())
```

### CrystalNN

=== "C++ Engine (Faster)"

    ```python
    from graph_id_cpp import CrystalNN, GraphIDGenerator

    gen = GraphIDGenerator(nn=CrystalNN())
    ```

=== "Python Engine"

    ```python
    from pymatgen.analysis.local_env import CrystalNN
    from graph_id.core.graph_id import GraphIDGenerator

    gen = GraphIDGenerator(nn=CrystalNN())
    ```

### Distance Clustering

Uses DBSCAN clustering on interatomic distances:

```python
from graph_id.core.distance_clustering_graph_id import DistanceClusteringGraphID

gen = DistanceClusteringGraphID(
    rank_k=3,    # Number of distance clusters to consider
    cutoff=6.0   # Maximum distance cutoff (Å)
)

graph_id = gen.get_id(structure)
```

## Traversal Depth Control

### Dynamic Depth (Default)

Depth is calculated from the graph diameter:

```
depth = diameter_factor × graph_diameter + additional_depth
```

```python
from graph_id_cpp import GraphIDGenerator

gen = GraphIDGenerator(
    diameter_factor=2,     # Multiplier for diameter
    additional_depth=1     # Extra depth added
)
```

### Fixed Depth

Override dynamic calculation with a fixed value:

```python
from graph_id import GraphIDMaker

# Use exactly depth=5 for all structures
maker = GraphIDMaker(depth=5)
```


## Topology-only Mode

Generate IDs based purely on connectivity, ignoring element types:

```python
from graph_id_cpp import GraphIDGenerator

gen = GraphIDGenerator(topology_only=True)
topo_id = gen.get_id(structure)
```

**Applications:**

- Finding isostructural compounds
- Structure type classification
- Comparing across element substitutions


## Wyckoff Position Mode

Include crystallographic symmetry information:

```python
from graph_id_cpp import GraphIDGenerator

gen = GraphIDGenerator(
    wyckoff=True,
    symmetry_tol=0.1  # Tolerance for symmetry detection
)

wyckoff_id = gen.get_id(structure)
```

The ID incorporates:
- Space group number
- Wyckoff letter for each site


## Site Sequence Reduction

For structures with repeated motifs, reduce the sequence:

```python
from graph_id import GraphIDMaker

maker = GraphIDMaker(reduce=True, depth=6)
```

This divides repeated site patterns.
It is useful for:

- Supercell comparisons
- Detecting structural relationships
