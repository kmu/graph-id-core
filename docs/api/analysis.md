# Analysis Module

The analysis module provides low-level utilities for structure graph construction and compositional sequence computation.

## StructureGraph

::: graph_id.analysis.graphs.StructureGraph
    options:
      show_root_heading: false
      members:
        - with_local_env_strategy
        - with_indivisual_state_comp_strategy
        - set_elemental_labels
        - set_wyckoffs
        - set_compositional_sequence_node_attr
        - get_loops
        - set_loops

Extended version of pymatgen's `StructureGraph` with additional methods for Graph ID generation.

### Import

```python
from graph_id.analysis.graphs import StructureGraph
```

### Class Methods

#### with_local_env_strategy

```python
@staticmethod
def with_local_env_strategy(structure, strategy, weights=False)
```

Constructor for StructureGraph using a neighbor detection strategy.

**Parameters:**

- `structure` (`Structure`): pymatgen Structure object
- `strategy` (`NearNeighbors`): A neighbor detection strategy
- `weights` (`bool`): If True, use weights from the strategy

**Returns:**

- `StructureGraph`: Constructed structure graph

**Example:**

```python
from graph_id.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN

sg = StructureGraph.with_local_env_strategy(structure, MinimumDistanceNN())
```

### Instance Methods

#### set_elemental_labels

```python
def set_elemental_labels(self)
```

Set elemental species strings as starting labels for compositional sequence computation.

#### set_wyckoffs

```python
def set_wyckoffs(self, symmetry_tol: float = 0.01)
```

Set Wyckoff position labels for each site.

**Parameters:**

- `symmetry_tol` (`float`): Tolerance for symmetry detection

#### set_compositional_sequence_node_attr

```python
def set_compositional_sequence_node_attr(
    self,
    hash_cs: bool = False,
    wyckoff: bool = False,
    additional_depth: int = 0,
    diameter_factor: int = 2,
    use_previous_cs: bool = False
)
```

Compute and set compositional sequences as node attributes.

**Parameters:**

- `hash_cs` (`bool`): Hash the compositional sequence during computation
- `wyckoff` (`bool`): Use Wyckoff-labeled sequences
- `additional_depth` (`int`): Extra traversal depth
- `diameter_factor` (`int`): Multiplier for graph diameter
- `use_previous_cs` (`bool`): Use previous CS as starting point

#### get_loops

```python
def get_loops(self, depth: int, index: int, shortest: bool = True)
```

Compute loops/rings starting from a given atom.

**Parameters:**

- `depth` (`int`): Maximum loop size to search
- `index` (`int`): Starting atom index
- `shortest` (`bool`): Stop when all theoretical shortest loops are found

**Returns:**

- `list`: List of loops, each as a list of (index, image) tuples

---

## CompositionalSequence

::: graph_id.analysis.compositional_sequence.CompositionalSequence
    options:
      show_root_heading: false
      members:
        - __init__
        - count_composition_for_neighbors
        - finalize_this_depth
        - get_current_starting_sites

Class for computing compositional sequences around an atom.

### Import

```python
from graph_id.analysis.compositional_sequence import CompositionalSequence
```

### Constructor

```python
CompositionalSequence(
    focused_site_i,
    starting_labels,
    hash_cs=False,
    use_previous_cs=False
)
```

**Parameters:**

- `focused_site_i` (`int`): Index of the central atom
- `starting_labels` (`list[str]`): Labels for each site
- `hash_cs` (`bool`): Hash sequences incrementally
- `use_previous_cs` (`bool`): Use previous sequence as labels

### Methods

#### count_composition_for_neighbors

```python
def count_composition_for_neighbors(self, nsites)
```

Count the composition of neighboring sites.

#### finalize_this_depth

```python
def finalize_this_depth(self)
```

Finalize counting for the current depth level.

### String Representation

The string representation gives the full compositional sequence:

```python
cs = CompositionalSequence(0, labels)
# ... compute neighbors ...
print(str(cs))  # "Na-Cl6-Na12-..."
```

---

## DistanceClusteringNN

::: graph_id.analysis.local_env.DistanceClusteringNN
    options:
      show_root_heading: false
      members:
        - __init__
        - get_nn_info
        - get_cutoff_cluster
        - structures_allowed

Neighbor detection based on DBSCAN clustering of interatomic distances.

### Import

```python
from graph_id.analysis.local_env import DistanceClusteringNN
```

### Constructor

```python
DistanceClusteringNN()
```

### Methods

#### get_nn_info

```python
def get_nn_info(
    self,
    structure: Structure,
    n: int,
    rank_k: int,
    cutoff: float = 6.0
) -> list[dict]
```

Get neighbor information for a specific site and distance cluster.

**Parameters:**

- `structure` (`Structure`): Input structure
- `n` (`int`): Site index
- `rank_k` (`int`): Cluster index (0-based)
- `cutoff` (`float`): Maximum distance cutoff

**Returns:**

- `list[dict]`: List of neighbor information dictionaries

#### get_cutoff_cluster

```python
def get_cutoff_cluster(
    self,
    structure: Structure,
    n: int,
    cutoff: float = 6.0
) -> list
```

Get distance cutoffs for each cluster using DBSCAN.

**Parameters:**

- `structure` (`Structure`): Input structure
- `n` (`int`): Site index
- `cutoff` (`float`): Maximum distance to consider

**Returns:**

- `list`: Sorted list of maximum distances for each cluster

### How DBSCAN Clustering Works

The algorithm:

1. Computes all pairwise distances within the cutoff
2. Runs DBSCAN with `eps=0.5` and `min_samples=2`
3. Groups distances into clusters
4. Returns cutoffs as the maximum distance in each cluster

This is useful for structures with distinct bond length populations.
