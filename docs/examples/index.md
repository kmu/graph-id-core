# Examples

This section provides practical examples of using Graph ID in various scenarios.

## Available Examples

<div class="grid cards" markdown>

- :material-brain:{ .lg .middle } **Machine Learning**

    ---

    Use Graph ID features for materials property prediction

    [:octicons-arrow-right-24: Machine Learning](machine-learning.md)

</div>

## Quick Examples

### Deduplicating a CIF Collection

```python
from pathlib import Path
from pymatgen.core import Structure
from graph_id import GraphIDMaker

maker = GraphIDMaker()
seen_ids = {}

for cif_file in Path("structures/").glob("*.cif"):
    structure = Structure.from_file(cif_file)
    graph_id = maker.get_id(structure)

    if graph_id not in seen_ids:
        seen_ids[graph_id] = cif_file
        print(f"New: {cif_file.name} -> {graph_id}")
    else:
        print(f"Duplicate: {cif_file.name} matches {seen_ids[graph_id].name}")
```

### Finding Isostructural Materials

```python
from graph_id_cpp import GraphIDGenerator

# Topology-only mode ignores element types
topo_gen = GraphIDGenerator(topology_only=True)

structures = load_structures()  # Your loading function

# Group by topology
topology_groups = {}
for name, struct in structures.items():
    topo_id = topo_gen.get_id(struct)
    if topo_id not in topology_groups:
        topology_groups[topo_id] = []
    topology_groups[topo_id].append(name)

# Print groups with multiple members
for topo_id, members in topology_groups.items():
    if len(members) > 1:
        print(f"Isostructural group ({topo_id[:8]}...): {members}")
```

### Comparing Database Entries

```python
from graph_id import GraphIDMaker
from graph_id_db import Finder

maker = GraphIDMaker()
finder = Finder()

# Your structure
my_structure = Structure.from_file("my_material.cif")
my_id = maker.get_id(my_structure)

# Search database
matches = finder.find(my_id)

if matches:
    print(f"Found in Materials Project: {matches}")
else:
    print("Novel structure - not in database!")
```

## Jupyter Notebooks

The repository includes Jupyter notebook examples:

- [`examples/Search_Structures_from_Database.ipynb`](https://github.com/kmu/graph-id-core/blob/main/examples/Search_Structures_from_Database.ipynb) - Database search tutorial

## Test Files

The test suite in `tests/` provides additional usage examples:

- `tests/cpp/test_graph_id.py` - Basic functionality tests
- `tests/cpp/test_near_neighbor.py` - Neighbor detection tests
- `tests/cpp/test_maker.py` - GraphIDMaker tests
