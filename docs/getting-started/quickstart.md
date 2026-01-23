# Quick Start

Get started with Graph ID in under a minute.

## Generate Your First ID

```python
from pymatgen.core import Structure, Lattice
from graph_id import GraphIDMaker

# Create a structure (NaCl rock salt)
structure = Structure.from_spacegroup(
    "Fm-3m",
    Lattice.cubic(5.692),
    ["Na", "Cl"],
    [[0, 0, 0], [0.5, 0.5, 0.5]]
)

# Generate Graph ID
maker = GraphIDMaker()
graph_id = maker.get_id(structure)
print(graph_id)  # NaCl-88c8e156db1b0fd9
```

## Compare Structures

```python
id1 = maker.get_id(structure1)
id2 = maker.get_id(structure2)

if id1 == id2:
    print("Structures are topologically identical!")
```

## Find Unique Structures

```python
from graph_id.core.graph_id import GraphIDGenerator

generator = GraphIDGenerator()
unique = generator.get_unique_structures(all_structures)
```

## Next Steps

- [Basic Usage](../user-guide/basic-usage.md) - Detailed guide
- [API Reference](../api/index.md) - Complete documentation
