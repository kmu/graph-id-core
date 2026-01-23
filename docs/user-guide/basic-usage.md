# Basic Usage

## GraphIDMaker vs GraphIDGenerator

| Class | Use Case |
|-------|----------|
| `GraphIDMaker` | Simple, high-level interface (recommended) |
| `GraphIDGenerator` | Full control, batch processing, deduplication |

## Loading Structures

```python
from pymatgen.core import Structure

# From CIF
structure = Structure.from_file("material.cif")

# From POSCAR/VASP
structure = Structure.from_file("POSCAR")

# From Materials Project
from mp_api.client import MPRester
with MPRester("YOUR_API_KEY") as mpr:
    docs = mpr.materials.search(formula="TiO2", fields=["structure"])
    structure = docs[0].structure

# From ASE
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
atoms = read("structure.xyz")
structure = AseAtomsAdaptor.get_structure(atoms)
```

## Generating IDs

### Single Structure

```python
from graph_id import GraphIDMaker

maker = GraphIDMaker()
graph_id = maker.get_id(structure)
# Output: "Fe2O3-a1b2c3d4e5f6g7h8"
```

### Batch Processing

```python
from graph_id.core.graph_id import GraphIDGenerator

generator = GraphIDGenerator()

# Sequential
ids = [generator.get_id_catch_error(s) for s in structures]

# Parallel (with progress bar)
ids = generator.get_many_ids(structures, parallel=True)
```

## Site-Level Information

```python
maker = GraphIDMaker()
site_ids = maker.get_site_ids(structure)

for idx, cs in sorted(site_ids.items()):
    element = structure[idx].species_string
    print(f"Site {idx} ({element}): {cs}")
```

## Deduplication

```python
from graph_id.core.graph_id import GraphIDGenerator

generator = GraphIDGenerator()
unique = generator.get_unique_structures(all_structures)
print(f"Reduced {len(all_structures)} → {len(unique)} unique")
```

## Output Format

| Setting | Output |
|---------|--------|
| Default | `NaCl-3D-88c8e156db1b0fd9` |
| `prepend_dimensionality=False` | `NaCl-88c8e156db1b0fd9` |
| `prepend_composition=False` | `3D-88c8e156db1b0fd9` |
| `topology_only=True` | `3D-88c8e156db1b0fd9` |

## Next Steps

- [Advanced Configuration](advanced-configuration.md) - Fine-tune parameters
- [API Reference](../api/index.md) - Complete documentation
