[![PyPI version](https://img.shields.io/pypi/v/graph-id-core.svg)](https://pypi.org/project/graph-id-core/)
[![Python versions](https://img.shields.io/pypi/pyversions/graph-id-core.svg)](https://pypi.org/project/graph-id-core/)
[![codecov](https://codecov.io/gh/kmu/graph-id-core/graph/badge.svg?token=AE2JIT3BAX)](https://codecov.io/gh/kmu/graph-id-core)

# Graph ID

Graph ID is a universal identifier system for atomistic structures including crystals and molecules. It generates unique, deterministic identifiers based on the topological and compositional properties of atomic structures, enabling efficient structure comparison, database indexing, and materials discovery.

## Overview

Graph ID works by:
1. Converting atomic structures into graph representations where atoms are nodes and bonds are edges
2. Analyzing the local chemical environment around each atom using compositional sequences
3. Computing a hash-based identifier that captures both topology and composition
4. Supporting various modes including topology-only comparisons and Wyckoff position analysis


## Features

- **Universal Structure Identification**: Generate unique IDs for any crystal or molecular structure
- **Topological Analysis**: Option to generate topology-only IDs for structure type comparison
- **Wyckoff Position Support**: Include crystallographic symmetry information in ID generation
- **Distance Clustering**: Advanced clustering-based analysis for complex structures
- **C++ Performance**: High-performance C++ backend with Python bindings
- **Multiple Neighbor Detection**: Support for various neighbor-finding algorithms (MinimumDistanceNN, CrystalNN, etc.)

## Installation

### From PyPI

```bash
pip install graph-id-core
pip install graph-id-db  # optional database component
```

### From Source

```bash
git clone https://github.com/kmu/graph-id-core.git
cd graph-id-core
git submodule update --init --recursive
pip install -e .
```

## Quick Start

### Basic Usage

```python
from pymatgen.core import Structure, Lattice
from graph_id import GraphIDMaker

# Create a structure (NaCl)
structure = Structure.from_spacegroup(
    "Fm-3m",
    Lattice.cubic(5.692),
    ["Na", "Cl"],
    [[0, 0, 0], [0.5, 0.5, 0.5]]
)

# Generate Graph ID
maker = GraphIDMaker()
graph_id = maker.get_id(structure)
print(graph_id)  # Output: NaCl-88c8e156db1b0fd9
```

### Loading from Files

```python
from pymatgen.core import Structure
from graph_id_cpp import GraphIDGenerator

# Load structure from file
structure = Structure.from_file("path/to/structure.cif")
generator = GraphIDGenerator()
graph_id = generator.get_id(structure)
```

### Advanced Configuration

```python
from graph_id_cpp import GraphIDGenerator
from pymatgen.analysis.local_env import CrystalNN

# Topology-only comparison (ignores composition)
topo_gen = GraphIDGenerator(topology_only=True)
topo_id = topo_gen.get_id(structure)

# Include Wyckoff positions
wyckoff_gen = GraphIDGenerator(wyckoff=True)
wyckoff_id = wyckoff_gen.get_id(structure)

# Use different neighbor detection
crystal_gen = GraphIDGenerator(nn=CrystalNN())  # Faster CrystalNN using C++ is also available
crystal_id = crystal_gen.get_id(structure)

```

### Search Structures from Database

Use `graph-id-db` to search structures in the Materials Project using precomputed Graph ID stored in `graph-id-db`

```python
# pip install graph-id-db
from graph_id_cpp import GraphIDGenerator

from pymatgen.core import Structure, Lattice

structure = Structure.from_spacegroup(
    "Fm-3m",
    Lattice.cubic(5.692),
    ["Na", "Cl"],
    [[0, 0, 0], [0.5, 0.5, 0.5]]
).get_primitive_structure()
gen = GraphIDGenerator()
graph_id = gen.get_id(structure)
print(f"Graph ID of NaCl is {graph_id}")

from graph_id_db import Finder

# Search for structures in graph-id-db using GraphID
finder = Finder()
finder.find(graph_id)
```

## Examples

More comprehensive examples can be found in the [`tests/`](tests/) and [`examples/`](examples/) directories.

## Applications

Graph ID is particularly useful for:

- **Materials Databases**: Efficient indexing and deduplication of structure databases
- **High-throughput Screening**: Rapid identification of unique structures in computational workflows
- **Polymorph Identification**: Distinguishing between different polymorphs of the same composition

## Web Service (experimental)

You can search materials using Graph ID at [matfinder.net](https://matfinder.net).

## Developer's notes

This repo is managed by `poetry`.

### Installation

1. Clone the repository:
```bash
git clone https://github.com/kmu/graph-id-core.git
cd graph-id-core
```

2. Initialize git submodules (required for the C++ build):
```bash
git submodule update --init --recursive
```

3. Install the package and dependencies using Poetry:
```bash
poetry install
```

4. Install `pre-commit`
```bash
pre-commit install
```

**Note:** The git submodules (`library/pybind11`, `library/eigen`, `library/gtl`) are required for building the C++ extension. Without them, the installation will fail during the CMake build step.

### Testing

```
poetry run pytest
```

If you have made changes to the C++ code, run `poetry run pip install -e --force-reinstall` to apply the changes before running the tests.

### Releasing

- Bump version in `pyproject.toml`.
- Create a new PR from `main` branch to `release` branch.
