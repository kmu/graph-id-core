# Graph ID
GraphID: universal graph-based identifiers of chemical structures for linking large material databases.

## Installation 

### From pypi

```
pip install graph-id-core
pip install graph-id-db # optional
```

### From GitHub

```
git clone git+https://github.com/kmu/graph-id-core.git
git submodule init
git submodule update
pip install -e .
```

## Usage

```
from pymatgen.core import Structure
from graph_id import GraphIDGenerator
from graph_id_cpp import GraphIDGenerator as CppGraphIDGenerator

structure = Structure.from_file(path/to/your/cif)
gen = GraphIDGenerator()
cpp_gen = CppGraphIDGenerator()

# GraphID implemented in Python
print(gen.get_id(structure))
# GraphID implemented in C++
print(cpp_gen.get_id(structure))