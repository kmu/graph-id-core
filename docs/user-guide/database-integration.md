# Database Integration

Graph ID is designed for efficient database operations. This guide covers integration with structure databases.

## Using graph-id-db

The `graph-id-db` package provides precomputed Graph IDs for several databases including the Materials Project.

### Installation

```bash
pip install graph-id-db
```

### Searching for Structures

```python
from graph_id_cpp import GraphIDGenerator
from graph_id_db import Finder
from pymatgen.core import Structure, Lattice

# Create a query structure
structure = Structure.from_spacegroup(
    "Fm-3m",
    Lattice.cubic(5.692),
    ["Na", "Cl"],
    [[0, 0, 0], [0.5, 0.5, 0.5]]
).get_primitive_structure()

# Generate its Graph ID
gen = GraphIDGenerator()
graph_id = gen.get_id(structure)
print(f"Query Graph ID: {graph_id}")

# Search in the database
finder = Finder()
results = finder.find(graph_id)

for result in results:
    print(f"Found: {result}")
```

### Understanding Results

The `Finder` returns matching Materials Project IDs:

```python
results = finder.find(graph_id)

if results:
    print(f"Found {len(results)} matching structures")
    for mp_id in results:
        print(f"  - {mp_id}")
else:
    print("No matches found - this may be a new structure!")
```

## Building Your Own Database

### SQLite Integration

```python
import sqlite3
from graph_id import GraphIDMaker
from pymatgen.core import Structure

# Initialize
maker = GraphIDMaker()
conn = sqlite3.connect("structures.db")
cursor = conn.cursor()

# Create table
cursor.execute("""
    CREATE TABLE IF NOT EXISTS structures (
        id INTEGER PRIMARY KEY,
        graph_id TEXT NOT NULL,
        formula TEXT,
        file_path TEXT,
        UNIQUE(graph_id)
    )
""")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_graph_id ON structures(graph_id)")

# Insert structure
def add_structure(file_path):
    structure = Structure.from_file(file_path)
    graph_id = maker.get_id(structure)
    formula = structure.composition.reduced_formula

    try:
        cursor.execute(
            "INSERT INTO structures (graph_id, formula, file_path) VALUES (?, ?, ?)",
            (graph_id, formula, file_path)
        )
        conn.commit()
        return True
    except sqlite3.IntegrityError:
        print(f"Duplicate structure: {graph_id}")
        return False

# Search for structure
def find_structure(query_structure):
    query_id = maker.get_id(query_structure)
    cursor.execute("SELECT * FROM structures WHERE graph_id = ?", (query_id,))
    return cursor.fetchall()
```

### MongoDB Integration

```python
from pymongo import MongoClient
from graph_id import GraphIDMaker

client = MongoClient("mongodb://localhost:27017/")
db = client["materials"]
collection = db["structures"]

# Create index on graph_id
collection.create_index("graph_id", unique=True)

maker = GraphIDMaker()

def add_structure(structure, metadata=None):
    graph_id = maker.get_id(structure)

    doc = {
        "graph_id": graph_id,
        "formula": structure.composition.reduced_formula,
        "structure": structure.as_dict(),
        **(metadata or {})
    }

    try:
        collection.insert_one(doc)
        return True
    except Exception:
        return False  # Duplicate

def find_by_graph_id(query_structure):
    query_id = maker.get_id(query_structure)
    return collection.find_one({"graph_id": query_id})

def find_by_formula(formula):
    return list(collection.find({"formula": formula}))
```
