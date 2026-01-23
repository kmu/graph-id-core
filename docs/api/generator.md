# GraphIDGenerator

Core Python implementation with full parameter control.

## API Reference

::: graph_id.core.graph_id.GraphIDGenerator
    options:
      show_root_heading: true
      members:
        - __init__
        - get_id
        - get_id_catch_error
        - get_many_ids
        - are_same
        - get_unique_structures
        - get_component_ids
        - prepare_structure_graph
        - version

## Quick Example

```python
from graph_id.core.graph_id import GraphIDGenerator

gen = GraphIDGenerator()
graph_id = gen.get_id(structure)
```

## Common Configurations

### Topology-Only Mode

```python
gen = GraphIDGenerator(topology_only=True)
# Output: "3D-a1b2c3d4e5f6g7h8" (no formula)
```

### With Wyckoff Positions

```python
gen = GraphIDGenerator(wyckoff=True, symmetry_tol=0.1)
```

### Minimal Output

```python
gen = GraphIDGenerator(
    prepend_composition=False,
    prepend_dimensionality=False
)
# Output: "a1b2c3d4e5f6g7h8" (hash only)
```

## Output Format

| Configuration | Format |
|--------------|--------|
| Default | `{formula}-{dim}D-{hash}` |
| `prepend_dimensionality=False` | `{formula}-{hash}` |
| `prepend_composition=False` | `{dim}D-{hash}` |
| Both `False` | `{hash}` |

!!! warning "Invalid Combinations"
    - `wyckoff=True` + `loop=True` → ValueError
    - `loop=True` + `topology_only=True` → ValueError

## See Also

- [GraphIDMaker](maker.md) - Simplified interface
- [DistanceClusteringGraphID](distance-clustering.md) - Clustering variant
