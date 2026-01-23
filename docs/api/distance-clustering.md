# DistanceClusteringGraphID

Specialized Graph ID generator using DBSCAN clustering on interatomic distances.

## API Reference

::: graph_id.core.distance_clustering_graph_id.DistanceClusteringGraphID
    options:
      show_root_heading: true
      members:
        - __init__
        - get_id
        - prepare_structure_graph

## Quick Example

```python
from graph_id.core.distance_clustering_graph_id import DistanceClusteringGraphID

gen = DistanceClusteringGraphID()
graph_id = gen.get_id(structure)
```

## When to Use

Use this variant when:

- Standard neighbor detection gives unexpected results
- Structure has multiple distinct bond length scales
- Working with MOFs, zeolites, or complex frameworks

## Configuration Examples

```python
# For complex structures (MOFs, zeolites)
gen = DistanceClusteringGraphID(
    rank_k=5,      # More distance clusters
    cutoff=10.0    # Larger search radius
)
```

!!! warning "Performance"
    Distance clustering is slower than standard `GraphIDGenerator`.
    Use only when standard methods don't work.

## See Also

- [GraphIDGenerator](generator.md) - Standard generator (faster)
- [Advanced Configuration](../user-guide/advanced-configuration.md)
