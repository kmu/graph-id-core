# GraphIDMaker

The recommended high-level interface for generating Graph IDs.

## API Reference

::: graph_id.app.maker.GraphIDMaker
    options:
      show_root_heading: true
      members:
        - __init__
        - get_id
        - get_site_ids
        - get_id_reducing_site_sequences

## Quick Example

```python
from graph_id import GraphIDMaker
from pymatgen.core import Structure

maker = GraphIDMaker()
structure = Structure.from_file("material.cif")
graph_id = maker.get_id(structure)
```

## Engine Selection

=== "C++ (Default)"

    ```python
    maker = GraphIDMaker(engine="c++")  # Fast, recommended
    ```

=== "Python"

    ```python
    maker = GraphIDMaker(engine="py")  # For debugging
    ```

!!! warning "Neighbor Compatibility"
    C++ engine requires C++ neighbor objects:
    ```python
    from graph_id_cpp import CrystalNN
    maker = GraphIDMaker(nn=CrystalNN(), engine="c++")
    ```

## See Also

- [GraphIDGenerator](generator.md) - Lower-level with more options
- [Basic Usage](../user-guide/basic-usage.md) - Usage guide
