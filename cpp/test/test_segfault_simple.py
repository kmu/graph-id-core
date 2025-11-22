"""Simple test to check if find_near_neighbors causes segfault"""
import sys
import os
import traceback

def test_segfault_simple():
    # Try to find test files
    test_paths = [
        'tests/py/test_files/mp-36.cif',
        '../py/test_files/mp-36.cif',
        'graph_id/tests/test_files/mp-36.cif',
        '../../graph_id/tests/test_files/mp-36.cif',
    ]

    test_file = None
    for path in test_paths:
        if os.path.exists(path):
            test_file = path
            break

    if not test_file:
        print("ERROR: Could not find test file", file=sys.stderr)
        sys.exit(1)


    import numpy as np
    import graph_id_cpp
    from pymatgen.core import Structure

    # Use a simple structure
    try:
        s = Structure.from_file(test_file)
        pbc = np.ascontiguousarray(s.lattice.pbc, dtype=int)
    except Exception as e:
        print(f"Structure error: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

    print("Testing find_near_neighbors...", flush=True)
    try:
        result = graph_id_cpp.find_near_neighbors(
            s.cart_coords,
            s.cart_coords,
            1.0,
            pbc,
            s.lattice.matrix,
            1e-8,
            1.0
        )
        print("SUCCESS: No segfault")
        print(f"Result type: {type(result)}")
        assert True        
    except SystemExit:
        raise
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        traceback.print_exc()
        # sys.exit(1)
        assert False
