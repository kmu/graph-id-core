# Installation

This guide covers different ways to install Graph ID on your system.

## Requirements

- Python 3.10 or higher
- A C++ compiler (for building from source)

## Install from PyPI

The simplest way to install Graph ID is from PyPI:

```bash
pip install graph-id-core
```

### Optional: Database Component

If you want to search structures in the Materials Project database:

```bash
pip install graph-id-db
```

## Install from Source

For development or to get the latest features:

### 1. Clone the Repository

```bash
git clone https://github.com/kmu/graph-id-core.git
cd graph-id-core
```

### 2. Initialize Git Submodules

The C++ build requires several submodules (pybind11, Eigen, GTL):

```bash
git submodule update --init --recursive
```

!!! warning "Required Step"
    Without initializing submodules, the installation will fail during the CMake build step.

### 3. Install the Package

=== "Using pip"

    ```bash
    pip install -e .
    ```

=== "Using Poetry"

    ```bash
    poetry install
    ```

## Next Steps

- [Quick Start Guide](quickstart.md) - Learn the basics
- [API Reference](../api/index.md) - Detailed API documentation
