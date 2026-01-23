# Contributing

Thank you for your interest in contributing to Graph ID!

## Development Setup

### Prerequisites

- Python 3.10+
- C++ compiler (for building the extension)
- Git

### Installation

1. Clone the repository:

```bash
git clone https://github.com/kmu/graph-id-core.git
cd graph-id-core
```

2. Initialize git submodules:

```bash
git submodule update --init --recursive
```

3. Install with Poetry:

```bash
poetry install
```

4. Install pre-commit hooks:

```bash
pre-commit install
```

## Running Tests

```bash
poetry run pytest
```

### After C++ Changes

If you modify C++ code, reinstall the package:

```bash
poetry run pip install -e . --force-reinstall
```

Then run tests again.

## Code Style

This project uses:

- **Ruff** for Python linting
- **Black** for Python formatting
- **isort** for import sorting

Pre-commit hooks will automatically check these.

## Pull Request Process

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes
4. Run tests: `poetry run pytest`
5. Commit with a descriptive message
6. Push to your fork
7. Open a Pull Request

## Reporting Issues

When reporting issues, please include:

- Python version
- Operating system
- Minimal reproducible example
- Full error traceback

## Architecture Overview

```
graph-id-core/
├── cpp/                    # C++ source code
│   ├── graph_id.cpp       # Main C++ implementation
│   └── ...
├── graph_id/              # Python package
│   ├── app/               # High-level interfaces
│   ├── core/              # Core implementations
│   └── analysis/          # Analysis utilities
├── library/               # Git submodules
│   ├── pybind11/          # Python-C++ bindings
│   ├── eigen/             # Linear algebra
│   └── gtl/               # Graph template library
└── tests/                 # Test suite
```

## Questions?

Feel free to open an issue for questions or discussion.
