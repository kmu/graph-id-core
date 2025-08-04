from . import graph_id_cpp as _cpp
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError  # for Python <3.8

try:
    __version__ = version("graph-id-core")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = getattr(_cpp, '__all__', [name for name in dir(_cpp) if not name.startswith('_')])
globals().update({name: getattr(_cpp, name) for name in __all__})

