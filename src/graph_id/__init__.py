import importlib
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError  # for Python <3.8

try:
    __version__ = version("graph-id-core")
except PackageNotFoundError:
    __version__ = "unknown"

def __getattr__(name):
    cpp = importlib.import_module(__name__ + ".graph_id_cpp")
    if hasattr(cpp, name):
        return getattr(cpp, name)
    raise AttributeError(f"module {__name__} has no attribute '{name}'")
