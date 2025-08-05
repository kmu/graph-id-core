import importlib
from importlib.metadata import version, PackageNotFoundError


try:
    __version__ = version("graph-id-core")
except PackageNotFoundError:
    __version__ = "unknown"

def __getattr__(name):
    cpp = importlib.import_module(__name__ + ".graph_id_cpp")
    if hasattr(cpp, name):
        return getattr(cpp, name)
    raise AttributeError(f"module {__name__} has no attribute '{name}'")
