from glob import glob
from pybind11.setup_helpers import Pybind11Extension, build_ext


def build(setup_kwargs):
    ext_modules = [
        Pybind11Extension(
            "graph_id_cpp",
            sorted(glob("cpp/*.cpp")),  # Sort source files for reproducibility
            cxx_std=20,
        ),
    ]
    setup_kwargs.update({
        "ext_modules": ext_modules,
        "cmd_class": {"build_ext": build_ext},
        "zip_safe": False,
    })
