import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension
from setuptools.command.build_ext import build_ext
import platform
# https://github.com/pybind/cmake_example/blob/master/setup.py

# Convert distutils platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
    "darwin": "Darwin",  # macOS
    "linux": "Linux",  # Linux
    "linux-x86_64": "Linux-x86_64",
    "linux-aarch64": "Linux-aarch64",
    "darwin-x86_64": "Darwin-x86_64",
    "darwin-arm64": "Darwin-arm64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve() only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # Pass in the version to C++.
        cmake_args += [f"-DVERSION_INFO={self.distribution.get_version()}"]

        if sys.platform.startswith("win32"):
            build_args += ["--config", cfg]
        elif sys.platform.startswith("darwin") or sys.platform.startswith("linux"):
            # Use Unix Makefiles for macOS and Linux
            cmake_args += ["-G", "Unix Makefiles"]
        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]
            else:
                # Auto-detect macOS architecture
                machine = platform.machine()
                if machine == "arm64":
                    cmake_args += ["-DCMAKE_OSX_ARCHITECTURES=arm64"]
                elif machine == "x86_64":
                    cmake_args += ["-DCMAKE_OSX_ARCHITECTURES=x86_64"]
        elif sys.platform.startswith("linux"):
            # Linux architecture detection
            machine = platform.machine()
            if machine == "aarch64":
                cmake_args += ["-DCMAKE_SYSTEM_PROCESSOR=aarch64"]
            elif machine == "x86_64":
                cmake_args += ["-DCMAKE_SYSTEM_PROCESSOR=x86_64"]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)
        env = {**os.environ}
        if "PYTHONPATH" in env:
            # Google Colab など CMake が pip でインストールされている環境で、import cmake できなくなり build が失敗する
            del env["PYTHONPATH"]
        subprocess.run(
            ["cmake", ext.sourcedir, *cmake_args],
            cwd=build_temp,
            check=True,
            env=env,
        )
        subprocess.run(
            ["cmake", "--build", ".", *build_args],
            cwd=build_temp,
            check=True,
            env=env,
        )


def build(setup_kwargs):
    ext_modules = [
        CMakeExtension("graph_id_cpp", sourcedir="."),
    ]
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": CMakeBuild},
            "zip_safe": False,
        }
    )
