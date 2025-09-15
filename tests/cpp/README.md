# Tests for C++ version

## Executing tests

**Building C++ version**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . -j
cd ..
```

**Testing C++ version**
```bash
poetry run pytest
```


## Benchmark

**Building C++ code**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
cd ..
```

**Executing benchmark**
```bash
python -m unittest discover -s cpp -p "test/benchmark*.py"
```


## Profiling

**Obtainig dependencies (Only works for Linux)**
```bash
apt-get install -y google-perftools libgoogle-perftools-dev kcachegrind
pip install yep
```

**Building C++ code**
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build . -j
cd ..
```

**Obtaining profiles**
```bash
python cpp/test/profile.py
google-pprof \
  ./build/graph_id_cpp$(python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))') \
  build/profile-small.prof --svg > build/profile-small.svg
google-pprof \
  ./build/graph_id_cpp$(python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))') \
  build/profile-medium.prof --svg > build/profile-medium.svg
google-pprof \
  ./build/graph_id_cpp$(python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))') \
  build/profile-large.prof --svg > build/profile-large.svg
```

