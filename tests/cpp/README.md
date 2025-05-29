# テスト方法

## テストの実行

**C++ コードのビルド**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . -j
cd ..
```

**テストの実行**
```bash
poetry run pytest
```


## ベンチマークの実行

**C++ コードのビルド**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
cd ..
```

**ベンチマークのの実行**
```bash
python -m unittest discover -s cpp -p "test/benchmark*.py"
```


## プロファイルの取得

**必要なライブラリの取得 (Linux のみ対応)**
```bash
apt-get install -y google-perftools libgoogle-perftools-dev kcachegrind
pip install yep
```

**C++ コードのビルド**
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build . -j
cd ..
```

**プロファイルの取得**
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

