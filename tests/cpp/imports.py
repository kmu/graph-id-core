import os
import sys

"""
テスト用に C++ 拡張をインポートする
C++ 拡張のビルド方法は
```
mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make
```
"""

path = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../build"))
sys.path.insert(0, path)
import graph_id  # noqa: E402,F401
