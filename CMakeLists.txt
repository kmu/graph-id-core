cmake_minimum_required(VERSION 3.5...3.26)
project(graph_id LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(library/pybind11)
pybind11_add_module(graph_id_cpp
        cpp/graph_id.cpp
        cpp/graph_id.h
        cpp/main.cpp
        cpp/near_neighbor.cpp
        cpp/near_neighbor.h
        cpp/structure_graph.cpp
        cpp/structure_graph.h
        cpp/core.cpp
        cpp/core.h
)
target_compile_definitions(graph_id_cpp
        PRIVATE VERSION_INFO=${VERSION_INFO})
target_include_directories(graph_id_cpp
        PRIVATE library/eigen library/gtl/include
)
target_compile_options(graph_id_cpp  PRIVATE
        # 警告をできるだけ出力する
        $<$<CXX_COMPILER_ID:MSVC>:/W4 /Wall>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

# リリースビルド時に最適化を有効にする
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(graph_id_cpp PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/O2 /Ob2 /Oi /Ot /Oy /GL /Gw /GF /GS- /Gy /fp:fast /arch:AVX2>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3 -ffast-math>
    )
endif ()
set(CMAKE_VERBOSE_MAKEFILE on) # TODO: 後で消す
