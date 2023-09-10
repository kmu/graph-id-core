cmake_minimum_required(VERSION 3.5...3.26)
project(graph_id LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(pybind11)
pybind11_add_module(graph_id_cpp cpp/main.cpp)
target_compile_definitions(graph_id_cpp
        PRIVATE VERSION_INFO=${VERSION_INFO})
