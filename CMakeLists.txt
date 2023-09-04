cmake_minimum_required(VERSION 3.5...3.26)
project(graph_id LANGUAGES CXX)

add_subdirectory(pybind11)
pybind11_add_module(graph_id_cpp cpp/main.cpp)
