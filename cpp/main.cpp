#include <pybind11/pybind11.h>
#include "iostream"

namespace py = pybind11;

PYBIND11_MODULE(graph_id_cpp, m) {
    m.attr("test") = "true";
}
