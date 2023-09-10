#include <pybind11/pybind11.h>
#include "iostream"

#define VALUE(string) #string
#define TO_LITERAL(string) VALUE(string)

namespace py = pybind11;

PYBIND11_MODULE(graph_id_cpp, m) {
#ifdef VERSION_INFO
    m.attr("__version__") = TO_LITERAL(VERSION_INFO);
#else
    m.attr("__version__") = "";
#endif

}
