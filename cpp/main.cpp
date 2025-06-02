#include <pybind11/pybind11.h>

#define VALUE(string) #string
#define TO_LITERAL(string) VALUE(string)

namespace py = pybind11;

void init_near_neighbor(py::module &m);
void init_core(py::module &m);
void init_structure_graph(py::module &m);
void init_graph_id(py::module &m);

PYBIND11_MODULE(graph_id_cpp, m) {
    // Initialize the module
    init_core(m);
    init_near_neighbor(m);
    init_structure_graph(m);
    init_graph_id(m);

    // Set the version defined in pyproject.toml to __version__
#ifdef VERSION_INFO
    m.attr("__version__") = TO_LITERAL(VERSION_INFO);
#else
    m.attr("__version__") = "";
#endif
}
