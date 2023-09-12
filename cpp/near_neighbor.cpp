#include "near_neighbor.h"
#include <vector>
#include <pybind11/stl.h>

std::vector <py::dict> MinimumDistanceNN::get_all_nn_info(py::object structure) {
    // TODO: Implement this method.
    return {};
}


void init_near_neighbor(pybind11::module &m) {
    py::class_<NearNeighbor>(m, "NearNeighbor")
            .def("get_all_nn_info", &NearNeighbor::get_all_nn_info);
    py::class_<MinimumDistanceNN, NearNeighbor>(m, "MinimumDistanceNN")
            .def(py::init<double, double, double>(), py::arg("tol") = 0.1, py::arg("cutoff") = 10.0, py::arg("get_all_sites") = false)
            .def("get_all_nn_info", &MinimumDistanceNN::get_all_nn_info);
}
