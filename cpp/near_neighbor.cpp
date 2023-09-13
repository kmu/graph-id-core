#include "near_neighbor.h"
#include <vector>
#include <pybind11/stl.h>
#include "wrapper.h"


py::list NearNeighbor::get_all_nn_info(py::object &structure) {
    auto s = structure.cast<PymatgenStructure>();
    const auto result = this->get_all_nn_info_cpp(Structure(s));
    py::list arr;
    for (const auto &infos: result) {
        py::list inner;
        for (const auto &info: infos) {
            py::dict d;
            d["site"] = s.sites()[info.site_index].obj;
            d["site_index"] = info.site_index;
            d["weight"] = info.weight;
            d["image"] = info.image;
            inner.append(d);
        }
        arr.append(inner);
    }
    return arr;
}


std::vector<std::vector<NearNeighborInfo>> MinimumDistanceNN::get_all_nn_info_cpp(const Structure &structure) {
    // TODO: implement this
    return std::vector<std::vector<NearNeighborInfo>>(structure.sites.size());
}


void init_near_neighbor(pybind11::module &m) {
    py::class_<NearNeighbor>(m, "NearNeighbor")
            .def_property_readonly("structures_allowed", &NearNeighbor::structures_allowed)
            .def_property_readonly("molecules_allowed", &NearNeighbor::molecules_allowed)
            .def("get_all_nn_info", &NearNeighbor::get_all_nn_info);

    py::class_<MinimumDistanceNN, NearNeighbor>(m, "MinimumDistanceNN")
            .def(py::init<double, double, double>(),
                 py::arg("tol") = 0.1,
                 py::arg("cutoff") = 10.0,
                 py::arg("get_all_sites") = false);
}

