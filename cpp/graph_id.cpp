#include "graph_id.h"
#include "structure_graph.h"

std::string GraphIDGenerator::get_id(const Structure &structure) const {
    const auto sg = StructureGraph::with_local_env_strategy(std::make_shared<Structure>(structure), *this->nn);
    return "";
}

std::string GraphIDGenerator::get_id_catch_error(const Structure &structure) const noexcept {
    try {
        return this->get_id(structure);
    } catch (...) {
        return "";
    }
}

std::vector<std::string> GraphIDGenerator::get_many_ids(const std::vector<Structure> &structures) const {
    // マルチスレッド化するときはここを変更する
    std::vector<std::string> ids;
    ids.reserve(structures.size());
    for (const auto &structure: structures) {
        this->get_id_catch_error(structure);
    }
    return ids;
}

std::vector<std::string> GraphIDGenerator::get_component_ids(const Structure &structure) const {
    // TODO
    return std::vector<std::string>(structure.count);
}

bool GraphIDGenerator::are_same(const Structure &structure1, const Structure &structure2) const {
    return this->get_id(structure1) == this->get_id(structure2);
}

StructureGraph GraphIDGenerator::prepare_structure_graph(std::shared_ptr<const Structure> &structure) const {
    return StructureGraph::with_local_env_strategy(structure, *this->nn);
}

void init_graph_id(pybind11::module &m) {
    py::class_<GraphIDGenerator>(m, "GraphIDGenerator")
            .def(py::init<std::shared_ptr<NearNeighbor>, bool, int, int, double, bool>(),
                 py::arg("nn") = nullptr,
                 py::arg("wyckoff") = false,
                 py::arg("depth_factor") = 2,
                 py::arg("additional_depth") = 1,
                 py::arg("symmetry_tol") = 0.1,
                 py::arg("topology_only") = false)
            .def("get_id", [](const GraphIDGenerator &gig, py::object &structure) {
                auto s = structure.cast<PymatgenStructure>();
                return gig.get_id(Structure(s));
            })
            .def("get_id_catch_error", [](const GraphIDGenerator &gig, py::object &structure) {
                auto s = structure.cast<PymatgenStructure>();
                return gig.get_id_catch_error(Structure(s));
            })
            .def("get_many_ids", [](const GraphIDGenerator &gig, py::list &structures) {
                std::vector<Structure> ss;
                ss.reserve(structures.size());
                for (const auto &structure: structures) {
                    ss.emplace_back(structure.cast<PymatgenStructure>());
                }
                return gig.get_many_ids(ss);
            })
            .def("get_component_ids", [](const GraphIDGenerator &gig, py::object &structure) {
                auto s = structure.cast<PymatgenStructure>();
                return gig.get_component_ids(Structure(s));
            })
            .def("are_same", [](const GraphIDGenerator &gig, py::object &structure1, py::object &structure2) {
                auto s1 = structure1.cast<PymatgenStructure>();
                auto s2 = structure2.cast<PymatgenStructure>();
                return gig.are_same(Structure(s1), Structure(s2));
            });
}
