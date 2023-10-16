#include "graph_id.h"
#include "structure_graph.h"

std::string GraphIDGenerator::get_id(const Structure &structure) const {
    auto s_ptr = std::shared_ptr<const Structure>(&structure, [](const Structure *) {});
    const auto sg = prepare_structure_graph(s_ptr);
    std::vector<std::string> cc_labels(sg.cc_nodes.size());
    for (size_t i = 0; i < sg.cc_nodes.size(); ++i) {
        std::vector<std::string> labels = sg.cc_cs[i];
        std::sort(labels.begin(), labels.end());
        cc_labels[i] = blake2b(join_string("-", labels));
    }
    std::sort(cc_labels.begin(), cc_labels.end());
    std::string gid = blake2b(join_string(":", cc_labels), 16);

    return elaborate_comp_dim(sg, gid);
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

std::string GraphIDGenerator::elaborate_comp_dim(const StructureGraph &sg, const std::string &gid) const {
    int dim = sg.get_dimensionality_larsen();
    if (!topology_only) {
        return sg.structure->py_structure.reduced_formula() + "-" + std::to_string(dim) + "D-" + gid;
    }
    return std::to_string(dim) + "D-" + gid;
}

std::vector<std::string> GraphIDGenerator::get_component_ids(const Structure &structure) const {
    // TODO
    return std::vector<std::string>(structure.count);
}

bool GraphIDGenerator::are_same(const Structure &structure1, const Structure &structure2) const {
    return this->get_id(structure1) == this->get_id(structure2);
}

StructureGraph GraphIDGenerator::prepare_structure_graph(std::shared_ptr<const Structure> &structure) const {
    auto sg = StructureGraph::with_local_env_strategy(structure, *this->nn);
    bool use_previous_cs = false;

    auto labels = structure->species_strings;
    auto prev_num_uniq = std::unique(labels.begin(), labels.end()) - labels.begin();

    if (wyckoff) {
        sg.set_wyckoffs_label(symmetry_tol);
    } else if (topology_only) {
        sg.labels = std::vector<std::string>(structure->count, "X");
    } else {
        sg.set_elemental_labels();
    }

    while (true) {
        sg.set_compositional_sequence_node_attr(
                true,
                wyckoff,
                additional_depth,
                depth_factor,
                use_previous_cs
        );

        labels.resize(0);
        for (const auto &v: sg.cc_cs)
            for (const auto &s: v)
                labels.emplace_back(s);
        auto new_num_uniq = std::unique(labels.begin(), labels.end()) - labels.begin();
        if (new_num_uniq == prev_num_uniq) {
            break;
        }
        prev_num_uniq = new_num_uniq;
    }

    return sg;
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
