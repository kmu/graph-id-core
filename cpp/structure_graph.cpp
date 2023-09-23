#include "structure_graph.h"

std::string blake2b(const std::string &s) {
    py::object hashlib = py::module_::import("hashlib");
    return hashlib("blake2b")(py::bytes(s)).attr("hexdigest")().cast<std::string>();
}

int graph_diameter(const std::vector<std::vector<NearNeighborInfo>> &graph) {
    // TODO
    return 0;
}

StructureGraph StructureGraph::with_local_env_strategy(
        const std::shared_ptr<const Structure> &structure,
        const NearNeighbor &strategy
) {
    auto sg = StructureGraph::with_empty_graph(structure);
    sg.graph = strategy.get_all_nn_info_cpp(*structure);
    sg.graph_diameter = ::graph_diameter(sg.graph);
    return sg;
}

StructureGraph StructureGraph::with_empty_graph(const std::shared_ptr<const Structure> &structure) {
    const auto n = structure->count;
    return StructureGraph{
            structure,
            std::vector<std::vector<NearNeighborInfo>>(n),
            std::vector<std::string>(n),
            0,
    };
}

void StructureGraph::set_elemental_labels() {
    this->labels = this->structure->species_strings;
}

void StructureGraph::set_wyckoffs_label(double symmetry_tol) {
    auto core = py::module_::import("pymatgen.core");
    auto symmetry = py::module_::import("pymatgen.symmetry.analyzer");
    auto Element = core.attr("Element");
    auto SpacegroupAnalyzer = symmetry.attr("SpacegroupAnalyzer");

    py::object siteless = this->structure->py_structure.copy().obj;
    for (int i = 0; i < this->structure->count; i++) {
        siteless.attr("replace")(i, Element("H"));
    }

    auto sga = SpacegroupAnalyzer(siteless);
    auto sym_dataset = sga.attr("get_symmetry_dataset")();
    if (sym_dataset.is_none()) {
        this->set_elemental_labels();
        return;
    }

    auto wyckoffs = sym_dataset["wyckoffs"];
    auto number = sym_dataset["number"];

    for (size_t site_i = 0; site_i < py::len(wyckoffs); site_i++) {
        auto wyckoff = wyckoffs[py::int_(site_i)];
        this->labels[site_i] = py::str("{}_{}_{}").format(this->structure->species_strings[site_i], wyckoff, number);
    }
}

void StructureGraph::set_compositional_sequence_node_attr(
        bool hash_cs,
        bool wyckoff,
        int additional_depth,
        int depth_factor
) {

}

void init_structure_graph(pybind11::module &m) {
    m.def("graph_diameter", &graph_diameter);

    py::class_<StructureGraph>(m, "StructureGraph")
            .def_static("with_local_env_strategy", [](PymatgenStructure &s, NearNeighbor &nn) {
                return StructureGraph::with_local_env_strategy(std::make_shared<Structure>(s), nn);
            })
            .def_static("with_empty_graph", [](PymatgenStructure &s) {
                return StructureGraph::with_empty_graph(std::make_shared<Structure>(s));
            })
            .def("set_elemental_labels", &StructureGraph::set_elemental_labels)
            .def("set_wyckoffs", &StructureGraph::set_wyckoffs_label, py::arg("symmetry_tol") = 0.1) // 互換性
            .def("set_wyckoffs_label", &StructureGraph::set_wyckoffs_label)
            .def("set_compositional_sequence_node_attr", &StructureGraph::set_compositional_sequence_node_attr)
            .def_property("labels", [](const StructureGraph &sg) { return sg.labels; },
                          [](StructureGraph &sg, const std::vector<std::string> &labels) { sg.labels = labels; });
}
