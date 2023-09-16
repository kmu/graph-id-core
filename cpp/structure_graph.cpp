#include "structure_graph.h"

void StructureGraph::set_elemental_labels() {

}

void StructureGraph::set_wyckoffs_label(double symmetry_tol) {

}

void StructureGraph::set_compositional_sequence_node_attr(
        bool hash_cs,
        bool wyckoff,
        int additional_depth,
        int depth_factor
) {

}

void init_structure_graph(pybind11::module &m) {
    py::class_<StructureGraph>(m, "StructureGraph")
            .def_static("with_local_env_strategy", &StructureGraph::with_local_env_strategy)
            .def_static("with_empty_graph", &StructureGraph::with_empty_graph)
            .def("set_elemental_labels", &StructureGraph::set_elemental_labels)
            .def("set_wyckoffs_label", &StructureGraph::set_wyckoffs_label)
            .def("set_compositional_sequence_node_attr", &StructureGraph::set_compositional_sequence_node_attr)
            .def_property("labels", [](const StructureGraph &sg) { return sg.labels; },
                          [](StructureGraph &sg, const std::vector<std::string> &labels) { sg.labels = labels; });
}
