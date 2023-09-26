#include "structure_graph.h"

std::string blake2b(const std::string &s) {
    py::object hashlib = py::module_::import("hashlib");
    return hashlib("blake2b")(py::bytes(s)).attr("hexdigest")().cast<std::string>();
}

StructureGraph StructureGraph::with_local_env_strategy(
        const std::shared_ptr<const Structure> &structure,
        const NearNeighbor &strategy
) {
    auto sg = StructureGraph::with_empty_graph(structure);
    const auto &nn = strategy.get_all_nn_info_cpp(*structure);
    assert(int(nn.size()) == structure->count);
    for (int from = 0; from < int(nn.size()); ++from) {
        for (size_t i = 0; i < nn[from].size(); ++i) {
            const auto &nni = nn[from][i];
            sg.add_edge(from, {0, 0, 0}, nni.site_index, nni.image, nni.weight);
            sg.add_edge(nni.site_index, nni.image, from, {0, 0, 0}, nni.weight);
        }
    }
    sg.set_cc_diameter();
    return sg;
}

StructureGraph StructureGraph::with_empty_graph(const std::shared_ptr<const Structure> &structure) {
    const auto n = structure->count;
    return StructureGraph{
            structure,
            std::vector<std::vector<NearNeighborInfo>>(n),
            {},
            std::vector<std::string>(n),
            {},
            {},
            {},
    };
}

void StructureGraph::add_edge(
        int from,
        std::array<int, 3> from_image,
        int to,
        std::array<int, 3> to_image,
        double weight
) {
    assert(0 <= from && from < int(this->graph.size()));
    assert(0 <= to && to < int(this->graph.size()));
    std::array<int, 3> image{};
    for (int i = 0; i < 3; ++i) image[i] = to_image[i] - from_image[i];

    // 自分自身への辺は無視する
    if (from == to && image == std::array<int, 3>{0, 0, 0}) {
        return;
    }

    // すでに追加されている辺は無視する
    if (this->graph_map.find(std::make_tuple(from, to, image)) != this->graph_map.end()) {
        return;
    }

    this->graph[from].emplace_back(NearNeighborInfo{to, weight, image});
    this->graph_map[std::make_tuple(from, to, image)] = int(this->graph[from].size() - 1);
}

// グラフの連結成分とその直径を計算する
void StructureGraph::set_cc_diameter() {
    const int n = int(graph.size());
    assert(n > 0);

    this->cc_nodes.clear();
    this->cc_diameter.clear();
    this->cc_cs.clear();

    std::vector<bool> visited(n, false);
    std::vector<int> queue;
    queue.reserve(n);
    size_t qi = 0;

    // 幅優先探索で連結成分を調べる
    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;
        visited[i] = true;
        queue.push_back(i);
        this->cc_nodes.emplace_back();
        this->cc_nodes.back().push_back(i);
        while (qi < queue.size()) {
            int u = queue[qi++]; // pop_front
            for (const auto &nni: graph[u]) {
                if (!visited[nni.site_index]) {
                    visited[nni.site_index] = true;
                    this->cc_nodes.back().push_back(nni.site_index);
                    queue.push_back(nni.site_index);
                }
            }
        }
    }

    for (auto &vec: this->cc_nodes) std::sort(vec.begin(), vec.end());

    Eigen::VectorXi d(n);
    for (const auto &nodes: this->cc_nodes) {
        // 連結成分ごとに幅優先探索を行う
        int d_max = 0;
        qi = 0;
        queue.resize(0);
        for (const int start: nodes) {
            queue.push_back(start);
            for (const int node: nodes) visited[node] = false;
            d[start] = 0;
            visited[start] = true;
            while (qi < int(queue.size())) {
                const int v = queue[qi++];
                for (const auto &nni: graph[v]) {
                    const int u = nni.site_index;
                    if (!visited[u]) {
                        visited[u] = true;
                        d[u] = d[v] + 1;
                        d_max = std::max(d_max, d[u]);
                        queue.push_back(u);
                    }
                }
            }
        }
        this->cc_diameter.push_back(d_max);
        this->cc_cs.emplace_back("");
    }
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

py::object StructureGraph::to_py() const {
    py::object PmgStructureGraph = py::module::import("pymatgen.analysis.graphs").attr("StructureGraph");
    py::object sg = PmgStructureGraph.attr("with_empty_graph")(this->structure->py_structure);
    for (size_t i = 0; i < this->graph.size(); i++) {
        for (const auto &nni: this->graph[i]) {
            sg.attr("add_edge")(
                    py::arg("from_index") = i,
                    py::arg("from_jimage") = py::make_tuple(0, 0, 0),
                    py::arg("to_index") = nni.site_index,
                    py::arg("to_jimage") = nni.image,
                    py::arg("weight") = nni.weight
            );
        }
    }
    return sg;
}

void init_structure_graph(pybind11::module &m) {
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
            .def("to_py", &StructureGraph::to_py)
            .def("get_connected_site_index", [](const StructureGraph &sg) {
                // テスト用
                py::list arr;
                for (size_t i = 0; i < sg.graph.size(); i++) {
                    for (const auto &nni: sg.graph[i]) {
                        arr.append(py::make_tuple(i, nni.site_index));
                    }
                }
                arr.attr("sort")();
                return arr;
            })
            .def_property("labels", [](const StructureGraph &sg) { return sg.labels; },
                          [](StructureGraph &sg, const std::vector<std::string> &labels) { sg.labels = labels; })
            .def_property_readonly("cc_nodes", [](const StructureGraph &sg) { return sg.cc_nodes; })
            .def_property_readonly("cc_diameter", [](const StructureGraph &sg) { return sg.cc_diameter; })
            .def_property_readonly("cc_cs", [](const StructureGraph &sg) { return sg.cc_cs; });
}
