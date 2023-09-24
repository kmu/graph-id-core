#include "structure_graph.h"

std::string blake2b(const std::string &s) {
    py::object hashlib = py::module_::import("hashlib");
    return hashlib("blake2b")(py::bytes(s)).attr("hexdigest")().cast<std::string>();
}

/// グラフの直径を計算する。グラフが連結ではない場合、各連結成分の最も大きい直径を返す。
int graph_diameter(const std::vector<std::vector<NearNeighborInfo>> &graph) {
    const int n = int(graph.size());
    if (graph.size() <= 1) return 0;

    int ret = 0;
    std::vector<int> queue;
    queue.reserve(n);
    Eigen::VectorXi d(n);
    for (int start = 0; start < n; start++) {
        d.setConstant(-1);
        queue.clear();
        queue.push_back(start);
        d[start] = 0;
        int qi = 0;
        while (qi < int(queue.size())) {
            const int v = queue[qi++];
            for (const auto &nni: graph[v]) {
                const int u = nni.site_index;
                if (d[u] == -1) {
                    d[u] = d[v] + 1;
                    queue.push_back(u);
                }
            }
        }
        ret = std::max(ret, d.maxCoeff());
    }

    return ret;
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
    sg.graph_diameter = ::graph_diameter(sg.graph);
    return sg;
}

StructureGraph StructureGraph::with_empty_graph(const std::shared_ptr<const Structure> &structure) {
    const auto n = structure->count;
    return StructureGraph{
            structure,
            std::vector<std::vector<NearNeighborInfo>>(n),
            {},
            std::vector<std::string>(n),
            0,
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
                          [](StructureGraph &sg, const std::vector<std::string> &labels) { sg.labels = labels; })
            .def_property_readonly("graph_diameter", [](const StructureGraph &sg) { return sg.graph_diameter; })
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
            });
}
