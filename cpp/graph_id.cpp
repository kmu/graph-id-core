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
    std::string gid = blake2b(join_string(":", cc_labels), digest_size);

    // return elaborate_comp_dim(sg, gid);
    return gid;
}

std::string GraphIDGenerator::get_id_with_structure_graph(py::object py_structure_graph) const {
    // PythonのStructureGraphから構造体を取得
    py::object py_structure = py_structure_graph.attr("structure");
    auto s = py_structure.cast<PymatgenStructure>();
    auto s_ptr = std::make_shared<const Structure>(s);
    
    const auto sg_from_py = StructureGraph::from_py(py_structure_graph);
    
    // 既存のStructureGraphに対してラベル設定とcompositional sequenceの処理を適用
    auto sg = prepare_structure_graph_from_existing(s_ptr,sg_from_py);
    
    std::vector<std::string> cc_labels(sg.cc_nodes.size());
    for (size_t i = 0; i < sg.cc_nodes.size(); ++i) {
        std::vector<std::string> labels = sg.cc_cs[i];
        std::sort(labels.begin(), labels.end());
        cc_labels[i] = blake2b(join_string("-", labels));
    }
    std::sort(cc_labels.begin(), cc_labels.end());
    std::string gid = blake2b(join_string(":", cc_labels), digest_size);

    return gid;
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

std::string GraphIDGenerator::get_distance_clustering_id(const Structure &structure) const {
    auto s_ptr = std::shared_ptr<const Structure>(&structure, [](const Structure *) {});
    std::vector<std::string> gids(this->rank_k);
    const auto _sg = prepare_minimum_distance_structure_graph(s_ptr);
    for (int idx = 0; idx < this->rank_k; idx++){
        // MinimumDistanceNNでStructureGraphを作成する
        std::vector<std::string> j_strs(structure.count);
        for (int j = 0; j < structure.count; j++){
            auto sg = _sg;
            // py::object py_structure_graph = _sg.to_py();
            // 原子jを含むedgeを削除する処理
            for (const auto& [key, value] : _sg.graph_map){
                if (std::get<0>(key) == j || std::get<1>(key) == j){
                    sg.break_edge(std::get<0>(key), std::get<1>(key), std::get<2>(key), true);
                    // py_structure_graph.attr("break_edge")(
                    //     py::arg("from_index") = std::get<0>(key),
                    //     py::arg("to_index") = std::get<1>(key),
                    //     py::arg("to_jimage") = std::get<2>(key),
                    //     py::arg("allow_reverse") = true
                    // );
                }
            }
            auto _s_ptr =  std::shared_ptr<const Structure>(&structure, [](const Structure *) {});
            // auto sg = StructureGraph::from_py(py_structure_graph, _s_ptr);
            auto _sg_ptr = std::shared_ptr<StructureGraph>(&sg, [](StructureGraph *) {});
            const auto sg_for_cc = prepare_disctance_clustering_structure_graph(j, _s_ptr, _sg_ptr, idx, this->cutoff);
            std::vector<std::string> cc_labels(sg_for_cc.cc_cs.size());
            for (size_t i = 0; i < sg_for_cc.cc_cs.size(); ++i) {
                std::vector<std::string> labels = sg_for_cc.cc_cs[i];
                std::sort(labels.begin(), labels.end());
                cc_labels[i] = blake2b(join_string("-", labels));
                // cc_labels[i] = blake2b(join_string("-", labels), 16);
            }
            std::sort(cc_labels.begin(), cc_labels.end());
            // std::string j_str = blake2b(join_string(":", cc_labels), 16);
            std::string j_str = join_string(":", cc_labels);
            j_strs.at(j) = j_str;
        }
        std::sort(j_strs.begin(), j_strs.end());
        std::string gid = blake2b(join_string(":", j_strs), digest_size);
        gids.at(idx) = gid;
    }
    // std::string all_gid = blake2b(join_string(":", gids), 16);
    std::string all_gid = join_string("", gids);
    // return elaborate_comp_dim(_sg, blake2b(all_gid, 16));
    return blake2b(all_gid, digest_size);
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
    } else if (loop) {
        sg.set_loops(depth_factor, additional_depth);
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

StructureGraph GraphIDGenerator::prepare_structure_graph_from_existing(std::shared_ptr<const Structure> &structure, const StructureGraph &sg_from_py) const {
    auto sg = sg_from_py;  // コピーを作成
    bool use_previous_cs = false;

    auto labels = structure->species_strings;
    auto prev_num_uniq = std::unique(labels.begin(), labels.end()) - labels.begin();

    if (wyckoff) {
        sg.set_wyckoffs_label(symmetry_tol);
    } else if (topology_only) {
        sg.labels = std::vector<std::string>(structure->count, "X");
    } else if (loop) {
        sg.set_loops(depth_factor, additional_depth);
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

// TODO 一つのクラスにまとめる。もしくは継承を使う。
StructureGraph GraphIDGenerator::prepare_minimum_distance_structure_graph(std::shared_ptr<const Structure> &structure) const {
    auto sg = StructureGraph::with_local_env_strategy(structure, MinimumDistanceNN());
    bool use_previous_cs = false;

    auto labels = structure->species_strings;
    auto prev_num_uniq = std::unique(labels.begin(), labels.end()) - labels.begin();

    if (wyckoff) {
        sg.set_wyckoffs_label(symmetry_tol);
    } else if (topology_only) {
        sg.labels = std::vector<std::string>(structure->count, "X");
    } else if (loop) {
        sg.set_loops(depth_factor, additional_depth);
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

StructureGraph GraphIDGenerator::prepare_disctance_clustering_structure_graph(int n, std::shared_ptr<const Structure> &structure, std::shared_ptr<StructureGraph> &_sg, int rank_k, double cutoff) const {
    // DistanceClusteringNNのインスタンスを渡す
    // auto sg = StructureGraph::with_individual_state_comp_strategy(structure, sg, DistanceClusteringNN, n, rank_k, cutoff);
    auto sg = StructureGraph::with_individual_state_comp_strategy(structure, *_sg, n, rank_k, cutoff);
    bool use_previous_cs = false;

    auto labels = structure->species_strings;
    auto prev_num_uniq = std::unique(labels.begin(), labels.end()) - labels.begin();

    if (wyckoff) {
        sg.set_wyckoffs_label(symmetry_tol);
    } else if (topology_only) {
        sg.labels = std::vector<std::string>(structure->count, "X");
    } else if (loop) {
        sg.set_loops(depth_factor, additional_depth);
    } else {
        sg.set_elemental_labels();
    }

    while (true) {
        sg.set_individual_compositional_sequence_node_attr(
                n,
                false, // hash_cs
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
            .def(py::init<std::shared_ptr<NearNeighbor>, bool, int, int, double, bool, bool, int, double, int>(),
                 py::arg("nn") = nullptr,
                 py::arg("wyckoff") = false,
                 py::arg("depth_factor") = 2,
                 py::arg("additional_depth") = 1,
                 py::arg("symmetry_tol") = 0.1,
                 py::arg("topology_only") = false,
                 py::arg("loop") = false,
                 py::arg("rank_k") = 3,
                 py::arg("cutoff") = 6.0,
                 py::arg("digest_size") = 8)
            .def("get_id", [](const GraphIDGenerator &gig, py::object &structure) {
                auto s = structure.cast<PymatgenStructure>();
                return gig.get_id(Structure(s));
            })
            .def("get_id_with_structure_graph", [](const GraphIDGenerator &gig, py::object &py_structure_graph) {
                return gig.get_id_with_structure_graph(py_structure_graph);
            })
            .def("get_distance_clustering_id", [](const GraphIDGenerator &gig, py::object &structure) {
                auto s = structure.cast<PymatgenStructure>();
                return gig.get_distance_clustering_id(Structure(s));
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
