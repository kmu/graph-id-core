#pragma once

#include <utility>
#include <vector>
#include <pybind11/pybind11.h>
#include "near_neighbor.h"
#include "core.h"

class StructureGraph {
public:
    std::shared_ptr<const Structure> structure;

    // 連結リスト形式のグラフ
    std::vector<std::vector<NearNeighborInfo>> graph;

    std::vector<std::string> labels;

    ~StructureGraph() = default;

    static StructureGraph
    with_local_env_strategy(const std::shared_ptr<const Structure> &structure, const NearNeighbor &strategy) {
        auto sg = with_empty_graph(structure);
        sg.graph = strategy.get_all_nn_info_cpp(*structure);
        return sg;
    }

    static StructureGraph with_empty_graph(const std::shared_ptr<const Structure> &structure) {
        const auto n = structure->count;
        return StructureGraph{
                structure,
                std::vector<std::vector<NearNeighborInfo>>(n),
                std::vector<std::string>(n),
        };
    }

    void set_elemental_labels();

    void set_wyckoffs_label(double symmetry_tol = 0.1);

    void set_compositional_sequence_node_attr(
            bool hash_cs,
            bool wyckoff,
            int additional_depth,
            int depth_factor
    );
};

void init_structure_graph(pybind11::module &m);
