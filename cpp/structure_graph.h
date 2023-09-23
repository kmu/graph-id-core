#pragma once

#include <utility>
#include <vector>
#include <pybind11/pybind11.h>
#include "near_neighbor.h"
#include "core.h"

std::string blake2b(const std::string &s);

// グラフの直径を計算する
int graph_diameter(const std::vector<std::vector<NearNeighborInfo>> &graph);

class StructureGraph {
public:
    std::shared_ptr<const Structure> structure;

    // 連結リスト形式のグラフ
    std::vector<std::vector<NearNeighborInfo>> graph;

    std::vector<std::string> labels;

    int graph_diameter;

    ~StructureGraph() = default;

    static StructureGraph with_local_env_strategy(
            const std::shared_ptr<const Structure> &structure,
            const NearNeighbor &strategy);

    static StructureGraph with_empty_graph(const std::shared_ptr<const Structure> &structure);

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
