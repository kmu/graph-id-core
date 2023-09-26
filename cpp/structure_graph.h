#pragma once

#include <map>
#include <utility>
#include <vector>
#include <pybind11/pybind11.h>
#include "near_neighbor.h"
#include "core.h"

std::string blake2b(const std::string &s);

class StructureGraph {
public:
    std::shared_ptr<const Structure> structure;

    // 連結リスト形式のグラフ
    std::vector<std::vector<NearNeighborInfo>> graph;

    // from, to, jiamge のタプルから graph[from] の NearNeighborInfo への index へのマップ
    // graph_map[from, to, image] = index のとき、graph[from][index] が to への NearNeighborInfo である。
    std::map<std::tuple<int, int, std::array<int, 3>>, int> graph_map;

    std::vector<std::string> labels;

    std::vector<std::string> cc_cs;
    std::vector<std::vector<int>> cc_nodes;
    std::vector<int> cc_diameter;

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

    py::object to_py() const;

private:
    void add_edge(
            int from,
            std::array<int, 3> from_image,
            int to,
            std::array<int, 3> to_image,
            double weight
    );

    void set_cc_diameter();
};

void init_structure_graph(pybind11::module &m);
