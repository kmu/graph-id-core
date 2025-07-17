#pragma once

#include <map>
#include <utility>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gtl/phmap.hpp>
#include "near_neighbor.h"
#include "core.h"

std::string blake2b(const std::string &s);

std::string blake2b(const std::string &s, int digest_size);

class StructureGraph {
public:
    std::shared_ptr<const Structure> structure;

    // A graph in a linked list format
    std::vector<std::vector<NearNeighborInfo>> graph;

    // A map from a tuple of from, to, and image to the index of NearNeighborInfo in graph[from]
    // When graph_map[from, to, image] = index, graph[from][index] is the NearNeighborInfo to to.
    std::map<std::tuple<int, int, std::array<int, 3>>, int> graph_map;

    std::vector<std::string> labels;

    std::vector<std::vector<std::string>> cc_cs;
    std::vector<std::vector<int>> cc_nodes;
    std::vector<int> cc_diameter;

    ~StructureGraph() = default;

    static StructureGraph with_local_env_strategy(
            const std::shared_ptr<const Structure> &structure,
            const NearNeighbor &strategy);

    static StructureGraph with_individual_state_comp_strategy(
            const std::shared_ptr<const Structure> &structure,
            StructureGraph &sg,
        //     const NearNeighbor &strategy,
            int n,
            int rank_k,
            double coutoff);

    static StructureGraph from_empty_graph(const std::shared_ptr<const Structure> &structure);

    static StructureGraph from_py(py::object py_sg);

    void set_elemental_labels();

    void set_wyckoffs_label(double symmetry_tol = 0.1);

    void set_loops(int depth_factor, int additional_depth);

    void set_compositional_sequence_node_attr(
            bool hash_cs,
            bool wyckoff,
            int additional_depth,
            int depth_factor,
            bool use_previous_cs
    );

    void set_individual_compositional_sequence_node_attr(
            int n,
            bool hash_cs,
            bool wyckoff,
            int additional_depth,
            int depth_factor,
            bool use_previous_cs
    );

    int get_dimensionality_larsen() const;

    py::object to_py() const;

public:
    void add_edge(
            int from,
            std::array<int, 3> from_image,
            int to,
            std::array<int, 3> to_image,
            double weight
    );

    void break_edge(
            int from,
            int to,
            std::array<int, 3> image,
            bool allow_reverse
    );

    void set_cc_diameter();

    static bool rank_increase(const gtl::flat_hash_set<std::array<int, 3>> &seen, const std::array<int, 3> &candidate);
    
    static int calculate_rank(const gtl::flat_hash_set<std::array<int, 3>> &vertices);
};

uint64_t connected_site_to_uint64(int site, std::array<int, 3> arr);

class CompositionalSequence {
public:
    CompositionalSequence() = default;

    bool hash_cs{false};
    std::string cs_for_hashing;
    std::vector<std::string> compositional_seq;
    int focused_site_i{0};
    std::vector<std::tuple<int, std::array<int, 3>>> new_sites;
    gtl::flat_hash_set<uint64_t> seen_sites;
    bool use_previous_sites{false};
    std::vector<std::string> const *labels;
    std::map<std::string, int> composition_counter;

    std::string string() const;

    std::vector<std::tuple<int, std::array<int, 3>>> get_current_starting_sites();

    void count_composition_for_neighbors(int site_i, std::array<int, 3> image);

    void finalize_this_depth();

    std::vector<std::string> get_sorted_composition_list_form() const;
};

void init_structure_graph(pybind11::module &m);
