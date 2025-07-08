#pragma once

#include <pybind11/pybind11.h>
#include "near_neighbor.h"
#include "structure_graph.h"

class GraphIDGenerator {
public:
    std::shared_ptr<const NearNeighbor> nn;
    bool wyckoff = false;
    int depth_factor = 2;
    int additional_depth = 1;
    double symmetry_tol = 0.1;
    bool topology_only = false;
    bool loop = false;
    int rank_k = 3;
    double cutoff = 6.0;
    int digest_size = 8;

    GraphIDGenerator(
            const std::shared_ptr<const NearNeighbor> &nn,
            bool wyckoff,
            int depth_factor,
            int additional_depth,
            double symmetry_tol,
            bool topology_only,
            bool loop,
            int rank_k,
            double cutoff,
            int digest_size
    ) : wyckoff(wyckoff), depth_factor(depth_factor), additional_depth(additional_depth),
        symmetry_tol(symmetry_tol), topology_only(topology_only), loop(loop),
        rank_k(rank_k), cutoff(cutoff), digest_size(digest_size) {
        if (nn) {
            this->nn = nn;
        } else {
            this->nn = std::make_shared<MinimumDistanceNN>();
        }
    }

    std::string get_id(const Structure &structure) const;
    std::string get_id_with_structure_graph(py::object py_structure_graph) const;
    std::string get_distance_clustering_id(const Structure &structure) const;

    std::string get_id_catch_error(const Structure &structure) const noexcept;

    std::vector<std::string> get_many_ids(const std::vector<Structure> &structures) const;

    std::vector<std::string> get_component_ids(const Structure &structure) const;

    std::string elaborate_comp_dim(const StructureGraph &sg, const std::string &gid) const;

    bool are_same(const Structure &structure1, const Structure &structure2) const;

private:
    StructureGraph prepare_structure_graph(std::shared_ptr<const Structure> &structure) const;
    StructureGraph prepare_structure_graph_from_existing(std::shared_ptr<const Structure> &structure, const StructureGraph &sg) const;
    StructureGraph prepare_minimum_distance_structure_graph(std::shared_ptr<const Structure> &structure) const;
    StructureGraph prepare_disctance_clustering_structure_graph(int n, std::shared_ptr<const Structure> &structure, std::shared_ptr<StructureGraph> &_sg, int rank_k, double cutoff) const;
};

void init_graph_id(pybind11::module &m);
