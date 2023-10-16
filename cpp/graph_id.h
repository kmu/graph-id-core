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

    GraphIDGenerator(
            const std::shared_ptr<const NearNeighbor> &nn,
            bool wyckoff,
            int depth_factor,
            int additional_depth,
            double symmetry_tol,
            bool topology_only
    ) : wyckoff(wyckoff), depth_factor(depth_factor), additional_depth(additional_depth),
        symmetry_tol(symmetry_tol), topology_only(topology_only) {
        if (nn) {
            this->nn = nn;
        } else {
            this->nn = std::make_shared<MinimumDistanceNN>();
        }
    }

    std::string get_id(const Structure &structure) const;

    std::string get_id_catch_error(const Structure &structure) const noexcept;

    std::vector<std::string> get_many_ids(const std::vector<Structure> &structures) const;

    std::vector<std::string> get_component_ids(const Structure &structure) const;

    std::string elaborate_comp_dim(const StructureGraph &sg, const std::string &gid) const;

    bool are_same(const Structure &structure1, const Structure &structure2) const;

private:
    StructureGraph prepare_structure_graph(std::shared_ptr<const Structure> &structure) const;
};

void init_graph_id(pybind11::module &m);
