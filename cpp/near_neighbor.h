#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "wrapper.h"

namespace py = pybind11;

struct NearNeighborInfo {
    int site_index;
    double weight;
    std::array<int, 3> image;
};

// Base class to determine near neighbors that typically include nearest
// neighbors and others that are within some tolerable distance.
//
// C++ implementation of pymatgen.analysis.local_env.NearNeighbor
class NearNeighbor {
public:
    virtual ~NearNeighbor() = default;

    virtual bool structures_allowed() { return false; };

    virtual bool molecules_allowed() { return false; };

    // Pymatgen の NearNeighbor の get_all_nn_info と同じ処理を行う。
    // サブクラスで実装される get_all_nn_info_cpp を呼び出し、その結果を list[list[dict]] に変換して返す。
    py::list get_all_nn_info(py::object &structure);

    virtual std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) = 0;
};

class MinimumDistanceNN : public NearNeighbor {
private:
    double tol;
    double cutoff;
    bool get_all_sites;
public:
    explicit MinimumDistanceNN(double tol = 0.1, double cutoff = 10.0, bool get_all_sites = false) {
        this->tol = tol;
        this->cutoff = cutoff;
        this->get_all_sites = get_all_sites;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) override;
};

void init_near_neighbor(py::module &m);
