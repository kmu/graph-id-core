#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace  py = pybind11;

// Base class to determine near neighbors that typically include nearest
// neighbors and others that are within some tolerable distance.
//
// C++ implementation of pymatgen.analysis.local_env.NearNeighbor
class NearNeighbor {
public:
    virtual ~NearNeighbor() = default;
    virtual bool structure_allowed() {return false;};
    virtual bool molecule_allowed() {return false;};
    virtual std::vector<py::dict> get_all_nn_info(py::object structure) = 0;
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
    std::vector<py::dict> get_all_nn_info(py::object structure) override;
};

void init_near_neighbor(py::module &m);
