#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gtl/phmap.hpp>
#include <utility>
#include "core.h"

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

    virtual std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const = 0;
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

    bool molecules_allowed() override { return true; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

class CutOffDictNN : public NearNeighbor {
private:
    double max_cut_off;
    gtl::flat_hash_map<std::pair<std::string, std::string>, double> cut_off_dict;
public:
    explicit CutOffDictNN(const std::optional<py::dict> &cut_off_dict) {
        max_cut_off = 0;
        if (!cut_off_dict || cut_off_dict->is_none()) return;
        for (auto item: *cut_off_dict) {
            auto key = item.first.cast<std::pair<std::string, std::string>>();
            auto value = item.second.cast<double>();
            this->cut_off_dict[key] = value;
            std::swap(key.first, key.second);
            this->cut_off_dict[key] = value;
            if (value > max_cut_off) {
                max_cut_off = value;
            }
        }
    };

    explicit CutOffDictNN(gtl::flat_hash_map<std::pair<std::string, std::string>, double> cut_off_dict) : cut_off_dict(
            std::move(cut_off_dict)) {
        max_cut_off = 0;
        for (const auto &item: this->cut_off_dict) {
            if (item.second > max_cut_off) {
                max_cut_off = item.second;
            }
        }
    };

    explicit CutOffDictNN() : max_cut_off(0) {};

    static CutOffDictNN from_preset(const std::string &preset) {
        py::module local_env = py::module::import("pymatgen.analysis.local_env");
        py::object nn = local_env.attr("CutOffDictNN").attr("from_preset")(preset);
        return CutOffDictNN(nn.attr("cut_off_dict"));
    }

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return true; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

struct FindNearNeighborsResult {
    int all_coords_idx;
    std::array<int, 3> image;
    double distances2;
};

std::vector<std::vector<FindNearNeighborsResult>> find_near_neighbors(
        const Eigen::Matrix3Xd &all_coords,
        const Eigen::Matrix3Xd &all_frac_coords,
        const Eigen::Matrix3Xd &center_coords,
        const Eigen::Matrix3Xd &center_frac_coords,
        double r,
        const Lattice &lattice,
        double min_r = 1,
        double tol = 1e-8
);

std::vector<std::vector<FindNearNeighborsResult>> find_near_neighbors(
        const Structure &structure,
        double r,
        double min_r = 1,
        double tol = 1e-8
);

Eigen::Matrix3Xd get_reciprocal_lattice(const Eigen::Matrix3d &lattice);

std::pair<Eigen::Vector3i, Eigen::Vector3i> get_bounds(
        const Eigen::Matrix3Xd &frac_coords,
        const Eigen::Vector3d &maxr,
        const std::array<bool, 3> &pbc
);

Eigen::VectorXi three_to_one(const Eigen::Matrix3Xi &label3d, const Eigen::Vector3i &n_cube);

int three_to_one1(const Eigen::Vector3i &label, const Eigen::Vector3i &n_cube);

Eigen::Matrix3Xi one_to_three(const Eigen::VectorXi &label1d, const Eigen::Vector3i &n_cube);

Eigen::Vector3i one_to_three1(int label, const Eigen::Vector3i &n_cube);

std::vector<std::vector<int>> get_cube_neighbors(const Eigen::Vector3i &n_cube);

void init_near_neighbor(py::module &m);
