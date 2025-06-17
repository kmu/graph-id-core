#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <gtl/phmap.hpp>
#include <utility>
#include "core.h"

namespace py = pybind11;


struct FindNearNeighborsResult {
    int all_coords_idx;
    std::array<int, 3> image;
    double distance;

    Eigen::Vector3d xyz(const Structure &s) const {
        return s.site_xyz.col(all_coords_idx) + s.lattice.matrix * Eigen::Vector3d(image[0], image[1], image[2]);
    }

    py::object to_periodic_neighbor(const Structure &s) const {
        auto PeriodicNeighbor = py::module::import("pymatgen.core.structure").attr("PeriodicNeighbor");
        return PeriodicNeighbor(
                py::arg("species") = s.py_structure.sites()[all_coords_idx].obj.attr("species"),
                py::arg("coords") = s.lattice.inv_matrix * s.site_xyz.col(all_coords_idx) +
                                    Eigen::Vector3d(image[0], image[1], image[2]),
                py::arg("lattice") = s.py_structure.lattice().obj,
                py::arg("properties") = s.py_structure.sites()[all_coords_idx].obj.attr("properties"),
                py::arg("nn_distance") = distance,
                py::arg("index") = all_coords_idx,
                py::arg("image") = py::make_tuple(image[0], image[1], image[2]),
                py::arg("label") = s.species_strings[all_coords_idx]
        );
    }
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

/// Helper method to calculate the solid angle of a set of coords from the center
double solid_angle(Eigen::Vector3d center, Eigen::Matrix3Xd coords);

/// Calculate the volume of a tetrahedron
double vol_tetra(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3);

double get_default_radius(py::object site);

double get_radius(py::object site);

double get_mean_fictive_ionic_radius(const Eigen::VectorXd &fictive_ionic_radii, double minimum_fir);

struct NearNeighborInfo {
    int site_index;
    double weight;
    std::array<int, 3> image;
    std::optional<py::dict> extra;

    Eigen::Vector3d xyz(const Structure &s) const {
        return s.site_xyz.col(site_index) + s.lattice.matrix * Eigen::Vector3d(image[0], image[1], image[2]);
    }
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

struct VoronoiPolyhedra {
    FindNearNeighborsResult site{};
    double solid_angle{NAN};
    double angle_normalized{NAN};
    double area{NAN};
    double face_dist{NAN};
    double volume{NAN};
    int n_verts{0};
    std::vector<int> verts;
    std::vector<int> adj_neighbors;
    Eigen::Vector3d normal;

    double operator[](const std::string &s) const {
        if (s == "solid_angle") return solid_angle;
        if (s == "angle_normalized") return angle_normalized;
        if (s == "area") return area;
        if (s == "face_dist") return face_dist;
        if (s == "volume") return volume;
        if (s == "n_verts") return n_verts;
        throw std::invalid_argument("invalid key");
    }

    py::dict to_dict(const Structure &s) const {
        py::dict ret;
        ret["site"] = site.to_periodic_neighbor(s);
        ret["solid_angle"] = solid_angle;
        if (!std::isnan(angle_normalized)) ret["angle_normalized"] = angle_normalized;
        ret["area"] = area;
        ret["face_dist"] = face_dist;
        ret["volume"] = volume;
        ret["n_verts"] = n_verts;
        if (!verts.empty()) ret["verts"] = verts;
        ret["adj_neighbors"] = adj_neighbors;
        ret["normal"] = normal;
        return ret;
    }
};

// A wrapper class for scipy.spatial.Voronoi
class Voronoi {
public:
    py::object obj;

    explicit Voronoi(
            const Eigen::Matrix3Xd &coords,
            bool furthest_site = false,
            bool incremental = false,
            const std::string &qhull_options = ""
    ) {
        py::object scipy_spatial = py::module::import("scipy.spatial");
        py::object Voronoi = scipy_spatial.attr("Voronoi");
        Eigen::MatrixX3d c = coords.transpose();
        py::object opt = qhull_options.empty() ? py::none() : py::object(py::str(qhull_options));
        obj = Voronoi(c, furthest_site, incremental, opt);
    }

    bool furthest_site() const {
        return obj.attr("furthest_site").cast<bool>();
    }

    py::array_t<double> max_bound() const {
        return obj.attr("max_bound").cast<py::array_t<double>>();
    }

    py::array_t<double> min_bound() const {
        return obj.attr("min_bound").cast<py::array_t<double>>();
    }

    int ndim() const {
        return obj.attr("ndim").cast<int>();
    }

    int npoints() const {
        return obj.attr("npoints").cast<int>();
    }

    // ndarray(npoints, 1)
    py::array_t<int> point_region() const {
        return obj.attr("point_region").cast<py::array_t<int>>();
    }

    // ndarray(npoints, ndim)
    py::array_t<double> points() const {
        return obj.attr("points").cast<py::array_t<double>>();
    }

    std::vector<std::vector<int>> regions() const {
        return obj.attr("regions").cast<std::vector<std::vector<int>>>();
    }

    // dict[(int, int), list[int]]
    py::dict ridge_dict() const {
        return obj.attr("ridge_dict").cast<py::dict>();
    }

    // ndarray(npoints, 2)
    py::array_t<int> ridge_points() const {
        return obj.attr("ridge_points").cast<py::array_t<int>>();
    }

    std::vector<std::vector<int>> ridge_vertices() const {
        return obj.attr("ridge_vertices").cast<std::vector<std::vector<int>>>();
    }

    // ndarray(npoints, ndim)
    py::array_t<double> vertices() const {
        return obj.attr("vertices").cast<py::array_t<double>>();
    }
};

class VoronoiNN : public NearNeighbor {
private:
    double tol;
    std::optional<std::vector<std::string>> targets;
    double cutoff;
    bool allow_pathological;
    std::string weight;
    bool extra_nn_info;
    bool compute_adj_neighbors;
public:
    explicit VoronoiNN(
            double tol = 0,
            const std::optional<std::vector<std::string>> &targets = std::nullopt,
            double cutoff = 13.0,
            bool allow_pathological = false,
            const std::string &weight = "solid_angle",
            bool extra_nn_info = true,
            bool compute_adj_neighbors = true
    ) {
        this->tol = tol;
        this->targets = targets;
        this->cutoff = cutoff;
        this->allow_pathological = allow_pathological;
        this->weight = weight;
        this->extra_nn_info = extra_nn_info;
        this->compute_adj_neighbors = compute_adj_neighbors;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;

    std::unordered_map<int, VoronoiPolyhedra> get_voronoi_polyhedra(const Structure &structure, int site_index) const;

    std::vector<std::unordered_map<int, VoronoiPolyhedra>> get_all_voronoi_polyhedra(const Structure &structure) const;

    std::unordered_map<int, VoronoiPolyhedra> extract_cell_info(
            int neighbor_index,
            const Structure &structure,
            const std::vector<FindNearNeighborsResult> &neighbors,
            const std::vector<std::string> &targets,
            const Voronoi &voro,
            bool compute_adj_neighbors = false) const;

    std::vector<NearNeighborInfo>
    extract_nn_info(const Structure &s, const std::unordered_map<int, VoronoiPolyhedra> &v) const;

    std::vector<NearNeighborInfo>
    extract_nn_info(const Structure &s, const std::unordered_map<int, VoronoiPolyhedra> &v,
                    const std::vector<std::string> &targets) const;

private:
    static double get_max_cutoff(const Structure &s);
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

class MinimumOKeeffeNN : public NearNeighbor {
private:
    double tol;
    double cutoff;
public:
    explicit MinimumOKeeffeNN(double tol = 0.1, double cutoff = 10.0) {
        this->tol = tol;
        this->cutoff = cutoff;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return true; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

class DistanceClusteringNN : public NearNeighbor{
private:
    double tol;
    int n;
    int rank_k;
    double cutoff;
public:
    explicit DistanceClusteringNN(double tol = 0.1, int n = 0, int rank_k = 3, double cutoff = 6.0) {
        this->tol = tol;
        this->n = n;
        this->rank_k = rank_k;
        this->cutoff = cutoff;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return true; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
    std::vector<double> get_cutoff_cluster(const Structure &structure, int n, double cutoff, const auto &nn) const;
};

class CrystalNN : public NearNeighbor {
private:
    bool weighted_cn;
    bool cation_anion;
    std::pair<double, double> distance_cutoffs;
    double x_diff_weight;
    bool porous_adjustment;
    double search_cutoff;
    int fingerprint_length;
public:
    explicit CrystalNN(
            bool weighted_cn = false,
            bool cation_anion = false,
            std::pair<double, double> distance_cutoffs = {0.5, 1},
            double x_diff_weight = 3.0,
            bool porous_adjustment = true,
            double search_cutoff = 7,
            int fingerprint_length = 0
    ) {
        this->weighted_cn = weighted_cn;
        this->cation_anion = cation_anion;
        this->distance_cutoffs = distance_cutoffs;
        this->x_diff_weight = x_diff_weight;
        this->porous_adjustment = porous_adjustment;
        this->search_cutoff = search_cutoff;
        this->fingerprint_length = fingerprint_length;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;

    struct NNData {
        std::vector<NearNeighborInfo> all_nninfo;
        std::unordered_map<int, double> cn_weights;
        std::unordered_map<int, std::vector<NearNeighborInfo>> cn_nninfo;
    };

    std::vector<NNData> get_all_nn_data(const Structure &structure, int length = 0) const;

private:
    static double semicircle_integral(const std::vector<double> &dist_bins, int idx);

    static void transform_to_length(NNData &nn_data, int length);
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

class BrunnerNN_reciprocal : public NearNeighbor {
private:
    double tol;
    double cutoff;
public:
    explicit BrunnerNN_reciprocal(double tol = 1e-4, double cutoff = 8.0) {
        this->tol = tol;
        this->cutoff = cutoff;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

class BrunnerNN_relative : public NearNeighbor {
private:
    double tol;
    double cutoff;
public:
    explicit BrunnerNN_relative(double tol = 1e-4, double cutoff = 8.0) {
        this->tol = tol;
        this->cutoff = cutoff;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

class BrunnerNN_real : public NearNeighbor {
private:
    double tol;
    double cutoff;
public:
    explicit BrunnerNN_real(double tol = 1e-4, double cutoff = 8.0) {
        this->tol = tol;
        this->cutoff = cutoff;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return false; };

    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

class EconNN : public NearNeighbor {
private:
    double tol;
    double cutoff;
    bool cation_anion;
    bool use_fictive_radius;
public:
    explicit EconNN(double tol = 0.2, double cutoff = 10.0, bool cation_anion = false,
                    bool use_fictive_radius = false) {
        this->tol = tol;
        this->cutoff = cutoff;
        this->cation_anion = cation_anion;
        this->use_fictive_radius = use_fictive_radius;
    };

    bool structures_allowed() override { return true; };

    bool molecules_allowed() override { return true; };


    std::vector<std::vector<NearNeighborInfo>> get_all_nn_info_cpp(const Structure &structure) const override;
};

void init_near_neighbor(py::module &m);
