#include "near_neighbor.h"
#include <vector>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "core.h"


py::list NearNeighbor::get_all_nn_info(py::object &structure) {
    auto s = structure.cast<PymatgenStructure>();
    auto pymatgen = py::module_::import("pymatgen.core");
    const bool is_structure = py::isinstance(structure, pymatgen.attr("Structure"));
    const bool is_molecule = py::isinstance(structure, pymatgen.attr("Molecule"));
    if (is_structure) {
        if (!this->structures_allowed()) {
            throw std::domain_error("This class does not support structures.");
        }
    } else if (is_molecule) {
        if (!this->molecules_allowed()) {
            throw std::domain_error("This class does not support molecules.");
        }
    } else {
        throw std::domain_error("argument must be pymatgen.core.Structure or pymatgen.core.Molecule.");
    }
    const auto result = this->get_all_nn_info_cpp(Structure(s));
    py::list arr;
    for (const auto &infos: result) {
        py::list inner;
        for (const auto &info: infos) {
            py::dict d;
            d["site_index"] = info.site_index;
            d["weight"] = info.weight;
            if (is_structure) {
                d["image"] = info.image;
            } else {
                d["image"] = py::none();
            }
            inner.append(d);
        }
        arr.append(inner);
    }
    return arr;
}


std::vector<std::vector<NearNeighborInfo>> MinimumDistanceNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto nn = find_near_neighbors(structure, this->cutoff);
    assert(int(nn.size()) == structure.count);
    if (this->get_all_sites) {
        std::vector<std::vector<NearNeighborInfo>> result(structure.count);
        for (int i = 0; i < structure.count; ++i) {
            if (nn[i].empty()) continue;
            Eigen::VectorXd d2(nn[i].size());
            for (int j = 0; j < int(nn[i].size()); ++j) {
                d2(j) = nn[i][j].distances2;
            }
            const Eigen::VectorXd d = d2.array().sqrt();
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (nn[i][j].distances2 < 1e-8) continue;
                result[i].emplace_back(NearNeighborInfo{
                        nn[i][j].all_coords_idx,
                        std::sqrt(nn[i][j].distances2),
                        nn[i][j].image
                });
            }
        }
        return result;
    } else {
        std::vector<std::vector<NearNeighborInfo>> result(structure.count);
        for (int i = 0; i < structure.count; ++i) {
            if (nn[i].empty()) continue;
            Eigen::VectorXd d2(nn[i].size());
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (nn[i][j].distances2 < 1e-8) {
                    d2(j) = 9999;
                } else {
                    d2(j) = nn[i][j].distances2;
                }
            }
            const Eigen::VectorXd d = d2.array().sqrt();
            const double min_distance = d.minCoeff();
            const double r = (1 + this->tol) * min_distance;
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (nn[i][j].distances2 > 1e-8 && d(j) < r) {
                    result[i].emplace_back(NearNeighborInfo{
                            nn[i][j].all_coords_idx,
                            min_distance / d(j),
                            nn[i][j].image
                    });
                }
            }
        }
        return result;
    }
}

std::vector<std::vector<NearNeighborInfo>> CutOffDictNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto nn = find_near_neighbors(structure, this->max_cut_off);
    assert(int(nn.size()) == structure.count);

    std::vector<std::vector<NearNeighborInfo>> result(structure.count);
    for (int i = 0; i < structure.count; ++i) {
        if (nn[i].empty()) continue;
        for (int j = 0; j < int(nn[i].size()); ++j) {
            if (nn[i][j].distances2 < 1e-8) continue;
            const auto key = std::make_pair(structure.species_strings[i],
                                            structure.species_strings[nn[i][j].all_coords_idx]);
            const auto it = this->cut_off_dict.find(key);
            if (it == this->cut_off_dict.end()) continue;
            double distance = std::sqrt(nn[i][j].distances2);
            if (distance < it->second) {
                result[i].emplace_back(NearNeighborInfo{nn[i][j].all_coords_idx, distance, nn[i][j].image});
            }
        }
    }

    return result;
}

/// pymatgen.optimization.neighbors.find_points_in_spheres と同じ処理を行う。
/// Python から利用する際はオーバーヘッドを考慮すると pymatgen で実装されたものを使ったほうが速いことが多い。
/// 原子の座標と中心点の座標を行列で与えると、中心点からの距離が r 以下の原子の情報を返す。
///
/// \param all_coords 原子のデカルト座標
/// \param all_frac_coords 原子の分数座標
/// \param center_coords 中心点のデカルト座標
/// \param center_frac_coords 中心点の分数座標
/// \param r 半径
/// \param lattice 格子
/// \param min_r
/// \param tol
/// \return 中心点の数と同じ要素数の配列を返す。i 番目のリストは i 番目の中心点から距離が r 以下の原子の情報を含む。
std::vector<std::vector<FindNearNeighborsResult>> find_near_neighbors(
        const Eigen::Matrix3Xd &all_coords,
        const Eigen::Matrix3Xd &all_frac_coords,
        const Eigen::Matrix3Xd &center_coords,
        const Eigen::Matrix3Xd &center_frac_coords,
        const double r,
        const Lattice &lattice,
        const double min_r,
        const double tol
) {
    if (all_coords.size() == 0 || center_coords.size() == 0) {
        return {};
    }
    assert(all_coords.size() == all_frac_coords.size());
    assert(center_coords.size() == center_frac_coords.size());
    if (r < min_r) {
        const auto res = find_near_neighbors(all_coords, all_frac_coords, center_coords, center_frac_coords,
                                             min_r + tol, lattice, min_r, tol);
        const double r2 = r * r;
        std::vector<std::vector<FindNearNeighborsResult>> result(res.size());
        for (size_t i = 0; i < res.size(); ++i) {
            for (const auto &x: res[i]) {
                if (x.distances2 <= r2) {
                    result[i].emplace_back(x);
                }
            }
        }
        return result;
    }


    const double r2 = r * r;
    const long n_center = center_coords.cols();
    const long n_total = all_coords.cols();
    const double ledge = std::max(0.1, r);
    std::vector<std::vector<FindNearNeighborsResult>> result(n_center);
    Eigen::Vector3d valid_max = center_coords.rowwise().maxCoeff();
    Eigen::Vector3d valid_min = center_coords.rowwise().minCoeff();
    valid_max.array() += (r + tol);
    valid_min.array() -= (r + tol);

    Eigen::Matrix3Xd reciprocal_lattice = get_reciprocal_lattice(lattice.matrix);
    Eigen::Vector3d max_r = (r) * reciprocal_lattice.colwise().norm() / (2 * pi);
    max_r = Eigen::ceil(max_r.array());

    Eigen::Vector3i min_bound, max_bound;
    std::tie(min_bound, max_bound) = get_bounds(center_frac_coords, max_r, lattice.pbc);

    // Process pbc
    Eigen::Matrix3Xd f_coords_in_cell = all_frac_coords;
    Eigen::Matrix3Xd offset_correction(3, n_total);
    for (int i = 0; i < 3; ++i) {
        if (lattice.pbc[i]) {
            offset_correction.row(i) = all_frac_coords.row(i).array().floor();
            f_coords_in_cell.row(i) -= offset_correction.row(i);
        } else {
            offset_correction.row(i).setZero();
            f_coords_in_cell.row(i) = all_frac_coords.row(i);
        }
    }
    Eigen::Matrix3Xd coords_in_cell = lattice.matrix * f_coords_in_cell;

    // Get translated images, coordinates and indices
    std::vector<std::tuple<int, std::array<int, 3>>> indices;
    std::vector<double> expanded_coords_vec;
    for (int i = min_bound.x(); i < max_bound.x(); ++i) {
        for (int j = min_bound.y(); j < max_bound.y(); ++j) {
            for (int k = min_bound.z(); k < max_bound.z(); ++k) {
                const Eigen::Vector3d tmp = lattice.matrix * Eigen::Vector3d(i, j, k);
                for (int l = 0; l < n_total; ++l) {
                    const Eigen::Vector3d v = tmp + coords_in_cell.col(l);
                    if ((v.array() < valid_max.array()).all() && (v.array() > valid_min.array()).all()) {
                        indices.emplace_back(l, std::array<int, 3>{i, j, k});
                        expanded_coords_vec.push_back(v.x());
                        expanded_coords_vec.push_back(v.y());
                        expanded_coords_vec.push_back(v.z());
                    }
                }
            }
        }
    }

    // if no valid neighbors were found return empty
    if (indices.empty()) {
        return result;
    }

    Eigen::Matrix3Xd expanded_coords = Eigen::Map<Eigen::Matrix3Xd>(
            expanded_coords_vec.data(), 3, int(expanded_coords_vec.size()) / 3);
    Eigen::Vector3i n_cube = Eigen::ceil((valid_max - valid_min).array() / ledge).cast<int>();
    int n_cube_all = n_cube.prod();

    // 分割数が少なすぎる時は分割すると逆に効率が悪くなるので分割しない
    if (n_cube_all >= 50) {
        Eigen::Matrix3Xi all_indices3 = Eigen::floor(
                ((expanded_coords.colwise() - valid_min).array() + 1e-8) / ledge).cast<int>();
        Eigen::VectorXi all_indices = three_to_one(all_indices3, n_cube);
        Eigen::Matrix3Xi center_indices3 = Eigen::floor(
                ((center_coords.colwise() - valid_min).array() + 1e-8) / ledge).cast<int>();
        Eigen::VectorXi center_indices = three_to_one(center_indices3, n_cube);

        // atom_indices[i] は i 番目のセルに含まれる原子の indices のリスト
        std::vector<std::vector<int>> atom_indices(n_cube_all);
        for (int i = 0; i < int(indices.size()); ++i) {
            atom_indices[all_indices(i)].push_back(i);
        }

        auto cube_neighbors = get_cube_neighbors(n_cube);

        for (int i = 0; i < n_center; ++i) {
            for (const int cube_index: cube_neighbors[center_indices(i)]) {
                for (const int j: atom_indices[cube_index]) {
                    const double d2 = (expanded_coords.col(j) - center_coords.col(i)).squaredNorm();
                    if (d2 < r2) {
                        const int all_coords_idx = std::get<0>(indices[j]);
                        auto offset = std::get<1>(indices[j]);
                        offset[0] -= int(offset_correction(0, all_coords_idx));
                        offset[1] -= int(offset_correction(1, all_coords_idx));
                        offset[2] -= int(offset_correction(2, all_coords_idx));
                        result[i].emplace_back(FindNearNeighborsResult{
                                all_coords_idx,
                                offset,
                                d2
                        });
                    }
                }
            }
        }
    } else {
        for (int i = 0; i < n_center; ++i) {
            for (int j = 0; j < expanded_coords.cols(); ++j) {
                const double d2 = (expanded_coords.col(j) - center_coords.col(i)).squaredNorm();
                if (d2 < r2) {
                    const int all_coords_idx = std::get<0>(indices[j]);
                    auto offset = std::get<1>(indices[j]);
                    offset[0] -= int(offset_correction(0, all_coords_idx));
                    offset[1] -= int(offset_correction(1, all_coords_idx));
                    offset[2] -= int(offset_correction(2, all_coords_idx));
                    result[i].emplace_back(FindNearNeighborsResult{
                            all_coords_idx,
                            offset,
                            d2
                    });
                }
            }
        }
    }


    return result;
}

std::vector<std::vector<FindNearNeighborsResult>> find_near_neighbors(
        const Structure &structure,
        double r,
        double min_r,
        double tol
) {
    const Eigen::Matrix3Xd frac = structure.lattice.inv_matrix * structure.site_xyz;
    return find_near_neighbors(structure.site_xyz, frac, structure.site_xyz, frac, r, structure.lattice, min_r, tol);
}

// Given the fractional coordinates and the number of repeation needed in each
// direction, maxr, compute the translational bounds in each dimension
std::pair<Eigen::Vector3i, Eigen::Vector3i> get_bounds(
        const Eigen::Matrix3Xd &frac_coords,
        const Eigen::Vector3d &maxr,
        const std::array<bool, 3> &pbc
) {
    Eigen::Vector3d max_fcoords = frac_coords.rowwise().maxCoeff();
    Eigen::Vector3d min_fcoords = frac_coords.rowwise().minCoeff();
    Eigen::Vector3i max_bounds = {1, 1, 1};
    Eigen::Vector3i min_bounds = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        if (pbc[i]) {
            min_bounds[i] = std::floor(min_fcoords[i] - maxr[i] - 1e-8);
            max_bounds[i] = std::ceil(max_fcoords[i] + maxr[i] + 1e-8);
        }
    }
    return {min_bounds, max_bounds};
}

// 逆格子ベクトルを求める
Eigen::Matrix3Xd get_reciprocal_lattice(const Eigen::Matrix3d &lattice) {
    Eigen::Matrix3d recp_lattice;
    for (int i = 0; i < 3; ++i) {
        Eigen::Vector3d cross = lattice.row((i + 1) % 3).cross(lattice.row((i + 2) % 3));
        double prod = lattice.row(i).dot(cross);
        recp_lattice.row(i) = 2.0 * pi * cross.array() / prod;
    }
    return recp_lattice;
}

Eigen::VectorXi three_to_one(const Eigen::Matrix3Xi &label3d, const Eigen::Vector3i &n_cube) {
    return label3d.row(0) * n_cube.y() * n_cube.z() + label3d.row(1) * n_cube.z() + label3d.row(2);
}

int three_to_one1(const Eigen::Vector3i &label, const Eigen::Vector3i &n_cube) {
    return three_to_one(label, n_cube)(0);
}

Eigen::Matrix3Xi one_to_three(const Eigen::VectorXi &label1d, const Eigen::Vector3i &n_cube) {
    int y = n_cube.y(), z = n_cube.z();
    Eigen::Matrix3Xi result(3, label1d.size());
    result.row(0) = label1d.array() / (y * z);
    result.row(1) = label1d.array() / z - label1d.array() / (y * z) * y;
    result.row(2) = label1d.array() - label1d.array() / z * z;
    return result;
}

Eigen::Vector3i one_to_three1(int label, const Eigen::Vector3i &n_cube) {
    return one_to_three(Eigen::VectorXi::Constant(1, label), n_cube).col(0);
}


// Get {cube_index: cube_neighbor_indices} map
std::vector<std::vector<int>> get_cube_neighbors(const Eigen::Vector3i &n_cube) {
    int n_cube_all = n_cube.prod();
    std::vector<std::vector<int>> result(n_cube_all);
    Eigen::Matrix<int, 3, 27> ovector;
    Eigen::Matrix3Xi cube_indices_3d(3, n_cube_all);

    int count = 0;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                ovector.col(count++) = Eigen::Vector3i(i, j, k);
            }
        }
    }

    count = 0;
    for (int i = 0; i < n_cube.x(); ++i) {
        for (int j = 0; j < n_cube.y(); ++j) {
            for (int k = 0; k < n_cube.z(); ++k) {
                cube_indices_3d.col(count++) = Eigen::Vector3i(i, j, k);
            }
        }
    }
    for (int i = 0; i < n_cube_all; ++i) {
        result[i].reserve(27);
        for (int j = 0; j < ovector.cols(); ++j) {
            Eigen::Vector3i index3 = ovector.col(j) + cube_indices_3d.col(i);
            if ((index3.array() >= 0).all() && (index3.array() < n_cube.array()).all()) {
                result[i].push_back(three_to_one1(index3, n_cube));
            }
        }
    }

    return result;
}

void init_near_neighbor(pybind11::module &m) {
    py::class_<NearNeighborInfo>(m, "NearNeighborInfo")
            .def(py::init<int, double, std::array<int, 3>>(),
                 py::arg("site_index"),
                 py::arg("weight"),
                 py::arg("image") = std::array<int, 3>{0, 0, 0})
            .def_property_readonly("site_index", [](const NearNeighborInfo &self) { return self.site_index; })
            .def_property_readonly("weight", [](const NearNeighborInfo &self) { return self.weight; })
            .def_property_readonly("image", [](const NearNeighborInfo &self) { return self.image; });

    py::class_<NearNeighbor, std::shared_ptr<NearNeighbor>>(m, "NearNeighbor")
            .def_property_readonly("structures_allowed", &NearNeighbor::structures_allowed)
            .def_property_readonly("molecules_allowed", &NearNeighbor::molecules_allowed)
            .def("get_all_nn_info", &NearNeighbor::get_all_nn_info);

    py::class_<MinimumDistanceNN, std::shared_ptr<MinimumDistanceNN>, NearNeighbor>(m, "MinimumDistanceNN")
            .def(py::init<double, double, double>(),
                 py::arg("tol") = 0.1,
                 py::arg("cutoff") = 10.0,
                 py::arg("get_all_sites") = false);

    py::class_<CutOffDictNN, std::shared_ptr<CutOffDictNN>, NearNeighbor>(m, "CutOffDictNN")
            .def(py::init<std::optional<py::dict>>(),
                 py::arg("cut_off_dict") = py::none())
            .def_static("from_preset", CutOffDictNN::from_preset);

    m.def("find_near_neighbors", [](
            const Eigen::MatrixX3d &all_coords,
            const Eigen::MatrixX3d &center_coords,
            const double r,
            const Eigen::Vector3i &pbc,
            const Eigen::Matrix3d &lattice,
            const double tol = 1e-8,
            const double min_r = 1.0) {
        const Eigen::Matrix3Xd A = all_coords.transpose();
        const Eigen::Matrix3Xd C = center_coords.transpose();
        const Eigen::Matrix3d L = lattice.transpose();
        const Eigen::Matrix3d L_inv = L.inverse();
        Lattice l;
        l.matrix = L;
        l.inv_matrix = L_inv;
        l.pbc = {pbc.x() != 0, pbc.y() != 0, pbc.z() != 0};
        const auto res = find_near_neighbors(A, L_inv * A, C, L_inv * C, r, l, min_r, tol);

        size_t total = 0;
        for (const auto &vec: res) total += vec.size();
        Eigen::VectorXi res1(total), res2(total);
        Eigen::MatrixX3d res_offset(total, 3);
        Eigen::VectorXd distances(total);
        int i = 0;
        for (int res_i = 0; res_i < int(res.size()); ++res_i) {
            for (const auto &x: res[res_i]) {
                res1(i) = res_i;
                res2(i) = x.all_coords_idx;
                res_offset.row(i) = Eigen::Vector3d(x.image[0], x.image[1], x.image[2]);
                distances(i) = std::sqrt(x.distances2);
                i++;
            }
        }

        return py::make_tuple(res1, res2, res_offset, distances);
    });
}

