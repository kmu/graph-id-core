#include "near_neighbor.h"
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/eval.h>
#include <gtl/phmap.hpp>
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
            if (info.extra) for (const auto &[k, v]: info.extra.value()) d[k] = v;
            inner.append(d);
        }
        arr.append(inner);
    }
    return arr;
}


std::vector<std::vector<NearNeighborInfo>> VoronoiNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto voro = this->get_all_voronoi_polyhedra(structure);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);
    for (int i = 0; i < structure.count; ++i) {
        result[i] = this->extract_nn_info(structure, voro[i]);
    }
    return result;
}

std::vector<NearNeighborInfo>
VoronoiNN::extract_nn_info(const Structure &s, const std::unordered_map<int, VoronoiPolyhedra> &voro) const {
    const auto &targets = this->targets ? this->targets.value() : s.species_strings;
    return extract_nn_info(s, voro, targets);
}

std::vector<NearNeighborInfo>
VoronoiNN::extract_nn_info(const Structure &s,
                           const std::unordered_map<int, VoronoiPolyhedra> &voro,
                           const std::vector<std::string> &targets
) const {
    gtl::flat_hash_set<std::string> target_set(targets.begin(), targets.end());
    std::vector<NearNeighborInfo> result;

    double max_weight = 0;
    for (const auto &[k, v]: voro) max_weight = std::max(max_weight, v[this->weight]);
    for (const auto &[k, v]: voro) {
        if (v[this->weight] > this->tol * max_weight && target_set.contains(s.species_strings[v.site.all_coords_idx])) {
            NearNeighborInfo info;
            info.site_index = v.site.all_coords_idx;
            info.image = v.site.image;
            info.weight = v[this->weight] / max_weight;
            if (this->extra_nn_info) {
                auto pi = v.to_dict(s);
                info.extra = py::dict(py::arg("poly_info") = pi);
                pi.attr("pop")("site");
            }
            result.push_back(std::move(info));
        }
    }

    return result;
}

std::unordered_map<int, VoronoiPolyhedra>
VoronoiNN::get_voronoi_polyhedra(const Structure &structure, int site_index) const {
    if (site_index < 0 || site_index >= structure.count) {
        throw py::index_error("site_index out of range."); // IndexError
    }

    const auto &targets = this->targets ? this->targets.value() : structure.species_strings;
    const Eigen::Matrix3Xd center = structure.site_xyz.col(site_index);
    double max_cutoff = get_max_cutoff(structure);
    double cutoff = this->cutoff;

    while (true) {
        try {
            auto neighbors = find_near_neighbors(
                    structure.site_xyz,
                    structure.lattice.inv_matrix * structure.site_xyz,
                    center,
                    structure.lattice.inv_matrix * center,
                    cutoff,
                    structure.lattice
            )[0];
            std::sort(neighbors.begin(), neighbors.end(), [](const auto &lhs, const auto &rhs) {
                return lhs.distance < rhs.distance;
            });

            Eigen::Matrix3Xd qvoronoi_input(3, neighbors.size());
            for (int i = 0; i < int(neighbors.size()); ++i) {
                qvoronoi_input.col(i) = neighbors[i].xyz(structure);
            }

            // Run the Voronoi tessellation
            auto voro = Voronoi(qvoronoi_input);  // can give seg fault if cutoff is too small

            return this->extract_cell_info(0, structure, neighbors, targets, voro, this->compute_adj_neighbors);
        } catch (py::error_already_set &eas) {
            if (!eas.matches(PyExc_RuntimeError)) {
                throw;
            }
        } catch (std::exception &e) {
            if (std::string(e.what()) !=
                "This structure is pathological, infinite vertex in the Voronoi construction") {
                throw;
            }
        }
        if (cutoff >= max_cutoff) {
            throw std::runtime_error("Error in Voronoi neighbor finding; max cutoff exceeded");
        }
        cutoff = std::min(cutoff * 2, max_cutoff + 0.001);
    }

    throw std::logic_error("unreachable");
}

std::vector<std::unordered_map<int, VoronoiPolyhedra>>
VoronoiNN::get_all_voronoi_polyhedra(const Structure &structure) const {
    if (structure.count == 1) {
        return {this->get_voronoi_polyhedra(structure, 0)};
    }
    assert(structure.count >= 2);

    const auto &targets = this->targets ? this->targets.value() : structure.species_strings;
    double max_cutoff = get_max_cutoff(structure);
    double cutoff = this->cutoff;

    while (true) {
        try {
            const auto neighbors = find_near_neighbors(structure, cutoff);
            assert(int(neighbors.size()) == structure.count);

            gtl::flat_hash_set<std::array<int, 4>, std::hash<std::array<int, 4>>> indices_set;
            std::vector<FindNearNeighborsResult> flat; // (site_index, image[0], image[1], image[2])
            std::vector<int> root_image_index(structure.count, -1);
            for (int i = 0; i < int(neighbors.size()); ++i) {
                for (const auto &nn: neighbors[i]) {
                    std::array<int, 4> arr = {nn.all_coords_idx, nn.image[0], nn.image[1], nn.image[2]};
                    if (indices_set.contains(arr)) continue;
                    if (nn.image[0] == 0 && nn.image[1] == 0 && nn.image[2] == 0) {
                        root_image_index[nn.all_coords_idx] = int(flat.size());
                    }
                    flat.push_back(nn);
                    indices_set.insert(arr);
                }
            }

            Eigen::Matrix3Xd qvoronoi_input(3, flat.size());
            int i = 0;
            for (const auto &nn: flat) {
                qvoronoi_input.col(i++) = nn.xyz(structure);
            }

            for (const int i: root_image_index) assert(0 <= i && i < int(flat.size()));

            auto voro = Voronoi(qvoronoi_input);

            std::vector<std::unordered_map<int, VoronoiPolyhedra>> result(structure.count);
            for (int i = 0; i < structure.count; ++i) {
                result[i] = this->extract_cell_info(root_image_index[i], structure, flat, targets, voro,
                                                    this->compute_adj_neighbors);
            }
            return result;
        } catch (py::error_already_set &eas) {
            if (!eas.matches(PyExc_RuntimeError)) {
                throw;
            }
        } catch (std::exception &e) {
            if (std::string(e.what()) !=
                "This structure is pathological, infinite vertex in the Voronoi construction") {
                throw;
            }
        }
        if (cutoff >= max_cutoff) {
            throw std::runtime_error("Error in Voronoi neighbor finding; max cutoff exceeded");
        }
        cutoff = std::min(cutoff * 2, max_cutoff + 0.001);
    }

    throw std::logic_error("unreachable");
}

std::unordered_map<int, VoronoiPolyhedra> VoronoiNN::extract_cell_info(
        int neighbor_index,
        const Structure &structure,
        const std::vector<FindNearNeighborsResult> &neighbors,
        const std::vector<std::string> &targets,
        const Voronoi &voro,
        bool compute_adj_neighbors
) const {
    assert(0 <= neighbor_index && neighbor_index < int(neighbors.size()));

    // Get the coordinates of every vertex
    const auto _all_vertices = voro.vertices();
    Eigen::Matrix3Xd all_vertices(3, _all_vertices.shape(0));
    for (int i = 0; i < _all_vertices.shape(0); ++i) {
        for (int j = 0; j < 3; ++j) {
            all_vertices(j, i) = _all_vertices.at(i, j);
        }
    }

    // Get the coordinates of the central site
    const Eigen::Vector3d center_coords = neighbors[neighbor_index].xyz(structure);

    // Iterate through all the faces in the tessellation
    std::unordered_map<int, VoronoiPolyhedra> results;
    for (const auto &[_nn, _vind]: voro.ridge_dict()) {
        auto nn = _nn.cast<py::tuple>();
        auto vind = _vind.cast<std::vector<int>>();
        if (nn.contains(neighbor_index)) {
            int other_neighbor_index = nn[0].cast<int>() == neighbor_index ? nn[1].cast<int>() : nn[0].cast<int>();
            assert(0 <= other_neighbor_index && other_neighbor_index < int(neighbors.size()));
            Eigen::Vector3d other_xyz = neighbors[other_neighbor_index].xyz(structure);
            if (std::find(vind.begin(), vind.end(), -1) != vind.end()) {
                if (this->allow_pathological) {
                    continue;
                } else {
                    throw std::runtime_error(
                            "This structure is pathological, infinite vertex in the Voronoi construction");
                }
            }

            // Get the solid angle of the face
            Eigen::Matrix3Xd facets(3, vind.size());
            for (int i = 0; i < int(vind.size()); ++i) {
                facets.col(i) << all_vertices.col(vind[i]);
            }
            double angle = solid_angle(center_coords, facets);

            // Compute the volume of associated with this face
            double volume = 0;
            // qvoronoi returns vertices in CCW order, so I can break
            // the face up in to segments (0,1,2), (0,2,3), ... to compute
            // its area where each number is a vertex size
            for (int i = 1; i < int(vind.size()) - 1; ++i) {
                int j = vind[i];
                int k = vind[i + 1];
                volume += vol_tetra(
                        center_coords,
                        all_vertices.col(vind[0]),
                        all_vertices.col(j),
                        all_vertices.col(k)
                );
            }

            // Compute the distance of the site to the face
            double face_dist = (center_coords - other_xyz).norm() / 2;

            // Compute the area of the face (knowing V=Ad/3)
            double face_area = 3 * volume / face_dist;

            // Compute the normal of the facet
            Eigen::Vector3d normal = other_xyz - center_coords;
            normal /= normal.norm();
            VoronoiPolyhedra v;
            v.site = neighbors[other_neighbor_index];
            v.normal = normal;
            v.solid_angle = angle;
            v.volume = volume;
            v.face_dist = face_dist;
            v.area = face_area;
            v.n_verts = int(vind.size());

            if (compute_adj_neighbors) v.verts = vind;

            results[other_neighbor_index] = v;
        }
    }

    // all sites should have at least two connected ridges in periodic system
    if (results.empty()) {
        throw py::value_error("No Voronoi neighbors found for site - try increasing cutoff");
    }

    // Get only target elements
    gtl::flat_hash_set<std::string> target_set(targets.begin(), targets.end());
    std::unordered_map<int, VoronoiPolyhedra> result_weighted;
    for (const auto &[nn_index, nn_stats]: results) {
        // Check if this is a target site
        py::object nn = structure.py_structure.sites()[nn_stats.site.all_coords_idx].obj;
        if (nn.attr("is_ordered").cast<bool>()) {
            if (target_set.contains(structure.species_strings[nn_stats.site.all_coords_idx])) {
                result_weighted[nn_index] = nn_stats;
            }
        } else { // if nn site is disordered
            for (auto disordered_sp: nn.attr("species")) {
                if (target_set.contains(disordered_sp.attr("formula").cast<std::string>())) {
                    result_weighted[nn_index] = nn_stats;
                }
            }
        }
    }

    // If desired, determine which neighbors are adjacent
    if (compute_adj_neighbors) {
        gtl::flat_hash_map<int, std::vector<int>> adj_neighbors;
        for (const auto &[a_ind, a_nn_info]: result_weighted) {
            gtl::flat_hash_set<int> a_verts(a_nn_info.verts.begin(), a_nn_info.verts.end());
            // Loop over all neighbors that have an index lower that this one
            // The goal here is to exploit the fact that neighbor adjacency is
            // symmetric (if A is adj to B, B is adj to A)
            for (auto &[b_ind, b_nn_info]: result_weighted) {
                if (b_ind > a_ind) continue;
                gtl::flat_hash_set<int> v;
                for (const int n: b_nn_info.verts) {
                    if (a_verts.contains(n)) v.insert(n);
                }
                if (v.size() == 2) {
                    adj_neighbors[a_ind].push_back(b_ind);
                    adj_neighbors[b_ind].push_back(a_ind);
                }
            }
        }

        for (const auto &[k, v]: adj_neighbors) {
            result_weighted[k].adj_neighbors = v;
        }
    }

    return result_weighted;
}

double VoronoiNN::get_max_cutoff(const Structure &structure) {
    // max cutoff is the longest diagonal of the cell + room for noise
    Eigen::Matrix3Xd corners(3, 4);
    corners << 1, 1, 1,
            -1, 1, 1,
            1, -1, 1,
            1, 1, -1;
    Eigen::VectorXd d_corners(4);

    for (auto i = 0; i < 4; ++i) {
        d_corners(i) = (structure.lattice.matrix * corners.col(i)).norm();
    }
    return d_corners.maxCoeff() + 0.01;
}

std::vector<std::vector<NearNeighborInfo>> MinimumDistanceNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto nn = find_near_neighbors(structure, this->cutoff);
    assert(int(nn.size()) == structure.count);
    if (this->get_all_sites) {
        std::vector<std::vector<NearNeighborInfo>> result(structure.count);
        for (int i = 0; i < structure.count; ++i) {
            if (nn[i].empty()) continue;
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (nn[i][j].distance < 1e-8) continue;
                result[i].emplace_back(NearNeighborInfo{
                        nn[i][j].all_coords_idx,
                        nn[i][j].distance,
                        nn[i][j].image,
                        py::dict()
                });
            }
        }
        return result;
    } else {
        std::vector<std::vector<NearNeighborInfo>> result(structure.count);
        for (int i = 0; i < structure.count; ++i) {
            if (nn[i].empty()) continue;
            Eigen::VectorXd d(nn[i].size());
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (nn[i][j].distance < 1e-8) {
                    d(j) = 9999;
                } else {
                    d(j) = nn[i][j].distance;
                }
            }
            const double min_distance = d.minCoeff();
            const double r = (1 + this->tol) * min_distance;
            for (int j = 0; j < int(nn[i].size()); ++j) {
                if (d(j) < r) {
                    result[i].emplace_back(NearNeighborInfo{
                            nn[i][j].all_coords_idx,
                            min_distance / d(j),
                            nn[i][j].image,
                            py::dict()
                    });
                }
            }
        }
        return result;
    }
}

std::vector<std::vector<NearNeighborInfo>> DistanceClusteringNN::get_all_nn_info_cpp(const Structure &structure) const {
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);

    const auto nn = find_near_neighbors(structure, this->cutoff);
    assert(int(nn.size()) == structure.count);

    std::vector<double> cutoff_cluster_list;
    cutoff_cluster_list = get_cutoff_cluster(structure, this->n, this->cutoff, nn);
    if (int(cutoff_cluster_list.size()) <= this->rank_k) {
        return result;
    }

    Eigen::VectorXd d(nn[this->n].size());
    for (int j = 0; j < int(nn[this->n].size()); ++j) {
        if (nn[this->n][j].distance < 1e-8) {
            d(j) = 9999;
        } else {
            d(j) = nn[this->n][j].distance;
        }
    }
    const double min_distance = d.minCoeff();
    for (int j = 0; j < int(nn[this->n].size()); ++j) {
        if (
            (this->rank_k > 0 
            && round(1000*d(j))/1000 <= round(1000*cutoff_cluster_list.at(this->rank_k))/1000
            && round(1000*d(j))/1000 > round(1000*cutoff_cluster_list.at(this->rank_k-1))/1000)
            || (this->rank_k == 0 
            && round(1000*d(j))/1000 <= round(1000*cutoff_cluster_list.at(this->rank_k))/1000)
        ) {
            result[this->n].emplace_back(NearNeighborInfo{
                    nn[this->n][j].all_coords_idx,
                    min_distance / d(j), // min_d / d(j) になぜしていたか分からないが、これ以降でdistの値を使わないので保留。
                    nn[this->n][j].image,
                    py::dict()
            });
        }
    }
    return result;
}

std::vector<double> DistanceClusteringNN::get_cutoff_cluster(const Structure &structure, int n, double cutoff, const auto &nn) const {
    // サイトごとに3つの閾値を決定する

    // std::vector<double> max_dist_list = {0.0, 1.2, 3.1, 5.2};
    assert(n < nn.size());

    int nn_size = int(nn[n].size());
    std::vector<std::vector<double>> distance_vec(nn_size, std::vector<double>(2));
    
    int count = 0;
    for (const auto &neighbor: nn[n]) {
        // const auto nn = neighs_dists
        double dist = neighbor.distance;
        distance_vec.at(count).at(0) = dist;
        distance_vec.at(count).at(1) = 0.0;
        count++;
    }

    std::vector<int> clustering_labels;
    
    //  py::scoped_interpreter guard{}; // Pythonインタープリタを開始

    py::object sklearn = py::module::import("sklearn.cluster");
    py::object DBSCAN = sklearn.attr("DBSCAN");

    // py::array_t<double> data;
    // Call DBSCAN
    // py::object dbscan = DBSCAN(0.5, 2);
    py::object dbscan = DBSCAN(py::arg("eps")=0.5, py::arg("min_samples")=2);
    dbscan.attr("fit")(distance_vec);
    py::list labels_py = dbscan.attr("labels_").cast<py::list>();
    std::vector<int> labels;
    for (const auto& item : labels_py) {
        labels.push_back(item.cast<int>());
    }

    int max_label = *std::max_element(begin(labels), end(labels));
    std::vector<double> max_dist_list(max_label+1);
    for (int label_number = 0; label_number <= max_label; label_number++) {
        double max_dist = 0.0;
        for (int i = 0; i < int(labels.size()); i++) {
            int label = labels.at(i);
            std::vector<double> distance = distance_vec.at(i);
            if (label == label_number) {
                max_dist = std::max(max_dist, distance.at(0));
            }
        }
        max_dist_list.at(label_number) = max_dist;
    }

    std::sort(max_dist_list.begin(), max_dist_list.end());
    return max_dist_list;
}


std::vector<std::vector<NearNeighborInfo>> MinimumOKeeffeNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto neighs_dists = find_near_neighbors(structure, this->cutoff);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);

    // get_okeeffe_distance_prediction をあらかじめ計算しておく
    py::object get_okeeffe_distance_prediction = py::module_::import("pymatgen.analysis.local_env").attr(
            "get_okeeffe_distance_prediction");
    gtl::flat_hash_map<std::pair<std::string, std::string>, double> okeeffe_distance_prediction;
    std::vector<std::string> elem(structure.count);
    for (int site_i = 0; site_i < structure.count; site_i++) {
        try {
            elem[site_i] = structure.py_structure.obj[py::int_(site_i)].attr("specie").attr(
                    "element").cast<std::string>();
        } catch (py::error_already_set &_) {
            elem[site_i] = structure.py_structure.obj[py::int_(site_i)].attr("species_string").cast<std::string>();
        }
    }

    std::vector<std::string> elem_uniq(elem.begin(), elem.end());
    std::sort(elem_uniq.begin(), elem_uniq.end());
    elem_uniq.erase(std::unique(elem_uniq.begin(), elem_uniq.end()), elem_uniq.end());
    for (int i = 0; i < int(elem_uniq.size()); i++) {
        for (int j = i; j < int(elem_uniq.size()); j++) {
            auto d = get_okeeffe_distance_prediction(elem_uniq[i], elem_uniq[j]).cast<double>();
            okeeffe_distance_prediction[std::make_pair(elem_uniq[i], elem_uniq[j])] = d;
            okeeffe_distance_prediction[std::make_pair(elem_uniq[j], elem_uniq[i])] = d;
        }
    }


    for (int site_i = 0; site_i < structure.count; site_i++) {
        if (neighs_dists[site_i].empty()) continue;

        std::vector<double> reldists_neighs;
        reldists_neighs.reserve(neighs_dists[site_i].size());
        for (const auto &nn: neighs_dists[site_i]) {
            double dist = nn.distance;
            if (nn.distance < 1e-8) dist = 1e9; // 距離が 0 のときは同じサイトの組なので無視するために十分大きな値にする
            reldists_neighs.push_back(
                    dist / okeeffe_distance_prediction[std::make_pair(elem[site_i], elem[nn.all_coords_idx])]);
        }

        double min_reldist = reldists_neighs[0];
        for (double d: reldists_neighs) if (min_reldist > d) min_reldist = d;
        for (int i = 0; i < int(neighs_dists[site_i].size()); i++) {
            if (reldists_neighs[i] < min_reldist * (1 + this->tol)) {
                result[site_i].emplace_back(NearNeighborInfo{
                        neighs_dists[site_i][i].all_coords_idx,
                        min_reldist / reldists_neighs[i],
                        neighs_dists[site_i][i].image,
                        py::dict()
                });
            }
        }
    }

    return result;
}

std::vector<std::vector<NearNeighborInfo>> CrystalNN::get_all_nn_info_cpp(const Structure &structure) const {
    auto all_nn_data = get_all_nn_data(structure);
    std::vector<std::vector<NearNeighborInfo>> result;
    result.reserve(structure.count);

    for (int i = 0; i < structure.count; i++) {
        auto &nn_data = all_nn_data[i];
        if (!this->weighted_cn) {
            int max_key = nn_data.cn_weights.begin()->first;
            for (const auto &[key, value]: nn_data.cn_weights) {
                if (nn_data.cn_weights[max_key] < value) {
                    max_key = key;
                }
            }
            auto &nn = nn_data.cn_nninfo[max_key];
            for (auto &entry: nn) entry.weight = 1;
            result.push_back(nn);
        } else {
            for (auto &entry: nn_data.all_nninfo) {
                double weight = 0;
                for (auto &[cn, cn_nninfo]: nn_data.cn_nninfo) {
                    for (auto &cn_entry: cn_nninfo) {
                        if (cn_entry.site_index == entry.site_index && cn_entry.image == entry.image) {
                            weight += nn_data.cn_weights[cn];
                        }
                    }
                }
                entry.weight = weight;
            }
            result.push_back(nn_data.all_nninfo);
        }
    }

    return result;
}

std::vector<CrystalNN::NNData> CrystalNN::get_all_nn_data(const Structure &structure, int length) const {
    std::vector<CrystalNN::NNData> result;
    result.reserve(structure.count);
    py::object py_structure = structure.py_structure.obj;

    if (length == 0) length = this->fingerprint_length;

    // get base VoronoiNN targets
    auto vnn = VoronoiNN(0, std::nullopt, this->search_cutoff, false, "solid_angle");
    auto all_voronoi = vnn.get_all_voronoi_polyhedra(structure);

    std::vector<double> site_radius(structure.count);
    std::vector<double> site_default_radius(structure.count);
    for (int site_i = 0; site_i < structure.count; site_i++) {
        site_radius[site_i] = get_radius(py_structure[py::int_(site_i)]);
        site_default_radius[site_i] = get_default_radius(py_structure[py::int_(site_i)]);
    }

    for (int site_i = 0; site_i < structure.count; site_i++) {
        auto voronoi = all_voronoi[site_i];
        std::vector<NearNeighborInfo> nn;

        if (this->cation_anion) {
            std::vector<std::string> targets;
            auto m_oxi = py_structure[py::int_(site_i)].attr("specie").attr("oxi_state");
            for (auto site: py_structure) {
                py::object oxi_state = py::getattr(site.attr("specie"), "oxi_state", py::none());
                if (!oxi_state.is_none() && oxi_state.cast<double>() * m_oxi.cast<double>() <= 0) {
                    targets.push_back(site.attr("species_string").cast<std::string>());
                }
            }
            if (targets.empty()) {
                throw py::value_error("No valid targets for site within cation_anion constraint!");
            }
            nn = vnn.extract_nn_info(structure, voronoi, targets);
        } else {
            nn = vnn.extract_nn_info(structure, voronoi);
        }


        // solid angle weights can be misleading in open / porous structures
        // adjust weights to correct for this behavior
        if (this->porous_adjustment) {
            for (auto &x: nn) {
                x.weight *= x.extra.value()["poly_info"]["solid_angle"].cast<double>() /
                            x.extra.value()["poly_info"]["area"].cast<double>();
            }
        }

        if (this->x_diff_weight > 0) {
            for (auto &entry: nn) {
                auto x1 = py_structure[py::int_(site_i)].attr("specie").attr("X").cast<double>();
                auto x2 = py_structure[py::int_(entry.site_index)].attr("specie").attr("X").cast<double>();
                double chemical_weight = 0;
                if (std::isnan(x1) || std::isnan(x2)) {
                    chemical_weight = 1;
                } else {
                    // note: 3.3 is max deltaX between 2 elements
                    chemical_weight = 1 + this->x_diff_weight * std::sqrt(std::abs(x1 - x2) / 3.3);
                }
                entry.weight *= chemical_weight;
            }
        }

        // sort nearest neighbors from highest to lowest weight
        std::sort(nn.begin(), nn.end(), [](const auto &lhs, const auto &rhs) {
            return lhs.weight > rhs.weight;
        });
        if (nn[0].weight == 0) {
            NNData data;
            data.cn_weights[0] = 1;
            data.cn_nninfo[0] = {};
            transform_to_length(data, length);
            result.push_back(std::move(data));
            continue;
        }

        // renormalize weights so the highest weight is 1.0
        double highest_weight = nn[0].weight;
        for (auto &entry: nn) {
            entry.weight /= highest_weight;
        }

        // adjust solid angle weights based on distance
        if (this->distance_cutoffs.first != 0 && this->distance_cutoffs.second != 0) {
            double r1 = site_radius[site_i];
            for (auto &entry: nn) {
                double r2 = site_radius[entry.site_index];
                double d = 0;
                if (r1 > 0 && r2 > 0) {
                    d = r1 + r2;
                } else {
                    warn("CrystalNN: cannot locate an appropriate radius, "
                         "covalent or atomic radii will be used, this can lead "
                         "to non-optimal results.");
                    d = site_default_radius[site_i] + site_default_radius[entry.site_index];
                }

                double dist = (structure.site_xyz.col(site_i) - entry.xyz(structure)).norm();
                double dist_weight = 0;

                double cutoff_low = d + this->distance_cutoffs.first;
                double cutoff_high = d + this->distance_cutoffs.second;
                if (dist <= cutoff_low) {
                    dist_weight = 1;
                } else if (dist < cutoff_high) {
                    dist_weight = (std::cos((dist - cutoff_low) / (cutoff_high - cutoff_low) * pi) + 1) * 0.5;
                }
                entry.weight *= dist_weight;
            }
        }

        // sort nearest neighbors from highest to lowest weight
        std::sort(nn.begin(), nn.end(), [](const auto &lhs, const auto &rhs) {
            return lhs.weight > rhs.weight;
        });
        if (nn[0].weight == 0) {
            NNData data;
            data.cn_weights[0] = 1;
            data.cn_nninfo[0] = {};
            transform_to_length(data, length);
            result.push_back(std::move(data));
            continue;
        }

        for (auto &entry: nn) {
            entry.weight = std::round(entry.weight * 1000) / 1000;
            if (entry.extra) entry.extra->attr("pop")("poly_info", py::none());
        }

        // remove entries with no weight
        nn.erase(std::remove_if(nn.begin(), nn.end(), [](const auto &entry) {
            return entry.weight == 0;
        }), nn.end());

        // get the transition distances, i.e. all distinct weights
        std::vector<double> dist_bins;
        for (const auto &entry: nn) {
            if (dist_bins.empty() || dist_bins[dist_bins.size() - 1] != entry.weight) {
                dist_bins.push_back(entry.weight);
            }
        }
        dist_bins.push_back(0);

        // main algorithm to determine fingerprint from bond weights
        NNData data;
        for (int idx = 0; idx < int(dist_bins.size()); idx++) {
            double val = dist_bins[idx];
            if (val != 0) {
                std::vector<NearNeighborInfo> nn_info;
                for (const auto &entry: nn) {
                    if (entry.weight >= val) {
                        nn_info.push_back(entry);
                    }
                }
                int cn = int(nn_info.size());
                data.cn_nninfo[cn] = nn_info;
                data.cn_weights[cn] = semicircle_integral(dist_bins, idx);
            }
        }

        // add zero coord
        double cn0_weight = 1;
        for (const auto &[_, v]: data.cn_weights) {
            cn0_weight -= v;
        }

        if (cn0_weight > 0) {
            data.cn_weights[0] = cn0_weight;
            data.cn_nninfo[0] = {};
        }

        data.all_nninfo = nn;
        transform_to_length(data, length);
        result.push_back(std::move(data));
    }

    return result;
}

double CrystalNN::semicircle_integral(const std::vector<double> &dist_bins, int idx) {
    double r = 1;
    double x1 = dist_bins[idx];
    double x2 = dist_bins[idx + 1];
    double area1, area2;

    if (dist_bins[idx] == 1) {
        area1 = 0.25 * pi * r * r;
    } else {
        area1 = 0.5 * ((x1 * std::sqrt(r * r - x1 * x1)) + (r * r * std::atan(x1 / std::sqrt(r * r - x1 * x1))));
    }

    area2 = 0.5 * ((x2 * std::sqrt(r * r - x2 * x2)) + (r * r * std::atan(x2 / std::sqrt(r * r - x2 * x2))));

    return (area1 - area2) / (0.25 * pi * r * r);
}

void CrystalNN::transform_to_length(CrystalNN::NNData &nn_data, int length) {
    if (length == 0) return;
    for (int cn = 0; cn < length; ++cn) {
        // if (!nn_data.cn_weights.contains(cn)) {
        if (nn_data.cn_weights.find(cn) == nn_data.cn_weights.end()){
            nn_data.cn_weights[cn] = 0;
            nn_data.cn_nninfo[cn] = {};
        }
    }
}

std::vector<std::vector<NearNeighborInfo>> CutOffDictNN::get_all_nn_info_cpp(const Structure &structure) const {
    const auto nn = find_near_neighbors(structure, this->max_cut_off);
    assert(int(nn.size()) == structure.count);

    std::vector<std::vector<NearNeighborInfo>> result(structure.count);
    for (int i = 0; i < structure.count; ++i) {
        if (nn[i].empty()) continue;
        for (int j = 0; j < int(nn[i].size()); ++j) {
            if (nn[i][j].distance < 1e-8) continue;
            const auto key = std::make_pair(structure.species_strings[i],
                                            structure.species_strings[nn[i][j].all_coords_idx]);
            const auto it = this->cut_off_dict.find(key);
            if (it == this->cut_off_dict.end()) continue;
            double distance = nn[i][j].distance;
            if (distance < it->second) {
                result[i].emplace_back(NearNeighborInfo{nn[i][j].all_coords_idx, distance, nn[i][j].image, py::dict()});
            }
        }
    }

    return result;
}

std::vector<std::vector<NearNeighborInfo>> BrunnerNN_reciprocal::get_all_nn_info_cpp(const Structure &structure) const {
    auto all_nn = find_near_neighbors(structure, this->cutoff);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);

    for (int site_i = 0; site_i < structure.count; site_i++) {
        auto &nn = all_nn[site_i];
        std::vector<double> ds;
        ds.reserve(nn.size());
        for (auto & i : nn) {
            if (i.distance < 1e-8) continue;
            ds.push_back(i.distance);
        }
        std::sort(ds.begin(), ds.end());
        if (ds.size() < 2) continue;

        std::vector<double> ns(ds.size() - 1);
        for (int i = 0; i < int(ns.size()); i++) ns[i] = 1 / ds[i] - 1 / ds[i + 1];

        size_t d_max_idx = std::max_element(ns.begin(), ns.end()) - ns.begin();
        double d_max = ds[d_max_idx];

        for (auto &nni: nn) {
            if (nni.distance < 1e-8) continue;
            if (nni.distance < d_max + this->tol) {
                result[site_i].emplace_back(NearNeighborInfo{
                        nni.all_coords_idx,
                        ds[0] / nni.distance,
                        nni.image,
                        py::dict()
                });
            }
        }
    }

    return result;
}

std::vector<std::vector<NearNeighborInfo>> BrunnerNN_relative::get_all_nn_info_cpp(const Structure &structure) const {
    auto all_nn = find_near_neighbors(structure, this->cutoff);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);

    for (int site_i = 0; site_i < structure.count; site_i++) {
        auto &nn = all_nn[site_i];
        std::vector<double> ds;
        ds.reserve(nn.size());
        for (auto & i : nn) {
            if (i.distance < 1e-8) continue;
            ds.push_back(i.distance);
        }
        std::sort(ds.begin(), ds.end());
        if (ds.size() < 2) continue;

        std::vector<double> ns(ds.size() - 1);
        for (int i = 0; i < int(ns.size()); i++) ns[i] = ds[i + 1] / ds[i];

        size_t d_max_idx = std::max_element(ns.begin(), ns.end()) - ns.begin();
        double d_max = ds[d_max_idx];

        for (auto &nni: nn) {
            if (nni.distance < 1e-8) continue;
            if (nni.distance < d_max + this->tol) {
                result[site_i].emplace_back(NearNeighborInfo{
                        nni.all_coords_idx,
                        ds[0] / nni.distance,
                        nni.image,
                        py::dict()
                });
            }
        }
    }

    return result;
}

std::vector<std::vector<NearNeighborInfo>> BrunnerNN_real::get_all_nn_info_cpp(const Structure &structure) const {
    auto all_nn = find_near_neighbors(structure, this->cutoff);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);

    for (int site_i = 0; site_i < structure.count; site_i++) {
        auto &nn = all_nn[site_i];
        std::vector<double> ds;
        ds.reserve(nn.size());
        for (auto & i : nn) {
            if (i.distance < 1e-8) continue;
            ds.push_back(i.distance);
        }
        std::sort(ds.begin(), ds.end());
        if (ds.size() < 2) continue;

        std::vector<double> ns(ds.size() - 1);
        for (int i = 0; i < int(ns.size()); i++) ns[i] = ds[i + 1] - ds[i];

        size_t d_max_idx = std::max_element(ns.begin(), ns.end()) - ns.begin();
        double d_max = ds[d_max_idx];

        for (auto &nni: nn) {
            if (nni.distance < 1e-8) continue;
            if (nni.distance < d_max + this->tol) {
                result[site_i].emplace_back(NearNeighborInfo{
                        nni.all_coords_idx,
                        ds[0] / nni.distance,
                        nni.image,
                        py::dict()
                });
            }
        }
    }

    return result;
}

std::vector<std::vector<NearNeighborInfo>> EconNN::get_all_nn_info_cpp(const Structure &structure) const {
    auto all_nn = find_near_neighbors(structure, this->cutoff);
    std::vector<std::vector<NearNeighborInfo>> result(structure.count);
    std::vector<double> oxi_state(structure.count);
    std::vector<bool> oxi_state_valid(structure.count);

    // oxi_state を取得
    for (int site_i = 0; site_i < structure.count; site_i++) {
        auto site = structure.py_structure.obj[py::int_(site_i)];
        auto oxi = py::getattr(site.attr("specie"), "oxi_state", py::none());
        if (!oxi.is_none()) {
            oxi_state[site_i] = oxi.cast<double>();
            oxi_state_valid[site_i] = true;
        }
    }

    // site_radius をあらかじめ計算
    std::vector<double> site_radius(structure.count);
    if (this->use_fictive_radius) {
        auto py_structure = structure.py_structure.obj;
        for (int site_i = 0; site_i < structure.count; site_i++) {
            double r = get_radius(py_structure[py::int_(site_i)]);
            if (r > 0) {
                site_radius[site_i] = r;
            } else {
                site_radius[site_i] = get_default_radius(py_structure[py::int_(site_i)]);
            };
        }
    }

    for (int site_i = 0; site_i < structure.count; site_i++) {
        if (all_nn[site_i].empty()) continue;

        auto &nn = all_nn[site_i];

        if (this->cation_anion && oxi_state_valid[site_i]) {
            if (oxi_state[site_i] >= 0) {
                nn.erase(std::remove_if(nn.begin(), nn.end(), [&](const auto &x) {
                    return oxi_state[x.all_coords_idx] > 0;
                }), nn.end());
            } else {
                nn.erase(std::remove_if(nn.begin(), nn.end(), [&](const auto &x) {
                    return oxi_state[x.all_coords_idx] <= 0;
                }), nn.end());
            }
        }

        nn.erase(std::remove_if(nn.begin(), nn.end(), [&](const auto &x) {
            return x.distance < 1e-8; // 距離が 0 のときは同じサイトの組なので無視する
        }), nn.end());


        Eigen::VectorXd fir(nn.size());
        if (this->use_fictive_radius) {
            for (int i = 0; i < int(nn.size()); i++) {
                // calculate fictive ionic radii
                fir[i] = nn[i].distance *
                         (site_radius[site_i] / (site_radius[site_i] + site_radius[nn[i].all_coords_idx]));
            }
        } else {
            for (int i = 0; i < int(nn.size()); i++) {
                fir[i] = nn[i].distance;
            }
        }

        // calculate mean fictive ionic radius
        double mefir = get_mean_fictive_ionic_radius(fir, fir.minCoeff());

        // iteratively solve MEFIR; follows equation 4 in Hoppe's EconN paper
        double prev_mefir = 1e100;
        while (std::abs(prev_mefir - mefir) > 1e-4) {
            prev_mefir = mefir;
            mefir = get_mean_fictive_ionic_radius(fir, mefir);
        }

        Eigen::VectorXd w = Eigen::exp(1 - (fir.array() / mefir).pow(6));
        for (int i = 0; i < int(nn.size()); i++) {
            if (nn[i].distance < this->cutoff && w[i] > this->tol) {
                result[site_i].emplace_back(NearNeighborInfo{
                        nn[i].all_coords_idx,
                        w[i],
                        nn[i].image,
                        py::dict()
                });
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
        std::vector<std::vector<FindNearNeighborsResult>> result(res.size());
        for (size_t i = 0; i < res.size(); ++i) {
            for (const auto &x: res[i]) {
                if (x.distance <= r) {
                    result[i].emplace_back(x);
                }
            }
        }
        return result;
    }


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
                    const double d = (expanded_coords.col(j) - center_coords.col(i)).norm();
                    if (d < r) {
                        const int all_coords_idx = std::get<0>(indices[j]);
                        auto offset = std::get<1>(indices[j]);
                        offset[0] -= int(offset_correction(0, all_coords_idx));
                        offset[1] -= int(offset_correction(1, all_coords_idx));
                        offset[2] -= int(offset_correction(2, all_coords_idx));
                        result[i].emplace_back(FindNearNeighborsResult{
                                all_coords_idx,
                                offset,
                                d
                        });
                    }
                }
            }
        }
    } else {
        for (int i = 0; i < n_center; ++i) {
            for (int j = 0; j < expanded_coords.cols(); ++j) {
                const double d = (expanded_coords.col(j) - center_coords.col(i)).norm();
                if (d < r) {
                    const int all_coords_idx = std::get<0>(indices[j]);
                    auto offset = std::get<1>(indices[j]);
                    offset[0] -= int(offset_correction(0, all_coords_idx));
                    offset[1] -= int(offset_correction(1, all_coords_idx));
                    offset[2] -= int(offset_correction(2, all_coords_idx));
                    result[i].emplace_back(FindNearNeighborsResult{
                            all_coords_idx,
                            offset,
                            d
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

double solid_angle(Eigen::Vector3d center, Eigen::Matrix3Xd coords) {
    Eigen::Matrix3Xd r = coords.colwise() - center;
    Eigen::VectorXd r_norm = r.colwise().norm();

    // Compute the solid angle for each tetrahedron that makes up the facet
    // Following: https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
    double angle = 0;
    for (int i = 1; i < r.cols() - 1; ++i) {
        int j = i + 1;
        double tp = std::abs(r.col(0).dot(r.col(i).cross(r.col(j))));
        double de = r_norm[0] * r_norm[i] * r_norm[j] +
                    r_norm[j] * r.col(0).dot(r.col(i)) +
                    r_norm[i] * r.col(0).dot(r.col(j)) +
                    r_norm[0] * r.col(i).dot(r.col(j));
        double _angle = de == 0 ? (tp > 0 ? 0.5 * pi : -0.5 * pi) : std::atan(tp / de);
        angle += 2 * (_angle > 0 ? _angle : _angle + pi);
    }
    return angle;
}

double vol_tetra(Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3) {
    return std::abs((v3 - v0).dot((v1 - v0).cross(v2 - v0))) / 6;
}

double get_default_radius(py::object site) {
    py::object CovalentRadius =
    py::module_::import("pymatgen.analysis.molecule_structure_comparator").attr("CovalentRadius");
    try {
        return CovalentRadius.attr("radius")[site.attr("specie").attr("symbol")].cast<double>();
    } catch (py::error_already_set &e) {
        return site.attr("specie").attr("atomic_radius").cast<double>();
    }
}

double get_radius(py::object site) {
    py::dict scope(py::arg("site") = site);
    py::exec(R"(
def _get_radius(site):
    """
    An internal method to get the expected radius for a site with
    oxidation state.

    Args:
        site: (Site)

    Returns:
        Oxidation-state dependent radius: ionic, covalent, or atomic.
        Returns 0 if no oxidation state or appropriate radius is found.
    """
    if hasattr(site.specie, "oxi_state"):
        el = site.specie.element
        oxi = site.specie.oxi_state

        if oxi == 0:
            return _get_default_radius(site)

        if oxi in el.ionic_radii:
            return el.ionic_radii[oxi]

        # e.g., oxi = 2.667, average together 2+ and 3+ radii
        if int(math.floor(oxi)) in el.ionic_radii and int(math.ceil(oxi)) in el.ionic_radii:
            oxi_low = el.ionic_radii[int(math.floor(oxi))]
            oxi_high = el.ionic_radii[int(math.ceil(oxi))]
            x = oxi - int(math.floor(oxi))
            return (1 - x) * oxi_low + x * oxi_high

        if oxi > 0 and el.average_cationic_radius > 0:
            return el.average_cationic_radius

        if el.average_anionic_radius > 0 > oxi:
            return el.average_anionic_radius

    else:
        import warnings
        warnings.warn(
            "No oxidation states specified on sites! For better results, set "
            "the site oxidation states in the structure."
        )
    return 0

ret = _get_radius(site)
    )", py::globals(), scope);
    return scope["ret"].cast<double>();
}

double get_mean_fictive_ionic_radius(const Eigen::VectorXd &f, double min_fir) {
    Eigen::VectorXd v = Eigen::exp(1 - (f.array() / min_fir).pow(6));
    return (f.array() * v.array()).sum() / v.sum();
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

    py::class_<VoronoiPolyhedra>(m, "VoronoiPolyhedra");

    py::class_<VoronoiNN, std::shared_ptr<VoronoiNN>, NearNeighbor>(m, "VoronoiNN")
            .def(py::init<double, std::optional<std::vector<std::string>>, double, bool, std::string, bool, bool>(),
                 py::arg("tol") = 0,
                 py::arg("targets") = py::none(),
                 py::arg("cutoff") = 13.0,
                 py::arg("allow_pathological") = false,
                 py::arg("weight") = "solid_angle",
                 py::arg("extra_nn_info") = true,
                 py::arg("compute_adj_neighbors") = true)
            .def("get_voronoi_polyhedra", [](VoronoiNN &self, const PymatgenStructure &s, int n) {
                Structure ss(s);
                py::dict ret;
                for (const auto &[key, val]: self.get_voronoi_polyhedra(ss, n)) {
                    ret[py::int_(key)] = val.to_dict(ss);
                }
                return ret;
            })
            .def("get_all_voronoi_polyhedra", [](VoronoiNN &self, const PymatgenStructure &s) {
                Structure ss(s);
                py::list ret;
                for (const auto &res: self.get_all_voronoi_polyhedra(Structure(s))) {
                    py::dict d;
                    for (const auto &[key, val]: res) {
                        d[py::int_(key)] = val.to_dict(ss);
                    }
                    ret.append(d);
                }
                return ret;
            });

    py::class_<MinimumDistanceNN, std::shared_ptr<MinimumDistanceNN>, NearNeighbor>(m, "MinimumDistanceNN")
            .def(py::init<double, double, double>(),
                 py::arg("tol") = 0.1,
                 py::arg("cutoff") = 10.0,
                 py::arg("get_all_sites") = false);

    py::class_<MinimumOKeeffeNN, std::shared_ptr<MinimumOKeeffeNN>, NearNeighbor>(m, "MinimumOKeeffeNN")
            .def(py::init<double, double>(),
                 py::arg("tol") = 0.1,
                 py::arg("cutoff") = 10.0);
    
    py::class_<DistanceClusteringNN, std::shared_ptr<DistanceClusteringNN>, NearNeighbor>(m, "DistanceClusteringNN")
            .def(py::init<double, int, int, double>(),
                 py::arg("tol") = 0.1,
                 py::arg("n") = 0,
                 py::arg("rank_k") = 3,
                 py::arg("cutoff") = 6.0);

    py::class_<CrystalNN, std::shared_ptr<CrystalNN>, NearNeighbor>(m, "CrystalNN")
            .def(py::init<bool, bool, std::pair<double, double>, double, bool, double, int>(),
                 py::arg("weighted_cn") = false,
                 py::arg("cation_anion") = false,
                 py::arg("distance_cutoffs") = std::pair<double, double>{0.5, 1},
                 py::arg("x_diff_weight") = 3.0,
                 py::arg("porous_adjustment") = true,
                 py::arg("search_cutoff") = 7,
                 py::arg("fingerprint_length") = 0);

    py::class_<CutOffDictNN, std::shared_ptr<CutOffDictNN>, NearNeighbor>(m, "CutOffDictNN")
            .def(py::init<std::optional<py::dict>>(),
                 py::arg("cut_off_dict") = py::none())
            .def_static("from_preset", CutOffDictNN::from_preset);

    py::class_<BrunnerNN_reciprocal, std::shared_ptr<BrunnerNN_reciprocal>, NearNeighbor>(m, "BrunnerNN_reciprocal")
            .def(py::init<double, double>(),
                 py::arg("tol") = 1e-4,
                 py::arg("cutoff") = 8.0);

    py::class_<BrunnerNN_relative, std::shared_ptr<BrunnerNN_relative>, NearNeighbor>(m, "BrunnerNN_relative")
            .def(py::init<double, double>(),
                 py::arg("tol") = 1e-4,
                 py::arg("cutoff") = 8.0);

    py::class_<BrunnerNN_real, std::shared_ptr<BrunnerNN_real>, NearNeighbor>(m, "BrunnerNN_real")
            .def(py::init<double, double>(),
                 py::arg("tol") = 1e-4,
                 py::arg("cutoff") = 8.0);

    py::class_<EconNN, std::shared_ptr<EconNN>, NearNeighbor>(m, "EconNN")
            .def(py::init<double, double, bool, bool>(),
                 py::arg("tol") = 0.2,
                 py::arg("cutoff") = 10.0,
                 py::arg("cation_anion") = false,
                 py::arg("use_fictive_radius") = false);

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
                distances(i) = x.distance;
                i++;
            }
        }

        return py::make_tuple(res1, res2, res_offset, distances);
    });
}

