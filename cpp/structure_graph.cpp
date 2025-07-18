#include "structure_graph.h"
#include <Eigen/Core>
#include <Eigen/LU>

std::string blake2b(const std::string &s) {
    py::object hashlib = py::module_::import("hashlib");

    return hashlib.attr("blake2b")(py::bytes(s)).attr("hexdigest")().cast<std::string>();
}

std::string blake2b(const std::string &s, int digest_size) {
    py::object hashlib = py::module_::import("hashlib");

    return hashlib.attr("blake2b")(py::bytes(s), py::arg("digest_size") = digest_size).attr(
            "hexdigest")().cast<std::string>();
}

StructureGraph StructureGraph::with_local_env_strategy(
        const std::shared_ptr<const Structure> &structure,
        const NearNeighbor &strategy
) {
    auto sg = StructureGraph::from_empty_graph(structure);
    const auto &nn = strategy.get_all_nn_info_cpp(*structure);
    assert(int(nn.size()) == structure->count);
    for (int from = 0; from < int(nn.size()); ++from) {
        for (size_t i = 0; i < nn[from].size(); ++i) {
            const auto &nni = nn[from][i];
            sg.add_edge(from, {0, 0, 0}, nni.site_index, nni.image, nni.weight);
            sg.add_edge(nni.site_index, nni.image, from, {0, 0, 0}, nni.weight);
        }
    }
    sg.set_cc_diameter();
    return sg;
}

StructureGraph StructureGraph::with_individual_state_comp_strategy(
        const std::shared_ptr<const Structure> &structure,
        // const std::shared_ptr<const StructureGraph> &sg,
        StructureGraph &sg,
        // const NearNeighbor &strategy,
        int n,
        int rank_k,
        double cutoff
) {
    // auto sg = StructureGraph::from_empty_graph(structure);
    // assert(strategy == DistanceClusteringNN);
    const auto &nn = DistanceClusteringNN(0.1, n, rank_k, cutoff).get_all_nn_info_cpp(*structure);
    assert(int(nn.size()) == structure->count);
    // for (int from = 0; from < int(nn.size()); ++from) {
    for (size_t i = 0; i < nn[n].size(); ++i) {
        const auto &nni = nn[n][i];
        sg.add_edge(n, {0, 0, 0}, nni.site_index, nni.image, nni.weight);
        sg.add_edge(nni.site_index, nni.image, n, {0, 0, 0}, nni.weight);
    }
    // }
    sg.set_cc_diameter();
    return sg;
}


StructureGraph StructureGraph::from_empty_graph(const std::shared_ptr<const Structure> &structure) {
    const auto n = structure->count;
    if (n == 0) throw py::value_error("Structure must have at least one site.");
    return StructureGraph{
            structure,
            std::vector<std::vector<NearNeighborInfo>>(n),
            {},
            std::vector<std::string>(n),
            {},
            {},
            {},
    };
}

StructureGraph StructureGraph::from_py(py::object py_sg) {
    py::object py_structure = py_sg.attr("structure");
    auto structure = py_structure.cast<PymatgenStructure>();
    auto s_ptr = std::make_shared<const Structure>(structure);
    auto sg = StructureGraph::from_empty_graph(s_ptr);

    // PythonのStructureGraphからエッジを取得
    py::list edges = py_sg.attr("graph").attr("edges");

    // エッジを追加
    for (auto& edge : edges) {
        // edgeの詳細を取得し、C++のStructureGraphに追加します。
        // ここではedgeがstd::pair<int, int>型と仮定しています。
        std::tuple<int, int, int> e = edge.cast<std::tuple<int, int, int>>();
        py::dict edges_property = py_sg.attr("graph").attr("edges")[edge];
        std::array<int, 3> to_jimage_array;
        auto to_jimage = edges_property["to_jimage"].cast<py::list>();
        for (size_t i = 0; i < to_jimage_array.size(); ++i) {
            to_jimage_array[i] = to_jimage[i].cast<int>();
        }
        
        // 双方向のエッジを追加（with_local_env_strategyと同じ方法）
        sg.add_edge(std::get<0>(e), {0, 0, 0}, std::get<1>(e), to_jimage_array, std::get<2>(e));
        
        // 逆向きのimageを計算
        std::array<int, 3> from_jimage_array;
        for (size_t i = 0; i < from_jimage_array.size(); ++i) {
            from_jimage_array[i] = -to_jimage_array[i];
        }
        
        sg.add_edge(std::get<1>(e), to_jimage_array, std::get<0>(e), {0, 0, 0}, std::get<2>(e));
    }
    sg.set_cc_diameter();
    
    return sg;
}




void StructureGraph::add_edge(
        int from,
        std::array<int, 3> from_image,
        int to,
        std::array<int, 3> to_image,
        double weight
) {
    assert(0 <= from && from < int(this->graph.size()));
    assert(0 <= to && to < int(this->graph.size()));
    std::array<int, 3> image{};
    for (int i = 0; i < 3; ++i) image[i] = to_image[i] - from_image[i];

    // 自分自身への辺は無視する
    if (from == to && image == std::array<int, 3>{0, 0, 0}) {
        return;
    }

    // すでに追加されている辺は無視する
    if (this->graph_map.find(std::make_tuple(from, to, image)) != this->graph_map.end()) {
        return;
    }

    this->graph[from].emplace_back(NearNeighborInfo{to, weight, image, std::nullopt});
    this->graph_map[std::make_tuple(from, to, image)] = int(this->graph[from].size() - 1);
}

void StructureGraph::break_edge(
        int from,
        int to,
        std::array<int, 3> image,
        bool allow_reverse
) {
    assert(0 <= from && from < int(this->graph.size()));
    assert(0 <= to && to < int(this->graph.size()));
    // std::array<int, 3> image{};

    // 自分自身への辺は無視する
    if (from == to && image == std::array<int, 3>{0, 0, 0}) {
        return;
    }
    // 辺が存在する場合はそのまま取り除く
    if (auto iter = graph_map.find(std::make_tuple(from, to, image)); iter != graph_map.end()) {
        // graph[from].erase(graph[from].begin() + iter);
        auto begin_it = graph[from].begin();
        std::advance(begin_it, iter->second);
        graph[from].erase(begin_it);
        graph_map.erase(std::make_tuple(from, to, image));
        // graph_mapでgraphの要素を削除したのでiterより値が大きいvalueを1減らす
        for (const auto& [key, value] : graph_map){
            if (value > iter->second && std::get<0>(key) == from){
                graph_map[key] = value - 1;
            }
        }
    }else{
        if(allow_reverse){
            // 逆向きの辺を削除する
            // 逆向きのimageを定義
            std::array<int, 3> jimage;
            for (int i = 0; i < 3; i++) jimage[i] = -image[i];
            // if(auto iter = graph_map.find(std::make_tuple(to, from, jimage)); iter != this->graph_map.end()){
            if(auto iter = graph_map.find(std::make_tuple(to, from, jimage)); iter != graph_map.end()){
                auto begin_it = graph[to].begin();
                std::advance(begin_it, iter->second);
                graph[to].erase(begin_it);
                graph_map.erase(std::make_tuple(to, from, jimage));
                // graph_mapでgraphの要素を削除したのでiterより値が大きいvalueを1減らす
                for (const auto& [key, value] : graph_map){
                    if (value > iter->second && std::get<0>(key) == to){
                        graph_map[key] = value - 1;
                    }
                }
            }else{
                throw std::logic_error(
                    "Edge cannot be broken between from_index and to_index; "
                    "no edge exists between those sites.");
            }
        }
    }
}

// グラフの連結成分とその直径を計算する
void StructureGraph::set_cc_diameter() {
    const int n = int(graph.size());
    assert(n > 0);

    this->cc_nodes.clear();
    this->cc_diameter.clear();
    this->cc_cs.clear();

    std::vector<bool> visited(n, false);
    std::vector<int> queue;
    queue.reserve(n);
    size_t qi = 0;

    // 幅優先探索で連結成分を調べる
    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;
        visited[i] = true;
        queue.push_back(i);
        this->cc_nodes.emplace_back();
        this->cc_nodes.back().push_back(i);
        while (qi < queue.size()) {
            int u = queue[qi++]; // pop_front
            for (const auto &nni: graph[u]) {
                if (!visited[nni.site_index]) {
                    visited[nni.site_index] = true;
                    this->cc_nodes.back().push_back(nni.site_index);
                    queue.push_back(nni.site_index);
                }
            }
        }
    }

    for (auto &vec: this->cc_nodes) std::sort(vec.begin(), vec.end());

    Eigen::VectorXi d(n);
    for (const auto &nodes: this->cc_nodes) {
        // 連結成分ごとに幅優先探索を行う
        int d_max = 0;
        qi = 0;
        queue.resize(0);
        for (const int start: nodes) {
            queue.push_back(start);
            for (const int node: nodes) visited[node] = false;
            d[start] = 0;
            visited[start] = true;
            while (qi < queue.size()) {
                const int v = queue[qi++];
                for (const auto &nni: graph[v]) {
                    const int u = nni.site_index;
                    if (!visited[u]) {
                        visited[u] = true;
                        d[u] = d[v] + 1;
                        d_max = std::max(d_max, d[u]);
                        queue.push_back(u);
                    }
                }
            }
        }
        this->cc_diameter.push_back(d_max);
    }
}

void StructureGraph::set_elemental_labels() {
    this->labels = this->structure->species_strings;
}

void StructureGraph::set_loops(int depth_factor, int additional_depth) {
    // Import the required Python module
    py::module GraphAnalysisModule = py::module::import("graph_id.analysis.graphs");

    // Convert the C++ StructureGraph to its Python equivalent
    py::object py_structure_graph = this->to_py();

    // Convert py_structure_graph from PmgStructureGraph to StructureGraph
    py::object DesiredStructureGraphClass = GraphAnalysisModule.attr("StructureGraph");
    py::object desired_structure_graph = DesiredStructureGraphClass.attr("from_pymatgen_structure_graph")(
            py_structure_graph);

    // Call the set_loops method on this Python object
    desired_structure_graph.attr("set_loops")(depth_factor, additional_depth);

    // Retrieve the labels
    py::list labels = desired_structure_graph.attr("starting_labels").cast<py::list>();

    // Assign the labels to the C++ object
    for (size_t site_i = 0; site_i < labels.size(); site_i++) {
        this->labels[site_i] = labels[site_i].cast<std::string>();
    }
}


void StructureGraph::set_wyckoffs_label(double symmetry_tol) {
    auto core = py::module_::import("pymatgen.core");
    auto symmetry = py::module_::import("pymatgen.symmetry.analyzer");
    auto Element = core.attr("Element");
    auto SpacegroupAnalyzer = symmetry.attr("SpacegroupAnalyzer");

    py::object siteless = this->structure->py_structure.copy().obj;
    for (int i = 0; i < this->structure->count; i++) {
        siteless.attr("replace")(i, Element("H"));
    }

    auto sga = SpacegroupAnalyzer(siteless);
    auto sym_dataset = sga.attr("get_symmetry_dataset")();
    if (sym_dataset.is_none()) {
        this->set_elemental_labels();
        return;
    }

    auto wyckoffs = sym_dataset["wyckoffs"];
    auto number = sym_dataset["number"];

    for (size_t site_i = 0; site_i < py::len(wyckoffs); site_i++) {
        auto wyckoff = wyckoffs[py::int_(site_i)];
        this->labels[site_i] = py::str("{}_{}_{}").format(this->structure->species_strings[site_i], wyckoff, number);
    }
}

uint64_t connected_site_to_uint64(int site, std::array<int, 3> arr) {
    unsigned long long int ret = 0;
    assert(0 <= site && site < (1 << 16));
    assert(-32768 <= arr[0] && arr[0] < 32768);
    assert(-32768 <= arr[1] && arr[1] < 32768);
    assert(-32768 <= arr[2] && arr[2] < 32768);
    ret |= uint16_t(site);
    ret <<= 16;
    ret |= uint16_t(arr[0] + 32768);
    ret <<= 16;
    ret |= uint16_t(arr[1] + 32768);
    ret <<= 16;
    ret |= uint16_t(arr[2] + 32768);
    return ret;
}

void StructureGraph::set_compositional_sequence_node_attr(
        bool hash_cs,
        bool wyckoff,
        int additional_depth,
        int depth_factor,
        bool use_previous_cs
) {
    cc_cs.resize(0);

    for (size_t cc_i = 0; cc_i < cc_nodes.size(); cc_i++) {
        std::vector<std::string> cs_list;
        cs_list.reserve(cc_nodes[cc_i].size());

        const int depth = cc_diameter[cc_i] * depth_factor + additional_depth;

        for (const int focused_site_i: cc_nodes[cc_i]) {
            if (PyErr_CheckSignals()) throw py::error_already_set();
            CompositionalSequence cs;
            cs.hash_cs = hash_cs;
            cs.focused_site_i = focused_site_i;
            cs.labels = &labels;
            cs.use_previous_sites = use_previous_cs || wyckoff;
            cs.new_sites = {{focused_site_i, {0, 0, 0}}};
            cs.seen_sites.insert(connected_site_to_uint64(focused_site_i, {0, 0, 0}));

            for (int di = 0; di < depth; ++di) {
                for (const auto &c_site: cs.get_current_starting_sites()) {
                    for (const auto &nni: graph[std::get<0>(c_site)]) {
                        auto &arr = std::get<1>(c_site);
                        cs.count_composition_for_neighbors(nni.site_index, {
                                nni.image[0] + arr[0],
                                nni.image[1] + arr[1],
                                nni.image[2] + arr[2]
                        });
                    }
                }
                cs.finalize_this_depth();
            }
            cs_list.emplace_back(cs.string());
        }
        cc_cs.emplace_back(std::move(cs_list));
    }
}

void StructureGraph::set_individual_compositional_sequence_node_attr(
        int n,
        bool hash_cs,
        bool wyckoff,
        int additional_depth,
        int depth_factor,
        bool use_previous_cs
) {
    cc_cs.resize(0);

    py::object distance_measures = py::module::import("networkx.algorithms.distance_measures");
    py::object diameter = distance_measures.attr("diameter");
    py::object networkx = py::module::import("networkx");
    py::object py_structure_graph = this->to_py();
    py::object ug = py_structure_graph.attr("graph").attr("to_undirected")();
    
    // for (size_t cc_i = 0; cc_i < cc_nodes.size(); cc_i++) {
    py::object cc = networkx.attr("connected_components")(ug);
    std::vector<std::set<int>> cpp_vec;

    for (const auto& item : cc) {
        std::set<int> inner_vec = py::cast<std::set<int>>(item);
        cpp_vec.push_back(inner_vec);
    }
    for (const auto &cc_vector: cpp_vec) {
        std::vector<std::string> cs_list;
        // cs_list.reserve(cc_nodes[cc_i].size());
        cs_list.reserve(cc_vector.size());

        // int d = cc_diameter[cc_i];
        // int d = diameter(ug.attr("subgraph")(cc_nodes[cc_i])).cast<int>();
        int d = diameter(ug.attr("subgraph")(cc_vector)).cast<int>();
        const int depth = d * depth_factor + additional_depth;
        // if (std::count(cc_nodes[cc_i].begin(), cc_nodes[cc_i].end(), n)) {
        if (std::count(cc_vector.begin(), cc_vector.end(), n)) {
            // for (const int focused_site_i: cc_nodes[cc_i]) {
            if (PyErr_CheckSignals()) throw py::error_already_set();
            CompositionalSequence cs;
            cs.hash_cs = hash_cs;
            cs.focused_site_i = n;
            cs.labels = &labels;
            cs.use_previous_sites = use_previous_cs || wyckoff;
            // cs.new_sites = {{focused_site_i, {0, 0, 0}}};
            cs.new_sites = {{n, {0, 0, 0}}};
            // cs.seen_sites.insert(connected_site_to_uint64(focused_site_i, {0, 0, 0}));
            cs.seen_sites.insert(connected_site_to_uint64(n, {0, 0, 0}));

            for (int di = 0; di < depth; ++di) {
                for (const auto &c_site: cs.get_current_starting_sites()) {
                    for (const auto &nni: graph[std::get<0>(c_site)]) {
                        auto &arr = std::get<1>(c_site);
                        cs.count_composition_for_neighbors(nni.site_index, {
                                nni.image[0] + arr[0],
                                nni.image[1] + arr[1],
                                nni.image[2] + arr[2]
                        });
                    }
                }
                cs.finalize_this_depth();
            }
            // cs_list.emplace_back(blake2b(cs.string(), 16));
            cs_list.emplace_back(cs.string());
        // }
            cc_cs.emplace_back(std::move(cs_list));
        }
    }
}

int StructureGraph::calculate_rank(const gtl::flat_hash_set<std::array<int, 3>> &vertices) {
    size_t n = vertices.size();
    if (n == 0) return -1;
    if (n == 1) return 0;
    
    // 最初の頂点を基準点として、他の頂点との相対位置を計算
    std::vector<std::array<int, 3>> vertices_vec(vertices.begin(), vertices.end());
    const auto& base_vertex = vertices_vec[0];
    Eigen::Matrix3Xd relative_positions(3, n - 1);
    
    for (size_t i = 1; i < vertices_vec.size(); ++i) {
        relative_positions.col(i - 1) << 
            vertices_vec[i][0] - base_vertex[0],
            vertices_vec[i][1] - base_vertex[1], 
            vertices_vec[i][2] - base_vertex[2];
    }
    
    return int(relative_positions.fullPivLu().rank());
}

bool
StructureGraph::rank_increase(const gtl::flat_hash_set<std::array<int, 3>> &seen, const std::array<int, 3> &candidate) {
    size_t n = seen.size();
    if (n == 0) return true;

    // 既存の頂点セットに新しい頂点を追加
    gtl::flat_hash_set<std::array<int, 3>> extended_vertices = seen;
    extended_vertices.insert(candidate);
    
    // 拡張された頂点セットのランクを計算
    int rank1 = calculate_rank(extended_vertices);
    int rank0 = n - 1;  // 既存の頂点セットのランク
    
    return rank1 > rank0;
}

/// Larsen らの方法で次元を計算する
/// 連結成分ごとに適当に一つ頂点を選び pymatgen.analysis.calculate_dimensionality_of_site と
/// 同じ方法で次元を計算する。
/// 連結成分の他の頂点もその頂点と同じ次元になるので、連結成分ごとに１頂点しか計算しない。
/// 連結成分の頂点のうち、一番次元が高いものを選ぶ。
int StructureGraph::get_dimensionality_larsen() const {
    int max_dim = 0;
    std::vector<gtl::flat_hash_set<std::array<int, 3>>> seen_comp_vertices(structure->count);
    for (const auto &nodes: this->cc_nodes) {
        assert(!nodes.empty());
        const int node = nodes[0];
        gtl::flat_hash_set<uint64_t> seen_vertices;
        std::deque<std::tuple<int, std::array<int, 3>>> queue;
        queue.emplace_back(node, std::array<int, 3>{0, 0, 0});
        while (!queue.empty()) {
            const auto t = queue.front();
            const auto [comp_i, image_i] = t;
            const auto ni = connected_site_to_uint64(comp_i, image_i);
            queue.pop_front();

            if (seen_vertices.find(ni) != seen_vertices.end()) continue;
            seen_vertices.insert(ni);

            if (!rank_increase(seen_comp_vertices[comp_i], image_i)) continue;

            seen_comp_vertices[comp_i].insert(image_i);

            for (const auto &nni: graph[comp_i]) {
                int comp_j = nni.site_index;
                std::array<int, 3> image_j{};
                for (int i = 0; i < 3; ++i) image_j[i] = nni.image[i] + image_i[i];
                if (seen_vertices.find(connected_site_to_uint64(comp_j, image_j)) != seen_vertices.end()) {
                    continue;
                }
                if (!rank_increase(seen_comp_vertices[comp_j], image_j)) continue;
                queue.emplace_back(comp_j, image_j);
            }
        }
        // Pythonの実装に合わせて、rank関数を使用して次元を計算
        int dim = calculate_rank(seen_comp_vertices[node]);

        max_dim = std::max(max_dim, dim);
        if (max_dim == 3) return 3;
    }
    return max_dim;
}

py::object StructureGraph::to_py() const {
    py::object PmgStructureGraph = py::module::import("pymatgen.analysis.graphs").attr("StructureGraph");
    py::object sg = PmgStructureGraph.attr("from_empty_graph")(this->structure->py_structure.obj);
    for (int i = 0; size_t(i) < this->graph.size(); i++) {
        for (const auto &nni: this->graph[i]) {
            if (i <= nni.site_index) {
                sg.attr("add_edge")(
                        py::arg("from_index") = i,
                        py::arg("from_jimage") = py::make_tuple(0, 0, 0),
                        py::arg("to_index") = nni.site_index,
                        py::arg("to_jimage") = nni.image,
                        py::arg("weight") = nni.weight,
                        py::arg("warn_duplicates") = false
                );
            }
        }
    }
    return sg;
}


std::string CompositionalSequence::string() const {
    if (hash_cs) {
        return (*labels)[focused_site_i] + "-" + cs_for_hashing;
    } else {
        bool empty_bool = true;
        for (const auto cs_str: compositional_seq) {
            if (cs_str != "") {
                empty_bool = false;
            }
        }
        if (empty_bool) {
            return (*labels)[focused_site_i] + "-";
        }
        else {
            return (*labels)[focused_site_i] + "-" + join_string("-", compositional_seq);
        }
        
    }
}


std::vector<std::tuple<int, std::array<int, 3>>> CompositionalSequence::get_current_starting_sites() {
    const auto ret = std::move(new_sites);
    new_sites = {};
    return ret;
}


void CompositionalSequence::count_composition_for_neighbors(int site_i, std::array<int, 3> image) {
    const std::tuple<int, std::array<int, 3>> t = std::make_tuple(site_i, image);
    auto n = connected_site_to_uint64(site_i, image);
    if (seen_sites.find(n) == seen_sites.end()) {
        seen_sites.insert(n);
        new_sites.emplace_back(t);
        this->composition_counter[(*labels)[site_i]] += 1;
    }
}


void CompositionalSequence::finalize_this_depth() {
    auto formula = get_sorted_composition_list_form();
    if (hash_cs) {
        cs_for_hashing = blake2b(cs_for_hashing + "-" + join_string("", formula));
    } else {
        compositional_seq.emplace_back(join_string("", formula));
    }
}


std::vector<std::string> CompositionalSequence::get_sorted_composition_list_form() const {
    std::vector<std::string> ret;
    ret.reserve(composition_counter.size());
    for (const auto &t: composition_counter) {
        if (t.second > 0) {
            ret.emplace_back(t.first + std::to_string(t.second));
        }
    }
    return ret;
}


void init_structure_graph(pybind11::module &m) {
    py::class_<StructureGraph>(m, "StructureGraph")
            .def_static("with_local_env_strategy", [](PymatgenStructure &s, NearNeighbor &nn) {
                return StructureGraph::with_local_env_strategy(std::make_shared<Structure>(s), nn);
            })
            .def_static("from_empty_graph", [](PymatgenStructure &s) {
                return StructureGraph::from_empty_graph(std::make_shared<Structure>(s));
            })
            .def("set_elemental_labels", &StructureGraph::set_elemental_labels)
            .def("set_wyckoffs", &StructureGraph::set_wyckoffs_label, py::arg("symmetry_tol") = 0.1) // 互換性
            .def("set_wyckoffs_label", &StructureGraph::set_wyckoffs_label)
            .def("set_compositional_sequence_node_attr",
                 &StructureGraph::set_compositional_sequence_node_attr,
                 py::arg("hash_cs") = false,
                 py::arg("wyckoff") = false,
                 py::arg("additional_depth") = 0,
                 py::arg("depth_factor") = 2,
                 py::arg("use_previous_cs") = false)
            .def("get_dimensionality_larsen", &StructureGraph::get_dimensionality_larsen)
            .def("to_py", &StructureGraph::to_py)
            .def_static("from_py", &StructureGraph::from_py)
            .def("get_connected_site_index", [](const StructureGraph &sg) {
                // テスト用
                py::list arr;
                for (size_t i = 0; i < sg.graph.size(); i++) {
                    for (const auto &nni: sg.graph[i]) {
                        arr.append(py::make_tuple(i, nni.site_index));
                    }
                }
                arr.attr("sort")();
                return arr;
            })
            .def_property("labels", [](const StructureGraph &sg) { return sg.labels; },
                          [](StructureGraph &sg, const std::vector<std::string> &labels) { sg.labels = labels; })
            .def_property_readonly("cc_nodes", [](const StructureGraph &sg) { return sg.cc_nodes; })
            .def_property_readonly("cc_diameter", [](const StructureGraph &sg) { return sg.cc_diameter; })
            .def_property_readonly("cc_cs_labels", [](const StructureGraph &sg) { return sg.cc_cs; })
            .def_property_readonly("cc_cs", [](const StructureGraph &sg) {
                // Python との互換性のため、cc_nodes, cc_cs_labels を使うと効率的
                py::list res;
                for (size_t i = 0; i < sg.cc_nodes.size(); i++) {
                    py::dict d;
                    d["site_i"] = py::set(py::list(py::cast(sg.cc_nodes[i])));
                    d["cs_list"] = sg.cc_cs[i];
                    res.append(d);
                }
                return res;
            });
}
