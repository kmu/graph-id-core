#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include <utility>

namespace py = pybind11;

const double pi = 3.1415926535897932384626433832795028841971;

// Typed wrapper of python list.
// py::object::cast<T>() is used to convert python object to C++ object.
template<typename T>
class List {
public:
    py::list obj;

    explicit List(py::list obj) : obj(std::move(obj)) {};

    class Iterator {
    public:
        List<T> *list;
        size_t i;

        Iterator(List<T> *list, size_t i) : list(list), i(i) {};

        Iterator operator++() {
            i++;
            return *this;
        }

        T operator*() {
            return (*list)[i];
        }

        bool operator!=(const Iterator &other) {
            return i != other.i;
        }

        bool operator==(const Iterator &other) {
            return i == other.i;
        }

        bool operator<(const Iterator &other) {
            return i < other.i;
        }
    };

    T operator[](int i) {
        return obj[i].template cast<T>();
    }

    Iterator begin() {
        return Iterator(this, 0);
    }

    Iterator end() {
        return Iterator(this, size());
    }

    size_t size() {
        return py::len(obj);
    }
};

// Wrapper of pymatgen.core.Site class
class PymatgenSite {
public:
    py::object obj;

    explicit PymatgenSite(py::object obj) : obj(std::move(obj)) {};

    double x() {
        return obj.attr("x").cast<double>();
    }

    double y() {
        return obj.attr("y").cast<double>();
    }

    double z() {
        return obj.attr("z").cast<double>();
    }

    std::string species_string() {
        return obj.attr("species_string").cast<std::string>();
    }
};

// Wrapper of pymatgen.core.PeriodicSite class
class PymatgenPeriodicSite : public PymatgenSite {
public:
    explicit PymatgenPeriodicSite(py::object obj) : PymatgenSite(std::move(obj)) {};

    double a() {
        return obj.attr("a").cast<double>();
    }

    double b() {
        return obj.attr("b").cast<double>();
    }

    double c() {
        return obj.attr("c").cast<double>();
    }
};

// Wrapper of pymatgen.core.Lattice class
class PymatgenLattice {
public:
    py::object obj;

    explicit PymatgenLattice(py::object obj) : obj(std::move(obj)) {};

    bool is_none() {
        return obj.is_none();
    }

    std::tuple<double, double, double> lengths() {
        return obj.attr("lengths").cast<std::tuple<double, double, double>>();
    }

    std::tuple<double, double, double> angles() {
        return obj.attr("angles").cast<std::tuple<double, double, double>>();
    }

    std::array<bool, 3> pbc() {
        return obj.attr("pbc").cast<std::array<bool, 3>>();
    }

    bool is_orthogonal() {
        return obj.attr("is_orthogonal").cast<bool>();
    }

    PymatgenLattice copy() {
        return PymatgenLattice(obj.attr("copy")());
    }


    py::array_t<double> matrix() {
        return obj.attr("matrix").cast<py::array_t<double>>();
    }
};

// Wrapper of pymatgen.core.IStructure class
class PymatgenStructure {
public:
    mutable py::object obj;

    PymatgenStructure() = default;

    explicit PymatgenStructure(py::object obj) : obj(std::move(obj)) {};

    PymatgenLattice lattice() const {
        if (py::hasattr(obj, "lattice")) {
            return PymatgenLattice(obj.attr("lattice"));
        } else {
            return PymatgenLattice(py::none());
        }
    }

    List<PymatgenPeriodicSite> sites() const {
        return List<PymatgenPeriodicSite>(obj.attr("sites"));
    }

    std::string formula() const {
        return obj.attr("composition").attr("formula").cast<std::string>();
    }

    std::string reduced_formula() const {
        return obj.attr("composition").attr("reduced_formula").cast<std::string>();
    }

    PymatgenStructure copy() const {
        return PymatgenStructure(obj.attr("copy")());
    }

    bool is_none() const {
        return obj.is_none();
    }
};


struct Lattice {
    Lattice() = default;

    explicit Lattice(PymatgenLattice l) {
        if (l.is_none()) {
            matrix << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;
            inv_matrix = matrix;
            pbc = {false, false, false};
        } else {
            auto m = l.matrix();
            matrix <<
                   m.at(0, 0), m.at(1, 0), m.at(2, 0),
                    m.at(0, 1), m.at(1, 1), m.at(2, 1),
                    m.at(0, 2), m.at(1, 2), m.at(2, 2);
            inv_matrix = matrix.inverse();
            pbc = l.pbc();
        }
    }

    // Eigen uses column-major ordering like Fortran arrays, but numpy uses row-major, so we store the Pymatgen matrix transposed.
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d inv_matrix;
    std::array<bool, 3> pbc{true, true, true};
};

// A copy of PymatgenStructure used to avoid the overhead of accessing Python objects repeatedly in C++.
// Can also handle Molecule transparently
struct Structure {
    Structure() = default;

    explicit Structure(PymatgenStructure s) : py_structure(s) {
        auto l = s.lattice();
        lattice = Lattice(l);
        count = int(s.sites().size());
        site_xyz.resize(3, count);
        species_strings.reserve(count);
        for (int i = 0; i < count; i++) {
            auto site = s.sites()[i];
            site_xyz.col(i) << site.x(), site.y(), site.z();
            species_strings.emplace_back(site.species_string());
        }
    }

    int count{0};
    Lattice lattice;
    Eigen::Matrix3Xd site_xyz;
    std::vector<std::string> species_strings;
    PymatgenStructure py_structure;
};

inline void warn(const std::string &s) {
    // Python で warnings.warn(s) するのと同じ。
    PyErr_WarnEx(PyExc_Warning, s.c_str(), 2);
}

inline void warn(py::object &s) {
    PyErr_WarnEx(PyExc_Warning, py::str(s).cast<std::string>().c_str(), 2);
}

std::string join_string(const std::string &delimiter, const std::vector<std::string> &strings);

template<class T>
size_t HashCombine(const size_t seed, const T &v) {
    return seed ^ (std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}

// Some primes between 2^63 and 2^64 for various uses.
// CityHash (http://code.google.com/p/cityhash/) の一部を参考にした。
static const uint64_t hash_seed0 = 0xc3a5c85c97cb3127ULL;
static const uint64_t hash_seed1 = 0xb492b66fbe98f273ULL;
static const uint64_t hash_seed2 = 0x9ae16a3b2f90404fULL;

// Hash function that combines two 64-bit integers x and y to create a 64-bit integer
// Based on parts of CityHash (http://code.google.com/p/cityhash/)
// Usage: hash_combine(hash_seed0, x) etc.
inline uint64_t hash_combine(uint64_t x, uint64_t y) {
    const uint64_t kMul = 0x9ddfea08eb382d69ULL;
    uint64_t a = (x ^ y) * kMul;
    a ^= (a >> 47);
    uint64_t b = (x ^ a) * kMul;
    b ^= (b >> 47);
    b *= kMul;
    return b;
}

template<>
struct std::hash<std::tuple<int, std::array<int, 3>>> {
    size_t operator()(const std::tuple<int, std::array<int, 3>> &t) const noexcept {
        return HashCombine(std::get<0>(t),
                           HashCombine(std::get<1>(t)[0], HashCombine(std::get<1>(t)[1], std::get<1>(t)[2])));
    }
};

template<>
struct std::hash<std::array<int, 3>> {
    size_t operator()(const std::array<int, 3> &t) const noexcept {
        return HashCombine(t[0], HashCombine(t[1], t[2]));
    }
};

template<>
struct std::hash<std::array<int, 4>> {
    size_t operator()(const std::array<int, 4> &t) const noexcept {
        return HashCombine(HashCombine(t[0], t[3]), HashCombine(t[1], t[2]));
    }
};

template<typename T, uint64_t seed = hash_seed0>
struct Hash {
    size_t operator()(const T &t) const noexcept {
        return hash_combine(seed, t);
    }
};

void init_core(py::module_ &m);
