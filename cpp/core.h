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
    py::object obj;

    explicit PymatgenStructure(py::object obj) : obj(std::move(obj)) {};

    PymatgenLattice lattice() {
        return PymatgenLattice(obj.attr("lattice"));
    }

    List<PymatgenPeriodicSite> sites() {
        return List<PymatgenPeriodicSite>(obj.attr("sites"));
    }
};


struct Lattice {
    Lattice() = default;

    explicit Lattice(PymatgenLattice l) {
        auto m = l.matrix();
        matrix <<
               m.at(0, 0), m.at(1, 0), m.at(2, 0),
                m.at(0, 1), m.at(1, 1), m.at(2, 1),
                m.at(0, 2), m.at(1, 2), m.at(2, 2);
        inv_matrix = matrix.inverse();
        pbc = l.pbc();
    }

    // Eigen は Fortran 配列のように列優先だが、numpy は行優先なので Pymatgen の行列を転地して格納する。
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d inv_matrix;
    std::array<bool, 3> pbc{true, true, true};
};

// C++ でいちいち python のオブジェクトにアクセスするオーバーヘッドを避けるために使う PymatgenStructure のコピー。
struct Structure {
    Structure() = default;

    explicit Structure(PymatgenStructure s) {
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
};


void init_core(py::module_ &m);
