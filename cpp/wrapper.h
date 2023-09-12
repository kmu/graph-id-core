#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <utility>

namespace py = pybind11;

// Typed wrapper of python list.
// py::object::cast<T>() is used to convert python object to C++ object.
template <typename T> class List {
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

inline void init_structure(py::module_& m) {
    // void f(PymatgenStructure& s) にたいして Python から f(Structure.from_file("...")) などと呼び出せるように
    // 暗黙の型変換を定義する。

    // PymatgenSite
    py::class_<PymatgenSite>(m, "PymatgenSite")
            .def(py::init<py::object>());
    py::implicitly_convertible<py::object, PymatgenSite>();

    // PymatgenPeriodicSite
    py::class_<PymatgenPeriodicSite, PymatgenSite>(m, "PymatgenPeriodicSite")
            .def(py::init<py::object>());
    py::implicitly_convertible<py::object, PymatgenPeriodicSite>();

    // PymatgenLattice
    py::class_<PymatgenLattice>(m, "PymatgenLattice")
            .def(py::init<py::object>());
    py::implicitly_convertible<py::object, PymatgenLattice>();

    // PymatgenStructure
    py::class_<PymatgenStructure>(m, "PymatgenStructure")
            .def(py::init<py::object>());
    py::implicitly_convertible<py::object, PymatgenStructure>();
}

