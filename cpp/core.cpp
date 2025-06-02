#include "core.h"

std::string join_string(const std::string &delimiter, const std::vector<std::string> &strings) {
    if (strings.empty()) return "";
    if (strings.size() == 1) return strings[0];
    size_t size = delimiter.size() * (strings.size() - 1);
    for (const auto &s: strings) {
        size += s.size();
    }
    std::string res;
    res.reserve(size);
    for (int i = 0; i < int(strings.size()) - 1; ++i) {
        res += strings[i];
        res += delimiter;
    }
    res += strings[strings.size() - 1];

    return res;
}

void init_core(py::module_ &m) {
    // Define implicit type conversions so that Python calls like f(Structure.from_file("..."))
    // can be made to void f(PymatgenStructure& s).

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

