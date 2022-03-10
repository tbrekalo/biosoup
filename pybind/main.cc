#include <limits>

#include "biosoup/nucleic_acid.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace py = pybind11;

PYBIND11_MODULE(biosouppy, m) {
  m.doc() = "biosoup python bindings";

  py::class_<biosoup::NucleicAcid>(m, "NucleicAcid")
      .def(py::init<const std::string&, const std::string&>())
      .def(py::init<const std::string&, const std::string&,
                    const std::string&>())
      .def("code", &biosoup::NucleicAcid::Code)
      .def("score", &biosoup::NucleicAcid::Score)
      .def("inflate_data", &biosoup::NucleicAcid::InflateData,
           py::arg("i") = 0U,
           py::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("inflate_quality", &biosoup::NucleicAcid::InflateQuality,
           py::arg("i") = 0U,
           py::arg("len") = std::numeric_limits<std::uint32_t>::max())
      .def("reverse_and_complement",
           &biosoup::NucleicAcid::ReverseAndComplement)
      .def_readwrite("id", &biosoup::NucleicAcid::id)
      .def_readwrite("name", &biosoup::NucleicAcid::name)
      .def_readonly("__len__", &biosoup::NucleicAcid::inflated_len);

  m.def("set_nucleic_acid_obj_cnt", [](const std::uint32_t val) -> void {
    biosoup::NucleicAcid::num_objects.store(val);
  });
}
