#include <limits>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/overlap.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(biosouppy, m) {
  m.doc() = "biosoup python bindings";

  py::class_<biosoup::NucleicAcid>(m, "NucleicAcid")
      .def(py::init<const std::string&, const std::string&>())
      .def(py::init<const std::string&, const std::string&,
                    const std::string&>())
      .def(
          "__deepcopy__",
          [](const biosoup::NucleicAcid& seq,
             py::dict) -> biosoup::NucleicAcid { return seq; },
          "memo"_a)
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

  py::class_<biosoup::Overlap>(m, "Overlap")
      .def(py::init<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                    std::uint32_t, std::uint32_t, std::uint32_t, bool>(),
           py::arg("lhs_id"), py::arg("lhs_begin"), py::arg("lhs_end"),
           py::arg("rhs_id"), py::arg("rhs_begin"), py::arg("rhs_end"),
           py::arg("score"), py::arg("strand") = true)
      .def(py::init<std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                    std::uint32_t, std::uint32_t, std::uint32_t,
                    const std::string&, bool>(),
           py::arg("lhs_id"), py::arg("lhs_begin"), py::arg("lhs_end"),
           py::arg("rhs_id"), py::arg("rhs_begin"), py::arg("rhs_end"),
           py::arg("score"), py::arg("alignment"), py::arg("strand") = true)
      .def(
          "__deepcopy__",
          [](const biosoup::Overlap& ovlp, py::dict) -> biosoup::Overlap {
            return ovlp;
          },
          "memo"_a)
      .def_readwrite("lhs_id", &biosoup::Overlap::lhs_id)
      .def_readwrite("lhs_begin", &biosoup::Overlap::lhs_begin)
      .def_readwrite("lhs_end", &biosoup::Overlap::lhs_end)
      .def_readwrite("rhs_id", &biosoup::Overlap::rhs_id)
      .def_readwrite("rhs_begin", &biosoup::Overlap::rhs_begin)
      .def_readwrite("rhs_end", &biosoup::Overlap::rhs_end)
      .def_readwrite("score", &biosoup::Overlap::score)
      .def_readwrite("strand", &biosoup::Overlap::strand)
      .def_readwrite("alignment", &biosoup::Overlap::alignment);
}
