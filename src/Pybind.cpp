#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ObjectPool.h"
#include "ForwardRayTracing.h"
#include "Utils.h"

namespace py = pybind11;

template<typename Real>
void define_sweep_result(pybind11::module_ &mod, const char *name) {
  using SweepR = SweepResult<Real>;
  py::class_<SweepR>(mod, name)
      .def_readonly("theta", &SweepR::theta)
      .def_readonly("phi", &SweepR::phi)
      .def_readonly("delta_theta", &SweepR::delta_theta)
      .def_readonly("delta_phi", &SweepR::delta_phi)
      .def_readonly("theta_roots", &SweepR::theta_roots)
      .def_readonly("phi_roots", &SweepR::phi_roots)
      .def_readonly("theta_roots_closest", &SweepR::theta_roots_closest);
}

template<typename Real, typename Complex>
void define_forward_ray_tracing(pybind11::module_ &mod, const char *name) {
  using RayTracing = ForwardRayTracing<Real, Complex>;
  py::class_<RayTracing, std::shared_ptr<RayTracing>>(mod, name)
      .def_readonly("a", &RayTracing::a)
      .def_readonly("rp", &RayTracing::rp)
      .def_readonly("rm", &RayTracing::rm)
      .def_readonly("r_s", &RayTracing::r_s)
      .def_readonly("theta_s", &RayTracing::theta_s)
      .def_readonly("r_o", &RayTracing::r_o)
      .def_readonly("r1", &RayTracing::r1)
      .def_readonly("r2", &RayTracing::r2)
      .def_readonly("r3", &RayTracing::r3)
      .def_readonly("r4", &RayTracing::r4)
      .def_readonly("r1_c", &RayTracing::r1_c)
      .def_readonly("r2_c", &RayTracing::r2_c)
      .def_readonly("r3_c", &RayTracing::r3_c)
      .def_readonly("r4_c", &RayTracing::r4_c)
      .def_readonly("t_f", &RayTracing::t_f)
      .def_readonly("theta_f", &RayTracing::theta_f)
      .def_readonly("phi_f", &RayTracing::phi_f)
      .def_readonly("m", &RayTracing::m)
      .def_readonly("n_half", &RayTracing::n_half)
      .def_readonly("ray_status", &RayTracing::ray_status);
}

PYBIND11_MODULE(py_forward_ray_tracing, mod) {
  py::enum_<RayStatus>(mod, "RayStatus")
      .value("NORMAL", RayStatus::NORMAL)
      .value("CONFINED", RayStatus::CONFINED)
      .value("ETA_OUT_OF_RANGE", RayStatus::ETA_OUT_OF_RANGE)
      .value("THETA_OUT_OF_RANGE", RayStatus::THETA_OUT_OF_RANGE)
      .value("UNKOWN_ERROR", RayStatus::UNKOWN_ERROR)
      .export_values();

  py::enum_<Sign>(mod, "Sign")
      .value("POSITIVE", Sign::POSITIVE)
      .value("NEGATIVE", Sign::NEGATIVE)
      .export_values();

  mod.def("calc_rc_d", &ForwardRayTracingUtils<double, std::complex<double>>::calc_rc_d,
          py::call_guard<py::gil_scoped_release>());
  mod.def("calc_lambda_q", &ForwardRayTracingUtils<double, std::complex<double>>::calc_lambda_q,
          py::call_guard<py::gil_scoped_release>());
  mod.def("sweep_rc_d", &ForwardRayTracingUtils<double, std::complex<double>>::sweep_rc_d,
          py::return_value_policy::move);
  mod.def("clean_cache", ForwardRayTracing<double, std::complex<double>>::clear_cache);

  define_forward_ray_tracing<double, std::complex<double>>(mod, "ForwardRayTracingFloat64");
  define_forward_ray_tracing<long double, std::complex<long double>>(mod, "ForwardRayTracingLongDouble");
  define_sweep_result<double>(mod, "SweepResultFloat64");
  define_sweep_result<long double>(mod, "SweepResultLongDouble");
#ifdef FLOAT128
  // define_forward_ray_tracing<boost::multiprecision::float128, boost::multiprecision::complex128>(mod, "ForwardRayTracingFloat128");
  // define_sweep_result<boost::multiprecision::float128>(mod, "SweepResultFloat128");
#endif
#ifdef BIGFLOAT
  // define_forward_ray_tracing<BigFloat, BigComplex>(mod, "ForwardRayTracingBigFloat");
  // define_sweep_result<BigFloat>(mod, "SweepResultBigFloat");
#endif
}
