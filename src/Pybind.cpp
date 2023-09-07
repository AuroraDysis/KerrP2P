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

template<typename Real>
void define_params(pybind11::module_ &mod, const char *name) {
  using Params = ForwardRayTracingParams<Real>;
  py::class_<Params>(mod, name)
      .def(py::init<>())
      .def_readwrite("a", &Params::a)
      .def_readwrite("r_s", &Params::r_s)
      .def_readwrite("theta_s", &Params::theta_s)
      .def_readwrite("r_o", &Params::r_o)
      .def_readwrite("nu_r", &Params::nu_r)
      .def_readwrite("nu_theta", &Params::nu_theta)
      .def_readwrite("rc", &Params::rc)
      .def_readwrite("lgd", &Params::lgd)
      .def_readwrite("lgd_sign", &Params::lgd_sign)
      .def_readwrite("lambda", &Params::lambda)
      .def_readwrite("q", &Params::q)
      .def_readwrite("calc_t_f", &Params::calc_t_f)
      .def("rc_d_to_lambda_q", &Params::rc_d_to_lambda_q);
}

template<typename Real, typename Complex>
void define_forward_ray_tracing_result(pybind11::module_ &mod, const char *name) {
  using ResultType = ForwardRayTracingResult<Real, Complex>;
  py::class_<ResultType>(mod, name)
      .def_readonly("a", &ResultType::a)
      .def_readonly("rp", &ResultType::rp)
      .def_readonly("rm", &ResultType::rm)
      .def_readonly("r_s", &ResultType::r_s)
      .def_readonly("theta_s", &ResultType::theta_s)
      .def_readonly("r_o", &ResultType::r_o)
      .def_readonly("r1", &ResultType::r1)
      .def_readonly("r2", &ResultType::r2)
      .def_readonly("r3", &ResultType::r3)
      .def_readonly("r4", &ResultType::r4)
      .def_readonly("r1_c", &ResultType::r1_c)
      .def_readonly("r2_c", &ResultType::r2_c)
      .def_readonly("r3_c", &ResultType::r3_c)
      .def_readonly("r4_c", &ResultType::r4_c)
      .def_readonly("t_f", &ResultType::t_f)
      .def_readonly("theta_f", &ResultType::theta_f)
      .def_readonly("phi_f", &ResultType::phi_f)
      .def_readonly("m", &ResultType::m)
      .def_readonly("n_half", &ResultType::n_half)
      .def_readonly("ray_status", &ResultType::ray_status);
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

  mod.def("calc_ray", &ForwardRayTracingUtils<double, std::complex<double>>::calc_ray,
          py::call_guard<py::gil_scoped_release>(), py::return_value_policy::move);
  mod.def("sweep_rc_d", &ForwardRayTracingUtils<double, std::complex<double>>::sweep_rc_d,
          py::return_value_policy::move);
  mod.def("clean_cache", ForwardRayTracing<double, std::complex<double>>::clear_cache);

  define_params<double>(mod, "ForwardRayTracingParamsFloat64");
  define_params<long double>(mod, "ForwardRayTracingParamsLongDouble");
  define_forward_ray_tracing_result<double, std::complex<double>>(mod, "ForwardRayTracingFloat64");
  define_forward_ray_tracing_result<long double, std::complex<long double>>(mod, "ForwardRayTracingLongDouble");
  define_sweep_result<double>(mod, "SweepResultFloat64");
  define_sweep_result<long double>(mod, "SweepResultLongDouble");
#ifdef FLOAT128
  // define_forward_ray_tracing_result<boost::multiprecision::float128, boost::multiprecision::complex128>(mod, "ForwardRayTracingFloat128");
  // define_sweep_result<boost::multiprecision::float128>(mod, "SweepResultFloat128");
#endif
#ifdef BIGFLOAT
  // define_forward_ray_tracing_result<BigFloat, BigComplex>(mod, "ForwardRayTracingBigFloat");
  // define_sweep_result<BigFloat>(mod, "SweepResultBigFloat");
#endif
}
