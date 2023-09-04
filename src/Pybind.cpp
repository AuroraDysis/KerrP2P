#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ForwardRayTracing.h"

namespace py = pybind11;

template <typename Real, typename Complex>
void define_forward_ray_tracing(pybind11::module_ &mod, const char *name) {
  using RayTracing = ForwardRayTracing<Real, Complex>;
  py::class_<RayTracing>(mod, name)
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
      .def_readonly("n_half", &RayTracing::n_half);
}

PYBIND11_MODULE(py_forward_ray_tracing, mod) {
  define_forward_ray_tracing<double, std::complex<double>>(mod, "ForwardRayTracingFloat64");
  define_forward_ray_tracing<long double, std::complex<long double>>(mod, "ForwardRayTracingLongDouble");
#ifdef FLOAT128
  define_forward_ray_tracing<boost::multiprecision::float128, boost::multiprecision::complex128>(mod,
                                                                                                 "ForwardRayTracingFloat128");
#endif
#ifdef BIGFLOAT
  py::class_<ForwardRayTracing<boost::multiprecision::mpfr_float_50, boost::multiprecision::mpc_complex_50>>(mod,
                                                                                                             "ForwardRayTracingBigFloat50");
#endif
}
