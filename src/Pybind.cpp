#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ObjectPool.h"
#include "ForwardRayTracing.h"
#include "Utils.h"

namespace py = pybind11;

template<typename Real, typename Complex>
void define_sweep_result(pybind11::module_ &mod, const char *name) {
    using SweepR = SweepResult<Real, Complex>;
    py::class_<SweepR>(mod, name)
            .def_readonly("theta", &SweepR::theta)
            .def_readonly("phi", &SweepR::phi)
            .def_readonly("lda", &SweepR::lambda)
            .def_readonly("eta", &SweepR::eta)
            .def_readonly("delta_theta", &SweepR::delta_theta)
            .def_readonly("delta_phi", &SweepR::delta_phi)
            .def_readonly("theta_roots", &SweepR::theta_roots)
            .def_readonly("phi_roots", &SweepR::phi_roots)
            .def_readonly("theta_roots_closest", &SweepR::theta_roots_closest)
            .def_readonly("results", &SweepR::results);
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
            .def_readwrite("lda", &Params::lambda)
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
            .def_readonly("eta", &ResultType::eta)
            .def_readonly("lda", &ResultType::lambda) // lambda
            .def_readonly("rc", &ResultType::rc)
            .def_readonly("lgd", &ResultType::lgd)
            .def_readonly("lgd_sign", &ResultType::lgd_sign)
            .def_readonly("ray_status", &ResultType::ray_status);
}

template<typename Real, typename Complex>
void define_find_root_result(pybind11::module_ &mod, const char *name) {
    using ResultType = FindRootResult<Real, Complex>;
    py::class_<ResultType>(mod, name)
            .def_readonly("success", &ResultType::success)
            .def_readonly("fail_reason", &ResultType::fail_reason)
            .def_readonly("root", &ResultType::root);
}

template<typename Real, typename Complex>
void define_methods(pybind11::module_& mod, const std::string &suffix) {
    mod.def(("calc_ray" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::calc_ray,
        py::call_guard<py::gil_scoped_release>(), py::return_value_policy::move);
    mod.def(("sweep_rc_d" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::sweep_rc_d,
        py::return_value_policy::move);
    mod.def(("find_root_period" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::find_root_period,
        py::call_guard<py::gil_scoped_release>(), py::return_value_policy::move);
    mod.def(("find_root" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::find_root,
        py::call_guard<py::gil_scoped_release>(), py::return_value_policy::move);
    mod.def(("clean_cache" + suffix).c_str(), ForwardRayTracing<Real, Complex>::clear_cache);
}

template <typename Real, typename Complex>
void define_all(pybind11::module_& mod, const std::string& suffix) {
    define_methods<Real, Complex>(mod, "_" + suffix);
    define_params<Real>(mod, ("ForwardRayTracingParams" + suffix).c_str());
    define_forward_ray_tracing_result<Real, Complex>(mod, ("ForwardRayTracing" + suffix).c_str());
    define_find_root_result<Real, Complex>(mod, ("FindRootResult" + suffix).c_str());
    define_sweep_result<Real, Complex>(mod, ("SweepResult" + suffix).c_str());
}

PYBIND11_MODULE(py_forward_ray_tracing, mod) {
    py::enum_<RayStatus>(mod, "RayStatus")
            .value("NORMAL", RayStatus::NORMAL)
            .value("CONFINED", RayStatus::CONFINED)
            .value("ETA_OUT_OF_RANGE", RayStatus::ETA_OUT_OF_RANGE)
            .value("THETA_OUT_OF_RANGE", RayStatus::THETA_OUT_OF_RANGE)
            .value("UNKOWN_ERROR", RayStatus::UNKOWN_ERROR)
            .value("ARGUMENT_ERROR", RayStatus::ARGUMENT_ERROR)
            .value("INTERNAL_ERROR", RayStatus::INTERNAL_ERROR)
            .export_values();

    py::enum_<Sign>(mod, "Sign")
            .value("POSITIVE", Sign::POSITIVE)
            .value("NEGATIVE", Sign::NEGATIVE)
            .export_values();

    define_all<double, std::complex<double>>(mod, "Float64");
    define_all<long double, std::complex<long double>>(mod, "LongDouble");

#ifdef FLOAT128
    // define_forward_ray_tracing_result<boost::multiprecision::float128, boost::multiprecision::complex128>(mod, "ForwardRayTracingFloat128");
    // define_sweep_result<boost::multiprecision::float128>(mod, "SweepResultFloat128");
#endif
    // define_forward_ray_tracing_result<BigFloat, BigComplex>(mod, "ForwardRayTracingBigFloat");
    // define_sweep_result<BigFloat>(mod, "SweepResultBigFloat");
}
