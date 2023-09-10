// #include <pybind11/numpy.h>
// #include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "ObjectPool.h"
#include "ForwardRayTracing.h"
#include "Utils.h"

namespace py = pybind11;

template<typename Real, typename Complex>
void define_sweep_result(pybind11::module_& mod, const char* name) {
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
void define_params(pybind11::module_& mod, const char* name) {
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
void define_forward_ray_tracing_result(pybind11::module_& mod, const char* name) {
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
void define_find_root_result(pybind11::module_& mod, const char* name) {
	using ResultType = FindRootResult<Real, Complex>;
	py::class_<ResultType>(mod, name)
		.def_readonly("success", &ResultType::success)
		.def_readonly("fail_reason", &ResultType::fail_reason)
		.def_readonly("root", &ResultType::root);
}

template<typename Real, typename Complex>
void define_methods(pybind11::module_& mod, const std::string& suffix) {
	mod.def(("calc_ray" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::calc_ray,
		py::call_guard<py::gil_scoped_release>(), py::return_value_policy::move);
	mod.def(("calc_ray_batch" + suffix).c_str(), &ForwardRayTracingUtils<Real, Complex>::calc_ray_batch,
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

template <typename T>
void define_numerical_type(pybind11::module_& mod, const char* name, bool is_complex) {
	py::class_<T> t(mod, name);
	t.def(py::init<>())
		.def(py::init<const std::string&>())
		.def(py::self + py::self)
		.def(py::self - py::self)
		.def(py::self * py::self)
		.def(py::self / py::self)
		.def(py::self += py::self)
		.def(py::self -= py::self)
		.def(py::self *= py::self)
		.def(py::self /= py::self)
		.def(py::self + float())
		.def(py::self - float())
		.def(py::self * float())
		.def(py::self / float())
		.def(py::self += float())
		.def(py::self -= float())
		.def(py::self *= float())
		.def(py::self /= float())
		.def(-py::self)
		.def(py::self == py::self)
		.def(py::self != py::self)
		.def("assign", [](T& self, const std::string& value) { self.assign(value); })
		.def("__str__", [](const T& self) { return self.str(std::numeric_limits<T>::max_digits10); })
		.def("__repr__", [](const T& self) { return self.str(std::numeric_limits<T>::max_digits10); });

	if (is_complex) {
		t.def_property("real", [](const T& self) { return self.real().str(std::numeric_limits<T>::max_digits10); },
			[](T& self, const std::string& value) {
				self.real().assign(value);
			})
			.def_property("imag", [](const T& self) { return self.imag().str(std::numeric_limits<T>::max_digits10); },
				[](T& self, const std::string& value) {
					self.imag().assign(value);
				});
	}
	else {
		t.def_property("value", [](const T& self) { return self.str(std::numeric_limits<T>::max_digits10); },
			[](T& self, const std::string& value) {
				self.assign(value);
			});
	}
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

	define_numerical_type<Float128>(mod, "Float128", false);
	define_numerical_type<Complex128>(mod, "Complex128", true);
	define_all<Float128, Complex128>(mod, "Float128");

	define_numerical_type<Float256>(mod, "Float256", false);
	define_numerical_type<Complex256>(mod, "Complex256", true);
	define_all<Float256, Complex256>(mod, "Float256");

	mod.attr("calc_ray") = mod.attr("calc_ray_Float64");
	mod.attr("calc_ray_batch") = mod.attr("calc_ray_batch_Float64");
	mod.attr("sweep_rc_d") = mod.attr("sweep_rc_d_Float64");
	mod.attr("find_root_period") = mod.attr("find_root_period_Float64");
	mod.attr("find_root") = mod.attr("find_root_Float64");
	mod.attr("clean_cache") = mod.attr("clean_cache_Float64");
	mod.attr("ForwardRayTracingParams") = mod.attr("ForwardRayTracingParamsFloat64");
	mod.attr("ForwardRayTracing") = mod.attr("ForwardRayTracingFloat64");
	mod.attr("FindRootResult") = mod.attr("FindRootResultFloat64");
	mod.attr("SweepResult") = mod.attr("SweepResultFloat64");
}
