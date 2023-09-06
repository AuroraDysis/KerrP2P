#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <oneapi/tbb.h>

#include "ObjectPool.h"
#include "ForwardRayTracing.h"

namespace py = pybind11;

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

template<typename Real, typename Complex>
struct PyForwardRayTracing {
  static std::shared_ptr<ForwardRayTracing<Real, Complex>>
  ray_tracing_rc_d(Real a, Real r_s, Real theta_s, Real r_o, Sign nu_r, Sign nu_theta, Real rc, Real d) {
    auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
    ray_tracing->calc_ray_by_rc_d(a, r_s, theta_s, r_o, nu_r, nu_theta, rc, d);
    return ray_tracing;
  }

  static std::shared_ptr<ForwardRayTracing<Real, Complex>> ray_tracing_lambda_q(
      Real a, Real r_s, Real theta_s, Real r_o, Sign nu_r, Sign nu_theta, Real lambda, Real q) {
    auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
    ray_tracing->calc_ray_by_lambda_q(a, r_s, theta_s, r_o, nu_r, nu_theta, lambda, q);
    return ray_tracing;
  }

  static std::tuple<py::array_t<Real>, py::array_t<Real>, py::array_t<Real>, py::array_t<Real>, py::array_t<Real>>
  ray_tracing_rc_d_area(Real a, Real r_s, Real theta_s, Real r_o, Sign nu_r, Sign nu_theta, Real theta_o, Real phi_o,
                        const std::vector<Real> &rc_list, const std::vector<Real> &d_list) {
    size_t rc_size = rc_list.size();
    size_t d_size = d_list.size();

    // d_size rows, rc_size cols
    py::array_t<Real> theta_mat({d_size, rc_size});
    py::array_t<Real> phi_mat({d_size, rc_size});

    auto theta_data = theta_mat.template mutable_unchecked<2>();
    auto phi_data = phi_mat.template mutable_unchecked<2>();

    // rc and d
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range2d<size_t>(0u, d_size, 0u, rc_size),
                              [&](const oneapi::tbb::blocked_range2d<size_t, size_t> &r) {
                                auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
                                Real two_pi = boost::math::constants::two_pi<Real>();
                                for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                                  for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                                    ray_tracing->calc_ray_by_rc_d(a, r_s, theta_s, r_o, nu_r, nu_theta, rc_list[j],
                                                                  d_list[i]);
                                    if (ray_tracing->ray_status == RayStatus::NORMAL) {
                                      theta_data(i, j) = ray_tracing->theta_f - theta_o;
                                      phi_data(i, j) = ray_tracing->phi_f - phi_o;
                                    } else {
                                      theta_data(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                      phi_data(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                    }
                                  }
                                }
                              });

    std::vector<std::pair<py::ssize_t, py::ssize_t>> theta_roots_index;
    std::vector<std::pair<py::ssize_t, py::ssize_t>> phi_roots_index;
    Real d_row, d_col;
    for (py::ssize_t i = 1; i < d_size; i++) {
      for (py::ssize_t j = 1; j < rc_size; j++) {
        d_row = theta_data(i, j) * theta_data(i, j - 1);
        d_col = theta_data(i, j) * theta_data(i - 1, j);
        if (!isnan(d_row) && !isnan(d_col) && (d_row <= 0 || d_col <= 0)) {
          theta_roots_index.emplace_back(i, j);
        }
        d_row = phi_data(i, j) * phi_data(i, j - 1);
        d_col = phi_data(i, j) * phi_data(i - 1, j);
        if (!isnan(d_row) && !isnan(d_col) && (d_row <= 0 || d_col <= 0)) {
          phi_roots_index.emplace_back(i, j);
        }
      }
    }

    std::array<py::ssize_t, 2> shape{static_cast<py::ssize_t>(theta_roots_index.size()), 2u};
    py::array_t<Real> theta_roots(shape);
    py::array_t<Real> theta_roots_closest(shape);
    shape[0] = static_cast<py::ssize_t>(phi_roots_index.size());
    py::array_t<Real> phi_roots(shape);

    auto theta_roots_data = theta_roots.template mutable_unchecked<2>();
    auto phi_roots_data = phi_roots.template mutable_unchecked<2>();

    for (py::ssize_t i = 0; i < theta_roots_index.size(); i++) {
      theta_roots_data(i, 0) = rc_list[theta_roots_index[i].second];
      theta_roots_data(i, 1) = d_list[theta_roots_index[i].first];
    }
    for (py::ssize_t i = 0; i < phi_roots_index.size(); i++) {
      phi_roots_data(i, 0) = rc_list[phi_roots_index[i].second];
      phi_roots_data(i, 1) = d_list[phi_roots_index[i].first];
    }

    auto theta_roots_closest_data = theta_roots_closest.template mutable_unchecked<2>();
    // find the closest point in phi_roots_index for each point in theta_roots_index
    for (py::ssize_t i = 0; i < theta_roots_index.size(); i++) {
      double min_dist = std::numeric_limits<double>::max();
      for (py::ssize_t j = 0; j < phi_roots_index.size(); j++) {
        double dist = square(theta_roots_index[i].first - phi_roots_index[j].first) +
                      square(theta_roots_index[i].second - phi_roots_index[j].second);
        if (dist < min_dist) {
          min_dist = dist;
          theta_roots_closest_data(i, 0) = phi_roots_data(j, 0);
          theta_roots_closest_data(i, 1) = phi_roots_data(j, 1);
        }
      }
    }

    return {theta_mat, phi_mat, theta_roots, phi_roots, theta_roots_closest};
  }
};

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

  mod.def("ray_tracing_rc_d", &PyForwardRayTracing<double, std::complex<double>>::ray_tracing_rc_d,
          py::call_guard<py::gil_scoped_release>());
  mod.def("ray_tracing_lambda_q", &PyForwardRayTracing<double, std::complex<double>>::ray_tracing_lambda_q,
          py::call_guard<py::gil_scoped_release>());
  mod.def("ray_tracing_rc_d_area", &PyForwardRayTracing<double, std::complex<double>>::ray_tracing_rc_d_area,
          py::return_value_policy::move);
  mod.def("clean_cache", ForwardRayTracing<double, std::complex<double>>::clear_cache);

  define_forward_ray_tracing<double, std::complex<double>>(mod, "ForwardRayTracingFloat64");
  define_forward_ray_tracing<long double, std::complex<long double>>(mod, "ForwardRayTracingLongDouble");
#ifdef FLOAT128
  // define_forward_ray_tracing<boost::multiprecision::float128, boost::multiprecision::complex128>(mod, "ForwardRayTracingFloat128");
#endif
#ifdef BIGFLOAT
  // define_forward_ray_tracing<BigFloat, BigComplex>(mod, "ForwardRayTracingBigFloat");
#endif
}
