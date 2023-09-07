#pragma once

#include "ForwardRayTracing.h"
#include "NelderMead.h"

#include <oneapi/tbb.h>
#include <Eigen/Dense>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<typename Real, typename Complex>
struct SweepResult {
  using PointVector = Eigen::Matrix<Real, Eigen::Dynamic, 2>;
  using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

  Matrix theta;
  Matrix phi;

  Matrix delta_theta;
  Matrix delta_phi;

  PointVector theta_roots;
  PointVector phi_roots;

  PointVector theta_roots_closest;

  std::vector<ForwardRayTracingResult<Real, Complex>> results;
};

template<typename Real>
inline Real mod_real(const Real &x, const Real &y) {
  using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
      boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;
  return x - y * RealToInt::convert(x / y);
}

template<typename Real, typename Complex>
class RootFunctor {
private:
  using Vector = Eigen::Matrix<Real, 2, 1>;

  ForwardRayTracingParams<Real> &params;
  const Real theta_o;
  const Real phi_o;
  const int period;
  const Real two_pi = boost::math::constants::two_pi<Real>();

public:
  std::shared_ptr<ForwardRayTracing<Real, Complex>> ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();

  RootFunctor(ForwardRayTracingParams<Real> &params_, int period_, Real theta_o_, Real phi_o_) : params(params_),
                                                                                                 period(period_),
                                                                                                 theta_o(std::move(
                                                                                                     theta_o_)),
                                                                                                 phi_o(
                                                                                                     std::move(
                                                                                                         phi_o_)) {
    ray_tracing->calc_t_f = false;
  }

  Real operator()(Vector x) {
    auto &rc = x[0];
    auto &lgd = x[1];
    params.rc = rc;
    params.lgd = lgd;
    params.rc_d_to_lambda_q();
    ray_tracing->calc_ray(params);

    if (ray_tracing->ray_status != RayStatus::NORMAL) {
      fmt::println("ray status: {}", ray_status_to_str(ray_tracing->ray_status));
      throw std::runtime_error("Ray tracing failed: ray status is not normal");
    }

    Vector residual;
    residual[0] = ray_tracing->theta_f - theta_o;
    residual[1] = ray_tracing->phi_f - phi_o - period * two_pi;

    if (isnan(residual[0]) || isnan(residual[1])) {
      throw std::runtime_error("Ray tracing failed: residual is NaN");
    }

#ifdef PRINT_DEBUG
    fmt::println("rc: {}, lgd: {}, theta_f: {}, phi_f: {}", x[0], x[1], ray_tracing->theta_f, ray_tracing->phi_f);
    fmt::println("residual: {}, {}", residual[0], residual[1]);
#endif
    return residual.norm();
  }
};

template<typename Real, typename Complex>
struct ForwardRayTracingUtils {
  static ForwardRayTracingResult<Real, Complex> calc_ray(const ForwardRayTracingParams<Real> &params) {
    auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
    ray_tracing->calc_ray(params);
    return ray_tracing->to_result();
  }

  static ForwardRayTracingResult<Real, Complex>
  find_result(ForwardRayTracingParams<Real> &params, int period, Real theta_o, Real phi_o) {
    // std::array<Real, 2> x = {params.rc, params.lgd};
    Eigen::Vector<Real, 2> x = Eigen::Vector<Real, 2>();
    x << params.rc, params.lgd;

    RootFunctor<Real, Complex> root_functor(params, period, std::move(theta_o), std::move(phi_o));

    NelderMeadOptimizerParams<Real> settings;
    auto xout = NelderMeadOptimizer(root_functor, x, settings);

    root_functor(xout);

    return root_functor.ray_tracing->to_result();
  }

  static ForwardRayTracingResult<Real, Complex> refine_result(ForwardRayTracingResult<Real, Complex> &res) {

  }

  static SweepResult<Real, Complex>
  sweep_rc_d(ForwardRayTracingParams<Real> &params, Real theta_o, Real phi_o, const std::vector<Real> &rc_list,
             const std::vector<Real> &lgd_list, size_t cutoff) {
    size_t rc_size = rc_list.size();
    size_t lgd_size = lgd_list.size();

    SweepResult<Real, Complex> sweep_result;

    auto &theta = sweep_result.theta;
    auto &phi = sweep_result.phi;

    theta.resize(lgd_size, rc_size);
    phi.resize(lgd_size, rc_size);

    // rc and d
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range2d<size_t>(0u, lgd_size, 0u, rc_size),
                              [&](const oneapi::tbb::blocked_range2d<size_t, size_t> &r) {
                                auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
                                Real two_pi = boost::math::constants::two_pi<Real>();
                                ForwardRayTracingParams<Real> local_params(params);
                                for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                                  for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                                    local_params.rc = rc_list[j];
                                    local_params.lgd = lgd_list[i];
                                    local_params.rc_d_to_lambda_q();
                                    ray_tracing->calc_ray(local_params);
                                    if (ray_tracing->ray_status == RayStatus::NORMAL) {
                                      theta(i, j) = ray_tracing->theta_f;
                                      phi(i, j) = ray_tracing->phi_f;
                                    } else {
                                      theta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                      phi(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                    }
                                  }
                                }
                              });

    auto &delta_theta = sweep_result.delta_theta;
    auto &delta_phi = sweep_result.delta_phi;

    delta_theta.resize(lgd_size, rc_size);
    delta_phi.resize(lgd_size, rc_size);

    delta_theta = theta.array() - theta_o;
    delta_phi = ((phi.array() - phi_o) * half<Real>()).sin();

    using Point = bg::model::point<int, 2, bg::cs::cartesian>;
    std::vector<Point> theta_roots_index;
    std::vector<Point> phi_roots_index;
    Real d_row, d_col;
    for (size_t i = 1; i < lgd_size; i++) {
      for (size_t j = 1; j < rc_size; j++) {
        d_row = delta_theta(i, j) * delta_theta(i, j - 1);
        d_col = delta_theta(i, j) * delta_theta(i - 1, j);
        if (!isnan(d_row) && !isnan(d_col) && (d_row <= 0 || d_col <= 0)) {
          theta_roots_index.emplace_back(i, j);
        }
        d_row = delta_phi(i, j) * delta_phi(i, j - 1);
        d_col = delta_phi(i, j) * delta_phi(i - 1, j);
        if (!isnan(d_row) && !isnan(d_col) && (d_row <= 0 || d_col <= 0)) {
          phi_roots_index.emplace_back(i, j);
        }
      }
    }

    auto &theta_roots = sweep_result.theta_roots;
    auto &phi_roots = sweep_result.phi_roots;

    theta_roots.resize(theta_roots_index.size(), 2);
    phi_roots.resize(phi_roots_index.size(), 2);

    for (size_t i = 0; i < theta_roots_index.size(); i++) {
      theta_roots(i, 0) = rc_list[theta_roots_index[i].template get<1>()];
      theta_roots(i, 1) = lgd_list[theta_roots_index[i].template get<0>()];
    }
    for (size_t i = 0; i < phi_roots_index.size(); i++) {
      phi_roots(i, 0) = rc_list[phi_roots_index[i].template get<1>()];
      phi_roots(i, 1) = lgd_list[phi_roots_index[i].template get<0>()];
    }

    std::vector<Point> theta_roots_closest_index;
    theta_roots_closest_index.reserve(std::max(theta_roots_index.size(), cutoff));
    std::vector<double> distances(theta_roots_index.size());
    bgi::rtree<Point, bgi::quadratic<16>> rtree(phi_roots_index);
    for (size_t i = 0; i < theta_roots_index.size(); i++) {
      Point p = theta_roots_index[i];
      rtree.query(bgi::nearest(p, 1), std::back_inserter(theta_roots_closest_index));
      distances[i] = bg::distance(p, theta_roots_closest_index.back());
    }

    // sort rows of theta_roots_closest_index by distances
    std::vector<size_t> indices(theta_roots_index.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&distances](size_t i1, size_t i2) { return distances[i1] < distances[i2]; });

    auto &theta_roots_closest = sweep_result.theta_roots_closest;
    theta_roots_closest.resize(theta_roots_index.size(), 2);
    for (size_t i = 0; i < theta_roots_index.size(); i++) {
      theta_roots_closest(i, 0) = rc_list[theta_roots_closest_index[indices[i]].template get<1>()];
      theta_roots_closest(i, 1) = lgd_list[theta_roots_closest_index[indices[i]].template get<0>()];
    }

    // find results
    auto &results = sweep_result.results;
    results.reserve(cutoff);
    tbb::parallel_for(tbb::blocked_range<size_t>(0u, cutoff),
                      [&](const tbb::blocked_range<size_t> &r) {
                        ForwardRayTracingParams<Real> local_params(params);
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                          size_t row = theta_roots_closest_index[indices[i]].template get<0>();
                          size_t col = theta_roots_closest_index[indices[i]].template get<1>();
                          local_params.rc = rc_list[col];
                          local_params.lgd = lgd_list[row];
                          local_params.rc_d_to_lambda_q();
                          using RealToInt = boost::numeric::converter<int, Real, boost::numeric::conversion_traits<int, Real>,
                              boost::numeric::def_overflow_handler, boost::numeric::Floor<Real>>;
                          int period = RealToInt::convert(phi(row, col) / boost::math::constants::two_pi<Real>());
                          try {
                            auto res = find_result(local_params, period, theta_o, phi_o);
                            res.rc = local_params.rc;
                            res.lgd = local_params.lgd;
                            results.push_back(std::move(res));
                          } catch (std::exception &e) {
                            fmt::print(stderr, "Error: {}\n", e.what());
                          }
                        }
                      });

    return sweep_result;
  }
};