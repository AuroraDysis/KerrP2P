#pragma once

#include "ForwardRayTracing.h"

#define OPTIM_ENABLE_EIGEN_WRAPPERS

#include "optim/optim.hpp"

#include <oneapi/tbb.h>
#include <Eigen/Dense>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<typename Real>
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
};


template<typename Real, typename Complex>
class RootFunctor {
private:
  const std::shared_ptr<ForwardRayTracingParams<Real>> params;
  const Real theta_o;
  const Real phi_o;

  const Real two_pi = boost::math::constants::two_pi<Real>();
  Real phi_tmp;

public:
  std::shared_ptr<ForwardRayTracing<Real, Complex>> ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();

  RootFunctor(ForwardRayTracingParams<Real> params_, Real theta_o_, Real phi_o_) : params(
      std::make_shared<ForwardRayTracingParams<Real>>(params_)), theta_o(std::move(theta_o_)), phi_o(std::move(phi_o_)) {
    ray_tracing->calc_t_f = false;
  }

  Eigen::Vector2d operator()(const Eigen::Vector2d &x, void *opt_data) {
    auto &rc = x[0];
    auto &lgd = x[1];
    params->rc = rc;
    params->lgd = lgd;
    params->rc_d_to_lambda_q();
    ray_tracing->calc_ray(*params);

    if (ray_tracing->ray_status != RayStatus::NORMAL) {
      fmt::println("ray status: {}", ray_status_to_str(ray_tracing->ray_status));
      throw std::runtime_error("Ray tracing failed: ray status is not normal");
    }

    Eigen::Vector2d residual;
    residual[0] = ray_tracing->theta_f - theta_o;
    phi_tmp = fmod(ray_tracing->phi_f, two_pi);
    residual[1] = (phi_tmp < 0 ? phi_tmp + two_pi : phi_tmp) - phi_o;

    if (isnan(residual[0]) || isnan(residual[1])) {
      throw std::runtime_error("Ray tracing failed: residual is NaN");
    }

#ifdef PRINT_DEBUG
    fmt::println("rc: {}, lgd: {}, theta_f: {}, phi_f: {}", x[0], x[1], ray_tracing->theta_f, ray_tracing->phi_f);
    fmt::println("residual: {}, {}", residual[0], residual[1]);
#endif
    return residual;
  }
};

template<typename Real, typename Complex>
struct ForwardRayTracingUtils {
  static ForwardRayTracingResult<Real, Complex> calc_ray(const ForwardRayTracingParams<Real> &params) {
    auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
    ray_tracing->calc_ray(params);
    return ray_tracing->to_result();
  }

  static ForwardRayTracingResult<Real, Complex> find_result(const ForwardRayTracingParams<Real> &params, Real theta_o, Real phi_o) {
    // std::array<Real, 2> x = {params.rc, params.lgd};
    Eigen::VectorXd x = Eigen::VectorXd::Zero(2);
    x << params.rc, params.lgd;

    RootFunctor<Real, Complex> root_functor(params, std::move(theta_o), std::move(phi_o));
    optim::algo_settings_t settings;
    settings.broyden_settings.par_rho = 0.5;
    optim::broyden_df(x, root_functor, nullptr, settings);
    root_functor(x, nullptr);

    return root_functor.ray_tracing->to_result();
  }

  static ForwardRayTracingResult<Real, Complex> refine_result(ForwardRayTracingResult<Real, Complex> &res) {

  }

  static SweepResult<Real>
  sweep_rc_d(Real a, Real r_s, Real theta_s, Real r_o, Sign nu_r, Sign nu_theta, Real theta_o, Real phi_o,
             const std::vector<Real> &rc_list, const std::vector<Real> &d_list) {
    size_t rc_size = rc_list.size();
    size_t d_size = d_list.size();

    SweepResult<Real> sweep_result;

    auto &theta = sweep_result.theta;
    auto &phi = sweep_result.phi;
    theta.resize(d_size, rc_size);
    phi.resize(d_size, rc_size);

    // rc and d
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range2d<size_t>(0u, d_size, 0u, rc_size),
                              [&](const oneapi::tbb::blocked_range2d<size_t, size_t> &r) {
                                auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
                                Real two_pi = boost::math::constants::two_pi<Real>();
                                Real phi_tmp;
                                for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                                  for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                                    ray_tracing->calc_ray_by_rc_d(a, r_s, theta_s, r_o, nu_r, nu_theta, rc_list[j],
                                                                  d_list[i]);
                                    if (ray_tracing->ray_status == RayStatus::NORMAL) {
                                      theta(i, j) = ray_tracing->theta_f;
                                      phi_tmp = fmod(ray_tracing->phi_f, two_pi);
                                      phi(i, j) = phi_tmp < 0 ? phi_tmp + two_pi : phi_tmp;
                                    } else {
                                      theta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                      theta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                    }
                                  }
                                }
                              });

    auto &delta_theta = sweep_result.delta_theta;
    auto &delta_phi = sweep_result.delta_phi;

    delta_theta.resize(d_size, rc_size);
    delta_phi.resize(d_size, rc_size);

    delta_theta = theta.array() - theta_o;
    delta_phi = phi.array() - phi_o;

    using Point = bg::model::point<int, 2, bg::cs::cartesian>;
    std::vector<Point> theta_roots_index;
    std::vector<Point> phi_roots_index;
    Real d_row, d_col;
    for (size_t i = 1; i < d_size; i++) {
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
      theta_roots(i, 1) = d_list[theta_roots_index[i].template get<0>()];
    }
    for (size_t i = 0; i < phi_roots_index.size(); i++) {
      phi_roots(i, 0) = rc_list[phi_roots_index[i].template get<1>()];
      phi_roots(i, 1) = d_list[phi_roots_index[i].template get<0>()];
    }

    std::vector<Point> theta_roots_closest_index;
    theta_roots_closest_index.reserve(theta_roots_index.size());
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
      theta_roots_closest(i, 1) = d_list[theta_roots_closest_index[indices[i]].template get<0>()];
    }

    return sweep_result;
  }
};