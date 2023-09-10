#pragma once

#include "ForwardRayTracing.h"

#include "optim/optim.hpp"

#include <optional>
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

    Matrix lambda;
    Matrix eta;

    Matrix delta_theta;
    Matrix delta_phi;

    PointVector theta_roots;
    PointVector phi_roots;

    PointVector theta_roots_closest;

    std::vector<ForwardRayTracingResult<Real, Complex>> results;
};

template<typename Real, typename Complex>
struct FindRootResult {
    bool success;
    std::string fail_reason;
    std::optional<ForwardRayTracingResult<Real, Complex>> root;
};

template<typename Real, typename Complex>
class RootFunctor {
private:
    using Vector = Eigen::Vector<Real, 2>;

    ForwardRayTracingParams<Real> &params;
    const Real theta_o;
    const Real phi_o;
    const int period;
    const bool fixed_period;
    const Real two_pi = boost::math::constants::two_pi<Real>();

public:
    std::shared_ptr<ForwardRayTracing<Real, Complex>> ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();

    RootFunctor(ForwardRayTracingParams<Real> &params_, Real theta_o_, Real phi_o_)
            : params(params_),
              period(std::numeric_limits<int>::max()),
              fixed_period(false),
              theta_o(std::move(
                      theta_o_)),
              phi_o(std::move(phi_o_)) {
        ray_tracing->calc_t_f = false;
    }

    RootFunctor(ForwardRayTracingParams<Real> &params_, int period_, Real theta_o_, Real phi_o_)
            : params(params_),
              period(period_),
              fixed_period(true),
              theta_o(std::move(
                      theta_o_)),
              phi_o(std::move(phi_o_)) {
        ray_tracing->calc_t_f = false;
    }

    Vector operator()(Vector x, void *opt_data) {
        auto &rc = x[0];
        auto &lgd = x[1];
        params.rc = rc;
        params.lgd = lgd;
        params.rc_d_to_lambda_q();
        ray_tracing->calc_ray(params);

        if (ray_tracing->ray_status != RayStatus::NORMAL) {
            fmt::println("ray status: {}", ray_status_to_str(ray_tracing->ray_status));
            return Vector::Constant(std::numeric_limits<Real>::quiet_NaN());
        }

        Vector residual;
        residual[0] = ray_tracing->theta_f - theta_o;
        if (fixed_period) {
            residual[1] = ray_tracing->phi_f - phi_o - period * two_pi;
        } else {
            residual[1] = sin((ray_tracing->phi_f - phi_o) * half<Real>());
        }

#ifdef PRINT_DEBUG
        fmt::println("rc: {}, lgd: {}, theta_f: {}, phi_f: {}", x[0], x[1], ray_tracing->theta_f, ray_tracing->phi_f);
        fmt::println("residual: {}, {}", residual[0], residual[1]);
#endif
        return residual;
    }
};

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename Real, typename Complex>
struct ForwardRayTracingUtils {
    static ForwardRayTracingResult<Real, Complex> calc_ray(const ForwardRayTracingParams<Real> &params) {
        auto ray_tracing = ForwardRayTracing<Real, Complex>::get_from_cache();
        ray_tracing->calc_ray(params);
        return ray_tracing->to_result();
    }

    static FindRootResult<Real, Complex>
    find_root_period(ForwardRayTracingParams<Real> &params, int period, Real theta_o, Real phi_o) {
        Eigen::Vector<Real, 2> x = Eigen::Vector<Real, 2>::Zero(2);
        x << params.rc, params.lgd;

        auto root_functor =
                period == std::numeric_limits<int>::max() ? RootFunctor<Real, Complex>(params, theta_o, phi_o)
                                                          : RootFunctor<Real, Complex>(params, period,
                                                                                       std::move(theta_o),
                                                                                       std::move(phi_o));
        optim::algo_settings_t<Real> settings;
        settings.rel_objfn_change_tol = ErrorLimit<Real>::Value;
        settings.rel_sol_change_tol = ErrorLimit<Real>::Value;
        optim::BroydenDF<Real, 2> solver;
        solver.broyden_df(x, root_functor, nullptr, settings);
        auto residual = root_functor(x, nullptr);

        FindRootResult<Real, Complex> result;

        if (root_functor.ray_tracing->ray_status != RayStatus::NORMAL) {
            result.success = false;
            result.fail_reason = fmt::format("ray status: {}", ray_status_to_str(root_functor.ray_tracing->ray_status));
            return result;
        }

        if (residual.norm() > ErrorLimit<Real>::Value * 1000) {
            result.success = false;
            result.fail_reason = fmt::format("residual: {}", residual.norm());
            return result;
        }

        result.success = true;
        result.root = root_functor.ray_tracing->to_result();
        return result;
    }

    static FindRootResult<Real, Complex>
    find_root(ForwardRayTracingParams<Real> &params, Real theta_o, Real phi_o) {
        return find_root_period(params, std::numeric_limits<int>::max(), theta_o, phi_o);
    }

    static ForwardRayTracingResult<Real, Complex> refine_result(ForwardRayTracingResult<Real, Complex> &res) {
    }

    static SweepResult<Real, Complex>
    sweep_rc_d(const ForwardRayTracingParams<Real> &params, Real theta_o, Real phi_o, const std::vector<Real> &rc_list,
               const std::vector<Real> &lgd_list, size_t cutoff) {
        size_t rc_size = rc_list.size();
        size_t lgd_size = lgd_list.size();

        SweepResult<Real, Complex> sweep_result;

        auto &theta = sweep_result.theta;
        auto &phi = sweep_result.phi;

        auto &delta_theta = sweep_result.delta_theta;
        auto &delta_phi = sweep_result.delta_phi;

        auto &lambda = sweep_result.lambda;
        auto &eta = sweep_result.eta;

        theta.resize(lgd_size, rc_size);
        phi.resize(lgd_size, rc_size);

        delta_theta.resize(lgd_size, rc_size);
        delta_phi.resize(lgd_size, rc_size);

        lambda.resize(lgd_size, rc_size);
        eta.resize(lgd_size, rc_size);

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
                                                  delta_theta(i, j) = theta(i, j) - theta_o;
                                                  delta_phi(i, j) = sin((phi(i, j) - phi_o) * half<Real>());
                                                  lambda(i, j) = ray_tracing->lambda;
                                                  eta(i, j) = ray_tracing->eta;
                                              } else {
                                                  theta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                                  phi(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                                  delta_theta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                                  delta_phi(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                                  lambda(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                                  eta(i, j) = std::numeric_limits<Real>::quiet_NaN();
                                              }
                                          }
                                      }
                                  });

        using Point = bg::model::point<int, 2, bg::cs::cartesian>;
        std::vector<Point> theta_roots_index;
        std::vector<Point> phi_roots_index;
        int d_row, d_col, d_row_lambda, d_col_lambda;
        for (size_t i = 1; i < lgd_size; i++) {
            for (size_t j = 1; j < rc_size; j++) {
                d_row = sgn(delta_theta(i, j)) * sgn(delta_theta(i, j - 1));
                d_col = sgn(delta_theta(i, j)) * sgn(delta_theta(i - 1, j));
                if (!isnan(delta_theta(i, j)) && !isnan(delta_theta(i, j - 1)) && !isnan(delta_theta(i - 1, j)) &&
                    (d_row <= 0 || d_col <= 0)) {
                    theta_roots_index.emplace_back(i, j);
                }
                d_row = sgn(delta_phi(i, j)) * sgn(delta_phi(i, j - 1));
                d_col = sgn(delta_phi(i, j)) * sgn(delta_phi(i - 1, j));
                d_row_lambda = sgn(lambda(i, j)) * sgn(lambda(i, j - 1));
                d_col_lambda = sgn(lambda(i, j)) * sgn(lambda(i - 1, j));
                if (!isnan(delta_phi(i, j)) && !isnan(delta_phi(i, j - 1)) && !isnan(delta_phi(i - 1, j)) &&
                    !isnan(lambda(i, j)) && !isnan(lambda(i, j - 1)) && !isnan(lambda(i - 1, j)) && d_row_lambda > 0 &&
                    d_col_lambda > 0 &&
                    (d_row <= 0 || d_col <= 0)) {
                    phi_roots_index.emplace_back(i, j);
                }
            }
        }

        if (theta_roots_index.empty() && phi_roots_index.empty()) {
            return sweep_result;
        }

        auto &theta_roots = sweep_result.theta_roots;
        theta_roots.resize(theta_roots_index.size(), 2);
        for (size_t i = 0; i < theta_roots_index.size(); i++) {
            theta_roots(i, 0) = rc_list[theta_roots_index[i].template get<1>()];
            theta_roots(i, 1) = lgd_list[theta_roots_index[i].template get<0>()];
        }

        if (phi_roots_index.empty()) {
            return sweep_result;
        }

        auto &phi_roots = sweep_result.phi_roots;
        phi_roots.resize(phi_roots_index.size(), 2);
        for (size_t i = 0; i < phi_roots_index.size(); i++) {
            phi_roots(i, 0) = rc_list[phi_roots_index[i].template get<1>()];
            phi_roots(i, 1) = lgd_list[phi_roots_index[i].template get<0>()];
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
            theta_roots_closest(i, 1) = lgd_list[theta_roots_closest_index[indices[i]].template get<0>()];
        }

        // find results
        auto &results = sweep_result.results;
        cutoff = std::min(cutoff, indices.size());
        results.reserve(cutoff);
        std::mutex results_mutex;
        tbb::parallel_for(tbb::blocked_range<size_t>(0u, cutoff),
                          [&](const tbb::blocked_range<size_t> &r) {
                              ForwardRayTracingParams<Real> local_params(params);
                              Real two_pi = boost::math::constants::two_pi<Real>();
                              for (size_t i = r.begin(); i != r.end(); ++i) {
                                  size_t row = theta_roots_closest_index[indices[i]].template get<0>();
                                  size_t col = theta_roots_closest_index[indices[i]].template get<1>();
                                  local_params.rc = rc_list[col];
                                  local_params.lgd = lgd_list[row];
                                  local_params.rc_d_to_lambda_q();
                                  int period = MY_FLOOR<Real>::convert(phi(row, col) / two_pi);
                                  auto root_res = find_root_period(local_params, period, theta_o, phi_o);
                                  if (root_res.success) {
                                      auto root = *std::move(root_res.root);
                                      root.rc = local_params.rc;
                                      root.lgd = local_params.lgd;
                                      std::lock_guard<std::mutex> lock(results_mutex);
                                      results.push_back(std::move(root));
                                  } else {
                                      fmt::println("find root failed: {}\n", root_res.fail_reason);
                                  }
                              }
                          });
        // remove duplicate by rc and lgd using std::set
        // abs(r1.rc - r2.rc) < 10000 * ErrorLimit<Real>::Value &&
        //                                   abs(r1.lgd - r2.lgd) < 10000 * ErrorLimit<Real>::Value
        std::vector<size_t> duplicated_index;
        for (size_t i = 0; i < results.size(); i++) {
            for (size_t j = i + 1; j < results.size(); j++) {
                if (abs(results[i].rc - results[j].rc) < 10000 * ErrorLimit<Real>::Value &&
                    abs(results[i].lgd - results[j].lgd) < 10000 * ErrorLimit<Real>::Value) {
                    duplicated_index.push_back(j);
                    break;
                }
            }
        }
        // remove duplicated results
        std::sort(duplicated_index.begin(), duplicated_index.end());
        for (size_t i = duplicated_index.size(); i > 0; i--) {
            results.erase(results.begin() + duplicated_index[i - 1]);
        }

        return sweep_result;
    }
};