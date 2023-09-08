#pragma once

// Modified from https://github.com/Enderdead/nelder-mead

#include <Eigen/Dense>
#include <vector>
#include <tuple>
#include <algorithm>
#include <functional>
#include <iostream>
#include <cstdio>
#include <utility>
#include <boost/math/constants/constants.hpp>

template <typename Real>
struct NelderMeadOptimizerParams {
  Real step = static_cast<Real>(1) / static_cast<Real>(10000);
  Real no_improve_thr = ErrorLimit<Real>::Value * 100;
  int no_improve_break = 10;
  int max_iter = 1000;
  Real alpha = 1;
  Real gamma = 2;
  Real rho = -boost::math::constants::half<Real>();
  Real sigma = boost::math::constants::half<Real>();
  bool log = false;
};

/**
* Nelder Mean function.
*
* This function compute the nelder_meand on your func as loss function. 
*
* @param f  function to minimize, must return a Real scalar.
* @param x_start  initial position.
* @param step  look-around radius in initial step.
* @param no_improv_thr  threshold on improve classification .
* @param no_improve_break  break after no_improve_break iterations without improvement.
* @param max_iter  break after exeed max_iter iterations.
* @param alpha  function to minimize, must return a Real scalar.
* @param gamma  function to minimize, must return a Real scalar.
* @param rho  function to minimize, must return a Real scalar.
* @param sigma  function to minimize, must return a Real scalar.

* @return best x find.
*/
template <int dim, typename Real, typename F>
Eigen::Matrix<Real, dim, 1> NelderMeadOptimizer(F &func, Eigen::Matrix<Real, dim, 1> x_start, const NelderMeadOptimizerParams<Real> &params) {
  using Tuple = std::tuple<Eigen::Matrix<Real, dim, 1>, Real>;
  using Vector = Eigen::Matrix<Real, dim, 1>;

  // unpack params
  const Real &step = params.step;
  const Real &no_improve_thr = params.no_improve_thr;
  int no_improve_break = params.no_improve_break;
  int max_iter = params.max_iter;
  const Real &alpha = params.alpha;
  const Real &gamma = params.gamma;
  const Real &rho = params.rho;
  const Real &sigma = params.sigma;
  const bool log = params.log;

  Real best, prev_best = func(x_start);
  int no_improve = 0;
  std::vector<Tuple> result;

  result.push_back(std::make_tuple(x_start, prev_best));

  for (int i = 0; i < dim; i++) {
    Vector x(x_start);
    x[i] += step;
    result.push_back({x, func(x)});
  }

  int iteration = 0;
  while (true) {

    // order
    std::sort(result.begin(), result.end(), [](const Tuple &a, const Tuple &b) -> bool {
      return (std::get<1>(a) < std::get<1>(b));
    });

    best = std::get<1>(result[0]);

    // break after max_iter
    if (max_iter && iteration >= max_iter) {
      return std::get<0>(result[0]);
    }

    iteration++;

    //break after no_improve_break iterations with no improvement
    if (log) std::cout << "... best so far:  " << best << std::endl;

    if (best < (prev_best - no_improve_thr)) {
      no_improve = 0;
      prev_best = best;
    } else {
      no_improve++;
    }

    if (no_improve >= no_improve_break) {
      return std::get<0>(result[0]);
    }

    //centroid
    Vector centroid_pt = Vector::Zero();
    for (auto it_pt = result.begin(); it_pt != (result.end() - 1); it_pt++) {
      centroid_pt += std::get<0>(*it_pt);
    }

    centroid_pt /= (result.size() - 1);

    // reflection
    Vector reflection_pt(centroid_pt);
    reflection_pt += alpha * (centroid_pt - std::get<0>(result[result.size() - 1]));
    Real reflection_score = func(reflection_pt);

    if ((std::get<1>(result[0]) <= reflection_score) && (reflection_score < std::get<1>(result[result.size() - 2]))) {
      result.pop_back();
      result.push_back(std::make_tuple(reflection_pt, reflection_score));
      continue;
    }

    // expansion
    if (reflection_score < std::get<1>(result[0])) {
      Vector expansion_pt(centroid_pt);
      expansion_pt += gamma * (centroid_pt - std::get<0>(result[result.size() - 1]));
      Real expansion_score = func(expansion_pt);
      if (expansion_score < reflection_score) {
        result.pop_back();
        result.push_back({expansion_pt, expansion_score});
        continue;

      } else {
        result.pop_back();
        result.push_back({reflection_pt, reflection_score});
        continue;

      }
    }
    // Contraction
    Vector contraction_pt(centroid_pt);
    contraction_pt += rho * (centroid_pt - std::get<0>(result[result.size() - 1]));
    Real contraction_score = func(contraction_pt);
    if (contraction_score < std::get<1>(result[result.size() - 1])) {
      result.pop_back();
      result.push_back({contraction_pt, contraction_score});
      continue;
    }

    // Reduction
    auto pt_1 = std::get<0>(result[0]);
    std::vector<Tuple> reduct_result;

    for (auto it_pt = result.begin(); it_pt != result.end(); it_pt++) {
      Vector new_pt(pt_1);
      new_pt += sigma * (std::get<0>(*it_pt) - pt_1);
      Real new_score = func(new_pt);
      reduct_result.push_back({new_pt, new_score});
    }

    result.clear();
    result.insert(result.end(), reduct_result.begin(), reduct_result.end());
  }
}
