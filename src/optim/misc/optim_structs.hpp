/*################################################################################
  ##
  ##   Copyright (C) 2016-2023 Keith O'Hara
  ##
  ##   This file is part of the OptimLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * Optimization control parameters
 */

#ifndef optim_structs_HPP
#define optim_structs_HPP

// Broyden

template<typename Real>
struct broyden_settings_t
{
    Real par_rho = 0.9;
    Real par_sigma_1 = 0.001;
    Real par_sigma_2 = 0.001;
};

// 
template<typename Real>
struct algo_settings_t {
    using ColVec_t = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    // print and convergence options

    int print_level = 0;
    int conv_failure_switch = 0;

    // error tolerance and maxiumum iterations

    size_t iter_max = 2000;

    Real grad_err_tol = 1E-08;
    Real rel_sol_change_tol = 1E-14;
    Real rel_objfn_change_tol = 1E-08;

    // bounds

    bool vals_bound = false;
    
    ColVec_t lower_bounds;
    ColVec_t upper_bounds;

    // values returned upon successful completion

    Real opt_fn_value;      // will be returned by the optimization algorithm
    ColVec_t opt_root_fn_values; // will be returned by the root-finding method

    size_t opt_iter;
    Real opt_error_value;

    // Broyden
    broyden_settings_t<Real> broyden_settings;
};

#endif
