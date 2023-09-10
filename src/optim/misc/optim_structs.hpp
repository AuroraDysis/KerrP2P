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

struct broyden_settings_t
{
    fp_t par_rho = 0.9;
    fp_t par_sigma_1 = 0.001;
    fp_t par_sigma_2 = 0.001;
};

// 

struct algo_settings_t
{
    // RNG seeding

    size_t rng_seed_value = std::random_device{}();

    // print and convergence options

    int print_level = 0;
    int conv_failure_switch = 0;

    // error tolerance and maxiumum iterations

    size_t iter_max = 2000;

    fp_t grad_err_tol = 1E-08;
    fp_t rel_sol_change_tol = 1E-14;
    fp_t rel_objfn_change_tol = 1E-08;

    // bounds

    bool vals_bound = false;
    
    ColVec_t lower_bounds;
    ColVec_t upper_bounds;

    // values returned upon successful completion

    fp_t opt_fn_value;      // will be returned by the optimization algorithm
    ColVec_t opt_root_fn_values; // will be returned by the root-finding method

    size_t opt_iter;
    fp_t opt_error_value;

    // Broyden
    broyden_settings_t broyden_settings;
};

#endif
