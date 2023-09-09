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
 * Conjugate Gradient method for non-linear optimization
 */

#ifndef _optim_cg_HPP
#define _optim_cg_HPP

/**
 * @brief The Nonlinear Conjugate Gradient (CG) Optimization Algorithm
 *
 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
 * @param opt_objfn the function to be minimized, taking three arguments:
 *   - \c vals_inp a vector of inputs;
 *   - \c grad_out a vector to store the gradient; and
 *   - \c opt_data additional data passed to the user-provided function.
 * @param opt_data additional data passed to the user-provided function
 *
 * @return a boolean value indicating successful completion of the optimization algorithm.
 */

bool 
cg(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data
);

/**
 * @brief The Nonlinear Conjugate Gradient (CG) Optimization Algorithm
 *
 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
 * @param opt_objfn the function to be minimized, taking three arguments:
 *   - \c vals_inp a vector of inputs;
 *   - \c grad_out a vector to store the gradient; and
 *   - \c opt_data additional data passed to the user-provided function.
 * @param opt_data additional data passed to the user-provided function.
 * @param settings parameters controlling the optimization routine.
 *
 * @return a boolean value indicating successful completion of the optimization algorithm.
 */

bool
cg(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data, 
    algo_settings_t& settings
);

//
// internal

namespace internal
{

bool 
cg_impl(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data, 
    algo_settings_t* settings_inp
);

// update function

inline
fp_t
cg_update(
    const ColVec_t& grad, 
    const ColVec_t& grad_p, 
    const ColVec_t& direc, 
    const uint_t iter, 
    const uint_t cg_method, 
    const fp_t cg_restart_threshold
)
{
    // threshold test
    fp_t ratio_value = std::abs( BMO_MATOPS_DOT_PROD(grad_p,grad) ) / BMO_MATOPS_DOT_PROD(grad_p,grad_p);

    if ( ratio_value > cg_restart_threshold ) {
        return 0.0;
    } else {
        fp_t beta = 1.0;

        switch (cg_method)
        {
            case 1: // Fletcher-Reeves (FR)
            {
                beta = BMO_MATOPS_DOT_PROD(grad_p,grad_p) / BMO_MATOPS_DOT_PROD(grad,grad);
                break;
            }

            case 2: // Polak-Ribiere (PR) + 
            {
                beta = BMO_MATOPS_DOT_PROD(grad_p, grad_p - grad) / BMO_MATOPS_DOT_PROD(grad,grad); // max(.,0.0) moved to end
                break;
            }

            case 3: // FR-PR hybrid, see eq. 5.48 in Nocedal and Wright
            {
                if (iter > 1) {
                    const fp_t beta_denom = BMO_MATOPS_DOT_PROD(grad, grad);
                    
                    const fp_t beta_FR = BMO_MATOPS_DOT_PROD(grad_p, grad_p) / beta_denom;
                    const fp_t beta_PR = BMO_MATOPS_DOT_PROD(grad_p, grad_p - grad) / beta_denom;
                    
                    if (beta_PR < - beta_FR) {
                        beta = -beta_FR;
                    } else if (std::abs(beta_PR) <= beta_FR) {
                        beta = beta_PR;
                    } else { // beta_PR > beta_FR
                        beta = beta_FR;
                    }
                } else {
                    // default to PR+
                    beta = BMO_MATOPS_DOT_PROD(grad_p,grad_p - grad) / BMO_MATOPS_DOT_PROD(grad,grad); // max(.,0.0) moved to end
                }
                break;
            }

            case 4: // Hestenes-Stiefel
            {
                beta = BMO_MATOPS_DOT_PROD(grad_p,grad_p - grad) / BMO_MATOPS_DOT_PROD(grad_p - grad,direc);
                break;
            }

            case 5: // Dai-Yuan
            {
                beta = BMO_MATOPS_DOT_PROD(grad_p,grad_p) / BMO_MATOPS_DOT_PROD(grad_p - grad,direc);
                break;
            }

            case 6: // Hager-Zhang
            {
                ColVec_t y = grad_p - grad;

                ColVec_t term_1 = y - 2*direc*(BMO_MATOPS_DOT_PROD(y,y) / BMO_MATOPS_DOT_PROD(y,direc));
                ColVec_t term_2 = grad_p / BMO_MATOPS_DOT_PROD(y,direc);

                beta = BMO_MATOPS_DOT_PROD(term_1,term_2);
                break;
            }
            
            default:
            {
                printf("error: unknown value for cg_method");
                break;
            }
        }

        //

        return std::max(beta, fp_t(0.0));
    }
}

}

//

inline
bool
internal::cg_impl(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data, 
    algo_settings_t* settings_inp
)
{
    // notation: 'p' stands for '+1'.
    
    bool success = false;
    
    const size_t n_vals = BMO_MATOPS_SIZE(init_out_vals);

    // CG settings

    algo_settings_t settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    const int print_level = settings.print_level;
    
    const uint_t conv_failure_switch = settings.conv_failure_switch;
    const size_t iter_max = settings.iter_max;
    const fp_t grad_err_tol = settings.grad_err_tol;
    fp_t rel_sol_change_tol = settings.rel_sol_change_tol;

    if (!settings.cg_settings.use_rel_sol_change_crit) {
        rel_sol_change_tol = -1.0;
    }

    const uint_t cg_method = settings.cg_settings.method; // cg update method
    const fp_t cg_restart_threshold = settings.cg_settings.restart_threshold;

    const fp_t wolfe_cons_1 = settings.cg_settings.wolfe_cons_1; // line search tuning parameter
    const fp_t wolfe_cons_2 = settings.cg_settings.wolfe_cons_2;

    const bool vals_bound = settings.vals_bound;
    
    const ColVec_t lower_bounds = settings.lower_bounds;
    const ColVec_t upper_bounds = settings.upper_bounds;

    const ColVecInt_t bounds_type = determine_bounds_type(vals_bound, n_vals, lower_bounds, upper_bounds);

    // lambda function for box constraints

    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* box_data)> box_objfn \
    = [opt_objfn, vals_bound, bounds_type, lower_bounds, upper_bounds] (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data) \
    -> fp_t 
    {
        if (vals_bound) {
            ColVec_t vals_inv_trans = inv_transform(vals_inp, bounds_type, lower_bounds, upper_bounds);
            fp_t ret;
            
            if (grad_out) {
                ColVec_t grad_obj = *grad_out;

                ret = opt_objfn(vals_inv_trans,&grad_obj,opt_data);

                ColVec_t jacob_vec = BMO_MATOPS_EXTRACT_DIAG( jacobian_adjust(vals_inp,bounds_type,lower_bounds,upper_bounds) );

                *grad_out = BMO_MATOPS_HADAMARD_PROD(jacob_vec, grad_obj);
            } else {
                ret = opt_objfn(vals_inv_trans, nullptr, opt_data);
            }

            return ret;
        } else {
            return opt_objfn(vals_inp, grad_out, opt_data);
        }
    };

    // initialization

    if (! BMO_MATOPS_IS_FINITE(init_out_vals) ) {
        printf("gd error: non-finite initial value(s).\n");
        return false;
    }

    ColVec_t x = init_out_vals;
    ColVec_t d = BMO_MATOPS_ZERO_COLVEC(n_vals);

    if (vals_bound) { // should we transform the parameters?
        x = transform(x, bounds_type, lower_bounds, upper_bounds);
    }

    ColVec_t grad(n_vals); // gradient
    box_objfn(x, &grad, opt_data);

    fp_t grad_err = BMO_MATOPS_L2NORM(grad);

    OPTIM_CG_TRACE(-1, grad_err, 0.0, x, d, grad, 0.0);

    if (grad_err <= grad_err_tol) {
        return true;
    }

    //

    fp_t t_init = 1.0; // initial value for line search

    d = - grad;
    ColVec_t x_p = x, grad_p = grad;

    fp_t t = line_search_mt(t_init, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);

    grad_err = BMO_MATOPS_L2NORM(grad_p);
    fp_t rel_sol_change = BMO_MATOPS_L1NORM( BMO_MATOPS_ARRAY_DIV_ARRAY( (x_p - x), (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x), OPTIM_FPN_SMALL_NUMBER)) ) );
    
    OPTIM_CG_TRACE(0, grad_err, rel_sol_change, x, d, grad, 0.0);

    if (grad_err <= grad_err_tol) {
        if (vals_bound) {
    	    init_out_vals = inv_transform(x_p, bounds_type, lower_bounds, upper_bounds);
    	} else {
            init_out_vals = x_p;
        }
        return true;
    }

    // begin loop

    size_t iter = 0;

    while (grad_err > grad_err_tol && iter < iter_max && rel_sol_change > rel_sol_change_tol) {
        ++iter;

        //

        fp_t beta = cg_update(grad, grad_p, d, iter, cg_method, cg_restart_threshold);

        ColVec_t d_p = - grad_p + beta*d;

        t_init = t * (BMO_MATOPS_DOT_PROD(grad,d) / BMO_MATOPS_DOT_PROD(grad_p,d_p));

        grad = grad_p;

        t = line_search_mt(t_init, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);

        //

        grad_err = BMO_MATOPS_L2NORM(grad_p);
        rel_sol_change = BMO_MATOPS_L1NORM( BMO_MATOPS_ARRAY_DIV_ARRAY( (x_p - x), (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x), OPTIM_FPN_SMALL_NUMBER)) ) );

        d = d_p;
        x = x_p;

        //

        OPTIM_CG_TRACE(iter, grad_err, rel_sol_change, x, d, grad, beta);
    }

    //

    if (vals_bound) {
        x_p = inv_transform(x_p, bounds_type, lower_bounds, upper_bounds);
    }

    error_reporting(init_out_vals, x_p, opt_objfn, opt_data, 
                    success, grad_err, grad_err_tol, iter, iter_max, 
                    conv_failure_switch, settings_inp);

    //

    return success;
}

inline
bool
cg(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data
)
{
    return internal::cg_impl(init_out_vals,opt_objfn,opt_data,nullptr);
}

inline
bool
cg(
    ColVec_t& init_out_vals, 
    std::function<fp_t (const ColVec_t& vals_inp, ColVec_t* grad_out, void* opt_data)> opt_objfn, 
    void* opt_data, 
    algo_settings_t& settings
)
{
    return internal::cg_impl(init_out_vals,opt_objfn,opt_data,&settings);
}

#endif
