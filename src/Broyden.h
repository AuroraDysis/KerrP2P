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
   * Li and Fukushima (2000) derivative-free variant of Broyden's method for solving systems of nonlinear equations
   */

#pragma once

#include "Common.h"

#include <iostream>
#include <Eigen/Dense>

#define BMO_MATOPS_EYE(n) Mat_t::Identity(n,n)
#define BMO_MATOPS_L1NORM(x) (x).array().abs().sum()
#define BMO_MATOPS_L2NORM(x) (x).norm()
#define BMO_MATOPS_SIZE(x) static_cast<size_t>((x).size())
#define BMO_MATOPS_TRANSPOSE(x) (x).transpose()
#define BMO_MATOPS_INV(x) (x).inverse()
#define BMO_MATOPS_ZERO_COLVEC(n) Vector::Zero(n)
#define BMO_MATOPS_DOT_PROD(x, y) (x).dot(y)
#define BMO_MATOPS_ABS(x) (x).cwiseAbs()
#define BMO_MATOPS_ARRAY_ADD_SCALAR(x, a) ((x).array() + (a)).matrix()
#define BMO_MATOPS_ARRAY_DIV_ARRAY(x, y)  ((x).array() / (y).array()).matrix()
#define BMO_MATOPS_COUT std::cout
#define BMO_MATOPS_TRANSPOSE_INPLACE(x) (x).transpose()
#define OPTIM_FPN_SMALL_NUMBER (ErrorLimit<Real>::Value * 100)

#ifndef OPTIM_NO_TRACE

#define OPTIM_BROYDEN_DF_TRACE(iter_in, rel_objfn_change_in, rel_sol_change_in, x_in, d_in, objfn_vec_in,   \
                               lambda_in, s_in, y_in, B_in)                                                 \
    if (print_level > 0) {                                                                                  \
        printf("\n");                                                                                       \
        std::cout << "  - Iteration:                          " << iter_in << "\n";                         \
        std::cout << "    Relative change in objective value: " << rel_objfn_change_in << "\n";             \
        std::cout << "    Relative change in solution:        " << rel_sol_change_in << "\n";               \
                                                                                                            \
        if (print_level >= 2) {                                                                             \
            printf("\n");                                                                                   \
            std::cout << "    Current values:\n";                                                           \
            BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE(x_in) << "\n";                                          \
        }                                                                                                   \
                                                                                                            \
        if (print_level >= 3) {                                                                             \
            printf("\n");                                                                                   \
            std::cout << "    Direction:\n";                                                                \
            BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE(d_in) << "\n";                                          \
            std::cout << "    f(x):\n";                                                                     \
            BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE(objfn_vec_in) << "\n";                                  \
        }                                                                                                   \
                                                                                                            \
        if (print_level >= 4) {                                                                             \
            printf("\n");                                                                                   \
            std::cout << "    lambda: " << lambda_in << "\n";                                               \
            std::cout << "    s:\n";                                                                        \
            BMO_MATOPS_COUT << s_in << "\n";                                                                \
            std::cout << "    y:\n";                                                                        \
            BMO_MATOPS_COUT << y_in << "\n";                                                                \
            std::cout << "    B:\n";                                                                        \
            BMO_MATOPS_COUT << B_in << "\n";                                                                \
        }                                                                                                   \
    }

#else

#define OPTIM_BROYDEN_DF_TRACE(iter_in, rel_objfn_change_in, rel_sol_change_in, x_in, d_in, objfn_vec_in,   \
                               lambda_in, s_in, y_in, B_in)

#endif


// Broyden
template<typename Real>
struct BroydenParams
{
	Real par_rho = 0.9;
	Real par_sigma_1 = 0.001;
	Real par_sigma_2 = 0.001;
};

template<typename Real, size_t dim>
struct AlgoParams {
	using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

	// print and convergence options

	int print_level = 0;
	int conv_failure_switch = 0;

	// error tolerance and maxiumum iterations

	size_t iter_max = 5000;

	Real grad_err_tol = 1E-08;
	//  tolerance value controlling convergence based on the relative change in optimal input values
	Real rel_sol_change_tol = ErrorLimit<Real>::Value;
	// tolerance value controlling convergence based on the relative change in objective function.
	Real rel_objfn_change_tol = 1E-08;

	// bounds

	bool vals_bound = false;

	Vector lower_bounds;
	Vector upper_bounds;

	// values returned upon successful completion

	Real opt_fn_value;      // will be returned by the optimization algorithm
	Vector opt_root_fn_values; // will be returned by the root-finding method

	size_t opt_iter;
	Real opt_error_value;

	// Broyden
	BroydenParams<Real> broyden_settings;
};

template<typename Real, size_t dim, typename RF>
struct BroydenDF {
private:
	using Mat_t = Eigen::Matrix<Real, dim, dim>;
	using Vector = Eigen::Matrix<Real, dim, 1>;

	inline
		Real
		df_eta(size_t k) {
		return 1.0 / (k * k);
	}

	inline
		Real
		df_proc_1(
			const Vector& x_vals,
			const Vector& direc,
			Real sigma_1,
			size_t k,
			RF &opt_objfn
		) {
		const Real beta = 0.9;
		const Real eta_k = df_eta(k);
		Real lambda = 1.0;

		// check: || F(x_k + lambda*d_k) || <= ||F(x_k)||*(1+eta_k) - sigma_1*||lambda*d_k||^2

		Real Fx = BMO_MATOPS_L2NORM(opt_objfn(x_vals));
		Real Fx_p = BMO_MATOPS_L2NORM(opt_objfn(x_vals + lambda * direc));
		Real direc_norm2 = BMO_MATOPS_DOT_PROD(direc, direc);

		Real term_2 = sigma_1 * (lambda * lambda) * direc_norm2;
		Real term_3 = eta_k * Fx;

		if (Fx_p <= Fx - term_2 + term_3) {
			return lambda;
		}

		// begin loop

		size_t iter = 0;
		size_t max_iter = 10000;

		while (iter < max_iter) {
			++iter;
			lambda *= beta; // lambda_i = beta^i;

			Fx_p = BMO_MATOPS_L2NORM(opt_objfn(x_vals + lambda * direc));
			term_2 = sigma_1 * (lambda * lambda) * direc_norm2;

			if (Fx_p <= Fx - term_2 + term_3) {
				break;
			}
		}

		//

		return lambda;
	}


	inline
		bool
		broyden_df_impl(
			Vector& init_out_vals,
			RF &opt_objfn,
			std::function<Mat_t(const Vector& vals_inp, void* jacob_data)> jacob_objfn,
			void* jacob_data,
			AlgoParams<Real, dim>* settings_inp
		) {
		// notation: 'p' stands for '+1'.

		bool success = false;

		const size_t n_vals = BMO_MATOPS_SIZE(init_out_vals);

		// Broyden settings

		AlgoParams<Real, dim> settings;

		if (settings_inp) {
			settings = *settings_inp;
		}

		const int print_level = settings.print_level;

		const size_t conv_failure_switch = settings.conv_failure_switch;
		const size_t iter_max = settings.iter_max;
		const Real rel_objfn_change_tol = settings.rel_objfn_change_tol;
		const Real rel_sol_change_tol = settings.rel_sol_change_tol;

		const Real rho = settings.broyden_settings.par_rho;
		const Real sigma_1 = settings.broyden_settings.par_sigma_1;
		const Real sigma_2 = settings.broyden_settings.par_sigma_2;

		// initialization

		Vector x = init_out_vals;
		Vector d = BMO_MATOPS_ZERO_COLVEC(n_vals);

		Mat_t B = BMO_MATOPS_INV(jacob_objfn(x, jacob_data)); // inverse Jacobian

		Vector objfn_vec = opt_objfn(x);

		Real rel_objfn_change = BMO_MATOPS_L2NORM(objfn_vec);

		OPTIM_BROYDEN_DF_TRACE(-1, rel_objfn_change, 0.0, x, d, objfn_vec, 0.0, d, d, B);

		if (rel_objfn_change <= rel_objfn_change_tol) {
			return true;
		}

		Real Fx = BMO_MATOPS_L2NORM(objfn_vec);

		//

		d = -B * objfn_vec; // step 1

		Vector objfn_vec_p = opt_objfn(x + d);

		Real Fx_p = BMO_MATOPS_L2NORM(objfn_vec_p);

		Real lambda;

		if (Fx_p <= rho * Fx - sigma_2 * BMO_MATOPS_DOT_PROD(d, d)) {
			// step 2
			lambda = 1.0;
		}
		else {
			// step 3
			lambda = df_proc_1(x, d, sigma_1, 0, opt_objfn);
		}

		Vector x_p = x + lambda * d; // step 4

		Vector s = x_p - x;
		Vector y = objfn_vec_p - objfn_vec;

		rel_objfn_change = BMO_MATOPS_L2NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(y, (BMO_MATOPS_ARRAY_ADD_SCALAR(
			BMO_MATOPS_ABS(objfn_vec), OPTIM_FPN_SMALL_NUMBER))));
		Real rel_sol_change = BMO_MATOPS_L1NORM(
			BMO_MATOPS_ARRAY_DIV_ARRAY(s, (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x), OPTIM_FPN_SMALL_NUMBER))));

		// B += (y - B*s) * s.t() / BMO_MATOPS_DOT_PROD(s,s); // step 5
		B += (s - B * y) * BMO_MATOPS_TRANSPOSE(y) / (BMO_MATOPS_DOT_PROD(y, y) + 1.0e-14); // update B

		OPTIM_BROYDEN_DF_TRACE(0, rel_objfn_change, rel_sol_change, x_p, d, objfn_vec_p, lambda, s, y, B);

		//

		x = x_p;
		objfn_vec = objfn_vec_p;
		Fx = Fx_p;

		// begin loop

		size_t iter = 0;

		while (rel_objfn_change > rel_objfn_change_tol && rel_sol_change > rel_sol_change_tol && iter < iter_max) {
			++iter;

			// d = arma::solve(B,-objfn_vec);
			d = -B * objfn_vec;
			objfn_vec_p = opt_objfn(x + d);

			//

			Fx_p = BMO_MATOPS_L2NORM(objfn_vec_p);

			if (Fx_p <= rho * Fx - sigma_2 * BMO_MATOPS_DOT_PROD(d, d)) {
				lambda = 1;
			}
			else {
				lambda = df_proc_1(x, d, sigma_1, iter, opt_objfn);
			}

			//

			x_p = x + lambda * d;

			s = x_p - x;
			y = objfn_vec_p - objfn_vec;

			if (iter % 5 == 0) {
				// B = jacob_objfn(x_p,jacob_data);
				B = BMO_MATOPS_INV(jacob_objfn(x_p, jacob_data));
			}
			else {
				// B += (y - B*s) * s.t() / BMO_MATOPS_DOT_PROD(s,s);
				B += (s - B * y) * BMO_MATOPS_TRANSPOSE(y) / (BMO_MATOPS_DOT_PROD(y, y) + 1.0e-14); // update B
			}

			rel_objfn_change = BMO_MATOPS_L2NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(y, (BMO_MATOPS_ARRAY_ADD_SCALAR(
				BMO_MATOPS_ABS(objfn_vec), OPTIM_FPN_SMALL_NUMBER))));
			rel_sol_change = BMO_MATOPS_L1NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(s, (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x),
				OPTIM_FPN_SMALL_NUMBER))));

			//

			x = x_p;
			objfn_vec = objfn_vec_p;
			Fx = Fx_p;

			//

			OPTIM_BROYDEN_DF_TRACE(iter, rel_objfn_change, rel_sol_change, x_p, d, objfn_vec_p, lambda, s, y, B);
		}

		//

		error_reporting(init_out_vals, x_p, opt_objfn,
			success, rel_objfn_change, rel_objfn_change_tol,
			iter, iter_max, conv_failure_switch, settings_inp);

		return success;
	}


	//
	// internal functions

	//

	inline
		bool
		broyden_df_impl(
			Vector& init_out_vals,
			RF &opt_objfn,
			AlgoParams<Real, dim>* settings_inp
		) {
		// notation: 'p' stands for '+1'.

		bool success = false;

		const size_t n_vals = BMO_MATOPS_SIZE(init_out_vals);

		// Broyden settings

		AlgoParams<Real, dim> settings;

		if (settings_inp) {
			settings = *settings_inp;
		}

		const int print_level = settings.print_level;
		const size_t conv_failure_switch = settings.conv_failure_switch;

		const size_t iter_max = settings.iter_max;
		const Real rel_objfn_change_tol = settings.rel_objfn_change_tol;
		const Real rel_sol_change_tol = settings.rel_sol_change_tol;

		const Real rho = settings.broyden_settings.par_rho;
		const Real sigma_1 = settings.broyden_settings.par_sigma_1;
		const Real sigma_2 = settings.broyden_settings.par_sigma_2;

		// initialization

		Vector x = init_out_vals;
		Vector d = BMO_MATOPS_ZERO_COLVEC(n_vals);

		Mat_t B = BMO_MATOPS_EYE(n_vals); // initial approx. to Jacobian

		Vector objfn_vec = opt_objfn(x);

		Real rel_objfn_change = BMO_MATOPS_L2NORM(objfn_vec);

		OPTIM_BROYDEN_DF_TRACE(-1, rel_objfn_change, 0.0, x, d, objfn_vec, 0.0, d, d, B);

		if (rel_objfn_change <= rel_objfn_change_tol) {
			return true;
		}

		Real Fx = BMO_MATOPS_L2NORM(objfn_vec);

		//

		d = -objfn_vec; // step 1

		Vector objfn_vec_p = opt_objfn(x + d);

		Real Fx_p = BMO_MATOPS_L2NORM(objfn_vec_p);

		Real lambda;

		if (Fx_p <= rho * Fx - sigma_2 * BMO_MATOPS_DOT_PROD(d, d)) {
			// step 2
			lambda = 1.0;
		}
		else {
			// step 3
			lambda = df_proc_1(x, d, sigma_1, 0, opt_objfn);
		}

		Vector x_p = x + lambda * d; // step 4

		Vector s = x_p - x;
		Vector y = objfn_vec_p - objfn_vec;

		rel_objfn_change = BMO_MATOPS_L2NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(y, (BMO_MATOPS_ARRAY_ADD_SCALAR(
			BMO_MATOPS_ABS(objfn_vec), OPTIM_FPN_SMALL_NUMBER))));
		Real rel_sol_change = BMO_MATOPS_L1NORM(
			BMO_MATOPS_ARRAY_DIV_ARRAY(s, (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x), OPTIM_FPN_SMALL_NUMBER))));

		// B += (y - B*s) * BMO_MATOPS_TRANSPOSE(s) / BMO_MATOPS_DOT_PROD(s,s); // step 5
		B += (s - B * y) * BMO_MATOPS_TRANSPOSE(y) / (BMO_MATOPS_DOT_PROD(y, y) + 1.0e-14);

		OPTIM_BROYDEN_DF_TRACE(0, rel_objfn_change, rel_sol_change, x_p, d, objfn_vec_p, lambda, s, y, B);

		if (rel_objfn_change <= rel_objfn_change_tol) {
			init_out_vals = x_p;
			return true;
		}

		//

		x = x_p;
		objfn_vec = objfn_vec_p;
		Fx = Fx_p;

		// begin loop

		size_t iter = 0;

		while (rel_objfn_change > rel_objfn_change_tol && rel_sol_change > rel_sol_change_tol && iter < iter_max) {
			++iter;

			// d = arma::solve(B,-objfn_vec);
			d = -B * objfn_vec;

			objfn_vec_p = opt_objfn(x + d);

			//

			Fx_p = BMO_MATOPS_L2NORM(objfn_vec_p);

			if (Fx_p <= rho * Fx - sigma_2 * BMO_MATOPS_DOT_PROD(d, d)) {
				lambda = 1.0;
			}
			else {
				lambda = df_proc_1(x, d, sigma_1, iter, opt_objfn);
			}

			//

			x_p = x + lambda * d;

			s = x_p - x;
			y = objfn_vec_p - objfn_vec;

			// B += (y - B*s) * s.t() / BMO_MATOPS_DOT_PROD(s,s);
			B += (s - B * y) * BMO_MATOPS_TRANSPOSE(y) / (BMO_MATOPS_DOT_PROD(y, y) + 1.0e-14); // update B

			//

			rel_objfn_change = BMO_MATOPS_L2NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(y, (BMO_MATOPS_ARRAY_ADD_SCALAR(
				BMO_MATOPS_ABS(objfn_vec), OPTIM_FPN_SMALL_NUMBER))));
			rel_sol_change = BMO_MATOPS_L1NORM(BMO_MATOPS_ARRAY_DIV_ARRAY(s, (BMO_MATOPS_ARRAY_ADD_SCALAR(BMO_MATOPS_ABS(x),
				OPTIM_FPN_SMALL_NUMBER))));

			x = x_p;
			objfn_vec = objfn_vec_p;
			Fx = Fx_p;

			//

			OPTIM_BROYDEN_DF_TRACE(iter, rel_objfn_change, rel_sol_change, x_p, d, objfn_vec_p, lambda, s, y, B);
		}

		//

		error_reporting(init_out_vals, x_p, opt_objfn, success, rel_objfn_change, rel_objfn_change_tol, iter,
			iter_max, conv_failure_switch, settings_inp);

		return success;
	}

	/*inline
		void
		error_reporting(
			Vector& out_vals,
			const Vector& x_p,
			std::function<Real(const Vector& vals_inp, Vector* grad_out, void* opt_data)> opt_objfn,
			void* opt_data,
			bool& success,
			const Real err,
			const Real err_tol,
			const size_t iter,
			const size_t iter_max,
			const int conv_failure_switch,
			AlgoParams<Real, dim>* settings_inp
		)
	{
		success = false;

		if (conv_failure_switch == 0) {
			out_vals = x_p;

			if (err <= err_tol && iter <= iter_max) {
				success = true;
			}
		}
		else if (conv_failure_switch == 1) {
			out_vals = x_p;

			if (err <= err_tol && iter <= iter_max) {
				success = true;
			}
			else if (err > err_tol && iter < iter_max) {
				printf("optim failure: iter_max not reached but algorithm stopped.\n");
				printf("optim failure: returned best guess.\n");

				std::cout << "iterations: " << iter << ". error: " << err << std::endl;
			}
			else {
				printf("optim failure: iter_max reached before convergence could be achieved.\n");
				printf("optim failure: returned best guess.\n");

				std::cout << "iterations: " << iter << ". error: " << err << std::endl;
			}
		}
		else if (conv_failure_switch == 2) {
			if (err <= err_tol && iter <= iter_max) {
				out_vals = x_p;
				success = true;
			}
			else {
				printf("optim failure: iter_max reached before convergence could be achieved.\n");
				printf("optim failure: best guess:\n");

				BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE_INPLACE(x_p) << "\n";
				std::cout << "iterations: " << iter << ". error: " << err << std::endl;
			}
		}
		else {
			printf("optim failure: unrecognized conv_failure_switch value.\n");
			success = false;
		}

		if (settings_inp) {
			settings_inp->opt_fn_value = opt_objfn(x_p, nullptr, opt_data);
			settings_inp->opt_iter = iter;
			settings_inp->opt_error_value = err;
		}
	}

	inline
		void
		error_reporting(
			Vector& out_vals,
			const Vector& x_p,
			std::function<Real(const Vector& vals_inp, Vector* grad_out, void* opt_data)> opt_objfn,
			void* opt_data,
			bool& success,
			const int conv_failure_switch,
			AlgoParams<Real, dim>* settings_inp
		)
	{
		if (conv_failure_switch == 0 || conv_failure_switch == 1) {
			out_vals = x_p;
		}
		else if (conv_failure_switch == 2) {
			if (success) {
				out_vals = x_p;
			}
		}
		else {
			printf("optim failure: unrecognized conv_failure_switch value.\n");
			success = false;
		}

		if (settings_inp) {
			settings_inp->opt_fn_value = opt_objfn(x_p, nullptr, opt_data);
		}
	}*/

	inline
		void
		error_reporting(
			Vector& out_vals,
			const Vector& x_p,
			RF &opt_objfn,
			bool& success,
			const Real err,
			const Real err_tol,
			const size_t iter,
			const size_t iter_max,
			const int conv_failure_switch,
			AlgoParams<Real, dim>* settings_inp
		)
	{
		success = false;

		if (conv_failure_switch == 0) {
			out_vals = x_p;

			if (err <= err_tol && iter <= iter_max) {
				success = true;
			}
		}
		else if (conv_failure_switch == 1) {
			out_vals = x_p;

			if (err <= err_tol && iter <= iter_max) {
				success = true;
			}
			else if (err > err_tol && iter < iter_max) {
				printf("optim failure: iter_max not reached but algorithm stopped.\n");
				printf("optim failure: returned best guess.\n");

				std::cout << "iterations: " << iter << ". error: " << err << std::endl;
			}
			else {
				printf("optim failure: iter_max reached before convergence could be achieved.\n");
				printf("optim failure: returned best guess.\n");

				std::cout << "error: " << err << std::endl;
			}
		}
		else if (conv_failure_switch == 2) {
			if (err <= err_tol && iter <= iter_max) {
				out_vals = x_p;
				success = true;
			}
			else if (err > err_tol && iter < iter_max) {
				printf("optim failure: iter_max not reached but algorithm stopped.\n");
				printf("optim failure: best guess:\n");

				BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE_INPLACE(x_p) << "\n";
				std::cout << "error: " << err << std::endl;
			}
			else {
				printf("optim failure: iter_max reached before convergence could be achieved.\n");
				printf("optim failure: best guess:\n");

				BMO_MATOPS_COUT << BMO_MATOPS_TRANSPOSE_INPLACE(x_p) << "\n";
				std::cout << "error: " << err << std::endl;
			}
		}
		else {
			printf("optim failure: unrecognized conv_failure_switch value.\n");
			success = false;
		}

		if (settings_inp) {
			settings_inp->opt_root_fn_values = opt_objfn(x_p);
			settings_inp->opt_iter = iter;
			settings_inp->opt_error_value = err;
		}
	}

	//

	//inline
	//	void
	//	error_reporting(
	//		Vector& out_vals,
	//		const Vector& x_p,
	//		std::function<Real(const Vector& vals_inp, Vector* grad_out, Mat_t* hess_out, void* opt_data)> opt_objfn,
	//		void* opt_data,
	//		bool& success,
	//		const Real err,
	//		const Real err_tol,
	//		const size_t iter,
	//		const size_t iter_max,
	//		const int conv_failure_switch,
	//		AlgoParams<Real, dim>* settings_inp
	//	)
	//{
	//	std::function<Real(const Vector& vals_inp, Vector* grad_out, void* opt_data)> lam_objfn \
	//		= [opt_objfn](const Vector& vals_inp, Vector* grad_out, void* opt_data)
	//		-> Real
	//		{
	//			return opt_objfn(vals_inp, grad_out, nullptr, opt_data);
	//		};

	//	//

	//	error_reporting(out_vals, x_p, lam_objfn, opt_data, success, err, err_tol, iter, iter_max, conv_failure_switch, settings_inp);
	//}
public:
	/**
	 * @brief Derivative-free variant of Broyden's method due to Li and Fukushima (2000)
	 *
	 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
	 * @param opt_objfn the function to be minimized, taking three arguments:
	 *   - \c vals_inp a vector of inputs; and
	 *   - \c opt_data additional data passed to the user-provided function.
	 * @param opt_data additional data passed to the user-provided function.
	 *
	 * @return a boolean value indicating successful completion of the optimization algorithm.
	 */
	inline
		bool
		broyden_df(
			Vector& init_out_vals,
			RF &opt_objfn
		) {
		return broyden_df_impl(init_out_vals, opt_objfn, opt_data, nullptr);
	}

	/**
	 * @brief Derivative-free variant of Broyden's method due to Li and Fukushima (2000)
	 *
	 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
	 * @param opt_objfn the function to be minimized, taking three arguments:
	 *   - \c vals_inp a vector of inputs; and
	 *   - \c opt_data additional data passed to the user-provided function.
	 * @param opt_data additional data passed to the user-provided function.
	 * @param settings parameters controlling the optimization routine.
	 *
	 * @return a boolean value indicating successful completion of the optimization algorithm.
	 */

	inline
		bool
		broyden_df(
			Vector& init_out_vals,
			RF &opt_objfn,
			AlgoParams<Real, dim>& settings
		) {
		return broyden_df_impl(init_out_vals, opt_objfn, &settings);
	}

	// derivative-free method with jacobian

	/**
	 * @brief Derivative-free variant of Broyden's method due to Li and Fukushima (2000)
	 *
	 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
	 * @param opt_objfn the function to be minimized, taking three arguments:
	 *   - \c vals_inp a vector of inputs; and
	 *   - \c opt_data additional data passed to the user-provided function.
	 * @param opt_data additional data passed to the user-provided function.
	 * @param jacob_objfn a function to calculate the Jacobian matrix, taking two arguments:
	 *   - \c vals_inp a vector of inputs; and
	 *   - \c jacob_data additional data passed to the Jacobian function.
	 * @param jacob_data additional data passed to the Jacobian function.
	 *
	 * @return a boolean value indicating successful completion of the optimization algorithm.
	 */

	 // derivative-free method with jacobian

	inline
		bool
		broyden_df(
			Vector& init_out_vals,
			RF &opt_objfn,
			std::function<Mat_t(const Vector& vals_inp, void* jacob_data)> jacob_objfn,
			void* jacob_data
		) {
		return broyden_df_impl(init_out_vals, opt_objfn, opt_data, jacob_objfn, jacob_data, nullptr);
	}

	/**
	 * @brief Derivative-free variant of Broyden's method due to Li and Fukushima (2000)
	 *
	 * @param init_out_vals a column vector of initial values, which will be replaced by the solution upon successful completion of the optimization algorithm.
	 * @param opt_objfn the function to be minimized, taking three arguments:
	 *   - \c vals_inp a vector of inputs; and
	 * @param jacob_objfn a function to calculate the Jacobian matrix, taking two arguments:
	 *   - \c vals_inp a vector of inputs; and
	 *   - \c jacob_data additional data passed to the Jacobian function.
	 * @param jacob_data additional data passed to the Jacobian function.
	 * @param settings parameters controlling the optimization routine.
	 *
	 * @return a boolean value indicating successful completion of the optimization algorithm.
	 */

	inline
		bool
		broyden_df(
			Vector& init_out_vals,
			RF &opt_objfn,
			std::function<Mat_t(const Vector& vals_inp, void* jacob_data)> jacob_objfn,
			void* jacob_data,
			AlgoParams<Real, dim>& settings
		) {
		return broyden_df_impl(init_out_vals, opt_objfn, jacob_objfn, jacob_data, &settings);
	}
};
