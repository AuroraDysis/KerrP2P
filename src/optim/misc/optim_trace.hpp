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
 * Optimization trace functions
 */

#ifndef OPTIM_TRACE_DEFS
#define OPTIM_TRACE_DEFS

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

#endif

