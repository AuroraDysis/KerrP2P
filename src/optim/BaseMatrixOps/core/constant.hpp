/*################################################################################
  ##
  ##   Copyright (C) 2016-2023 Keith O'Hara
  ##
  ##   This file is part of the BaseMatrixOps C++ library.
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
 * Vector/matrix of user-defined values
 */

#ifndef BMO_MATOPS_CONSTANT_COLVEC

//

#ifdef BMO_ENABLE_ARMA_WRAPPERS
    #define BMO_MATOPS_CONSTANT_COLVEC(n,a) bmo::ColVec_t(n,arma::fill::value(a))
    #define BMO_MATOPS_CONSTANT_ROWVEC(n,a) bmo::RowVec_t(n,arma::fill::value(a))
    #define BMO_MATOPS_CONSTANT_MAT(n,k,a) bmo::Mat_t(n,k,arma::fill::value(a))
#endif

#ifdef BMO_ENABLE_EIGEN_WRAPPERS
    #define BMO_MATOPS_CONSTANT_COLVEC(n,a) bmo::ColVec_t::Constant(n,a)
    #define BMO_MATOPS_CONSTANT_ROWVEC(n,a) bmo::RowVec_t::Constant(n,a)
    #define BMO_MATOPS_CONSTANT_MAT(n,k,a) bmo::Mat_t::Constant(n,k,a)
#endif

//

#endif
