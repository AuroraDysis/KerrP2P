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
 * Vector/matrix of ones
 */

#ifndef BMO_MATOPS_ONE_COLVEC

//

#ifdef BMO_ENABLE_ARMA_WRAPPERS
    #define BMO_MATOPS_ONE_COLVEC(n) bmo::ColVec_t(n,arma::fill::ones)
    #define BMO_MATOPS_ONE_ROWVEC(n) bmo::RowVec_t(n,arma::fill::ones)
    #define BMO_MATOPS_ONE_MAT(n,k) bmo::Mat_t(n,k,arma::fill::ones)
#endif

#ifdef BMO_ENABLE_EIGEN_WRAPPERS
    #define BMO_MATOPS_ONE_COLVEC(n) bmo::ColVec_t::Ones(n)
    #define BMO_MATOPS_ONE_ROWVEC(n) bmo::RowVec_t::Ones(n)
    #define BMO_MATOPS_ONE_MAT(n,k) bmo::Mat_t::Ones(n,k)
#endif

//

#endif
