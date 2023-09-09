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

#ifndef BMO_MATOPS_RANDN_VEC

//

#ifdef BMO_ENABLE_ARMA_WRAPPERS
    #define BMO_MATOPS_RANDN_VEC(j) arma::randn(j,1)
    #define BMO_MATOPS_RANDN_ROWVEC(j) arma::randn(1,j)
    #define BMO_MATOPS_RANDN_MAT(j,k) arma::randn(j,k)
#endif

#ifdef BMO_ENABLE_EIGEN_WRAPPERS
    // inline
    // ColVec_t
    // bmo_eigen_randn_colvec(size_t nr)
    // {
    //     static std::mt19937 gen{ std::random_device{}() };
    //     static std::normal_distribution<> dist;

    //     // return Eigen::VectorXd{ nr }.unaryExpr([&](double x) { (void)(x); return dist(gen); });
    //     return ColVec_t{ nr }.unaryExpr([&](double x) { (void)(x); return dist(gen); });
    // }

    // inline
    // Mat_t
    // bmo_eigen_randn_mat(size_t nr, size_t nc)
    // {
    //     static std::mt19937 gen{ std::random_device{}() };
    //     static std::normal_distribution<> dist;

    //     // return Eigen::MatrixXd{ nr, nc }.unaryExpr([&](double x) { (void)(x); return dist(gen); });
    //     return Mat_t{ nr, nc }.unaryExpr([&](double x) { (void)(x); return dist(gen); });
    // }

    // #define BMO_MATOPS_RANDN_VEC(j) bmo_eigen_randn_colvec(j)
    // #define BMO_MATOPS_RANDN_ROWVEC(j) (bmo_eigen_randn_colvec(j)).transpose()
    // #define BMO_MATOPS_RANDN_MAT(j,k) bmo_eigen_randn_mat(j,k)
    #define BMO_MATOPS_RANDN_VEC(j) bmo::stats::rnorm_vec<fp_t>(j)
    #define BMO_MATOPS_RANDN_ROWVEC(j) (bmo::stats::rnorm_vec<fp_t>(j)).transpose()
    #define BMO_MATOPS_RANDN_MAT(j,k) bmo::stats::rsnorm_mat<fp_t>(j,k)
#endif

//

#endif
