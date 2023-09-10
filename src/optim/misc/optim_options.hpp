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

#pragma once

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#ifndef optimlib_inline
    #define optimlib_inline
#endif

// floating point number type

#define OPTIM_FPN_TYPE double
#undef OPTIM_FPN_SMALL_NUMBER
#define OPTIM_FPN_SMALL_NUMBER fp_t(1e-08)

//

namespace optim
{
    using uint_t = unsigned int;
    using fp_t = OPTIM_FPN_TYPE;

    using rand_engine_t = std::mt19937_64;
}

#include <iostream>
#include <Eigen/Dense>

// template<typename eT, int iTr, int iTc>
// using EigenMat = Eigen::Matrix<eT,iTr,iTc>;

namespace optim
{
  using Mat_t = Eigen::Matrix<fp_t, Eigen::Dynamic, Eigen::Dynamic>;

  using ColVec_t = Eigen::Matrix<fp_t, Eigen::Dynamic, 1>;
  using RowVec_t = Eigen::Matrix<fp_t, 1, Eigen::Dynamic>;

  using ColVecInt_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  using RowVecInt_t = Eigen::Matrix<int, 1, Eigen::Dynamic>;

  using ColVecUInt_t = Eigen::Matrix<size_t, Eigen::Dynamic, 1>;
  using RowVecUInt_t = Eigen::Matrix<size_t, 1, Eigen::Dynamic>;
}

#ifndef BMO_RNG_ENGINE_TYPE
    #define BMO_RNG_ENGINE_TYPE optim::rand_engine_t
#endif

#ifndef BMO_CORE_TYPES
    #define BMO_CORE_TYPES

    namespace bmo
    {
        using fp_t = OPTIM_FPN_TYPE;

        using ColVec_t = optim::ColVec_t;
        using RowVec_t = optim::RowVec_t;
        using ColVecInt_t = optim::ColVecInt_t;
        using RowVecInt_t = optim::RowVecInt_t;
        using ColVecUInt_t = optim::ColVecUInt_t;

        using Mat_t = optim::Mat_t;
    }
#endif

//#include "BaseMatrixOps/BaseMatrixOps.hpp"
