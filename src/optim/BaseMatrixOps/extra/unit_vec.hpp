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
 * unit vector (a zero vector with j-th element equal to one)
 */

#ifndef BMO_EXTRA_UNIT_VEC
#define BMO_EXTRA_UNIT_VEC

inline
ColVec_t
unit_vec(const size_t j, 
         const size_t n)
{
    ColVec_t ret = BMO_MATOPS_ZERO_COLVEC(n);
    ret(j) = 1;

    return ret;
}

#endif
