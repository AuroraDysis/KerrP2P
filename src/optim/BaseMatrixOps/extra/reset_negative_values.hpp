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

#ifndef BMO_EXTRA_RESET_NEG_VALUES
#define BMO_EXTRA_RESET_NEG_VALUES

inline
void
reset_negative_values(const ColVec_t& vec_in, ColVec_t& vec_out)
{
    const size_t n = BMO_MATOPS_SIZE(vec_in);

    for (size_t i = 0; i < n; ++i) {
        if (vec_in(i) <= 0.0) {
            vec_out(i) = 0.0;
        }
    }
}

inline
void
reset_negative_rows(const ColVec_t& vec_in, Mat_t& mat_out)
{
    const size_t nr = BMO_MATOPS_SIZE(vec_in);
    const size_t nc = BMO_MATOPS_NCOL(mat_out);

    for (size_t i = 0; i < nr; ++i) {
        if (vec_in(i) <= 0.0) {
            for (size_t j = 0; j < nc; ++j) {
                mat_out(i,j) = 0.0;
            }
        }
    }
}

#endif
