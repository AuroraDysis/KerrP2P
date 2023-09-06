/*################################################################################
  ##
  ##   Copyright (C) 2011-2022 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
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
 * Random integer
 */

//
// scalar output

namespace internal
{

template<typename T>
inline
size_t
rind_compute(const T a_par, const T b_par, rand_engine_t& engine)
{
    return static_cast<size_t>( runif_compute(a_par, b_par + T(1), engine, true) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
inline
size_t
rind_type_check(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return rind_compute(static_cast<TC>(a_par),static_cast<TC>(b_par),engine);
}

}

template<typename T1, typename T2>
inline
size_t
rind(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return internal::rind_type_check(a_par,b_par,engine);
}

template<typename T1, typename T2>
inline
common_return_t<T1,T2> 
rind(const T1 a_par, const T2 b_par, const size_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rind(a_par,b_par,engine);
}

//

template<typename T1, typename T2, typename vT = ColVec_t>
inline
vT
rind_vec(size_t n_vals, const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    vT ret_vec(n_vals);

    for (size_t i=0; i < n_vals; ++i) {
        ret_vec(i) = rind(a_par, b_par, engine);
    }

    return ret_vec;
}

//

template<typename T1, typename T2>
inline
Mat_t
rind_mat(size_t nr, size_t nc, const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    Mat_t ret_mat(nr,nc);

    for (size_t j=0; j < nc; ++j) {
        for (size_t i=0; i < nr; ++i) {
            ret_mat(i,j) = rind(a_par, b_par, engine);
        }
    }

    return ret_mat;
}
