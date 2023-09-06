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
 * Sample from a uniform distribution
 */

//
// scalar output

namespace internal
{

template<typename T>
constexpr
bool
unif_sanity_check(const T a_par, const T b_par)
noexcept
{
    return( any_nan(a_par,b_par) ? \
                false :
            //
            a_par >= b_par ? \
                false :
            //
                true );
}

template<typename T>
inline
T
runif_compute(const T a_par, const T b_par, rand_engine_t& engine, 
              const bool skip_sanity_check = false)
{
    if (!skip_sanity_check) {
        if (!unif_sanity_check(a_par,b_par)) {
            return std::numeric_limits<T>::quiet_NaN();
        }
    }
    
    // convert from [a,b) to (a,b)

    T a_par_adj = std::nextafter(a_par, b_par);
    std::uniform_real_distribution<T> unif_dist(a_par_adj, b_par);

    return unif_dist(engine);
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
inline
TC
runif_type_check(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return runif_compute(static_cast<TC>(a_par),static_cast<TC>(b_par),engine);
}

}

template<typename T1, typename T2>
inline
common_return_t<T1,T2> 
runif(const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    return internal::runif_type_check(a_par,b_par,engine);
}

template<typename T1, typename T2>
inline
common_return_t<T1,T2> 
runif(const T1 a_par, const T2 b_par, const size_t seed_val)
{
    rand_engine_t engine(seed_val);
    return runif(a_par,b_par,engine);
}

template<typename T>
inline
T
runif(rand_engine_t& engine)
{
    return runif(T(0),T(1),engine);
}

template<typename T>
inline
T
runif()
{
    return runif(T(0),T(1));
}

//

namespace internal
{

template<typename T1, typename T2, typename vT>
inline
void
runif_vec_inplace(size_t n_vals, const T1 a_par, const T2 b_par, rand_engine_t& engine, vT& vec_out)
{
    for (size_t i=0; i < n_vals; ++i) {
        vec_out(i) = runif(a_par, b_par, engine);
    }
}

template<typename fT, typename vT>
inline
void
runif_vec_inplace(size_t n_vals, rand_engine_t& engine, vT& vec_out)
{
    for (size_t i=0; i < n_vals; ++i) {
        vec_out(i) = runif(fT(0), fT(1), engine);
    }
}

}

template<typename T1, typename T2, typename vT = ColVec_t>
inline
vT
runif_vec(size_t n_vals, const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    vT ret_vec(n_vals);

    internal::runif_vec_inplace(n_vals, a_par, b_par, engine, ret_vec);

    return ret_vec;
}

template<typename T, typename vT = ColVec_t>
inline
vT
runif_vec(size_t n_vals, rand_engine_t& engine)
{
    return runif_vec(n_vals, T(0), T(1), engine);
}

template<typename T, typename vT = ColVec_t>
inline
vT
runif_vec(size_t n_vals)
{
    rand_engine_t engine(std::random_device{}());
    
    return runif_vec(n_vals, T(0), T(1), engine);
}

//

template<typename T1, typename T2>
inline
Mat_t
runif_mat(size_t nr, size_t nc, const T1 a_par, const T2 b_par, rand_engine_t& engine)
{
    Mat_t ret_mat(nr,nc);

    for (size_t j=0; j < nc; ++j) {
        for (size_t i=0; i < nr; ++i) {
            ret_mat(i,j) = runif(a_par, b_par, engine);
        }
    }

    return ret_mat;
}

template<typename T>
inline
Mat_t
runif_mat(size_t nr, size_t nc, rand_engine_t& engine)
{
    return runif_mat(nr, nc, T(0), T(1), engine);
}

template<typename T>
inline
Mat_t
runif_mat(size_t nr, size_t nc)
{
    rand_engine_t engine(std::random_device{}());

    return runif_mat(nr, nc, T(0), T(1), engine);
}
