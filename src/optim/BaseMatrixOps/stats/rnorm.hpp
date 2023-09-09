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
 * Sample from a normal distribution
 */

//
// single draw

namespace internal
{

template<typename T>
constexpr
bool
norm_sanity_check(const T mu_par, const T sigma_par)
noexcept
{
    return( any_nan(mu_par,sigma_par) ? \
                false :
            //
            sigma_par < T(0) ? \
                false :
            //
                true );
}

template<typename T>
inline
T
rnorm_compute(const T mu_par, const T sigma_par, rand_engine_t& engine)
{
    if (!norm_sanity_check(mu_par,sigma_par)) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    //

    std::normal_distribution<T> norm_dist(T(0),T(1));

    return mu_par + sigma_par*norm_dist(engine);
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
inline
TC
rnorm_type_check(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    return rnorm_compute(static_cast<TC>(mu_par),static_cast<TC>(sigma_par),engine);
}

}

template<typename T1, typename T2>
inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    return internal::rnorm_type_check(mu_par,sigma_par,engine);
}

template<typename T1, typename T2>
inline
common_return_t<T1,T2>
rnorm(const T1 mu_par, const T2 sigma_par, const size_t seed_val)
{
    rand_engine_t engine(seed_val);
    return rnorm(mu_par,sigma_par,engine);
}

template<typename T>
inline
T
rnorm(rand_engine_t& engine)
{
    return rnorm(T(0),T(1),engine);
}

template<typename T>
inline
T
rnorm()
{
    return rnorm(T(0),T(1));
}

//

namespace internal
{

template<typename T1, typename T2, typename vT>
inline
void
rnorm_vec_inplace(size_t n_vals, const T1 mu_par, const T2 sigma_par, rand_engine_t& engine, vT& vec_out)
{
    for (size_t i=0; i < n_vals; ++i) {
        vec_out(i) = rnorm(mu_par, sigma_par, engine);
    }
}

template<typename fT, typename vT>
inline
void
rnorm_vec_inplace(size_t n_vals, rand_engine_t& engine, vT& vec_out)
{
    for (size_t i=0; i < n_vals; ++i) {
        vec_out(i) = rnorm(fT(0), fT(1), engine);
    }
}

}

template<typename T1, typename T2, typename vT = ColVec_t>
inline
vT
rnorm_vec(size_t n_vals, const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    vT ret_vec(n_vals);

    internal::rnorm_vec_inplace(n_vals, mu_par, sigma_par, engine, ret_vec);

    return ret_vec;
}

template<typename T, typename vT = ColVec_t>
inline
vT
rnorm_vec(size_t n_vals, rand_engine_t& engine)
{
    return rnorm_vec(n_vals, T(0), T(1), engine);
}

template<typename T, typename vT = ColVec_t>
inline
vT
rnorm_vec(size_t n_vals)
{
    rand_engine_t engine(std::random_device{}());
    
    return rnorm_vec(n_vals, T(0), T(1), engine);
}

//

template<typename T1, typename T2>
inline
Mat_t
rnorm_mat(size_t nr, size_t nc, const T1 mu_par, const T2 sigma_par, rand_engine_t& engine)
{
    Mat_t ret_mat(nr,nc);

    for (size_t j=0; j < nc; ++j) {
        for (size_t i=0; i < nr; ++i) {
            ret_mat(i,j) = rnorm(mu_par, sigma_par, engine);
        }
    }

    return ret_mat;
}

template<typename T>
inline
Mat_t
rsnorm_mat(size_t nr, size_t nc, rand_engine_t& engine)
{
    return rnorm_mat(nr, nc, T(0), T(1), engine);
}

template<typename T>
inline
Mat_t
rsnorm_mat(size_t nr, size_t nc)
{
    rand_engine_t engine(std::random_device{}());
    
    return rnorm_mat(nr, nc, T(0), T(1), engine);
}
