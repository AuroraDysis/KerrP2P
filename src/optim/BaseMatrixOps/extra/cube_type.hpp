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

#if defined(BMO_ENABLE_ARMA_WRAPPERS) || defined(BMO_ENABLE_EIGEN_WRAPPERS)

template<typename T>
class Cube_t
{
    public:
#ifdef BMO_ENABLE_EIGEN_WRAPPERS
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatData_t;
#else
        typedef arma::Mat<T> MatData_t;
#endif
        typedef std::vector< Cube_t<T>::MatData_t > CubeData_t;

        const size_t n_row;  // number of rows in each matrix
        const size_t n_col;  // number of columns in each matrix
        const size_t n_mat;  // number of matrices in each cube

        ~Cube_t() = default;

        Cube_t(const Cube_t<T>& cube_inp);
        Cube_t(Cube_t<T>&& cube_inp);

        explicit Cube_t();
        explicit Cube_t(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp);

        Cube_t<T>& operator=(const Cube_t<T>& cube_inp);
        Cube_t<T>& operator=(Cube_t<T>&& cube_inp);

        T& operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind);
        const T& operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind) const;

        Cube_t<T>::MatData_t& mat(const size_t mat_ind);
        const Cube_t<T>::MatData_t& mat(const size_t mat_ind) const;

        void setZero();
        void setZero(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp);

        CubeData_t& get_raw_data();
        const CubeData_t& get_raw_data() const;

        void reset_raw_data();
        void reset_dims();

    private:
        CubeData_t raw_data_;

        void set_dims(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp);
};

//

template<typename T>
inline
Cube_t<T>::Cube_t()
    : n_row(0), n_col(0), n_mat(0)
{
    raw_data_ = Cube_t<T>::CubeData_t(0, Cube_t<T>::MatData_t(0, 0));
}

template<typename T>
inline
Cube_t<T>::Cube_t(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp)
    : n_row(n_row_inp), n_col(n_col_inp), n_mat(n_mat_inp)
{
    raw_data_ = Cube_t<T>::CubeData_t(n_mat_inp, Cube_t<T>::MatData_t(n_row_inp, n_col_inp));
}

template<typename T>
inline
Cube_t<T>::Cube_t(const Cube_t<T>& cube_inp)
    : n_row(cube_inp.n_row), n_col(cube_inp.n_col), n_mat(cube_inp.n_mat)
{
    raw_data_ = cube_inp.get_raw_data();
}

template<typename T>
inline
Cube_t<T>::Cube_t(Cube_t<T>&& cube_inp)
    : n_row(cube_inp.n_row), n_col(cube_inp.n_col), n_mat(cube_inp.n_mat)
{
    raw_data_ = std::move(cube_inp.get_raw_data());
    cube_inp.reset_dims();
}

//

template<typename T>
inline
Cube_t<T>&
Cube_t<T>::operator=(const Cube_t<T>& cube_inp)
{
    raw_data_ = cube_inp.get_raw_data();

    this->set_dims(cube_inp.n_row, cube_inp.n_col, cube_inp.n_mat);

    return *this;
}

template<typename T>
inline
Cube_t<T>&
Cube_t<T>::operator=(Cube_t<T>&& cube_inp)
{
    raw_data_ = std::move(cube_inp.get_raw_data());

    this->set_dims(cube_inp.n_row, cube_inp.n_col, cube_inp.n_mat);

    cube_inp.reset_dims();

    return *this;
}

// access

template<typename T>
inline
T&
Cube_t<T>::operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind)
{
    return raw_data_[mat_ind](row_ind, col_ind);
}

template<typename T>
inline
const T&
Cube_t<T>::operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind)
const
{
    return raw_data_[mat_ind](row_ind, col_ind);
}

template<typename T>
inline
typename Cube_t<T>::MatData_t&
Cube_t<T>::mat(const size_t mat_ind)
{
    return raw_data_[mat_ind];
}

template<typename T>
inline
const typename Cube_t<T>::MatData_t&
Cube_t<T>::mat(const size_t mat_ind)
const
{
    return raw_data_[mat_ind];
}

// set

template<typename T>
inline
void
Cube_t<T>::setZero()
{
    if (n_mat > size_t(0)) {
        for (size_t mat_ind = 0; mat_ind < n_mat; ++mat_ind) {
            BMO_MATOPS_SET_ZERO(raw_data_[mat_ind]);
        }
    }
}

template<typename T>
inline
void
Cube_t<T>::setZero(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp)
{
    raw_data_.resize(n_mat_inp);

    if (n_mat_inp > size_t(0)) {
        for (size_t mat_ind = 0; mat_ind < n_mat_inp; ++mat_ind) {
            BMO_MATOPS_SET_ZERO_MAT(raw_data_[mat_ind], n_row_inp, n_col_inp);
        }
    }

    set_dims(n_row_inp, n_col_inp, n_mat_inp);
}

//

template<typename T>
inline
typename Cube_t<T>::CubeData_t&
Cube_t<T>::get_raw_data()
{
    return raw_data_;
}

template<typename T>
inline
const typename Cube_t<T>::CubeData_t&
Cube_t<T>::get_raw_data()
const
{
    return raw_data_;
}

template<typename T>
inline
void
Cube_t<T>::reset_raw_data()
{
    raw_data_.clear();
}

template<typename T>
inline
void
Cube_t<T>::reset_dims()
{
    const_cast<size_t&>(n_row) = 0;
    const_cast<size_t&>(n_col) = 0;
    const_cast<size_t&>(n_mat) = 0;
}

template<typename T>
inline
void
Cube_t<T>::set_dims(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp)
{
    const_cast<size_t&>(n_row) = n_row_inp;
    const_cast<size_t&>(n_col) = n_col_inp;
    const_cast<size_t&>(n_mat) = n_mat_inp;
}

#endif
