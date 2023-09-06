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
class Tess_t
{
    public:
        typedef std::vector< Cube_t<T> > TessData_t;

        const size_t n_row;  // number of rows in each matrix
        const size_t n_col;  // number of columns in each matrix
        const size_t n_mat;  // number of matrices in each cube
        const size_t n_cube; // number of cubes in the tesseract

        ~Tess_t() = default;

        Tess_t(const Tess_t<T>& tess_inp);
        Tess_t(Tess_t<T>&& tess_inp);

        explicit Tess_t();
        explicit Tess_t(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp);

        Tess_t<T>& operator=(const Tess_t<T>& tess_inp);
        Tess_t<T>& operator=(Tess_t<T>&& tess_inp);

        T& operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind, const size_t cube_ind);
        const T& operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind, const size_t cube_ind) const;

        Cube_t<T>& cube(const size_t cube_ind);
        const Cube_t<T>& cube(const size_t cube_ind) const;

        void setZero();
        void setZero(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp);

        TessData_t& get_raw_data();
        const TessData_t& get_raw_data() const;

        void reset_raw_data();
        void reset_dims();

    private:
        TessData_t raw_data_;

        void set_dims(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp);
};

//

template<typename T>
inline
Tess_t<T>::Tess_t()
    : n_row(0), n_col(0), n_mat(0), n_cube(0)
{
    raw_data_ = Tess_t<T>::TessData_t(0, Cube_t<T>(0, 0, 0));
}

template<typename T>
inline
Tess_t<T>::Tess_t(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp)
    : n_row(n_row_inp), n_col(n_col_inp), n_mat(n_mat_inp), n_cube(n_cube_inp)
{
    raw_data_ = Tess_t<T>::TessData_t(n_cube_inp, Cube_t<T>(n_row_inp, n_col_inp, n_mat_inp));
}

template<typename T>
inline
Tess_t<T>::Tess_t(const Tess_t<T>& tess_inp)
    : n_row(tess_inp.n_row), n_col(tess_inp.n_col), n_mat(tess_inp.n_mat), n_cube(tess_inp.n_cube)
{
    raw_data_ = tess_inp.get_raw_data();
}

template<typename T>
inline
Tess_t<T>::Tess_t(Tess_t<T>&& tess_inp)
    : n_row(tess_inp.n_row), n_col(tess_inp.n_col), n_mat(tess_inp.n_mat), n_cube(tess_inp.n_cube)
{
    raw_data_ = std::move(tess_inp.get_raw_data());
    tess_inp.reset_dims();
}

//

template<typename T>
inline
Tess_t<T>&
Tess_t<T>::operator=(const Tess_t<T>& tess_inp)
{
    raw_data_ = tess_inp.get_raw_data();

    this->set_dims(tess_inp.n_row, tess_inp.n_col, tess_inp.n_mat, tess_inp.n_cube);

    return *this;
}

template<typename T>
inline
Tess_t<T>&
Tess_t<T>::operator=(Tess_t<T>&& tess_inp)
{
    raw_data_ = std::move(tess_inp.get_raw_data());

    this->set_dims(tess_inp.n_row, tess_inp.n_col, tess_inp.n_mat, tess_inp.n_cube);

    tess_inp.reset_dims();

    return *this;
}

// access

template<typename T>
inline
T&
Tess_t<T>::operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind, const size_t cube_ind)
{
    return raw_data_[cube_ind](row_ind, col_ind, mat_ind);
}

template<typename T>
inline
const T&
Tess_t<T>::operator()(const size_t row_ind, const size_t col_ind, const size_t mat_ind, const size_t cube_ind)
const
{
    return &raw_data_[cube_ind](row_ind, col_ind, mat_ind);
}

template<typename T>
inline
Cube_t<T>&
Tess_t<T>::cube(const size_t cube_ind)
{
    return raw_data_[cube_ind];
}

template<typename T>
inline
const Cube_t<T>&
Tess_t<T>::cube(const size_t cube_ind)
const
{
    return raw_data_[cube_ind];
}

// set

template<typename T>
inline
void
Tess_t<T>::setZero()
{
    if (n_cube > size_t(0)) {
        for (size_t cube_ind = 0; cube_ind < n_cube; ++cube_ind) {
            raw_data_[cube_ind].setZero();
        }
    }
}

template<typename T>
inline
void
Tess_t<T>::setZero(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp)
{
    raw_data_.resize(n_cube_inp);

    if (n_cube_inp > size_t(0)) {
        for (size_t cube_ind = 0; cube_ind < n_cube_inp; ++cube_ind) {
            raw_data_[cube_ind].setZero(n_row_inp, n_col_inp, n_mat_inp);
        }
    }

    set_dims(n_row_inp, n_col_inp, n_mat_inp, n_cube_inp);
}

//

template<typename T>
inline
typename Tess_t<T>::TessData_t&
Tess_t<T>::get_raw_data()
{
    return raw_data_;
}

template<typename T>
inline
const typename Tess_t<T>::TessData_t&
Tess_t<T>::get_raw_data()
const
{
    return raw_data_;
}

template<typename T>
inline
void
Tess_t<T>::reset_raw_data()
{
    raw_data_.clear();
}

template<typename T>
inline
void
Tess_t<T>::reset_dims()
{
    const_cast<size_t&>(n_row)  = 0;
    const_cast<size_t&>(n_col)  = 0;
    const_cast<size_t&>(n_mat)  = 0;
    const_cast<size_t&>(n_cube) = 0;
}

template<typename T>
inline
void
Tess_t<T>::set_dims(const size_t n_row_inp, const size_t n_col_inp, const size_t n_mat_inp, const size_t n_cube_inp)
{
    const_cast<size_t&>(n_row)  = n_row_inp;
    const_cast<size_t&>(n_col)  = n_col_inp;
    const_cast<size_t&>(n_mat)  = n_mat_inp;
    const_cast<size_t&>(n_cube) = n_cube_inp;
}

#endif
