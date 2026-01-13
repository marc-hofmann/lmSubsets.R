// Copyright  2009-2020  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs.cc/>).
//
// 'mcs' is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'mcs' is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'mcs'.  If not, see <http://www.gnu.org/licenses/>.



#ifndef MCS_CORE_DETAIL_MATRIX_HH
#define MCS_CORE_DETAIL_MATRIX_HH



#include <memory>  // std::unique_ptr



#include "mcs/core/detail/matrix_impl.hh"
#include "mcs/core/detail/vector.hh"



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
class matrix
{

    friend void
    swap(matrix& a, matrix& b) noexcept
    {
        a.swap(b);
    }

    friend matrix<Scalar&>;

    friend matrix<const Scalar>;

    friend matrix<const Scalar&>;



private:

    std::unique_ptr<Scalar[]> data_;

    matrix_impl<Scalar> impl_;



public:

    matrix() = delete;

    matrix(int nrow, int ncol) noexcept;

    matrix(const matrix& other) noexcept;

    matrix(matrix&& other) noexcept;

    matrix(const matrix<Scalar&>& other) noexcept;

    matrix(matrix<Scalar&>&& other) noexcept;

    matrix(const matrix<const Scalar>& other) noexcept;

    matrix(matrix<const Scalar>&& other) noexcept;

    matrix(const matrix<const Scalar&>& other) noexcept;

    matrix(matrix<const Scalar&>&& other) noexcept;

    ~matrix() noexcept = default;



public:

    matrix&
    operator =(Scalar s) noexcept;

    matrix&
    operator =(const matrix& other) noexcept;

    matrix&
    operator =(matrix&& other) noexcept;

    matrix&
    operator =(const matrix<Scalar&>& other) noexcept;

    matrix&
    operator =(matrix<Scalar&>&& other) noexcept;

    matrix&
    operator =(const matrix<const Scalar>& other) noexcept;

    matrix&
    operator =(matrix<const Scalar>&& other) noexcept;

    matrix&
    operator =(const matrix<const Scalar&>& other) noexcept;

    matrix&
    operator =(matrix<const Scalar&>&& other) noexcept;



public:

    void
    swap(matrix& other) noexcept;



public:

    Scalar&
    operator ()(int i, int j) noexcept;

    const Scalar&
    operator ()(int i, int j) const noexcept;

    vector<Scalar&>
    operator ()(int i, subscript jj) noexcept;

    vector<const Scalar&>
    operator ()(int i, subscript jj) const noexcept;

    vector<Scalar&>
    operator ()(subscript ii, int j) noexcept;

    vector<const Scalar&>
    operator ()(subscript ii, int j) const noexcept;

    matrix<Scalar&>
    operator ()(subscript ii, subscript jj) noexcept;

    matrix<const Scalar&>
    operator ()(subscript ii, subscript jj) const noexcept;



public:

    int
    nrow() const noexcept;

    int
    ncol() const noexcept;

    int
    ldim() const noexcept;

    Scalar*
    base() noexcept;

    const Scalar*
    base() const noexcept;

    Scalar*
    ptr(int i, int j) noexcept;

    const Scalar*
    ptr(int i, int j) const noexcept;

    Scalar&
    elem(int i, int j) noexcept;

    const Scalar&
    elem(int i, int j) const noexcept;

    vector<Scalar&>
    elem(int i, subscript jj) noexcept;

    vector<const Scalar&>
    elem(int i, subscript jj) const noexcept;

    vector<Scalar&>
    elem(subscript ii, int j) noexcept;

    vector<const Scalar&>
    elem(subscript ii, int j) const noexcept;

    matrix<Scalar&>
    elem(subscript ii, subscript jj) noexcept;

    matrix<const Scalar&>
    elem(subscript ii, subscript jj) const noexcept;

    vector<Scalar&>
    row(int i) noexcept;

    vector<const Scalar&>
    row(int i) const noexcept;

    matrix<Scalar&>
    row(subscript ii) noexcept;

    matrix<const Scalar&>
    row(subscript ii) const noexcept;

    vector<Scalar&>
    col(int j) noexcept;

    vector<const Scalar&>
    col(int j) const noexcept;

    matrix<Scalar&>
    col(subscript jj) noexcept;

    matrix<const Scalar&>
    col(subscript jj) const noexcept;

};



template<typename Scalar>
class matrix<Scalar&>
{

    friend void
    swap(matrix& a, matrix& b) noexcept
    {
        a.swap(b);
    }



    friend matrix<Scalar>;

    friend matrix<const Scalar>;

    friend matrix<const Scalar&>;



private:

    matrix_impl<Scalar> impl_;



public:

    matrix() = delete;

    matrix(int nrow, int ncol, Scalar* base) noexcept;

    matrix(int nrow, int ncol, int ldim, Scalar* base) noexcept;

    matrix(matrix_impl<Scalar>&& impl) noexcept;

    matrix(matrix<Scalar>& other) noexcept;

    matrix(const matrix<Scalar>& other) = delete;

    matrix(matrix<Scalar>&& other) = delete;

    matrix(matrix& other) noexcept;

    matrix(const matrix& other) = delete;

    matrix(matrix&& other) noexcept;

    matrix(const matrix<const Scalar>& other) = delete;

    matrix(matrix<const Scalar>&& other) = delete;

    matrix(const matrix<const Scalar&>& other) = delete;

    matrix(matrix<const Scalar&>&& other) = delete;

    ~matrix() noexcept = default;



public:

    matrix&
    operator =(Scalar s) noexcept;

    matrix&
    operator =(const matrix<Scalar>& other) noexcept;

    matrix&
    operator =(matrix<Scalar>&& other) noexcept;

    matrix&
    operator =(const matrix& other) noexcept;

    matrix&
    operator =(matrix&& other) noexcept;

    matrix&
    operator =(const matrix<const Scalar>& other) noexcept;

    matrix&
    operator =(matrix<const Scalar>&& other) noexcept;

    matrix&
    operator =(const matrix<const Scalar&>& other) noexcept;

    matrix&
    operator =(matrix<const Scalar&>&& other) noexcept;



public:

    void
    swap(matrix& other) noexcept;



public:

    Scalar&
    operator ()(int i, int j) noexcept;

    const Scalar&
    operator ()(int i, int j) const noexcept;

    vector<Scalar&>
    operator ()(int i , subscript jj) noexcept;

    vector<const Scalar&>
    operator ()(int i , subscript jj) const noexcept;

    vector<Scalar&>
    operator ()(subscript ii, int j) noexcept;

    vector<const Scalar&>
    operator ()(subscript ii, int j) const noexcept;

    matrix
    operator ()(subscript ii, subscript jj) noexcept;

    matrix<const Scalar&>
    operator ()(subscript ii, subscript jj) const noexcept;



public:

    int
    nrow() const noexcept;

    int
    ncol() const noexcept;

    int
    ldim() const noexcept;

    Scalar*
    base() noexcept;

    const Scalar*
    base() const noexcept;

    Scalar*
    ptr(int i, int j) noexcept;

    const Scalar*
    ptr(int i, int j) const noexcept;

    Scalar&
    elem(int i, int j) noexcept;

    const Scalar&
    elem(int i, int j) const noexcept;

    vector<Scalar&>
    elem(int i, subscript jj) noexcept;

    vector<const Scalar&>
    elem(int i , subscript jj) const noexcept;

    vector<Scalar&>
    elem(subscript ii, int j) noexcept;

    vector<const Scalar&>
    elem(subscript ii, int j) const noexcept;

    matrix
    elem(subscript ii, subscript jj) noexcept;

    matrix<const Scalar&>
    elem(subscript ii, subscript jj) const noexcept;

    vector<Scalar&>
    row(int i) noexcept;

    vector<const Scalar&>
    row(int i) const noexcept;

    matrix
    row(subscript ii) noexcept;

    matrix<const Scalar&>
    row(subscript ii) const noexcept;

    vector<Scalar&>
    col(int j) noexcept;

    vector<const Scalar&>
    col(int j) const noexcept;

    matrix
    col(subscript jj) noexcept;

    matrix<const Scalar&>
    col(subscript jj) const noexcept;

};



template<typename Scalar>
class matrix<const Scalar>
{

    friend matrix<Scalar>;

    friend matrix<Scalar&>;

    friend matrix<const Scalar&>;



private:

    std::unique_ptr<Scalar[]> data_;

    matrix_impl<Scalar> impl_;



public:

    matrix() = delete;

    matrix(int nrow, int ncol) noexcept;

    matrix(const matrix<Scalar>& other) noexcept;

    matrix(matrix<Scalar>&& other) noexcept;

    matrix(const matrix<Scalar&>& other) noexcept;

    matrix(matrix<Scalar&>&& other) noexcept;

    matrix(const matrix& other) noexcept;

    matrix(matrix&& other) noexcept;

    matrix(const matrix<const Scalar&>& other) noexcept;

    matrix(matrix<const Scalar&>&& other) noexcept;

    ~matrix() noexcept = default;



public:

    matrix&
    operator =(const matrix<Scalar>& other) = delete;

    matrix&
    operator =(matrix<Scalar>&& other) = delete;

    matrix&
    operator =(const matrix<Scalar&>& other) = delete;

    matrix&
    operator =(matrix<Scalar&>&& other) = delete;

    matrix&
    operator =(const matrix& other) = delete;

    matrix&
    operator =(matrix&& other) = delete;

    matrix&
    operator =(const matrix<const Scalar&>& other) = delete;

    matrix&
    operator =(matrix<const Scalar&>&& other) = delete;



public:

    void
    swap(matrix& other) = delete;



public:

    Scalar const&
    operator ()(int i, int j) const noexcept;

    vector<const Scalar&>
    operator ()(int i, subscript jj) const noexcept;

    vector<const Scalar&>
    operator ()(subscript ii, int j) const noexcept;

    matrix<const Scalar&>
    operator ()(subscript ii, subscript jj) const noexcept;



public:

    int
    nrow() const noexcept;

    int
    ncol() const noexcept;

    int
    ldim() const noexcept;

    const Scalar*
    base() const noexcept;

    const Scalar*
    ptr(int i, int j) const noexcept;

    const Scalar&
    elem(int i, int j) const noexcept;

    vector<const Scalar&>
    elem(int i, subscript jj) const noexcept;

    vector<const Scalar&>
    elem(subscript ii, int j) const noexcept;

    matrix<const Scalar&>
    elem(subscript ii, subscript jj) const noexcept;

    vector<const Scalar&>
    row(int i) const noexcept;

    matrix<const Scalar&>
    row(subscript ii) const noexcept;

    vector<const Scalar&>
    col(int j) const noexcept;

    matrix<const Scalar&>
    col(subscript jj) const noexcept;

};



template<typename Scalar>
class matrix<const Scalar&>
{

    friend matrix<Scalar>;

    friend matrix<Scalar&>;

    friend matrix<const Scalar>;



private:

    matrix_impl<Scalar> impl_;

    std::unique_ptr<Scalar[]> keep_alive_;



public:

    matrix() = delete;

    matrix(int nrow, int ncol, Scalar* base) noexcept;

    matrix(int nrow, int ncol, int ldim, Scalar* base) noexcept;

    matrix(matrix_impl<Scalar>&& impl) noexcept;

    matrix(const matrix<Scalar>& other) noexcept;

    matrix(matrix<Scalar>&& other) noexcept;

    matrix(const matrix<Scalar&>& other) noexcept;

    matrix(matrix<Scalar&>&& other) noexcept;

    matrix(const matrix<const Scalar>& other) noexcept;

    matrix(matrix<const Scalar>&& other) noexcept;

    matrix(const matrix& other) noexcept;

    matrix(matrix&& other) noexcept;

    ~matrix() noexcept = default;



public:

    matrix&
    operator =(const matrix<Scalar>& other) = delete;

    matrix&
    operator =(matrix<Scalar>&& other) = delete;

    matrix&
    operator =(const matrix<Scalar&>& other) = delete;

    matrix&
    operator =(matrix<Scalar&>&& other) = delete;

    matrix&
    operator =(const matrix<const Scalar>& other) = delete;

    matrix&
    operator =(matrix<const Scalar>&& other) = delete;

    matrix&
    operator =(const matrix& other) = delete;

    matrix&
    operator =(matrix&& other) = delete;



public:

    void
    swap(matrix& other) = delete;



public:

    Scalar const&
    operator ()(int i, int j) const noexcept;

    vector<const Scalar&>
    operator ()(int i, subscript jj) const noexcept;

    vector<const Scalar&>
    operator ()(subscript ii, int j) const noexcept;

    matrix
    operator ()(subscript ii, subscript jj) const noexcept;



public:

    int
    nrow() const noexcept;

    int
    ncol() const noexcept;

    int
    ldim() const noexcept;

    const Scalar*
    base() const noexcept;

    const Scalar*
    ptr(int i, int j) const noexcept;

    const Scalar&
    elem(int i, int j) const noexcept;

    vector<const Scalar&>
    elem(int i, subscript jj) const noexcept;

    vector<const Scalar&>
    elem(subscript ii, int j) const noexcept;

    matrix
    elem(subscript ii, subscript jj) const noexcept;

    vector<const Scalar&>
    row(int i) const noexcept;

    matrix
    row(subscript ii) const noexcept;

    vector<const Scalar&>
    col(int j) const noexcept;

    matrix
    col(subscript jj) const noexcept;

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#include "mcs/core/detail/matrix.inc"
#endif
