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
#error Do not include this file.
#endif



#include <memory>  // std::make_unique
#include <utility>  // std::move



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
inline
matrix<Scalar>::matrix(
    const int nrow,
    const int ncol
) noexcept :
    data_(std::make_unique<Scalar[]>(nrow * ncol)),
    impl_(nrow, ncol, data_.get())
{
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(const matrix& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(matrix&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(const matrix<Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(matrix<Scalar&>&& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(const matrix<const Scalar>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(matrix<const Scalar>&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(const matrix<const Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<Scalar>::matrix(matrix<const Scalar&>&& other) noexcept :
    data_(std::move(other.keep_alive_)),
    impl_(std::move(other.impl_))
{
    if (!data_)
    {
        data_ = std::make_unique<Scalar[]>(other.nrow() * other.ncol());
        impl_.reset(other.nrow(), other.ncol(), data_.get());

        impl_.copy(other.impl_);
    }
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(Scalar s) noexcept
{
    impl_.fill(s);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(const matrix& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(matrix&& other) noexcept
{
    impl_.swap(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(const matrix<Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(matrix<Scalar&>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(const matrix<const Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(matrix<const Scalar>&& other) noexcept
{
    impl_.swap(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(const matrix<const Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar>&
matrix<Scalar>::operator =(matrix<const Scalar&>&& other) noexcept
{
    if (other.keep_alive_)
    {
        data_.swap(other.keep_alive_);
        impl_.swap(other.impl_);
    }
    else
    {
        impl_.copy(other.impl_);
    }

    return *this;
}



template<typename Scalar>
inline void
matrix<Scalar>::swap(matrix& other) noexcept
{
    data_.swap(other.data_);
    impl_.swap(other.impl_);
}



template<typename Scalar>
inline Scalar&
matrix<Scalar>::operator ()(
    const int i,
    const int j
) noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<Scalar>::operator ()(
    const int i,
    const int j
) const noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::operator ()(
    const int i,
    const subscript jj
) noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::operator ()(
    const int i,
    const subscript jj
) const noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::operator ()(
    const subscript ii,
    const int j
) noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::operator ()(
    const subscript ii,
    const int j
) const noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar>::operator ()(
    const subscript ii,
    const subscript jj
) noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar>::operator ()(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline int
matrix<Scalar>::nrow() const noexcept
{
    return impl_.nrow();
}



template<typename Scalar>
inline int
matrix<Scalar>::ncol() const noexcept
{
    return impl_.ncol();
}



template<typename Scalar>
inline int
matrix<Scalar>::ldim() const noexcept
{
    return impl_.ldim();
}



template<typename Scalar>
inline Scalar*
matrix<Scalar>::base() noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
matrix<Scalar>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline Scalar*
matrix<Scalar>::ptr(
    const int i,
    const int j
) noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline const Scalar*
matrix<Scalar>::ptr(
    const int i,
    const int j
) const noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline Scalar&
matrix<Scalar>::elem(
    const int i,
    const int j
) noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<Scalar>::elem(
    const int i,
    const int j
) const noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::elem(
    const int i,
    const subscript jj) noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::elem(
    const int i,
    const subscript jj
) const noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::elem(
    const subscript ii,
    const int j
) noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::elem(
    const subscript ii,
    const int j
) const noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar>::elem(
    const subscript ii,
    const subscript jj
) noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar>::elem(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::row(const int i) noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::row(const int i) const noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar>::row(const subscript ii) noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar>::row(const subscript ii) const noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar>::col(const int j) noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar>::col(const int j) const noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar>::col(const subscript jj) noexcept
{
    return impl_.col(jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar>::col(const subscript jj) const noexcept
{
    return impl_.col(jj);
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(
    const int nrow,
    const int ncol,
    Scalar* const base
) noexcept :
    impl_(nrow, ncol, base)
{
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(
    const int nrow,
    const int ncol,
    const int ldim,
    Scalar* const base
) noexcept :
    impl_(nrow, ncol, ldim, base)
{
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(matrix_impl<Scalar>&& impl) noexcept :
    impl_(std::move(impl))
{
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(matrix<Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(matrix& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<Scalar&>::matrix(matrix&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(Scalar s) noexcept
{
    impl_.fill(s);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(const matrix<Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(matrix<Scalar>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(const matrix& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(matrix&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(const matrix<const Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(matrix<const Scalar>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(const matrix<const Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline matrix<Scalar&>&
matrix<Scalar&>::operator =(matrix<const Scalar&>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline void
matrix<Scalar&>::swap(matrix& other) noexcept
{
    impl_.swap(other.impl_);
}



template<typename Scalar>
inline Scalar&
matrix<Scalar&>::operator ()(
    const int i,
    const int j
) noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<Scalar&>::operator ()(
    const int i,
    const int j
) const noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::operator ()(
    const int i,
    const subscript jj
) noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::operator ()(
    const int i,
    const subscript jj
) const noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::operator ()(
    const subscript ii,
    const int j
) noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::operator ()(
    const subscript ii,
    const int j
) const noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar&>::operator ()(
    const subscript ii,
    const subscript jj
) noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar&>::operator ()(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline int
matrix<Scalar&>::nrow() const noexcept
{
    return impl_.nrow();
}



template<typename Scalar>
inline int
matrix<Scalar&>::ncol() const noexcept
{
    return impl_.ncol();
}



template<typename Scalar>
inline int
matrix<Scalar&>::ldim() const noexcept
{
    return impl_.ldim();
}



template<typename Scalar>
inline Scalar*
matrix<Scalar&>::base() noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
matrix<Scalar&>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline Scalar*
matrix<Scalar&>::ptr(
    const int i,
    const int j
) noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline const Scalar*
matrix<Scalar&>::ptr(
    const int i,
    const int j
) const noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline Scalar&
matrix<Scalar&>::elem(
    const int i,
    const int j
) noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<Scalar&>::elem(
    const int i,
    const int j
) const noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::elem(
    const int i,
    const subscript jj
) noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::elem(
    const int i,
    const subscript jj
) const noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::elem(
    const subscript ii,
    const int j
) noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::elem(
    const subscript ii,
    const int j
) const noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar&>::elem(
    const subscript ii,
    const subscript jj
) noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar&>::elem(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::row(const int i) noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::row(const int i) const noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar&>::row(const subscript ii) noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar&>::row(const subscript ii) const noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline vector<Scalar&>
matrix<Scalar&>::col(const int j) noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<Scalar&>::col(const int j) const noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline matrix<Scalar&>
matrix<Scalar&>::col(const subscript jj) noexcept
{
    return impl_.col(jj);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<Scalar&>::col(const subscript jj) const noexcept
{
    return impl_.col(jj);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(
    const int nrow,
    const int ncol
) noexcept :
    data_(std::make_unique<Scalar[]>(nrow * ncol)),
    impl_(nrow, ncol, data_.get())
{
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(const matrix<Scalar>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(matrix<Scalar>&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(const matrix<Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(matrix<Scalar&>&& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(const matrix& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(matrix&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(const matrix<const Scalar&>& other)
    noexcept :
    data_(std::make_unique<Scalar[]>(other.nrow() * other.ncol())),
    impl_(other.nrow(), other.ncol(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
matrix<const Scalar>::matrix(matrix<const Scalar&>&& other) noexcept :
    data_(std::move(other.keep_alive_)),
    impl_(std::move(other.impl_))
{
    if (!data_)
    {
        data_ = std::make_unique<Scalar[]>(other.nrow() * other.ncol());
        impl_.reset(other.nrow(), other.ncol(), data_.get());

        impl_.copy(other.impl_);
    }
}



template<typename Scalar>
inline Scalar const&
matrix<const Scalar>::operator ()(
    const int i,
    const int j
) const noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::operator ()(
    const int i,
    const subscript jj
) const noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::operator ()(
    const subscript ii,
    const int j
) const noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar>::operator ()(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline int
matrix<const Scalar>::nrow() const noexcept
{
    return impl_.nrow();
}



template<typename Scalar>
inline int
matrix<const Scalar>::ncol() const noexcept
{
    return impl_.ncol();
}



template<typename Scalar>
inline int
matrix<const Scalar>::ldim() const noexcept
{
    return impl_.ldim();
}



template<typename Scalar>
inline const Scalar*
matrix<const Scalar>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
matrix<const Scalar>::ptr(
    const int i,
    const int j
) const noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<const Scalar>::elem(
    const int i,
    const int j
) const noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::elem(
    const int i,
    const subscript jj
) const noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::elem(
    const subscript ii,
    const int j
) const noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar>::elem(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::row(const int i) const noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar>::row(const subscript ii) const noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar>::col(const int j) const noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar>::col(const subscript jj) const noexcept
{
    return impl_.col(jj);
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(
    const int nrow,
    const int ncol,
    Scalar* const base
) noexcept :
    impl_(nrow, ncol, base)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(
    const int nrow,
    const int ncol,
    const int ldim,
    Scalar* const base
) noexcept :
    impl_(nrow, ncol, ldim, base)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(matrix_impl<Scalar>&& impl) noexcept :
    impl_(std::move(impl))
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(const matrix<Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(matrix<Scalar>&& other) noexcept :
    impl_(other.impl_),
    keep_alive_(std::move(other.data_))
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(const matrix<Scalar&>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(matrix<Scalar&>&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(const matrix<const Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(matrix<const Scalar>&& other) noexcept :
    impl_(other.impl_),
    keep_alive_(std::move(other.data_))
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(const matrix& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
matrix<const Scalar&>::matrix(matrix&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline Scalar const&
matrix<const Scalar&>::operator ()(
    const int i,
    const int j
) const noexcept
{
    return elem(i, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::operator ()(
    const int i,
    const subscript jj
) const noexcept
{
    return elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::operator ()(
    const subscript ii,
    const int j
) const noexcept
{
    return elem(ii, j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar&>::operator ()(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return elem(ii, jj);
}



template<typename Scalar>
inline int
matrix<const Scalar&>::nrow() const noexcept
{
    return impl_.nrow();
}



template<typename Scalar>
inline int
matrix<const Scalar&>::ncol() const noexcept
{
    return impl_.ncol();
}



template<typename Scalar>
inline int
matrix<const Scalar&>::ldim() const noexcept
{
    return impl_.ldim();
}



template<typename Scalar>
inline const Scalar*
matrix<const Scalar&>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
matrix<const Scalar&>::ptr(
    const int i,
    const int j
) const noexcept
{
    return impl_.ptr(i, j);
}



template<typename Scalar>
inline const Scalar&
matrix<const Scalar&>::elem(
    const int i,
    const int j
) const noexcept
{
    return impl_.elem(i, j);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::elem(
    const int i,
    const subscript jj
) const noexcept
{
    return impl_.elem(i, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::elem(
    const subscript ii,
    const int j
) const noexcept
{
    return impl_.elem(ii, j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar&>::elem(
    const subscript ii,
    const subscript jj
) const noexcept
{
    return impl_.elem(ii, jj);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::row(const int i) const noexcept
{
    return impl_.row(i);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar&>::row(const subscript ii) const noexcept
{
    return impl_.row(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
matrix<const Scalar&>::col(const int j) const noexcept
{
    return impl_.col(j);
}



template<typename Scalar>
inline matrix<const Scalar&>
matrix<const Scalar&>::col(const subscript jj) const noexcept
{
    return impl_.col(jj);
}



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs
