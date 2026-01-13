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



#ifndef MCS_CORE_DETAIL_VECTOR_HH
#error Do not include this file.
#endif



#include <memory>  // std::make_unique
#include <utility>  // std::move



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
inline
vector<Scalar>::vector(const int len) noexcept :
    data_(std::make_unique<Scalar[]>(len)),
    impl_(len, data_.get())
{
}



template<typename Scalar>
inline
vector<Scalar>::vector(const vector& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<Scalar>::vector(vector&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
vector<Scalar>::vector(const vector<Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<Scalar>::vector(vector<Scalar&>&& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<Scalar>::vector(const vector<const Scalar>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<Scalar>::vector(vector<const Scalar>&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
vector<Scalar>::vector(const vector<const Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<Scalar>::vector(vector<const Scalar&>&& other) noexcept :
    data_(std::move(other.keep_alive_)),
    impl_(std::move(other.impl_))
{
    if (!data_)
    {
        data_ = std::make_unique<Scalar[]>(other.len());
        impl_.reset(other.len(), data_.get());

        impl_.copy(other.impl_);
    }
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(Scalar s) noexcept
{
    impl_.fill(s);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(const vector& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(vector&& other) noexcept
{
    impl_.swap(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(const vector<Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(vector<Scalar&>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(const vector<const Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(vector<const Scalar>&& other) noexcept
{
    impl_.swap(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(const vector<const Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar>&
vector<Scalar>::operator =(vector<const Scalar&>&& other) noexcept
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
vector<Scalar>::swap(vector& other) noexcept
{
    data_.swap(other.data_);
    impl_.swap(other.impl_);
}



template<typename Scalar>
inline Scalar&
vector<Scalar>::operator ()(int const i) noexcept
{
    return elem(i);
}



template<typename Scalar>
inline const Scalar&
vector<Scalar>::operator ()(const int i) const noexcept
{
    return elem(i);
}



template<typename Scalar>
inline vector<Scalar&>
vector<Scalar>::operator ()(const subscript ii) noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<Scalar>::operator ()(const subscript ii) const noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline int
vector<Scalar>::len() const noexcept
{
    return impl_.len();
}



template<typename Scalar>
inline int
vector<Scalar>::inc() const noexcept
{
    return impl_.inc();
}



template<typename Scalar>
inline Scalar*
vector<Scalar>::base() noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
vector<Scalar>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline Scalar*
vector<Scalar>::ptr(const int i) noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline const Scalar*
vector<Scalar>::ptr(const int i) const noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline Scalar&
vector<Scalar>::elem(const int i) noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline const Scalar&
vector<Scalar>::elem(const int i) const noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline vector<Scalar&>
vector<Scalar>::elem(const subscript ii) noexcept
{
    return impl_.elem(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<Scalar>::elem(const subscript ii) const noexcept
{
    return impl_.elem(ii);
}



template<typename Scalar>
inline
vector<Scalar&>::vector(
    const int len,
    Scalar* const base
) noexcept :
    impl_(len, base)
{
}



template<typename Scalar>
inline
vector<Scalar&>::vector(
    const int len,
    const int inc,
    Scalar* const base
) noexcept :
    impl_(len, inc, base)
{
}



template<typename Scalar>
inline
vector<Scalar&>::vector(vector_impl<Scalar>&& impl) noexcept :
    impl_(std::move(impl))
{
}



template<typename Scalar>
inline
vector<Scalar&>::vector(vector<Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<Scalar&>::vector(vector& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<Scalar&>::vector(vector&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(Scalar s) noexcept
{
    impl_.fill(s);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(const vector<Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(vector<Scalar>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(const vector& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(vector&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(const vector<const Scalar>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(vector<const Scalar>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(const vector<const Scalar&>& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline vector<Scalar&>&
vector<Scalar&>::operator =(vector<const Scalar&>&& other) noexcept
{
    impl_.copy(other.impl_);

    return *this;
}



template<typename Scalar>
inline void
vector<Scalar&>::swap(vector& other) noexcept
{
    impl_.swap(other.impl_);
}



template<typename Scalar>
inline Scalar&
vector<Scalar&>::operator ()(const int i) noexcept
{
    return elem(i);
}



template<typename Scalar>
inline const Scalar&
vector<Scalar&>::operator ()(const int i) const noexcept
{
    return elem(i);
}



template<typename Scalar>
inline vector<Scalar&>
vector<Scalar&>::operator ()(const subscript ii) noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<Scalar&>::operator ()(const subscript ii) const noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline int
vector<Scalar&>::len() const noexcept
{
    return impl_.len();
}



template<typename Scalar>
inline int
vector<Scalar&>::inc() const noexcept
{
    return impl_.inc();
}



template<typename Scalar>
inline Scalar*
vector<Scalar&>::base() noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
vector<Scalar&>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline Scalar*
vector<Scalar&>::ptr(const int i) noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline const Scalar*
vector<Scalar&>::ptr(const int i) const noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline Scalar&
vector<Scalar&>::elem(const int i) noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline const Scalar&
vector<Scalar&>::elem(const int i) const noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline vector<Scalar&>
vector<Scalar&>::elem(const subscript ii) noexcept
{
    return impl_.elem(ii);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<Scalar&>::elem(subscript const ii) const noexcept
{
    return impl_.elem(ii);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(const int len) noexcept :
    data_(std::make_unique<Scalar[]>(len)),
    impl_(len, data_.get())
{
}



template<typename Scalar>
inline
vector<const Scalar>::vector(const vector<Scalar>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(vector<Scalar>&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
vector<const Scalar>::vector(const vector<Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(vector<Scalar&>&& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(const vector<const Scalar>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(vector<const Scalar>&& other) noexcept :
    data_(std::move(other.data_)),
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
vector<const Scalar>::vector(const vector<const Scalar&>& other) noexcept :
    data_(std::make_unique<Scalar[]>(other.len())),
    impl_(other.len(), data_.get())
{
    impl_.copy(other.impl_);
}



template<typename Scalar>
inline
vector<const Scalar>::vector(vector<const Scalar&>&& other) noexcept :
    data_(std::move(other.keep_alive_)),
    impl_(std::move(other.impl_))
{
    if (!data_)
    {
        data_ = std::make_unique<Scalar[]>(other.len());
        impl_.reset(other.len(), data_.get());

        impl_.copy(other.impl_);
    }
}



template<typename Scalar>
inline const Scalar&
vector<const Scalar>::operator ()(const int i) const noexcept
{
    return elem(i);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<const Scalar>::operator ()(const subscript ii) const noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline int
vector<const Scalar>::len() const noexcept
{
    return impl_.len();
}



template<typename Scalar>
inline int
vector<const Scalar>::inc() const noexcept
{
    return impl_.inc();
}



template<typename Scalar>
inline const Scalar*
vector<const Scalar>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
vector<const Scalar>::ptr(const int i) const noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline const Scalar&
vector<const Scalar>::elem(const int i) const noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<const Scalar>::elem(const subscript ii) const noexcept
{
    return impl_.elem(ii);
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(
    const int len,
    const Scalar* const base
) noexcept :
    impl_(len, base)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(
    const int len,
    const int inc,
    const Scalar* const base
) noexcept :
    impl_(len, inc, base)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(vector_impl<Scalar>&& impl) noexcept :
    impl_(std::move(impl))
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(const vector<Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(vector<Scalar>&& other) noexcept :
    impl_(other.impl_),
    keep_alive_(std::move(other.data_))
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(const vector<Scalar&>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(vector<Scalar&>&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(const vector<const Scalar>& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(vector<const Scalar>&& other) noexcept :
    impl_(other.impl_),
    keep_alive_(std::move(other.data_))
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(const vector& other) noexcept :
    impl_(other.impl_)
{
}



template<typename Scalar>
inline
vector<const Scalar&>::vector(vector&& other) noexcept :
    impl_(std::move(other.impl_))
{
}



template<typename Scalar>
inline const Scalar&
vector<const Scalar&>::operator ()(const int i) const noexcept
{
    return elem(i);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<const Scalar&>::operator ()(const subscript ii) const noexcept
{
    return elem(ii);
}



template<typename Scalar>
inline int
vector<const Scalar&>::len() const noexcept
{
    return impl_.len();
}



template<typename Scalar>
inline int
vector<const Scalar&>::inc() const noexcept
{
    return impl_.inc();
}



template<typename Scalar>
inline const Scalar*
vector<const Scalar&>::base() const noexcept
{
    return impl_.base();
}



template<typename Scalar>
inline const Scalar*
vector<const Scalar&>::ptr(const int i) const noexcept
{
    return impl_.ptr(i);
}



template<typename Scalar>
inline const Scalar&
vector<const Scalar&>::elem(const int i) const noexcept
{
    return impl_.elem(i);
}



template<typename Scalar>
inline vector<const Scalar&>
vector<const Scalar&>::elem(const subscript ii) const noexcept
{
    return impl_.elem(ii);
}



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs
