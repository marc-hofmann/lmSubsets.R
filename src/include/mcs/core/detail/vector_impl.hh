// Copyright 2018  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs/>).
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



#ifndef MCS_CORE_DETAIL_VECTOR_IMPL_HH
#define MCS_CORE_DETAIL_VECTOR_IMPL_HH



#include <utility>  // std::swap



#include "mcs/core/detail/subscript.hh"



namespace mcs    {
namespace core   {
namespace detail {




template<typename Scalar>
class vector_impl
{

private:

    int len_;

    int inc_;

    Scalar* base_;



public:

    vector_impl() noexcept :
        vector_impl(0, 0, nullptr)
    {
    }



    vector_impl(
        const int len,
        const Scalar* const base
    ) noexcept :
        vector_impl(len, 1, base)
    {
    }



    vector_impl(
        const int len,
        const int inc,
        const Scalar* const base
    ) noexcept :
        len_(len),
        inc_(inc),
        base_(const_cast<Scalar*>(base))
    {
    }



    vector_impl(const vector_impl& other) noexcept :
        vector_impl(other.len_, other.inc_, other.base_)
    {
    }



    vector_impl(vector_impl&& other) noexcept :
        vector_impl(other)
    {
        other.len_ = 0;
        other.inc_ = 0;
        other.base_ = nullptr;
    }



    ~vector_impl() noexcept = default;



public:

    vector_impl&
    operator =(const vector_impl& other) = delete;

    vector_impl&
    operator =(vector_impl&& other) = delete;



public:

    void
    reset(
        const int len,
        const Scalar* const base
    ) noexcept
    {
        reset(len, 1, base);
    }



    void
    reset(
        const int len,
        const int inc,
        const Scalar* const base
    ) noexcept
    {
        len_ = len;
        inc_ = inc;
        base_ = const_cast<Scalar*>(base);
    }



    void
    fill(Scalar s) noexcept
    {
        Scalar* dst = base_;

        for (int i = 0; i < len_; ++i)
        {
            *dst = s;

            dst += inc_;
        }
    }



    void
    copy(const vector_impl& other) noexcept
    {
        const Scalar* src = other.base_;
        Scalar* dst = base_;

        for (int i = 0; i < len_; ++i)
        {
            *dst = *src;

            src += other.inc_;
            dst += inc_;
        }
    }



    void
    swap(vector_impl& other) noexcept
    {
        using std::swap;

        swap(len_, other.len_);
        swap(inc_, other.inc_);
        swap(base_, other.base_);
    }



    void
    ref(const vector_impl& other) noexcept
    {
        len_ = other.len_;
        inc_ = other.inc_;
        base_ = other.base_;
    }



public:

    int
    len() const noexcept
    {
        return len_;
    }



    int
    inc() const noexcept
    {
        return inc_;
    }



    Scalar*
    base() const noexcept
    {
        return base_;
    }



    Scalar*
    ptr(int i) const noexcept
    {
        return base_ + i * inc_;
    }


    Scalar&
    elem(int i) const noexcept
    {
        return *ptr(i);
    }


    vector_impl
    elem(subscript ii) const noexcept
    {
        return {ii.len, inc_, ptr(ii.off)};
    }

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#endif
