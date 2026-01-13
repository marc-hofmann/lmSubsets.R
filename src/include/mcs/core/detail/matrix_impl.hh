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



#ifndef MCS_CORE_DETAIL_MATRIX_IMPL_HH
#define MCS_CORE_DETAIL_MATRIX_IMPL_HH



#include <algorithm>  // std::copy_n, std::fill_n
#include <utility>  // std::swap



#include "mcs/core/detail/subscript.hh"
#include "mcs/core/detail/vector_impl.hh"



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
class matrix_impl
{

    using vector_impl = detail::vector_impl<Scalar>;



private:

    int nrow_;

    int ncol_;

    int ldim_;

    Scalar* base_;



public:

    matrix_impl() noexcept :
        matrix_impl(0, 0, 0, nullptr)
    {
    }



    matrix_impl(
        const int nrow,
        const int ncol,
        const Scalar* const base
    ) noexcept :
        matrix_impl(nrow, ncol, nrow, base)
    {
    }




    matrix_impl(
        const int nrow,
        const int ncol,
        const int ldim,
        const Scalar* const base
    ) noexcept :
        nrow_(nrow),
        ncol_(ncol),
        ldim_(ldim),
        base_(const_cast<Scalar*>(base))
    {
    }



    matrix_impl(const matrix_impl& other) noexcept :
        matrix_impl(other.nrow_, other.ncol_, other.ldim_, other.base_)
    {
    }



    matrix_impl(matrix_impl&& other) noexcept :
        matrix_impl(other)
    {
        other.nrow_ = 0;
        other.ncol_ = 0;
        other.ldim_ = 0;
        other.base_ = nullptr;
    }



    ~matrix_impl() noexcept = default;



public:

    matrix_impl&
    operator =(const matrix_impl& other) = delete;

    matrix_impl&
    operator =(matrix_impl&& other) = delete;



public:

    void
    reset(
        const int nrow,
        const int ncol,
        const Scalar* const base
    ) noexcept
    {
        reset(nrow, ncol, nrow, base);
    }



    void
    reset(
        const int nrow,
        const int ncol,
        const int ldim,
        const Scalar* const base
    ) noexcept
    {
        nrow_ = nrow;
        ncol_ = ncol;
        ldim_ = ldim;
        base_ = const_cast<Scalar*>(base);
    }



    void
    fill(Scalar s) noexcept
    {
        Scalar* dst = base_;

        for (int j = 0; j < ncol_; ++j)
        {
            std::fill_n(dst, nrow_, s);

            dst += ldim_;
        }
    }



    void
    copy(const matrix_impl& other) noexcept
    {
        const Scalar* src = other.base_;
        Scalar* dst = base_;

        const int src_ldim = other.ldim_;
        const int dst_ldim = ldim_;

        for (int j = 0; j < ncol_; ++j)
        {
            std::copy_n(src, nrow_, dst);

            src += src_ldim;
            dst += dst_ldim;
        }
    }



    void
    swap(matrix_impl& other) noexcept
    {
        using std::swap;

        swap(nrow_, other.nrow_);
        swap(ncol_, other.ncol_);
        swap(ldim_, other.ldim_);
        swap(base_, other.base_);
    }



    void
    ref(const matrix_impl& other) noexcept
    {
        nrow_ = other.nrow_;
        ncol_ = other.ncol_;
        ldim_ = other.ldim_;
        base_ = other.base_;
    }



public:

    int
    nrow() const noexcept
    {
        return nrow_;
    }



    int
    ncol() const noexcept
    {
        return ncol_;
    }



    int
    ldim() const noexcept
    {
        return ldim_;
    }



    Scalar*
    base() const noexcept
    {
        return base_;
    }



    Scalar*
    ptr(int i, int j) const noexcept
    {
        return base_ + (j * ldim_) + i;
    }



    Scalar&
    elem(int i, int j) const noexcept
    {
        return *ptr(i, j);
    }



    vector_impl
    elem(int i, subscript jj) const noexcept
    {
        return {jj.len, ldim_, ptr(i, jj.off)};
    }



    vector_impl
    elem(subscript ii, int j) const noexcept
    {
        return {ii.len, 1, ptr(ii.off, j)};
    }



    matrix_impl
    elem(subscript ii, subscript jj) const noexcept
    {
        return {ii.len, jj.len, ldim_, ptr(ii.off, jj.off)};
    }



    vector_impl
    row(int i) const noexcept
    {
        return {ncol_, ldim_, base_ + i};
    }



    matrix_impl
    row(subscript ii) const noexcept
    {
        return {ii.len, ncol_, ldim_, base_ + ii.off};
    }



    vector_impl
    col(int j) const noexcept
    {
        return {nrow_, 1, base_ + (j * ldim_)};
    }



    matrix_impl
    col(subscript jj) const noexcept
    {
        return {nrow_, jj.len, ldim_, base_ + (jj.off * ldim_)};
    }

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#endif
