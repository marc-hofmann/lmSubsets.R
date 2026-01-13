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
#define MCS_CORE_DETAIL_VECTOR_HH



#include <memory>  // std::unique_ptr



#include "mcs/core/detail/vector_impl.hh"



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
class vector
{

    friend void
    swap(vector& a, vector& b) noexcept
    {
        a.swap(b);
    }



    friend vector<Scalar&>;

    friend vector<const Scalar>;

    friend vector<const Scalar&>;



private:

    std::unique_ptr<Scalar[]> data_;

    vector_impl<Scalar> impl_;



public:

    vector() = delete;

    vector(int len) noexcept;

    vector(const vector& other) noexcept;

    vector(vector&& other) noexcept;

    vector(const vector<Scalar&>& other) noexcept;

    vector(vector<Scalar&>&& other) noexcept;

    vector(const vector<const Scalar>& other) noexcept;

    vector(vector<const Scalar>&& other) noexcept;

    vector(const vector<const Scalar&>& other) noexcept;

    vector(vector<const Scalar&>&& other) noexcept;

    ~vector() noexcept = default;



public:

    vector&
    operator =(Scalar s) noexcept;

    vector&
    operator =(const vector& other) noexcept;

    vector&
    operator =(vector&& other) noexcept;

    vector&
    operator =(const vector<Scalar&>& other) noexcept;

    vector&
    operator =(vector<Scalar&>&& other) noexcept;

    vector&
    operator =(const vector<const Scalar>& other) noexcept;

    vector&
    operator =(vector<const Scalar>&& other) noexcept;

    vector&
    operator =(const vector<const Scalar&>& other) noexcept;

    vector&
    operator =(vector<const Scalar&>&& other) noexcept;



public:

    void
    swap(vector& other) noexcept;



public:

    Scalar&
    operator ()(int i) noexcept;

    const Scalar&
    operator ()(int i) const noexcept;

    vector<Scalar&>
    operator ()(subscript ii) noexcept;

    vector<const Scalar&>
    operator ()(subscript ii) const noexcept;



public:

    int
    len() const noexcept;

    int
    inc() const noexcept;

    Scalar*
    base() noexcept;

    const Scalar*
    base() const noexcept;

    Scalar*
    ptr(int i) noexcept;

    const Scalar*
    ptr(int i) const noexcept;

    Scalar&
    elem(int i) noexcept;

    const Scalar&
    elem(int i) const noexcept;

    vector<Scalar&>
    elem(subscript ii) noexcept;

    vector<const Scalar&>
    elem(subscript ii) const noexcept;

};



template<typename Scalar>
class vector<Scalar&>
{

    friend void
    swap(vector& a, vector& b) noexcept
    {
        a.swap(b);
    }



    friend vector<Scalar>;

    friend vector<const Scalar>;

    friend vector<const Scalar&>;



private:

    vector_impl<Scalar> impl_;



public:

    vector() = delete;

    vector(int len, Scalar* base) noexcept;

    vector(int len, int inc, Scalar* base) noexcept;

    vector(vector_impl<Scalar>&& impl) noexcept;

    vector(vector<Scalar>& other) noexcept;

    vector(const vector<Scalar>& other) = delete;

    vector(vector<Scalar>&& other) = delete;

    vector(vector& other) noexcept;

    vector(const vector& other) = delete;

    vector(vector&& other) noexcept;

    vector(const vector<const Scalar>& other) = delete;

    vector(vector<const Scalar>&& other) = delete;

    vector(const vector<const Scalar&>& other) = delete;

    vector(vector<const Scalar&>&& other) = delete;

    ~vector() noexcept = default;



public:

    vector&
    operator =(Scalar s) noexcept;

    vector&
    operator =(const vector<Scalar>& other) noexcept;

    vector&
    operator =(vector<Scalar>&& other) noexcept;

    vector&
    operator =(const vector& other) noexcept;

    vector&
    operator =(vector&& other) noexcept;

    vector&
    operator =(const vector<const Scalar>& other) noexcept;

    vector&
    operator =(vector<const Scalar>&& other) noexcept;

    vector&
    operator =(const vector<const Scalar&>& other) noexcept;

    vector&
    operator =(vector<const Scalar&>&& other) noexcept;



public:

    void
    swap(vector& other) noexcept;



public:

    Scalar&
    operator ()(int i) noexcept;

    const Scalar&
    operator ()(int i) const noexcept;

    vector
    operator ()(subscript ii) noexcept;

    vector<const Scalar&>
    operator ()(subscript ii) const noexcept;



public:

    int
    len() const noexcept;

    int
    inc() const noexcept;

    Scalar*
    base() noexcept;

    Scalar const*
    base() const noexcept;

    Scalar*
    ptr(int i) noexcept;

    const Scalar*
    ptr(int i) const noexcept;

    Scalar&
    elem(int i) noexcept;

    const Scalar&
    elem(int i) const noexcept;

    vector
    elem(subscript ii) noexcept;

    vector<const Scalar&>
    elem(subscript ii) const noexcept;

};



template<typename Scalar>
class vector<const Scalar>
{

    friend vector<Scalar>;

    friend vector<Scalar&>;

    friend vector<const Scalar&>;



private:

    std::unique_ptr<Scalar[]> data_;

    vector_impl<Scalar> impl_;



public:

    vector() = delete;

    vector(int len) noexcept;

    vector(const vector<Scalar>& other) noexcept;

    vector(vector<Scalar>&& other) noexcept;

    vector(const vector<Scalar&>& other) noexcept;

    vector(vector<Scalar&>&& other) noexcept;

    vector(const vector& other) noexcept;

    vector(vector&& other) noexcept;

    vector(const vector<const Scalar&>& other) noexcept;

    vector(vector<const Scalar&>&& other) noexcept;

    ~vector() noexcept = default;



public:

    vector&
    operator =(const vector<Scalar>& other) = delete;

    vector&
    operator =(vector<Scalar>&& other) = delete;

    vector&
    operator =(const vector<Scalar&>& other) = delete;

    vector&
    operator =(vector<Scalar&>&& other) = delete;

    vector&
    operator =(const vector& other) = delete;

    vector&
    operator =(vector&& other) = delete;

    vector&
    operator =(const vector<const Scalar&>& other) = delete;

    vector&
    operator =(vector<const Scalar&>&& other) = delete;



public:

    void
    swap(vector& other) = delete;



public:

    const Scalar&
    operator ()(int i) const noexcept;

    vector<const Scalar&>
    operator ()(subscript ii) const noexcept;



public:

    int
    len() const noexcept;

    int
    inc() const noexcept;

    const Scalar*
    base() const noexcept;

    const Scalar*
    ptr(int i) const noexcept;

    const Scalar&
    elem(int i) const noexcept;

    vector<const Scalar&>
    elem(subscript ii) const noexcept;

};



template<typename Scalar>
class vector<const Scalar&>
{

    friend vector<Scalar>;

    friend vector<Scalar&>;

    friend vector<const Scalar>;



private:

    vector_impl<Scalar> impl_;

    std::unique_ptr<Scalar[]> keep_alive_;



public:

    vector() = delete;

    vector(int len, const Scalar* base) noexcept;

    vector(int len, int inc, const Scalar* base) noexcept;

    vector(vector_impl<Scalar>&& impl) noexcept;

    vector(const vector<Scalar>& other) noexcept;

    vector(vector<Scalar>&& other) noexcept;

    vector(const vector<Scalar&>& other) noexcept;

    vector(vector<Scalar&>&& other) noexcept;

    vector(const vector<const Scalar>& other) noexcept;

    vector(vector<const Scalar>&& other) noexcept;

    vector(const vector& other) noexcept;

    vector(vector&& other) noexcept;

    ~vector() noexcept = default;



public:

    vector&
    operator =(const vector<Scalar>& other) = delete;

    vector&
    operator =(vector<Scalar>&& other) = delete;

    vector&
    operator =(const vector<Scalar&>& other) = delete;

    vector&
    operator =(vector<Scalar&>&& other) = delete;

    vector&
    operator =(const vector<const Scalar>& other) = delete;

    vector&
    operator =(vector<const Scalar>&& other) = delete;

    vector&
    operator =(const vector& other) = delete;

    vector&
    operator =(vector&& other) = delete;



public:

    void
    swap(vector& other) = delete;



public:

    const Scalar&
    operator ()(int i) const noexcept;

    vector
    operator ()(subscript ii) const noexcept;



public:

    int
    len() const noexcept;

    int
    inc() const noexcept;

    const Scalar*
    base() const noexcept;

    const Scalar*
    ptr(int i) const noexcept;

    const Scalar&
    elem(int i) const noexcept;

    vector
    elem(subscript ii) const noexcept;

};




}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#include "mcs/core/detail/vector.inc"
#endif
