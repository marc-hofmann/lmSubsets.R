// Copyright 2018  Marc Hofmann
// Copyright 2008-2018  Free Software Foundation, Inc.
//
// This file is part of the 'lmSubsets' R extension.  It is based on
// the GNU ISO C++ Library.
//
// 'lmSubsets' is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'lmSubsets' is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



#if __cplusplus < 201402L



#include <cstddef>  // std::size_t
#include <memory>  // std::unique_ptr
#include <type_traits>  // std::remove_extent
#include <utility>  // std::forward



// type_traits
namespace std
{



template<typename _Tp>
using remove_extent_t = typename remove_extent<_Tp>::type;



}  // namespace std



// iterator
namespace std
{



template<typename _Container>
inline constexpr auto
cbegin(const _Container& __cont) noexcept(noexcept(std::begin(__cont)))
    -> decltype(std::begin(__cont))
{
    return std::begin(__cont);
}



template<typename _Container>
inline constexpr auto
cend(const _Container& __cont) noexcept(noexcept(std::end(__cont)))
    -> decltype(std::end(__cont))
{
    return std::end(__cont);
}



template<typename _Iterator>
inline reverse_iterator<_Iterator>
make_reverse_iterator(_Iterator __i)
{
    return reverse_iterator<_Iterator>(__i);
}



}  // namespace std



// memory
namespace std
{



template<typename _Tp>
struct _MakeUniq
{
    typedef unique_ptr<_Tp> __single_object;
};



template<typename _Tp>
struct _MakeUniq<_Tp[]>
{
    typedef unique_ptr<_Tp[]> __array;
};



template<typename _Tp, size_t _Bound>
struct _MakeUniq<_Tp[_Bound]>
{
    struct __invalid_type { };
};



template<typename _Tp, typename... _Args>
inline typename _MakeUniq<_Tp>::__single_object
make_unique(_Args&&... __args)
{
    return unique_ptr<_Tp>(new _Tp(std::forward<_Args>(__args)...));
}



template<typename _Tp>
inline typename _MakeUniq<_Tp>::__array
make_unique(size_t __num)
{
    return unique_ptr<_Tp>(new remove_extent_t<_Tp>[__num]());
}



template<typename _Tp, typename... _Args>
inline typename _MakeUniq<_Tp>::__invalid_type
make_unique(_Args&&...) = delete;



}  // namespace std



#endif
