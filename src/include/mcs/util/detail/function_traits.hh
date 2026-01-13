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



#ifndef MCS_UTIL_DETAIL_FUNCTION_TRAITS_HH
#define MCS_UTIL_DETAIL_FUNCTION_TRAITS_HH



#include <iostream>
#include <tuple>
#include <utility>



namespace mcs    {
namespace util   {
namespace detail {



template<typename T>
struct identity
{
    using type = T;
};



template<typename T>
using identity_t = typename identity<T>::type;



// see also:  https://stackoverflow.com/a/7943765
template<typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{
};



template<typename Class,
         typename Result,
         typename... Args>
struct function_traits<Result(Class::*)(Args...) const> :
    public function_traits<Result(Args...) const>
{
};



template<typename Class,
         typename Result,
         typename... Args>
struct function_traits<Result(Class::*)(Args...)> :
    public function_traits<Result(Args...)>
{
};



template<typename Result,
         typename... Args>
struct function_traits<Result(*)(Args...)> :
    public function_traits<Result(Args...)>
{
};



template<typename Result,
         typename... Args>
struct function_traits<Result(Args...) const> :
    public function_traits<Result(Args...)>
{
};



template<typename Result,
         typename... Args>
struct function_traits<Result(Args...)>
{

    enum
    {
        arity = sizeof...(Args)
    };




    using result_type = Result;

    using signature = identity_t<Result(Args...)>;



    template <std::size_t i>
    struct arg
    {
        using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
    };

};



}  // end namespace detail
}  // end namespace util
}  // end namespace mcs



#endif
