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



#ifndef MCS_UTIL_DETAIL_ALGO_HH
#define MCS_UTIL_DETAIL_ALGO_HH



#include <algorithm>  // std::for_each, std::generate, std::reverse, std::sort,
                      // std::sort_heap, std::transform
#include <iterator>  // std::back_inserter, std::begin, std::cbegin, std::cend,
                     // std::end
#include <numeric>  // std::iota
#include <vector>



namespace mcs    {
namespace util   {
namespace detail {



template<typename RandomIt,
         typename InputIt,
         typename OutputIt>
void
arrange(
    RandomIt first,
    RandomIt last,
    InputIt pos,
    OutputIt result
) noexcept
{
    arrange_n(first, last - first, pos, result);
}



template<typename RandomIt,
         typename InputIt,
         typename OutputIt>
void
arrange_n(
    RandomIt first,
    const int count,
    InputIt pos,
    OutputIt result
) noexcept
{
    std::generate_n(result, count, [first, pos]() mutable {
            return first[*(pos++)];
        });
}



template<typename Container,
         typename UnaryFunction>
void
for_each(
    const Container& input,
    const UnaryFunction& func
) noexcept
{
    using std::cbegin;
    using std::cend;

    std::for_each(cbegin(input), cend(input), func);
}



template<typename T>
auto
iota(
    const T& value,
    const int count
) noexcept
{
    std::vector<int> ret(count);
    std::iota(ret.begin(), ret.end(), value);
    return ret;
}



template<typename T>
auto
repeat(
    const T& value,
    const int count
) noexcept
{
    return std::vector<T>(count, value);
}



template<typename Container>
auto
reverse(
    const Container& input
) noexcept
{
    using value_type = typename Container::value_type;

    using std::cbegin;
    using std::cend;

    std::vector<value_type> ret(cbegin(input), cend(input));
    std::reverse(ret.begin(), ret.end());
    return ret;
}



template<typename ContainerA,
         typename ContainerB>
auto
concat(
    const ContainerA& a,
    const ContainerB& b
) noexcept
{
    using value_type = typename ContainerA::value_type;

    using std::cbegin;
    using std::cend;

    std::vector<value_type> ret;
    std::copy(cbegin(a), cend(a), std::back_inserter(ret));
    std::copy(cbegin(b), cend(b), std::back_inserter(ret));
    return ret;
}



template<typename Container,
         typename UnaryFunction>
auto
transform(
    const Container& input,
    const UnaryFunction& func
) noexcept
{
    using value_type = typename Container::value_type;
    using result_type = decltype(func(std::declval<value_type>()));

    using std::cbegin;
    using std::cend;

    std::vector<result_type> ret;
    std::transform(cbegin(input), cend(input), std::back_inserter(ret), func);
    return ret;
}



template<typename Container,
         typename Compare>
auto
sort(
    const Container& input,
    const Compare& comp
) noexcept
{
    using value_type = typename Container::value_type;

    using std::cbegin;
    using std::cend;

    std::vector<value_type> ret(cbegin(input), cend(input));
    std::sort(ret.begin(), ret.end(), comp);
    return ret;
}



template<typename Container,
         typename Compare>
auto
sort_heap(
    const Container& input,
    const Compare& comp
) noexcept
{
    using value_type = typename Container::value_type;

    using std::cbegin;
    using std::cend;

    std::vector<value_type> ret(cbegin(input), cend(input));
    std::sort_heap(ret.begin(), ret.end(), comp);
    return ret;
}



template<typename UnaryFunction>
auto
map(const UnaryFunction& func) noexcept
{
    return [f = func](const auto& _) {
        return transform(_, f);
    };
}



template<typename T>
auto
plus(const T& value) noexcept
{
    return [value = value](const T& x) {
        return value + x;
    };
}



}  // end namespace detail
}  // end namespace util
}  // end namespace mcs



#endif
