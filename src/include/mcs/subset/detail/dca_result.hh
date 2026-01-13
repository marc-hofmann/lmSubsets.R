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



#ifndef MCS_SUBSET_DETAIL_DCA_RESULT_HH
#define MCS_SUBSET_DETAIL_DCA_RESULT_HH



#include <limits>  // std::numeric_limits
#include <utility>  // std::swap
#include <vector>



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_result
{

public:

    using value_type = int;

    using iterator = typename std::vector<int>::iterator;

    using const_iterator = typename std::vector<int>::const_iterator;



    friend void
    swap(dca_result& a, dca_result& b) noexcept
    {
        a.swap(b);
    }



    friend const_iterator
    cbegin(const dca_result& r) noexcept
    {
        return r.cbegin();
    }



    friend const_iterator
    cend(const dca_result& r) noexcept
    {
        return r.cend();
    }



private:

    std::vector<int> subset_;

    Scalar key_;



public:

    dca_result() noexcept :
        key_(std::numeric_limits<Scalar>::quiet_NaN())
    {
    }



    dca_result(
        const std::vector<int>& subset,
        const Scalar key
    ) noexcept :
        subset_(subset),
        key_(key)
    {
        if (subset.size() == 0)
        {
            key_ = std::numeric_limits<Scalar>::quiet_NaN();
        }
    }



public:

    void
    swap(dca_result& other) noexcept
    {
        subset_.swap(other.subset_);
        std::swap(key_, other.key_);
    }



    const_iterator
    cbegin() const noexcept
    {
        return subset_.cbegin();
    }



    const_iterator
    cend() const noexcept
    {
        return subset_.cend();
    }



public:

    operator bool() const noexcept
    {
        return size() > 0;
    }



    int&
    operator [](const int i) noexcept
    {
        return subset_[i];
    }



    const int&
    operator [](const int i) const noexcept
    {
        return subset_[i];
    }



public:

    int
    size() const noexcept
    {
        return subset_.size();
    }



    const std::vector<int>&
    subset() const noexcept
    {
        return subset_;
    }



    Scalar
    key() const noexcept
    {
        return key_;
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
