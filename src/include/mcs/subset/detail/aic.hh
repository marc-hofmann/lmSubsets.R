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



#ifndef MCS_SUBSET_DETAIL_AIC_HH
#define MCS_SUBSET_DETAIL_AIC_HH



#include "mcs/subset/detail/log_lik.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class aic
{

    using log_lik = detail::log_lik<Scalar>;



private:

    Scalar k_;

    log_lik ll_;



public:

    constexpr
    aic(
        const Scalar k,
        const int nobs
    ) noexcept :
        k_(k),
        ll_(nobs)
    {
    }


    constexpr Scalar
    operator ()(
        const int size,
        const Scalar rss
    ) const noexcept
    {
        // size + 1  to account for sd as estimated parameter
        // const int npar = size + 1;

        return Scalar(-2.0) * ll_(rss) + k_ * Scalar(size + 1);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs


#endif
