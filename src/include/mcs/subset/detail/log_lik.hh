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



#ifndef MCS_SUBSET_DETAIL_LOG_LIK_HH
#define MCS_SUBSET_DETAIL_LOG_LIK_HH



#define _USE_MATH_DEFINES  // M_PI
#include <cmath>  // std::log, M_PI



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class log_lik
{

private:

    //
    // NOTE:  'std::log' is not 'constexpr'
    // static constexpr Scalar LOG_2PI_ = std::log(Scalar(2.0) * Scalar(M_PI));
    //
    // from 'Rmath.h':
    static constexpr Scalar LOG_2PI_ = 1.837877066409345483560659472811;


private:

    Scalar nobs_half_;

    Scalar log_nobs_;



public:

    constexpr
    log_lik(const int nobs) noexcept :
        nobs_half_(Scalar(0.5) * Scalar(nobs)),
        log_nobs_(std::log(nobs))
    {
    }



    constexpr Scalar
    operator ()(const Scalar rss) const noexcept
    {
        const Scalar log_rss = std::log(rss);

        return -nobs_half_ * (LOG_2PI_ - log_nobs_ + log_rss + Scalar(1.0));
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
