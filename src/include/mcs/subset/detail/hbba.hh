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



#ifndef MCS_SUBSET_DETAIL_HBBA_HH
#define MCS_SUBSET_DETAIL_HBBA_HH



#include <vector>



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename DcaState>
int
hbba_best(
    DcaState& state,
    const Scalar tau
) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        const Scalar min_cost = state.min_cost();

        for (int j = k; j < n - 1; ++j)
        {
            if (tau * state.cost_bound(j) >= min_cost)
            {
                break;
            }

            state.drop_column(j);
        }

        ++node_cnt;
    }

    return node_cnt;
}



template<typename Scalar,
         typename DcaState>
int
hbba_all(
    DcaState& state,
    const std::vector<Scalar>& tau
) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        for (int j_sup = n - 1; j_sup > k; --j_sup)
        {
            if (tau[j_sup - 1] * state.rss_bound() < state.min_rss(j_sup))
            {
                for (int j = k; j < j_sup; ++j)
                {
                    state.drop_column(j);
                }

                break;
            }
        }

        ++node_cnt;
    }

    return node_cnt;
}



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
