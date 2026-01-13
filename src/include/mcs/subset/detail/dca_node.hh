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



#ifndef MCS_SUBSET_DETAIL_DCA_NODE_HH
#define MCS_SUBSET_DETAIL_DCA_NODE_HH



#include <algorithm>  // std::max_element, std::sort
#include <cmath>  // std::pow
#include <iterator>  // std::distance
#include <numeric>  // std::iota
#include <utility>  // std::move, std::swap
#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/core/matrix.hh"

#include "mcs/subset/detail/dca_qrz.hh"
#include "mcs/subset/detail/dca_subset.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_node
{

    friend void
    swap(dca_node& a, dca_node& b) noexcept
    {
        a.swap(b);
    }

    using matrix = mcs::core::matrix<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    using dca_qrz = detail::dca_qrz<Scalar>;



private:

    std::vector<int> subset_;

    int mark_;

    matrix rz_mat_;



public:

    dca_node(const int root_size) noexcept :
        rz_mat_(root_size + 1, root_size + 1)
    {
        subset_.reserve(root_size);
    }



public:

    void
    swap(dca_node& other) noexcept
    {
        subset_.swap(other.subset_);
        std::swap(mark_, other.mark_);
        rz_mat_.swap(other.rz_mat_);
    }



public:

    void
    root(matrix_cspan rz_mat) noexcept
    {
        const int n = rz_mat.ncol() - 1;

        for (int j = 0; j < n; ++j)  subset_.push_back(j);
        mark_ = 0;
        rz_mat_ = rz_mat;
    }



    const std::vector<int>&
    subset() const noexcept
    {
        return subset_;
    }



    int
    size() const noexcept
    {
        return subset_.size();
    }



    int
    mark() const noexcept
    {
        return mark_;
    }



    int
    rank() const noexcept
    {
        return size() - mark();
    }



    Scalar
    rss() const noexcept
    {
        const int n = size();

        return std::pow(rz_mat_(n, n), 2);
    }



    template<typename Function>
    void
    for_each(Function f) const noexcept
    {
        const int n = size();
        const int k = mark_;

        gsl::span<const int> s = subset_;

        const Scalar* z_ptr = rz_mat_.ptr(n, n);
        Scalar rss = 0;

        for (int j = n; j > k; --j, --z_ptr)
        {
            rss += std::pow(*z_ptr, 2);

            f(s.first(j), rss);
        }
    }



    void
    drop_column(
        const int mark,
        dca_node& result,
        const dca_qrz& qrz
    ) const noexcept
    {
        const int n = size();
        const int k = mark;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        dca_subset::drop_column(subset_, k, result.subset_);
        result.mark_ = k;
        qrz.drop_column(rz_span, k, result.rz_mat_);
    }



    void
    preorder_complete(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1,
        std::vector<int>& aux_2
    ) const noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        std::iota(aux_2.begin(), aux_2.begin() + p, 0);
        std::sort(aux_2.begin(), aux_2.begin() + p,
                  [&aux_1](const int i, const int j) -> bool {
                      return aux_1[i] > aux_1[j];
                  });

        dca_subset::permute_complete(subset_, k, aux_2, result.subset_);
        result.mark_ = k;
        qrz.permute_complete(rz_span, k, aux_2, result.rz_mat_);
    }



    void
    preorder_partial_1(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1
    ) noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        const auto max = std::max_element(aux_1.begin(), aux_1.begin() + p);
        const int j = std::distance(aux_1.begin(), max);

        dca_subset::permute_partial_1(subset_, k, j);
        qrz.permute_partial_1(rz_span, k, j);

        swap(result);
    }



    void
    preorder_partial_2(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1
    ) noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        const auto max = std::max_element(aux_1.begin(), aux_1.begin() + p);
        const int j = std::distance(aux_1.begin(), max);

        dca_subset::permute_partial_2(subset_, k, j);
        qrz.permute_partial_2(rz_span, k, j);

        swap(result);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
