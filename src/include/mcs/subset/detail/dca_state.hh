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



#ifndef MCS_SUBSET_DETAIL_DCA_STATE_HH
#define MCS_SUBSET_DETAIL_DCA_STATE_HH



#include <type_traits>  // std::is_same
#include <utility>  // std::declval
#include <vector>



#include "mcs/core/matrix.hh"

#include "mcs/subset/detail/dca_node.hh"
#include "mcs/subset/detail/dca_partial.hh"
#include "mcs/subset/detail/dca_qrz.hh"

#include "mcs/util/algo.hh"  // algo::concat, algo::iota, algo::map, algo::plus,
                             // algo::repeat, algo::transform
#include "mcs/util/function_traits.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename NodeXfer>
class dca_state_base
{

    using dca_node = detail::dca_node<Scalar>;

    using dca_qrz = detail::dca_qrz<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    using node_xfer_inst = decltype(
        std::declval<NodeXfer>().
        make(std::declval<dca_state_base>())
    );



public:

    std::vector<dca_node> node_stk_;

    typename decltype(node_stk_)::iterator cur_node_;

    typename decltype(node_stk_)::iterator nxt_node_;

    node_xfer_inst node_xfer_;

    dca_qrz qrz_;

    int root_size_;

    int root_mark_;

    int root_rank_;

    Scalar root_rss_;



public:

    dca_state_base(
        matrix_cspan ay_mat,
        const int mark,
        const NodeXfer& node_xfer
    ) noexcept :
        qrz_(ay_mat.ncol() - 1),
        root_size_(ay_mat.ncol() - 1),
        root_mark_(mark),
        root_rank_(root_size_ - root_mark_)
    {
        const int n = root_size_;
        const int k = root_mark_;
        const int p = root_rank_;

        node_stk_.reserve(p);
        for (int i = 0; i < p; ++i)
        {
            node_stk_.emplace_back(p);
        }

        cur_node_ = node_stk_.begin();

        nxt_node_ = cur_node_ + 1;
        nxt_node_->root(qrz_.rz(ay_mat)({k, p+1}, {k, p+1}));

        root_rss_ = nxt_node_->rss();

        node_xfer_ = node_xfer.make(*this);
    }



public:

    bool
    is_final() const noexcept
    {
        return cur_node_ == nxt_node_;
    }



    void
    next_node() noexcept
    {
        node_xfer_(*nxt_node_, *cur_node_);
        --nxt_node_;
    }



    int
    node_size() const  noexcept
    {
        return root_mark_ + cur_node_->size();
    }



    int
    node_mark() const noexcept
    {
        return root_mark_ + cur_node_->mark();
    }



    Scalar
    node_rss() const noexcept
    {
        return cur_node_->rss();
    }



    void
    drop_column(const int mark) noexcept
    {
        ++nxt_node_;
        cur_node_->drop_column(mark - root_mark_, *nxt_node_, qrz_);
    }



    int
    root_size() const noexcept
    {
        return root_size_;
    }



    int
    root_mark() const noexcept
    {
        return root_mark_;
    }



    int
    root_rank() const noexcept
    {
        return root_rank_;
    }



    Scalar
    root_rss() const noexcept
    {
        return root_rss_;
    }



    const dca_qrz&
    qrz() const noexcept
    {
        return qrz_;
    }

};



template<typename Scalar,
         typename NodeXfer>
class dca_state_all : private dca_state_base<Scalar, NodeXfer>
{

    using base = dca_state_base<Scalar, NodeXfer>;

    using dca_partial = detail::dca_partial_all<Scalar>;

    using dca_result = detail::dca_result<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;



private:

    dca_partial partial_;

    int nbest_;



public:

    dca_state_all(
        matrix_cspan ay_mat,
        const int mark,
        const int nbest,
        const NodeXfer& node_xfer
    ) noexcept :
        base(ay_mat, mark, node_xfer),
        partial_(base::root_rank(), nbest),
        nbest_(nbest)
    {
    }



public:

    using base::is_final;

    using base::node_size;

    using base::node_mark;

    using base::node_rss;

    using base::drop_column;

    using base::root_rss;



    Scalar
    rss_inf() const noexcept
    {
        return root_rss();
    }



    void
    next_node() noexcept
    {
        base::next_node();
        partial_.update(*base::cur_node_);
    }



    Scalar
    rss_bound() const noexcept
    {
        return base::cur_node_->rss();
    }



    Scalar
    min_rss(const int size) const noexcept
    {
        return partial_.min_rss(size - base::root_mark());
    }



    std::vector<std::vector<dca_result>>
    table() const noexcept
    {
        const int root_mark = this->root_mark();

        const auto prefix = util::iota(0, root_mark);

        const auto xform = [&prefix, &root_mark](
            const dca_result& r
        ) -> dca_result {
            if (!r)  return {};

            return {
                util::concat(
                    prefix,
                    util::transform(r, util::plus(root_mark))
                ),
                r.key()
            };
        };

        return util::concat(
            util::repeat(util::repeat(dca_result(), nbest_), root_mark),
            util::transform(partial_.results(), util::map(xform))
        );
    }

};



template<typename Scalar,
         typename CostFunc,
         typename NodeXfer>
class dca_state_best : private dca_state_base<Scalar, NodeXfer>
{

    using cost_func_traits = mcs::util::function_traits<CostFunc>;

    static_assert(
        std::is_same<
            typename cost_func_traits::signature,
            double(int,double)
        >::value,
        "cost function must be 'double(int,double)'"
    );



    using base = dca_state_base<Scalar, NodeXfer>;

    using dca_partial = detail::dca_partial_best<Scalar>;

    using dca_result = detail::dca_result<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;



private:

    dca_partial partial_;

    CostFunc cost_func_;

    Scalar cost_inf_;



public:

    dca_state_best(
        matrix_cspan ay_mat,
        const int mark,
        const CostFunc& cost_func,
        const int nbest,
        const NodeXfer& node_xfer
    ) noexcept :
        base(ay_mat, mark, node_xfer),
        partial_(base::root_rank(), nbest),
        cost_func_(cost_func)
    {
        cost_inf_ = cost_func_(mark + 1, base::root_rss());
    }



public:

    using base::is_final;

    using base::node_size;

    using base::node_mark;

    using base::node_rss;

    using base::drop_column;



    Scalar
    cost_inf() const noexcept
    {
        return cost_inf_;
    }



    void
    next_node() noexcept
    {
        base::next_node();
        partial_.update(
            *base::cur_node_,
            [this, root_mark = base::root_mark()](
                const int size,
                const Scalar rss
            ) -> Scalar {
                return cost_func_(root_mark + size, rss);
            }
        );
    }



    Scalar
    cost_bound(const int mark) const noexcept
    {
        return cost_func_(mark + 1, base::cur_node_->rss());
    }



    Scalar
    min_cost() const noexcept
    {
        return partial_.min_cost();
    }



    std::vector<dca_result>
    table() const noexcept
    {
        const int root_mark = this->root_mark();

        const auto prefix = util::iota(0, root_mark);

        const auto xform = [&prefix, &root_mark](
            const dca_result& r
        ) -> dca_result {
            if (!r)  return {};

            return {
                util::concat(
                    prefix,
                    util::transform(r, util::plus(root_mark))
                ),
                r.key()
            };
        };

        return util::transform(partial_.results(), xform);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
