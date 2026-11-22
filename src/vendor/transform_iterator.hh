// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3 of the License, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef TRANSFORM_ITERATOR_HH
#define TRANSFORM_ITERATOR_HH

#include <iterator>
#include <boost/iterator/transform_iterator.hpp>

template <class Predicate, class Iterator>
class transform_random_access_iterator:
    public boost::transform_iterator<Predicate, Iterator>
{
public:
    using iterator_concept = std::random_access_iterator_tag;

    typedef Iterator iter_t;
    typedef boost::transform_iterator<Predicate, Iterator> base_t;
    using typename base_t::difference_type;

    transform_random_access_iterator() = default;
    transform_random_access_iterator(const base_t& iter): base_t(iter) {}
    transform_random_access_iterator(const Iterator& iter, const Predicate& pred = Predicate())
        : base_t(iter, pred) {}

    const base_t& base() const { return *this; }

    transform_random_access_iterator& operator++() { base_t::operator++(); return *this; }
    transform_random_access_iterator  operator++(int) { auto t = *this; ++*this; return t; }
    transform_random_access_iterator& operator--() { base_t::operator--(); return *this; }
    transform_random_access_iterator  operator--(int) { auto t = *this; --*this; return t; }

    transform_random_access_iterator& operator+=(difference_type n) { base_t::operator+=(n); return *this; }
    transform_random_access_iterator& operator-=(difference_type n) { base_t::operator-=(n); return *this; }

    decltype(auto) operator[](difference_type n) const { return *(*this + n); }

    friend transform_random_access_iterator operator+(transform_random_access_iterator i, difference_type n) { return i += n; }
    friend transform_random_access_iterator operator+(difference_type n, transform_random_access_iterator i) { return i += n; }
    friend transform_random_access_iterator operator-(transform_random_access_iterator i, difference_type n) { return i -= n; }
    friend difference_type operator-(const transform_random_access_iterator& a,
                                     const transform_random_access_iterator& b)
    {
        return a.base() - b.base();
    }
};

namespace std
{
template <class Predicate, class Iterator>
struct iterator_traits<transform_random_access_iterator<Predicate, Iterator> >
{
    typedef transform_random_access_iterator<Predicate, Iterator> titer_t;
    typedef typename titer_t::base_t base_t;
    typedef typename iterator_traits<base_t>::difference_type difference_type;
    typedef typename iterator_traits<base_t>::value_type value_type;
    typedef typename iterator_traits<base_t>::value_type reference;
    typedef typename iterator_traits<base_t>::pointer pointer;
    typedef typename iterator_traits<typename titer_t::iter_t>::iterator_category iterator_category;
};
}

#endif // TRANSFORM_ITERATOR_HH
