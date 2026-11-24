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

#ifndef GRAPH_ADJACENCY_HH
#define GRAPH_ADJACENCY_HH

#include <vector>
#include <deque>
#include <utility>
#include <boost/iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/range/irange.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/iterator_facade.hpp>

#ifdef __clang__
#include <boost/algorithm/minmax_element.hpp>
#endif

#include "transform_iterator.hh"

#include "graph_adaptor.hh"
#include "openmp.hh"
#include "hash_map_wrap.hh"

namespace boost
{

// ========================================================================
// Forward declarations
// ========================================================================

template <class Vertex>
class adj_list;

// forward declaration of manipulation functions
template <class Vertex>
std::pair<typename adj_list<Vertex>::vertex_iterator,
          typename adj_list<Vertex>::vertex_iterator>
vertices(const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_iterator,
          typename adj_list<Vertex>::edge_iterator>
edges(const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
edge(Vertex s, Vertex t, const adj_list<Vertex>& g);

template <class Vertex, class F>
void edge_range_iter(Vertex s, Vertex t, const adj_list<Vertex>& g, F&& f);

template <class Vertex>
size_t out_degree(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
size_t in_degree(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
size_t degree(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
out_edges(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
in_edges(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
_all_edges_out(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
_all_edges_in(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::all_edge_iterator,
          typename adj_list<Vertex>::all_edge_iterator>
all_edges(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
adjacent_vertices(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
out_neighbors(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
in_neighbors(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
all_neighbors(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
size_t num_vertices(const adj_list<Vertex>& g);

template <class Vertex>
size_t num_edges(const adj_list<Vertex>& g);

template <class Vertex>
Vertex add_vertex(adj_list<Vertex>& g, size_t n = 1);

template <class Vertex>
void clear_vertex(Vertex v, adj_list<Vertex>& g);

template <class Vertex, class Pred>
void clear_vertex(Vertex v, adj_list<Vertex>& g, Pred&& pred);

template <class Vertex>
void remove_vertex(Vertex v, adj_list<Vertex>& g);

template <class Vertex>
void remove_vertex_fast(Vertex v, adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
add_edge(Vertex s, Vertex t, adj_list<Vertex>& g);

template <class Vertex>
void remove_edge(Vertex s, Vertex t, adj_list<Vertex>& g);

template <class Vertex>
void remove_edge(const typename adj_list<Vertex>::edge_descriptor& e,
                 adj_list<Vertex>& g);

// ========================================================================
// adj_list<Vertex>
// ========================================================================
//
// adj_list is a very simple adjacency list implementation for bidirectional
// graphs based on std::vector, meant to be reasonably efficient both
// memory-wise and computationally. It maintains a list of in and out-edges for
// each vertex, and each edge has a built-in index (which is replicated in both
// lists). For each edge, a total of 4 integers is necessary: the source and
// target vertices, in the in_edges and out_edges lists, respectively, and the
// (same) edge index in both lists. The integer type is given by the Vertex
// template parameter. It achieves about half as much memory as
// boost::adjacency_list with an edge index property map and the same integer
// type.

// The complexity guarantees and iterator invalidation rules are the same as
// boost::adjacency_list with vector storage selectors for both vertex and edge
// lists.

namespace detail
{

template <class V>
void soft_resize(V&& v, size_t size, auto&&... init)
{
    if (size > v.capacity())
        v.reserve(std::max(2 * v.capacity(), size));
    v.resize(size, init...);
}

template <class Vertex>
struct adj_edge_descriptor
{
    adj_edge_descriptor()
        : s(std::numeric_limits<Vertex>::max()),
          t(std::numeric_limits<Vertex>::max()),
          idx(std::numeric_limits<Vertex>::max()) {};

    adj_edge_descriptor(Vertex s_, Vertex t_, Vertex idx_)
        : s(s_), t(t_), idx(idx_) {}

    bool operator==(const adj_edge_descriptor& other) const
    {
        return idx == other.idx;
    }

    bool operator!=(const adj_edge_descriptor& other) const
    {
        return idx != other.idx;
    }

    bool operator<(const adj_edge_descriptor& other) const
    {
        return idx < other.idx;
    }

    Vertex s, t, idx;
};
} // namespace detail

template <class Vertex = size_t>
class adj_list
{
public:
    struct graph_tag {};
    typedef Vertex vertex_t;

    typedef detail::adj_edge_descriptor<Vertex> edge_descriptor;

    typedef std::vector<std::pair<vertex_t, vertex_t>> edge_list_t;
    typedef std::vector<std::pair<size_t, edge_list_t>> vertex_list_t;
    typedef typename integer_range<Vertex>::iterator vertex_iterator;

    adj_list(size_t n = 0):
        _edges(n), _n_edges(0), _edge_idx_range(0),
        _keep_epos(false), _keep_ehash(false) {}

    struct get_vertex
    {
        typedef const vertex_t& result_type;
        const vertex_t& operator()(const std::pair<vertex_t, vertex_t>& v) const
        { return v.first; }
    };

    typedef transform_random_access_iterator<get_vertex,
                                             typename edge_list_t::const_iterator>
        adjacency_iterator;

    typedef adjacency_iterator in_adjacency_iterator;

    static_assert(std::input_or_output_iterator<adjacency_iterator>);
    static_assert(std::random_access_iterator<adjacency_iterator>);

    template <class Base>
    struct edge_iter_facade:
        public boost::iterator_facade<Base,
                                      edge_descriptor,
                                      std::random_access_iterator_tag,
                                      edge_descriptor>
    {};

    template <class Dereference>
    struct base_edge_iterator:
        public edge_iter_facade<base_edge_iterator<Dereference>>
    {
        base_edge_iterator() {}
        base_edge_iterator(vertex_t v, const typename edge_list_t::const_iterator& iter)
            : _v(v), _iter(iter)
        {}

    private:
        friend class boost::iterator_core_access;
        void increment() { ++_iter; }

        void decrement() { --_iter; }

        template <class Distance>
        void advance(Distance n) { _iter += n; }

        auto distance_to(base_edge_iterator const& other) const
        {
            return other._iter - _iter;
        }

        bool equal(base_edge_iterator const& other) const
        {
            return _iter == other._iter;
        }

        edge_descriptor dereference() const
        {
            return Dereference::def(_v, *_iter, *this);
        }

    public:
        edge_descriptor operator*() const
        {
            return edge_iter_facade<base_edge_iterator>::operator*();
        }

        base_edge_iterator& operator+=(auto n) { _iter += n; return *this; }
        base_edge_iterator& operator-=(auto n) { _iter -= n; return *this; }

        decltype(auto) operator[](auto n) const { return *(*this + n); }

        friend base_edge_iterator operator+(base_edge_iterator i, auto n) { return i += n; }
        friend base_edge_iterator operator+(auto n, base_edge_iterator i) { return i += n; }
        friend base_edge_iterator operator-(base_edge_iterator i, auto n) { return i -= n; }
        friend auto operator-(const base_edge_iterator& a,  const base_edge_iterator& b)
        {
            return a._iter - b._iter;
        }

    protected:
        vertex_t _v;
        typename edge_list_t::const_iterator _iter;
    };

    struct make_out_edge
    {
        template <class Iter>
        static edge_descriptor def(vertex_t src,
                                   const std::pair<vertex_t, vertex_t>& v,
                                   Iter&&)
        { return edge_descriptor(src, v.first, v.second); }

        static edge_descriptor def(vertex_t src,
                                   const std::pair<vertex_t, vertex_t>& v)
        { return def(src, v, nullptr); }
    };

    struct make_in_edge
    {
        template <class Iter>
        static edge_descriptor def(vertex_t tgt,
                                   const std::pair<vertex_t, vertex_t>& v,
                                   Iter&&)
        { return edge_descriptor(v.first, tgt, v.second); }

        static edge_descriptor def(vertex_t tgt,
                                   const std::pair<vertex_t, vertex_t>& v)
        { return def(tgt, v, nullptr); }
    };

    typedef base_edge_iterator<make_out_edge> out_edge_iterator;
    typedef base_edge_iterator<make_in_edge> in_edge_iterator;

    static_assert(std::input_or_output_iterator<out_edge_iterator>);
    static_assert(std::input_or_output_iterator<in_edge_iterator>);
    static_assert(std::random_access_iterator<out_edge_iterator>);
    static_assert(std::random_access_iterator<in_edge_iterator>);

    template <class Iter>
    struct make_in_or_out_edge
    {
        template <class I>
        static edge_descriptor def(vertex_t u,
                                   const std::pair<vertex_t, vertex_t>& v,
                                   const I& i)
        {
            const Iter& iter = reinterpret_cast<const Iter&>(i);
            if (iter._iter < iter._pos)
                return edge_descriptor(u, v.first, v.second);
            else
                return edge_descriptor(v.first, u, v.second);
        }
    };

    struct all_edge_iterator:
        public base_edge_iterator<make_in_or_out_edge<all_edge_iterator>>
    {
        typedef base_edge_iterator<make_in_or_out_edge<all_edge_iterator>> base_t;
        all_edge_iterator() {}
        all_edge_iterator(vertex_t v,
                          const typename edge_list_t::const_iterator& iter,
                          const typename edge_list_t::const_iterator& pos)
            : base_t(v, iter),
              _pos(pos)
        {}

        const base_t& base() const { return *this; }

        all_edge_iterator& operator++() { base_t::operator++(); return *this; }
        all_edge_iterator  operator++(int) { auto t = *this; ++*this; return t; }
        all_edge_iterator& operator--() { base_t::operator--(); return *this; }
        all_edge_iterator  operator--(int) { auto t = *this; --*this; return t; }
        all_edge_iterator& operator+=(auto n) { base_t::operator+=(n); return *this; }
        all_edge_iterator& operator-=(auto n) { base_t::operator-=(n); return *this; }

        decltype(auto) operator[](auto n) const { return *(*this + n); }

        friend all_edge_iterator operator+(all_edge_iterator i, auto n) { return i += n; }
        friend all_edge_iterator operator+(auto n, all_edge_iterator i) { return i += n; }
        friend all_edge_iterator operator-(all_edge_iterator i, auto n) { return i -= n; }
        friend auto operator-(const all_edge_iterator& a,
                              const all_edge_iterator& b)
        {
            return a.base() - b.base();
        }

    private:
        friend struct make_in_or_out_edge<all_edge_iterator>;
        typename edge_list_t::const_iterator _pos;
    };

    static_assert(std::input_or_output_iterator<all_edge_iterator>);
    static_assert(std::random_access_iterator<all_edge_iterator>);

    class edge_iterator:
        public boost::iterator_facade<edge_iterator,
                                      edge_descriptor,
                                      boost::forward_traversal_tag,
                                      edge_descriptor>
    {
    public:
        edge_iterator() {}
        explicit edge_iterator(const typename vertex_list_t::const_iterator& vi_begin,
                               const typename vertex_list_t::const_iterator& vi_end,
                               const typename vertex_list_t::const_iterator& vi,
                               const typename edge_list_t::const_iterator& ei)
            : _vi_begin(vi_begin), _vi_end(vi_end), _vi(vi), _ei(ei)
        {
            // move position to first edge
            skip();
        }

    private:
        friend class boost::iterator_core_access;

        void skip()
        {
            //skip empty vertices
            while (_vi != _vi_end &&
                   _ei == _vi->second.begin() + _vi->first)
            {
                ++_vi;
                if (_vi != _vi_end)
                    _ei = _vi->second.begin();
            }
        }

        void increment()
        {
            ++_ei;
            skip();
        }

        bool equal(edge_iterator const& other) const
        {
            if (_vi_begin == _vi_end)
                return _vi == other._vi;
            return _vi == other._vi && _ei == other._ei;
        }

        edge_descriptor dereference() const
        {
            return edge_descriptor(vertex_t(_vi - _vi_begin),
                                   _ei->first, _ei->second);
        }

        typename vertex_list_t::const_iterator _vi_begin;
        typename vertex_list_t::const_iterator _vi_end;
        typename vertex_list_t::const_iterator _vi;
        typename edge_list_t::const_iterator _ei;
    };

    static_assert(std::input_or_output_iterator<edge_iterator>);

    void reindex_edges()
    {
        _free_idx.clear();
        _edge_idx_range = 0;
        for (auto& es : _edges)
            es.second.resize(es.first);
        for (size_t i = 0; i < _edges.size(); ++i)
        {
            auto pos = _edges[i].first;
            auto& es = _edges[i].second;
            for (size_t j = 0; j < pos; ++j)
            {
                auto& oe = es[j];
                Vertex v = oe.first;
                oe.second = _edge_idx_range;
                _edges[v].second.emplace_back(i, _edge_idx_range);
                _edge_idx_range++;
            }
        }

        if (_concurrent)
        {
            _concurrent = false;
            set_concurrent(true);
        }

        if (_keep_epos)
            rebuild_epos();

        if (_keep_ehash)
            rebuild_ehash();
    }

    void set_concurrent(bool concurrent)
    {
        if (concurrent == _concurrent)
            return;

        _concurrent = concurrent;
        if (_concurrent)
        {
            size_t nt = get_num_threads();
            _free_idx_m.resize(nt);
            size_t pos = 0;
            while (!_free_idx.empty())
            {
                _free_idx_m[pos++ % nt].push_back(_free_idx.back());
                _free_idx.pop_back();
            }
            _free_idx.shrink_to_fit();
        }
        else
        {
            for (auto& f : _free_idx_m)
            {
                _free_idx.insert(_free_idx.end(), f.begin(), f.end());
                f.clear();
                f.shrink_to_fit();
            }
            _free_idx_m.clear();
            _free_idx_m.shrink_to_fit();
        }
    }

    bool get_concurrent()
    {
        return _concurrent;
    }

    void set_keep_epos(bool keep)
    {
        if (keep)
        {
            if (!_keep_epos)
                rebuild_epos();
        }
        else
        {
            std::decay_t<decltype(_epos)> tmp;
            _epos.swap(tmp);
        }
        _keep_epos = keep;
    }

    bool get_keep_epos()
    {
        return _keep_epos;
    }

    void reserve_epos(size_t capacity)
    {
        if (_keep_epos && _epos.size() < capacity)
            _epos.resize(capacity);
        if (_keep_ehash && _ehpos.size() < capacity)
            _ehpos.resize(capacity);
    }

    void set_keep_ehash(bool keep)
    {
        if (keep)
        {
            if (!_keep_ehash)
                rebuild_ehash();
        }
        else
        {
            std::decay_t<decltype(_ehash)> htmp;
            _ehash.swap(htmp);
            std::decay_t<decltype(_ehpos)> ptmp;
            _ehpos.swap(ptmp);
        }
        _keep_ehash = keep;
    }

    bool get_keep_ehash()
    {
        return _keep_ehash;
    }

    size_t get_edge_index_range() const { return _edge_idx_range; }

    static constexpr Vertex null_vertex() { return std::numeric_limits<Vertex>::max(); }

    void shrink_to_fit()
    {
        _edge_idx_range = 0;
        #pragma omp parallel for schedule(runtime) reduction(max:_edge_idx_range)
        for (size_t i = 0; i < _edges.size(); ++i)
        {
            _edges[i].second.shrink_to_fit();
            for (auto& x : _edges[i].second)
                _edge_idx_range = std::max(_edge_idx_range, x.second + 1);
            if (_keep_ehash)
            {
                auto& h =_ehash[i];
                std::decay_t<decltype(h)> h_new(h);
                h.swap(h_new);
            }
        }
        _edges.shrink_to_fit();

        bool concurrent = get_concurrent();
        set_concurrent(false);
        auto iter_idx = std::remove_if(_free_idx.begin(),
                                       _free_idx.end(),
                                       [&](auto idx) -> bool
                                       {return idx >= _edge_idx_range;});
        _free_idx.erase(iter_idx, _free_idx.end());
        _free_idx.shrink_to_fit();
        set_concurrent(concurrent);

        if (_keep_epos)
            _epos.resize(_edge_idx_range);
        _epos.shrink_to_fit();
        _ehash.shrink_to_fit();
        if (_keep_ehash)
            _ehpos.resize(_edge_idx_range);
        _ehpos.shrink_to_fit();
    }

    void reverse_edge(edge_descriptor& e) const
    {
        auto& elist = _edges[e.s];
        auto pos = elist.first;
        auto& es = elist.second;
        if (_keep_epos)
        {
            auto& epos = _epos[e.idx];
            if (epos.first >= pos || es[epos.first].second != e.idx)
                std::swap(e.s, e.t);
        }
        else
        {
            bool found = false;
            for (size_t i = 0; i < pos; ++i)
            {
                if (es[i].second == e.idx)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
                std::swap(e.s, e.t);
        }
    }

private:
    vertex_list_t _edges;        // adjacency lists
    size_t _n_edges;             // total number of edges
    size_t _edge_idx_range;      // next edge index in range

    // indices of deleted edges to be used up for new edges to avoid very large
    // indices, and unnecessary property map memory use
    std::vector<size_t> _free_idx;
    std::vector<std::vector<size_t>> _free_idx_m; // concurrent version

    // optional concurrent version of the above
    bool _concurrent = false;

    // optional edge positions inside adjacency lists
    bool _keep_epos = false;
    std::vector<std::pair<uint32_t, uint32_t>> _epos; // out, in

    // optional hash table with adjacencies
    bool _keep_ehash = false;
    std::vector<gt_hash_map<Vertex, std::vector<Vertex>>> _ehash;
    std::vector<uint64_t> _ehpos;

    Vertex get_free_idx()
    {
        Vertex idx;
        if (_concurrent)
        {
            auto& free_idx = _free_idx_m[get_thread_num()];
            if (free_idx.empty())
            {
                #pragma omp atomic capture
                idx = _edge_idx_range++;
            }
            else
            {
                idx = free_idx.back();
                free_idx.pop_back();
            }
        }
        else
        {
            if (_free_idx.empty())
            {
                idx = _edge_idx_range++;
            }
            else
            {
                idx = _free_idx.back();
                _free_idx.pop_back();
            }
        }
        return idx;
    }

    void put_free_index(Vertex idx)
    {
        if (_concurrent)
            _free_idx_m[get_thread_num()].push_back(idx);
        else
            _free_idx.push_back(idx);
    }

    auto& get_epos(Vertex idx)
    {
        return _epos[idx];
    }

    void set_epos(Vertex idx, const std::pair<size_t, size_t>& epos)
    {
        if (idx >= _epos.size())
            boost::detail::soft_resize(_epos, idx + 1);
        _epos[idx] = epos;
    }

    void rebuild_epos()
    {
        _epos.resize(_edge_idx_range);

        #pragma omp parallel for schedule(runtime) if (_edges.size() > get_openmp_min_thresh())
        for (size_t v = 0; v < _edges.size(); ++v)
        {
            auto& pes = _edges[v];
            auto pos = pes.first;
            auto& es = pes.second;
            for (size_t j = 0; j < es.size(); ++j)
            {
                size_t idx = es[j].second;
                if (j < pos)
                    _epos[idx].first = j;
                else
                    _epos[idx].second = j;
            }
        }
        //check_epos();
    }

    void check_epos()
    {
#ifndef NDEBUG
        for (auto& pes : _edges)
        {
            auto pos = pes.first;
            auto& es = pes.second;
            for (size_t j = 0; j < es.size(); ++j)
            {
                assert(es[j].second < _epos.size());
                if (j < pos)
                    assert(_epos[es[j].second].first == j);
                else
                    assert(_epos[es[j].second].second == j);
            }
        }
#endif
    }

    void rebuild_ehash()
    {
        _ehash.clear();
        _ehash.resize(_edges.size());
        _ehpos.resize(_edge_idx_range);
        #pragma omp parallel for schedule(runtime) if (_edges.size() > get_openmp_min_thresh())
        for (size_t v = 0; v < _edges.size(); ++v)
        {
            auto [iter, eend] = out_edges(v, *this);
            for (; iter != eend; ++iter)
                add_ehash(*iter);
        }
    }

    void add_ehash(const edge_descriptor& e)
    {
        auto& es = _ehash[e.s][e.t];
        es.push_back(e.idx);
        if (_ehpos.size() <= e.idx)
            boost::detail::soft_resize(_ehpos, e.idx + 1);
        _ehpos[e.idx] = es.size() - 1;
    }

    void remove_ehash(const edge_descriptor& e)
    {
        auto& es = _ehash[e.s][e.t];
        auto pos = _ehpos[e.idx];
        _ehpos[es.back()] = pos;
        es[pos] = es.back();
        es.pop_back();
        if (es.empty())
            _ehash[e.s].erase(e.t);
    }

    // manipulation functions
    friend std::pair<vertex_iterator, vertex_iterator>
    vertices<>(const adj_list<Vertex>& g);

    friend std::pair<edge_iterator, edge_iterator>
    edges<>(const adj_list<Vertex>& g);

    friend std::pair<edge_descriptor, bool>
    edge<>(Vertex s, Vertex t, const adj_list<Vertex>& g);

    template <class V, class F>
    friend void edge_range_iter(V s, V t, const adj_list<V>& g, F&& f);

    friend size_t out_degree<>(Vertex v, const adj_list<Vertex>& g);

    friend size_t in_degree<>(Vertex v, const adj_list<Vertex>& g);

    friend size_t degree<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<out_edge_iterator, out_edge_iterator>
    out_edges<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<in_edge_iterator, in_edge_iterator>
    in_edges<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<out_edge_iterator, out_edge_iterator>
    _all_edges_out<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<in_edge_iterator, in_edge_iterator>
    _all_edges_in<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<all_edge_iterator, all_edge_iterator>
    all_edges<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    adjacent_vertices<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    out_neighbors<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    in_neighbors<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    all_neighbors<>(Vertex v, const adj_list<Vertex>& g);

    friend size_t num_vertices<>(const adj_list<Vertex>& g);

    friend size_t num_edges<>(const adj_list<Vertex>& g);

    friend Vertex add_vertex<>(adj_list<Vertex>& g, size_t n);

    template <class V, class Pred>
    friend void clear_vertex(V v, adj_list<V>& g, Pred&& pred);

    friend void remove_vertex<>(Vertex v, adj_list<Vertex>& g);

    friend void remove_vertex_fast<>(Vertex v, adj_list<Vertex>& g);

    friend std::pair<edge_descriptor, bool>
    add_edge<>(Vertex s, Vertex t, adj_list<Vertex>& g);

    friend void remove_edge<>(Vertex s, Vertex t, adj_list<Vertex>& g);

    friend void remove_edge<>(const edge_descriptor& e, adj_list<Vertex>& g);
};

//========================================================================
// Undirected version
//========================================================================

template <class Vertex = size_t>
using uadj_list = undirected_adaptor<adj_list<Vertex>>;

//========================================================================
// Graph traits and BGL scaffolding
//========================================================================

struct adj_list_traversal_tag
    : public vertex_list_graph_tag,
      public edge_list_graph_tag,
      public adjacency_graph_tag,
      public bidirectional_graph_tag,
      public adjacency_matrix_tag {};

template <class Vertex>
struct graph_traits<adj_list<Vertex> >
{
    typedef Vertex vertex_descriptor;
    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    typedef typename adj_list<Vertex>::edge_iterator edge_iterator;
    typedef typename adj_list<Vertex>::adjacency_iterator adjacency_iterator;
    typedef typename adj_list<Vertex>::in_adjacency_iterator in_adjacency_iterator;

    typedef typename adj_list<Vertex>::out_edge_iterator out_edge_iterator;
    typedef typename adj_list<Vertex>::in_edge_iterator in_edge_iterator;

    typedef typename adj_list<Vertex>::vertex_iterator vertex_iterator;

    typedef bidirectional_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef adj_list_traversal_tag traversal_category;

    typedef Vertex vertices_size_type;
    typedef Vertex edges_size_type;
    typedef size_t degree_size_type;

    static constexpr Vertex null_vertex() { return adj_list<Vertex>::null_vertex(); }

private:
    BOOST_STATIC_ASSERT((std::is_convertible<typename std::iterator_traits<out_edge_iterator>::iterator_category,
                                             std::random_access_iterator_tag>::value));
    BOOST_STATIC_ASSERT((std::is_convertible<typename std::iterator_traits<in_edge_iterator>::iterator_category,
                                             std::random_access_iterator_tag>::value));
    BOOST_STATIC_ASSERT((std::is_convertible<typename std::iterator_traits<adjacency_iterator>::iterator_category,
                                             std::random_access_iterator_tag>::value));
};

template <class Vertex>
struct graph_traits<const adj_list<Vertex> >
    : public graph_traits<adj_list<Vertex> >
{
};

template <class Vertex>
struct edge_property_type<adj_list<Vertex> >
{
    typedef void type;
};

template <class Vertex>
struct vertex_property_type<adj_list<Vertex> >
{
    typedef void type;
};

template <class Vertex>
struct graph_property_type<adj_list<Vertex> >
{
    typedef void type;
};

//========================================================================
// Graph access and manipulation functions
//========================================================================

template <class Vertex>
std::pair<typename adj_list<Vertex>::vertex_iterator,
          typename adj_list<Vertex>::vertex_iterator>
vertices(const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::vertex_iterator vi_t;
    return {vi_t(0), vi_t(g._edges.size())};
}


template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_iterator,
          typename adj_list<Vertex>::edge_iterator>
edges(const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::edge_list_t::const_iterator ei_t;
    typedef typename adj_list<Vertex>::vertex_list_t::const_iterator vi_t;
    ei_t ei_begin, ei_end;
    vi_t last_vi;
    if (g._edges.empty())
    {
        last_vi = g._edges.end();
    }
    else
    {
        ei_begin = g._edges[0].second.begin();
        last_vi = g._edges.end() - 1;
        ei_end = last_vi->second.begin() + last_vi->first;
    }
    typename adj_list<Vertex>::edge_iterator ebegin(g._edges.begin(),
                                                    g._edges.end(),
                                                    g._edges.begin(),
                                                    ei_begin);
    typename adj_list<Vertex>::edge_iterator eend(g._edges.begin(),
                                                  g._edges.end(),
                                                  last_vi,
                                                  ei_end);
    return {ebegin, eend};
}

template <class Vertex>
Vertex vertex(size_t i, const adj_list<Vertex>&)
{
    return i;
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
edge(Vertex s, Vertex t, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    if (g._keep_ehash)
    {
        auto& es = g._ehash[s];
        auto iter = es.find(t);
        if (iter != es.end())
            return {edge_descriptor(s, t, iter->second.front()), true};
        return {edge_descriptor(), false};
    }

    if (out_degree(s, g) < in_degree(t, g))
    {
        const auto& pes = g._edges[s];
        auto pos = pes.first;
        const auto& es = pes.second;
        auto end = es.begin() + pos;
        auto iter = std::find_if(es.begin(), end,
                                 [&](const auto& e) { return e.first == t; });
        if (iter != end)
            return {edge_descriptor(s, t, iter->second), true};
        return {edge_descriptor(), false};
    }
    else
    {
        const auto& pes = g._edges[t];
        auto pos = pes.first;
        const auto& es = pes.second;
        auto iter = std::find_if(es.begin() + pos, es.end(),
                                 [&](const auto& e) { return e.first == s; });
        if (iter != es.end())
            return {edge_descriptor(s, t, iter->second), true};
        return {edge_descriptor(), false};
    }
}

template <class Vertex, class F>
void edge_range_iter(Vertex s, Vertex t, const adj_list<Vertex>& g, F&& f)
{
    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    typedef decltype(f(edge_descriptor())) ret_t;
    if (g._keep_ehash)
    {
        auto& es = g._ehash[s];
        auto iter = es.find(t);
        if (iter != es.end())
        {
            for (auto idx : iter->second)
            {
                if constexpr (std::is_same_v<ret_t, bool>)
                {
                    if (!f(edge_descriptor(s, t, idx)))
                        break;
                }
                else
                {
                    f(edge_descriptor(s, t, idx));
                }
            }
        }
    }
    else
    {
        if (out_degree(s, g) < in_degree(t, g))
        {
            auto [iter, eend] = out_edges(s, g);
            for (; iter != eend; ++iter)
            {
                if (target(*iter, g) == t)
                {
                    if constexpr (std::is_same_v<ret_t, bool>)
                    {
                        if (!f(*iter))
                            break;
                    }
                    else
                    {
                        f(*iter);
                    }
                }
            }
        }
        else
        {
            auto [iter, eend] = in_edges(t, g);
            for (; iter != eend; ++iter)
            {
                if (source(*iter, g) == s)
                {
                    if constexpr (std::is_same_v<ret_t, bool>)
                    {
                        if (!f(*iter))
                            break;
                    }
                    else
                    {
                        f(*iter);
                    }
                }
            }
        }
    }
}

template <class Vertex>
size_t out_degree(Vertex v, const adj_list<Vertex>& g)
{
    const auto& pes = g._edges[v];
    return pes.first;
}

template <class Vertex>
size_t in_degree(Vertex v, const adj_list<Vertex>& g)
{
    const auto& pes = g._edges[v];
    auto pos = pes.first;
    auto& es = pes.second;
    return es.size() - pos;
}

template <class Vertex>
size_t degree(Vertex v, const adj_list<Vertex>& g)
{
    return g._edges[v].second.size();
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
out_edges(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::out_edge_iterator ei_t;
    const auto& pes = g._edges[v];
    auto pos = pes.first;
    auto& es = pes.second;
    return {ei_t(v, es.begin()), ei_t(v, es.begin() + pos)};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
in_edges(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::in_edge_iterator ei_t;
    const auto& pes = g._edges[v];
    auto pos = pes.first;
    auto& es = pes.second;
    return {ei_t(v, es.begin() + pos), ei_t(v, es.end())};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
_all_edges_out(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::out_edge_iterator ei_t;
    const auto& pes = g._edges[v];
    auto& es = pes.second;
    return {ei_t(v, es.begin()), ei_t(v, es.end())};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
_all_edges_in(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::in_edge_iterator ei_t;
    const auto& pes = g._edges[v];
    auto& es = pes.second;
    return {ei_t(v, es.begin()), ei_t(v, es.end())};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::all_edge_iterator,
          typename adj_list<Vertex>::all_edge_iterator>
all_edges(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::all_edge_iterator ei_t;
    const auto& pes = g._edges[v];
    auto& es = pes.second;
    auto pos = es.begin() + pes.first;
    return {ei_t(v, es.begin(), pos), ei_t(v, es.end(), pos)};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
out_neighbors(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::adjacency_iterator ai_t;
    const auto& pes = g._edges[v];
    auto pos = pes.first;
    auto& es = pes.second;
    return {ai_t(es.begin()), ai_t(es.begin() + pos)};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
in_neighbors(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::adjacency_iterator ai_t;
    const auto& pes = g._edges[v];
    auto pos = pes.first;
    auto& es = pes.second;
    return {ai_t(es.begin() + pos), ai_t(es.end())};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
all_neighbors(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::adjacency_iterator ai_t;
    const auto& pes = g._edges[v];
    auto& es = pes.second;
    return {ai_t(es.begin()), ai_t(es.end())};
}

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
adjacent_vertices(Vertex v, const adj_list<Vertex>& g)
{
    return out_neighbors(v, g);
}


template <class Vertex>
size_t num_vertices(const adj_list<Vertex>& g)
{
    return g._edges.size();
}

template <class Vertex>
size_t num_edges(const adj_list<Vertex>& g)
{
    return g._n_edges;
}


template <class Vertex>
typename std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
add_edge(Vertex s, Vertex t, adj_list<Vertex>& g)
{
    // get index from free list, if available
    Vertex idx = g.get_free_idx();

    // put target on back of source's out-list (middle of total list)
    auto& s_pes = g._edges[s];
    auto& s_pos = s_pes.first;
    auto& s_es = s_pes.second;

    if (s_pos < s_es.size())
    {
        //in-list is not empty: push first element to the back
        s_es.push_back(s_es[s_pos]);
        s_es[s_pos] = {t, idx};
        if (g._keep_epos)
            g.get_epos(s_es.back().second).second = s_es.size() - 1;
    }
    else
    {
        s_es.emplace_back(t, idx);
    }
    s_pos++;

    // put source on back of target's in-list
    auto& t_es = g._edges[t].second;
    t_es.emplace_back(s, idx);

    #pragma omp atomic
    g._n_edges++;

    if (g._keep_epos)
    {
        g.set_epos(idx, {s_pos - 1, t_es.size() - 1}); // out, in
        //g.check_epos();
    }

    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;

    edge_descriptor e(s, t, idx);

    if (g._keep_ehash)
        g.add_ehash(e);

    return {e, true};
}

template <class Vertex>
void remove_edge(Vertex s, Vertex t, adj_list<Vertex>& g)
{
    remove_edge(edge(s, t, g).first, g);
}

template <class Vertex>
void remove_edge(const typename adj_list<Vertex>::edge_descriptor& e,
                 adj_list<Vertex>& g)
{
    auto s = e.s;
    auto t = e.t;
    auto idx = e.idx;
    auto& s_pes = g._edges[s];
    auto& s_pos = s_pes.first;
    auto& s_es  = s_pes.second;
    auto& t_pes = g._edges[t];
    auto& t_pos = t_pes.first;
    auto& t_es  = t_pes.second;

    if (!g._keep_epos) // O(k_s + k_t)
    {
        // remove and shift
        auto remove_e = [&] (auto& elist, auto&& begin, auto&& end)
            {
                auto iter = std::find_if(begin, end,
                                         [&] (const auto& ei) -> bool
                                         { return idx == ei.second; });
                assert(iter != end);
                elist.erase(iter);
            };

        remove_e(s_es, s_es.begin(), s_es.begin() + s_pos);
        s_pos--;
        remove_e(t_es, t_es.begin() + t_pos, t_es.end());
    }
    else // O(1)
    {
        //g.check_epos();
        assert (idx < g._epos.size());

        // swap with back, and pop back
        auto remove_e = [&] (auto& elist, auto&& begin, auto&& end,
                             auto&& get_pos, bool swap_back)
        {
            auto back_iter = begin + ((end - begin) - 1);
            auto& back = *back_iter;
            auto j = get_pos(idx);
            assert(j < elist.size());
            assert(elist[j].second == idx);
            elist[j] = back;
            get_pos(back.second) = j;
            if (swap_back && end != elist.end())
            {
                // in-edge list; swap middle with very back
                back = elist.back();
                g._epos[back.second].second = back_iter - elist.begin();
            }
            elist.pop_back();
        };

        remove_e(s_es, s_es.begin(), s_es.begin() + s_pos,
                 [&](size_t i) -> auto& {return g.get_epos(i).first;}, true);
        s_pos--;
        remove_e(t_es, t_es.begin() + t_pos, t_es.end(),
                 [&](size_t i) -> auto& {return g.get_epos(i).second;}, false);

        //g.check_epos();
    }

    g.put_free_index(idx);

    #pragma omp atomic
    g._n_edges--;

    if (g._keep_ehash)
        g.remove_ehash(e);
}

template <class Vertex>
Vertex add_vertex(adj_list<Vertex>& g, size_t n)
{
    if (n == 1)
    {
        g._edges.emplace_back();
        if (g._keep_ehash)
            g._ehash.emplace_back();
        return g._edges.size() - 1;
    }
    else
    {
        auto v = g._edges.size();
        g._edges.resize(v + n);
        if (g._keep_ehash)
            g._ehash.resize(v + n);
        return v;
    }
}

template <class Vertex, class Pred>
void clear_vertex(Vertex v, adj_list<Vertex>& g, Pred&& pred)
{
    typename adj_list<Vertex>::make_out_edge mk_out_edge;
    typename adj_list<Vertex>::make_in_edge mk_in_edge;

    if (!g._keep_epos)
    {
        if (g._keep_ehash)
        {
            auto [iter, eend] = out_edges(v, g);
            for (; iter != eend; ++iter)
                g.remove_ehash(*iter);
        }

        auto& pos = g._edges[v].first;
        auto& es  = g._edges[v].second;
        for (size_t i = 0; i < es.size(); ++i)
        {
            auto u = es[i].first;
            auto idx = es[i].second;
            if (u == v)
            {
                if (i < pos && pred(mk_out_edge.def(v, es[i])))
                    g._free_idx.push_back(idx);
                continue;
            }
            auto& u_pos = g._edges[u].first;
            auto& u_es = g._edges[u].second;
            if (i < pos)
            {
                if (!pred(mk_out_edge.def(v, es[i])))
                    continue;
                auto iter =
                    std::remove_if(u_es.begin() + u_pos,
                                   u_es.end(),
                                   [&](auto& e)
                                   {
                                       return e.second == idx;
                                   });
                u_es.erase(iter, u_es.end());
            }
            else
            {
                if (!pred(mk_in_edge.def(v, es[i])))
                    continue;
                auto iter = std::remove_if(u_es.begin(),
                                           u_es.begin() + u_pos,
                                           [&](auto& e)
                                           {
                                               return e.second == idx;
                                           });
                u_es.erase(iter, u_es.begin() + u_pos);
                u_pos = iter - u_es.begin();
            }
            g._free_idx.push_back(idx);
        }
        auto iter =
            std::remove_if(es.begin() + pos, es.end(),
                           [&](auto& e)
                           {
                               return pred(mk_in_edge.def(v, e));
                           });
        size_t k = es.end() - iter;
        es.erase(iter, es.end());
        iter =
            std::remove_if(es.begin(), es.begin() + pos,
                           [&](auto& e)
                           {
                               return pred(mk_out_edge.def(v, e));
                           });
        k += std::count_if(iter, es.begin() + pos,
                           [&](auto& e){ return e.first != v; });
        es.erase(iter, es.begin() + pos);
        pos = iter - es.begin();
        g._n_edges -= k;
    }
    else
    {
        auto& pos = g._edges[v].first;
        auto& es = g._edges[v].second;
        std::vector<typename adj_list<Vertex>::edge_descriptor> res;
        res.reserve(es.size());
        for (size_t j = 0; j < es.size(); ++j)
        {
            auto& e = es[j];
            auto ed = (j < pos) ?
                mk_out_edge.def(v, e) :
                mk_in_edge.def(v, e);
            if (!pred(ed) || (j >= pos && e.first == v))
                continue;
            res.push_back(ed);
        }
        for (auto& ed : res)
            remove_edge(ed, g);
        //g.check_epos();
    }
}

template <class Vertex>
void clear_vertex(Vertex v, adj_list<Vertex>& g)
{
    clear_vertex(v, g, [](auto&&){ return true; });
}

// O(V + E)
template <class Vertex>
void remove_vertex(Vertex v, adj_list<Vertex>& g)
{
    clear_vertex(v, g);
    g._edges.erase(g._edges.begin() + v);

    size_t N = g._edges.size();

    if (v == N)
        return;

    #pragma omp parallel for schedule(runtime) if (N > get_openmp_min_thresh())
    for (size_t i = 0; i < N; ++i)
    {
        for (auto& e : g._edges[i].second)
        {
            if (e.first > v)
                e.first--;
        }
    }

    if (g._keep_ehash)
        g.rebuild_ehash();
}

// O(k + k_last)
template <class Vertex>
void remove_vertex_fast(Vertex v, adj_list<Vertex>& g)
{
    Vertex back = g._edges.size() - 1;

    clear_vertex(v, g);
    if (v < back)
    {
        if (g._keep_ehash)
        {
            auto [iter, eend] = out_edges(back, g);
            for (; iter != eend; ++iter)
                g.remove_ehash(*iter);
        }

        g._edges[v] = g._edges[back];

        auto pos = g._edges[v].first;
        auto& es = g._edges[v].second;
        for (size_t i = 0; i < es.size(); ++i)
        {
            auto& eu = es[i];
            Vertex u = eu.first;
            if (u == back)
            {
                eu.first = v; // self-loop
            }
            else
            {
                // change label in adjacent lists
                if (!g._keep_epos)
                {
                    auto u_pos = g._edges[u].first;
                    auto& u_es = g._edges[u].second;
                    size_t begin = (i < pos) ? u_pos : 0;
                    size_t end = (i < pos) ? u_es.size() : u_pos;
                    for (size_t j = begin; j < end; ++j)
                    {
                        auto& e = u_es[j];
                        if (e.first == back)
                            e.first = v;
                    }
                }
                else
                {
                    size_t idx = eu.second;
                    auto u_pos = (i < pos) ?
                        g._epos[idx].second :
                        g._epos[idx].first;
                    assert(g._edges[u].second[u_pos].first == back);
                    g._edges[u].second[u_pos].first = v;
                }
            }
        }

        if (g._keep_ehash)
        {
            auto [iter, eend] = out_edges(v, g);
            for (; iter != eend; ++iter)
                g.add_ehash(*iter);
        }
    }
    g._edges.pop_back();
    if (g._keep_ehash)
        g._ehash.pop_back();
}

template <class Vertex>
Vertex source(const typename adj_list<Vertex>::edge_descriptor& e,
              const adj_list<Vertex>&)
{
    return e.s;
}

template <class Vertex>
Vertex target(const typename adj_list<Vertex>::edge_descriptor& e,
              const adj_list<Vertex>&)
{
    return e.t;
}

//========================================================================
// Vertex and edge index property maps
//========================================================================

template <class Vertex>
struct property_map<adj_list<Vertex>, vertex_index_t>
{
    typedef identity_property_map type;
    typedef type const_type;
};

template <class Vertex>
struct property_map<const adj_list<Vertex>, vertex_index_t>
{
    typedef identity_property_map type;
    typedef type const_type;
};

template <class Vertex>
identity_property_map
get(vertex_index_t, adj_list<Vertex>&)
{
    return identity_property_map();
}

template <class Vertex>
identity_property_map
get(vertex_index_t, const adj_list<Vertex>&)
{
    return identity_property_map();
}

template<class Vertex>
class adj_edge_index_property_map:
    public put_get_helper<Vertex, adj_edge_index_property_map<Vertex> >
{
public:
    typedef typename adj_list<Vertex>::edge_descriptor key_type;
    typedef Vertex reference;
    typedef Vertex value_type;
    typedef boost::readable_property_map_tag category;

    reference operator[](const key_type& k) const {return k.idx;}
};

template <class Vertex>
struct property_map<adj_list<Vertex>, edge_index_t>
{
    typedef adj_edge_index_property_map<Vertex> type;
    typedef type const_type;

};

template <class Vertex>
adj_edge_index_property_map<Vertex>
get(edge_index_t, const adj_list<Vertex>&)
{
    return adj_edge_index_property_map<Vertex>();
}

// hashing of edge descriptors

namespace detail
{
template <class Vertex>
size_t hash_value(const adj_edge_descriptor<Vertex>& e)
{
    return e.idx;
}
} // namespace detail
} // namespace boost

namespace std
{
template <class Vertex>
struct hash<boost::detail::adj_edge_descriptor<Vertex>>
{
    template <class Edge>
    std::size_t operator()(Edge const& e) const
    {
        return boost::detail::hash_value(e);
    }
};
}  // namespace std

#endif //GRAPH_ADJACENCY_HH
