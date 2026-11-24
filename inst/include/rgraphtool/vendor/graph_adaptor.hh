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

#ifndef GRAPH_ADAPTOR_HH
#define GRAPH_ADAPTOR_HH

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>

#include "transform_iterator.hh"

namespace boost {

//==============================================================================
// undirected_adaptor
// This class encapsulates a directed graph with parallel edges and provides a
// view of the graph as undirected with parallel edges.
//==============================================================================
template <class Graph>
class undirected_adaptor
    : private Graph
{
public:
    template <class... Ts>
    undirected_adaptor(Ts&&... args) : Graph(args...) {}

    operator Graph() const = delete;

    typedef Graph original_graph_t;
    typedef typename graph_traits<Graph>::out_edge_iterator all_edge_iterator;

    Graph& original_graph() { return reinterpret_cast<Graph&>(*this); }
    const Graph& original_graph() const { return reinterpret_cast<const Graph&>(*this); }
};

template <class Graph>
auto& get_undirected_adaptor(Graph& g)
{
    if constexpr (std::is_convertible_v<typename boost::graph_traits<Graph>::directed_category,
                                        boost::directed_tag>)
        return reinterpret_cast<undirected_adaptor<Graph>&>(g);
    else
        return g;
}

template <class Graph>
auto& get_undirected_adaptor(const Graph& g)
{
    if constexpr (std::is_convertible_v<typename boost::graph_traits<Graph>::directed_category,
                                        boost::directed_tag>)
        return reinterpret_cast<const undirected_adaptor<Graph>&>(g);
    else
        return g;
}

//==============================================================================
// graph_traits<undirected_adaptor>
// this defines all the necessary types associated with undirected_adaptor
//==============================================================================
template <class Graph>
struct graph_traits<undirected_adaptor<Graph> >
    : public graph_traits<Graph>
{
    typedef undirected_tag directed_category;
};

template <class Graph>
struct graph_traits<const undirected_adaptor<Graph>>
    : public graph_traits<undirected_adaptor<Graph>> {};

//==============================================================================
// Nonmember functions
// these provide manipulation of the graph
//==============================================================================

//==============================================================================
// source(e,g)
//==============================================================================
template <class Graph>
[[gnu::pure]]
auto
source(const typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor& e,
       const undirected_adaptor<Graph>& g)
{
    return source(e, g.original_graph());
}

//==============================================================================
// target(e,g)
//==============================================================================
template <class Graph>
[[gnu::pure]]
auto
target(const typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor& e,
       const undirected_adaptor<Graph>& g)
{
    return target(e, g.original_graph());
}

//==============================================================================
// vertex(n,g)
//==============================================================================
template <class Graph>
[[gnu::pure]]
auto
vertex(typename graph_traits<undirected_adaptor<Graph> >::vertices_size_type n,
       const undirected_adaptor<Graph>& g)
{
    return vertex(n, g.original_graph());
}

//==============================================================================
// vertices(g)
//==============================================================================
template <class Graph>
auto
vertices(const undirected_adaptor<Graph>& g)
{
    return vertices(g.original_graph());
}

//==============================================================================
// edges(g)
//==============================================================================
template <class Graph>
auto
edges(const undirected_adaptor<Graph>& g)
{
    return edges(g.original_graph());
}

//==============================================================================
// edge(u, v, g)
//==============================================================================
template <class Graph>
auto
edge(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
     typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
     const undirected_adaptor<Graph>& g)
{
    auto res = edge(u, v, g.original_graph());

    if (!res.second)
    {
        res = edge(v, u, g.original_graph());
        std::swap(res.first.s, res.first.t);
    }

    return res;
}

template <class Vertex, class Graph, class F>
void edge_range_iter(Vertex s, Vertex t, const undirected_adaptor<Graph>& g, F&& f)
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef decltype(f(edge_t())) ret_t;

    bool ret = true;
    edge_range_iter(s, t, g.original_graph(),
                    [&](const auto& e)
                    {
                        if constexpr (std::is_same_v<ret_t, bool>)
                            return ret = f(e);
                        else
                            f(e);
                    });
    if constexpr (std::is_same_v<ret_t, bool>)
    {
        if (!ret)
            return;
    }
    edge_range_iter(t, s, g.original_graph(),
                    [&](auto e)
                    {
                        std::swap(e.s, e.t);
                        if constexpr (std::is_same_v<ret_t, bool>)
                            return f(e);
                        else
                            f(e);
                    });
}

//==============================================================================
// out_edges(u,g)
//==============================================================================
template <class Graph>
auto
out_edges(typename graph_traits<undirected_adaptor<Graph>>::vertex_descriptor u,
          const undirected_adaptor<Graph>& g)
{
    return _all_edges_out(u, g.original_graph());
}

template <class Graph>
auto
_all_edges_out(typename graph_traits<undirected_adaptor<Graph>>::vertex_descriptor u,
               const undirected_adaptor<Graph>& g)
{
    return out_edges(u, g);
}

//==============================================================================
// in_edges(u,g)
//==============================================================================
template <class Graph>
auto
in_edges(typename graph_traits<undirected_adaptor<Graph>>::vertex_descriptor,
         const undirected_adaptor<Graph>&)
{
    typedef typename graph_traits<undirected_adaptor<Graph>>::in_edge_iterator
        iter_t;
    return std::make_pair(iter_t(), iter_t());
}

template <class Graph>
auto
_all_edges_in(typename graph_traits<undirected_adaptor<Graph>>::vertex_descriptor u,
               const undirected_adaptor<Graph>& g)
{
    return _all_edges_in(u, g.original_graph());
}

template <class Graph>
auto
all_edges(typename graph_traits<undirected_adaptor<Graph>>::vertex_descriptor u,
          const undirected_adaptor<Graph>& g)
{
    return out_edges(u, g);
}

//==============================================================================
// out_neighbors(u, g)
//==============================================================================
template <class Graph>
auto
out_neighbors(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
               const undirected_adaptor<Graph>& g)
{
    return all_neighbors(u, g.original_graph());
}

//==============================================================================
// in_neighbors(u, g)
//==============================================================================
template <class Graph>
auto
in_neighbors(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
              const undirected_adaptor<Graph>& g)
{
    return out_neighbors(u, g);
}

//==============================================================================
// all_neighbors(u, g)
//==============================================================================
template <class Graph>
auto
all_neighbors(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
               const undirected_adaptor<Graph>& g)
{
    return out_neighbors(u, g);
}

//==============================================================================
// adjacent_vertices(u,g)
//==============================================================================
template <class Graph>
auto
adjacent_vertices
    (typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
     const undirected_adaptor<Graph>& g)
{
    return out_neighbors(u, g);
}

//==============================================================================
// num_vertices(g)
//==============================================================================
template <class Graph>
auto
num_vertices(const undirected_adaptor<Graph>& g)
{
    return num_vertices(g.original_graph());
}

//==============================================================================
// num_edges(g)
//==============================================================================
template <class Graph>
auto
num_edges(const undirected_adaptor<Graph>& g)
{
    return num_edges(g.original_graph());
}

//==============================================================================
// out_degree(u,g)
//==============================================================================
template <class Graph>
auto
out_degree(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
           const undirected_adaptor<Graph>& g)
{
    return degree(u, g.original_graph());
}

//==============================================================================
// in_degree(u,g)
//==============================================================================
template <class Graph>
auto
in_degree(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor,
          const undirected_adaptor<Graph>&)
{
    return 0;
}

//==============================================================================
// degree(u,g)
//==============================================================================
template <class Graph>
auto
degree(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
       const undirected_adaptor<Graph>& g)
{
    return out_degree(u, g);
}


//==============================================================================
// add_vertex(g)
//==============================================================================
template <class Graph>
auto
add_vertex(undirected_adaptor<Graph>& g, size_t n = 1)
{
    return add_vertex(g.original_graph(), n);
}

//==============================================================================
// add_vertex(vp,g)
//==============================================================================
template <class Graph, class VertexProperties>
auto
add_vertex(const VertexProperties& p, undirected_adaptor<Graph>& g)
{
    return add_vertex(p, g.original_graph());
}

//==============================================================================
// clear_vertex(u,g)
//==============================================================================
template <class Graph>
void clear_vertex(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
                  undirected_adaptor<Graph>& g)
{
    clear_vertex(u, g.original_graph());
}

//==============================================================================
// clear_vertex(u,g,pred)
//==============================================================================
template <class Graph, class Pred>
void clear_vertex(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
                  undirected_adaptor<Graph>& g, Pred&& pred)
{
    clear_vertex(u, g.original_graph(), pred);
}

//==============================================================================
// remove_vertex(u,g)
//==============================================================================
template <class Graph>
void remove_vertex(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
                   undirected_adaptor<Graph>& g)
{
    remove_vertex(u, g.original_graph());
}

//==============================================================================
// remove_vertex_fast(u,g)
//==============================================================================
template <class Graph>
void remove_vertex_fast(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
                        undirected_adaptor<Graph>& g)
{
    remove_vertex_fast(u, g.original_graph());
}

//==============================================================================
// add_edge(u,v,g)
//==============================================================================
template <class Graph>
std::pair<typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
         undirected_adaptor<Graph>& g)
{
    return add_edge(u, v, g.original_graph());
}

//==============================================================================
// add_edge(u,v,ep,g)
//==============================================================================
template <class Graph, class EdgeProperties>
auto
add_edge(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
         const EdgeProperties& ep, undirected_adaptor<Graph>& g)
{
    return add_edge(u, v, ep, g.original_graph());
}

//==============================================================================
// remove_edge(u,v,g)
//==============================================================================
template <class Graph>
void remove_edge(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor u,
                 typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
                 undirected_adaptor<Graph>& g)
{
    auto e = edge(u, v, g);
    if (e.second)
        remove_edge(e.first, g);
}

//==============================================================================
// remove_edge(e,g)
//==============================================================================
template <class Graph>
void remove_edge(typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor e,
                 undirected_adaptor<Graph>& g)
{
    auto& u = g.original_graph();
    u.reverse_edge(e);
    remove_edge(e, u);
}

//==============================================================================
// remove_edge(e_iter,g)
//==============================================================================
template <class Graph>
void remove_edge(const typename graph_traits<undirected_adaptor<Graph> >::out_edge_iterator& iter,
                 undirected_adaptor<Graph>& g)
{
    remove_edge(*iter, g);
}

//==============================================================================
// remove_out_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
void remove_out_edge_if(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
                        Predicate predicate, undirected_adaptor<Graph>& g)
{
    std::vector<typename undirected_adaptor<Graph>::EdgeDescriptor> removed_edges;
    auto edge_range = out_edges(v,g);
    for (auto iter = edge_range.first; iter != edge_range.second; ++iter)
        if (predicate(*iter))
            removed_edges.push_back(*iter);
    for (auto& e : removed_edges)
        remove_edge(e, g);
}

//==============================================================================
// remove_in_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
void remove_in_edge_if(typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
                       Predicate predicate, undirected_adaptor<Graph>& g)
{
    std::vector<typename undirected_adaptor<Graph>::EdgeDescriptor> removed_edges;
    auto edge_range = in_edges(v,g);
    for (auto iter = edge_range.first; iter != edge_range.second; ++iter)
        if (predicate(*iter))
            removed_edges.push_back(*iter);
    for (auto& e : removed_edges)
        remove_edge(e, g);
}


//==============================================================================
// Property maps
//==============================================================================

//==============================================================================
// vertex_property<undirected_adaptor>
//==============================================================================
template <class Graph>
class vertex_property<undirected_adaptor<Graph> >
{
public:
    typedef typename vertex_property<Graph>::type type;
};

//==============================================================================
// vertex_property_type<undirected_adaptor>
//==============================================================================
template <class Graph>
class vertex_property_type<undirected_adaptor<Graph> >
{
public:
    typedef typename vertex_property_type<Graph>::type type;
};

//==============================================================================
// edge_property<undirected_adaptor>
//==============================================================================
template <class Graph>
class edge_property<undirected_adaptor<Graph> >
{
public:
    typedef typename edge_property<Graph>::type type;
};

//==============================================================================
// edge_property_type<undirected_adaptor>
//==============================================================================
template <class Graph>
class edge_property_type<undirected_adaptor<Graph> >
{
public:
    typedef typename edge_property_type<Graph>::type type;
};

//==============================================================================
// property_map<UndirecterdAdaptor, PropertyTag>
//==============================================================================
template <class Graph, class PropertyTag>
class property_map<undirected_adaptor<Graph>, PropertyTag>
{
public:
    typedef typename property_map<Graph, PropertyTag>::type type;
    typedef typename property_map<Graph, PropertyTag>::const_type const_type;
};

//==============================================================================
// property_map<undirected_adaptor, T Bundle::*>
//==============================================================================
template <typename Graph, typename T, typename Bundle>
class property_map<undirected_adaptor<Graph>, T Bundle::*>
{
public:
    typedef typename property_map<Graph, T Bundle::*>::type type;
    typedef typename property_map<Graph, T Bundle::*>::const_type const_type;
};


//==============================================================================
// get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
typename property_map<undirected_adaptor<Graph>, PropertyTag>::type
get(PropertyTag tag, undirected_adaptor<Graph>& g)
{
    return get(tag, g.original_graph());
}

//==============================================================================
// const get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
typename property_map<undirected_adaptor<Graph>, PropertyTag>::const_type
get(PropertyTag tag, const undirected_adaptor<Graph>& g)
{
    return get(tag, g.original_graph());
}

//==============================================================================
// get(tag,g,v)
//==============================================================================
template <class PropertyTag, class Graph>
typename property_traits
    <typename property_map<undirected_adaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const undirected_adaptor<Graph>& g,
    typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v)
{
    return get(tag, g.original_graph(), v);
}

//==============================================================================
// get(tag,g,e)
//==============================================================================
template <class PropertyTag, class Graph>
typename property_traits
    <typename property_map<undirected_adaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const undirected_adaptor<Graph>& g,
    const typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor& e)
{
    return get(tag, g.original_graph(), e);
}

//==============================================================================
// put(tag, g, v, value)
//==============================================================================
template <class Graph, class PropertyTag, class Value>
void put(PropertyTag tag, undirected_adaptor<Graph>& g,
         typename graph_traits<undirected_adaptor<Graph> >::vertex_descriptor v,
         const Value& value)
{
    put(tag, g.original_graph(), v, value);
}

//==============================================================================
// put(tag, g, e, value)
//==============================================================================
template <class Graph, class PropertyTag, class X, class Value>
void put(PropertyTag tag, const undirected_adaptor<Graph>& g,
         const typename graph_traits<undirected_adaptor<Graph> >::edge_descriptor& e,
         const Value &value)
{
    put(tag, g.original_graph(), e, value);
}

//==============================================================================
// get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(undirected_adaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.original_graph(), tag);
}

//==============================================================================
// const get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
const typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(const undirected_adaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.original_graph(), tag);
}

} // namespace boost


#endif // GRAPH_ADAPTOR_HH
