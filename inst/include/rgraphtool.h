// rgraphtool: Public C++ extension interface for R/C++ user extensions and
// downstream packages (LinkingTo: rgraphtool, Rcpp, BH).
//
// Covers zero-copy graph access, range-based vertex/edge/neighbor iteration,
// vertex and edge property map views, and checked vs. unchecked access
// (analogous to graph-tool's PropertyMap and get_unchecked() C++ API).
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool adjacency & design)
// Copyright (C) 2026 Diogo Melo (R/C++ port and extension API)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#ifndef RGRAPHTOOL_PUBLIC_H
#define RGRAPHTOOL_PUBLIC_H

#include <Rcpp.h>
#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rgraphtool/vendor/graph_adjacency.hh"

namespace rgraphtool {

using Adjacency = boost::adj_list<size_t>;
using Vertex = boost::graph_traits<Adjacency>::vertex_descriptor;
using Edge = boost::graph_traits<Adjacency>::edge_descriptor;

struct Graph {
    Adjacency adjacency;
    bool directed;
    std::map<std::string, std::vector<double>> edge_props;
    explicit Graph(bool directed_) : directed(directed_) {}
};

// Extract Graph& from either an rgraphtool_graph R list (e.g. g) or its
// external pointer (g$pointer).
inline Graph& get_graph(SEXP obj) {
    SEXP ptr = obj;
    if (TYPEOF(obj) == VECSXP) {
        Rcpp::List lst(obj);
        if (!lst.containsElementNamed("pointer")) {
            Rcpp::stop("Expected an rgraphtool_graph object or its external pointer");
        }
        ptr = lst["pointer"];
    }
    if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrAddr(ptr) == nullptr) {
        Rcpp::stop("Graph pointer is invalid; graphs cannot be restored from serialized R objects");
    }
    return *Rcpp::XPtr<Graph>(ptr);
}

// ---------------------------------------------------------------------------
// Graph inspection helpers
// ---------------------------------------------------------------------------

inline size_t num_vertices(const Graph& g) {
    return boost::num_vertices(g.adjacency);
}

inline size_t num_edges(const Graph& g) {
    return boost::num_edges(g.adjacency);
}

inline bool is_directed(const Graph& g) {
    return g.directed;
}

inline Vertex source(const Edge& e, const Graph& g) {
    return boost::source(e, g.adjacency);
}

inline Vertex target(const Edge& e, const Graph& g) {
    return boost::target(e, g.adjacency);
}

inline size_t edge_index(const Edge& e) {
    return e.idx;
}

inline Vertex checked_vertex(size_t v, const Graph& g) {
    if (v >= num_vertices(g)) {
        throw std::out_of_range("Vertex index out of bounds: " + std::to_string(v) +
                                " (num_vertices = " + std::to_string(num_vertices(g)) + ")");
    }
    return v;
}

inline size_t out_degree(Vertex v, const Graph& g) {
    if (g.directed) {
        return boost::out_degree(v, g.adjacency);
    }
    return boost::out_degree(v, g.adjacency) + boost::in_degree(v, g.adjacency);
}

inline size_t in_degree(Vertex v, const Graph& g) {
    if (g.directed) {
        return boost::in_degree(v, g.adjacency);
    }
    return boost::out_degree(v, g.adjacency) + boost::in_degree(v, g.adjacency);
}

inline size_t total_degree(Vertex v, const Graph& g) {
    return boost::out_degree(v, g.adjacency) + boost::in_degree(v, g.adjacency);
}

// ---------------------------------------------------------------------------
// Range-based graph iteration helpers
// ---------------------------------------------------------------------------

template <typename Iterator>
class IteratorRange {
public:
    IteratorRange(Iterator first, Iterator last) : first_(first), last_(last) {}
    explicit IteratorRange(std::pair<Iterator, Iterator> p)
        : first_(p.first), last_(p.second) {}

    Iterator begin() const { return first_; }
    Iterator end() const { return last_; }

private:
    Iterator first_;
    Iterator last_;
};

inline IteratorRange<Adjacency::vertex_iterator>
vertices_range(const Graph& g) {
    return IteratorRange<Adjacency::vertex_iterator>(boost::vertices(g.adjacency));
}

inline IteratorRange<Adjacency::edge_iterator>
edges_range(const Graph& g) {
    return IteratorRange<Adjacency::edge_iterator>(boost::edges(g.adjacency));
}

// For directed graphs, iterates outgoing edges (source == v).
// For undirected graphs, iterates all incident edges oriented with source == v.
inline IteratorRange<Adjacency::out_edge_iterator>
out_edges_range(Vertex v, const Graph& g) {
    if (g.directed) {
        return IteratorRange<Adjacency::out_edge_iterator>(boost::out_edges(v, g.adjacency));
    }
    return IteratorRange<Adjacency::out_edge_iterator>(boost::_all_edges_out(v, g.adjacency));
}

// For directed graphs, iterates incoming edges (target == v).
// For undirected graphs, iterates all incident edges oriented with target == v.
inline IteratorRange<Adjacency::in_edge_iterator>
in_edges_range(Vertex v, const Graph& g) {
    if (g.directed) {
        return IteratorRange<Adjacency::in_edge_iterator>(boost::in_edges(v, g.adjacency));
    }
    return IteratorRange<Adjacency::in_edge_iterator>(boost::_all_edges_in(v, g.adjacency));
}

inline IteratorRange<Adjacency::all_edge_iterator>
all_edges_range(Vertex v, const Graph& g) {
    return IteratorRange<Adjacency::all_edge_iterator>(boost::all_edges(v, g.adjacency));
}

inline IteratorRange<Adjacency::adjacency_iterator>
out_neighbors_range(Vertex v, const Graph& g) {
    if (g.directed) {
        return IteratorRange<Adjacency::adjacency_iterator>(boost::out_neighbors(v, g.adjacency));
    }
    return IteratorRange<Adjacency::adjacency_iterator>(boost::all_neighbors(v, g.adjacency));
}

inline IteratorRange<Adjacency::adjacency_iterator>
in_neighbors_range(Vertex v, const Graph& g) {
    if (g.directed) {
        return IteratorRange<Adjacency::adjacency_iterator>(boost::in_neighbors(v, g.adjacency));
    }
    return IteratorRange<Adjacency::adjacency_iterator>(boost::all_neighbors(v, g.adjacency));
}

inline IteratorRange<Adjacency::adjacency_iterator>
all_neighbors_range(Vertex v, const Graph& g) {
    return IteratorRange<Adjacency::adjacency_iterator>(boost::all_neighbors(v, g.adjacency));
}

// ---------------------------------------------------------------------------
// Vertex Property Maps: Checked (bounds-checked) and Unchecked (raw pointer)
// ---------------------------------------------------------------------------

template <typename T>
class UncheckedVertexPropertyMap {
public:
    UncheckedVertexPropertyMap(T* data, size_t size) : data_(data), size_(size) {}

    inline T& operator[](Vertex v) noexcept { return data_[v]; }
    inline const T& operator[](Vertex v) const noexcept { return data_[v]; }

    inline size_t size() const noexcept { return size_; }
    inline T* data() noexcept { return data_; }
    inline const T* data() const noexcept { return data_; }

private:
    T* data_;
    size_t size_;
};

template <typename T>
class VertexPropertyMap {
public:
    // Non-owning view over external buffer (e.g. Rcpp vector memory)
    VertexPropertyMap(T* data, size_t size, size_t expected_vertices)
        : owned_(nullptr), data_(data), size_(size), expected_(expected_vertices) {}

    // Owning vector-backed property map
    explicit VertexPropertyMap(size_t n, T init = T{})
        : owned_(std::make_shared<std::vector<T>>(n, init)),
          data_(owned_->data()), size_(n), expected_(n) {}

    inline T& operator[](Vertex v) {
        if (v >= size_) {
            throw std::out_of_range("VertexPropertyMap index out of bounds: " +
                                    std::to_string(v) + " >= " + std::to_string(size_));
        }
        return data_[v];
    }

    inline const T& operator[](Vertex v) const {
        if (v >= size_) {
            throw std::out_of_range("VertexPropertyMap index out of bounds: " +
                                    std::to_string(v) + " >= " + std::to_string(size_));
        }
        return data_[v];
    }

    inline T& at(Vertex v) { return (*this)[v]; }
    inline const T& at(Vertex v) const { return (*this)[v]; }

    // Validate capacity once, then return zero-overhead unchecked accessor.
    UncheckedVertexPropertyMap<T> get_unchecked() {
        if (size_ < expected_) {
            throw std::out_of_range(
                "Cannot create UncheckedVertexPropertyMap: buffer length (" +
                std::to_string(size_) + ") is smaller than num_vertices (" +
                std::to_string(expected_) + ")");
        }
        return UncheckedVertexPropertyMap<T>(data_, size_);
    }

    inline size_t size() const noexcept { return size_; }
    inline T* data() noexcept { return data_; }
    inline const T* data() const noexcept { return data_; }

private:
    std::shared_ptr<std::vector<T>> owned_;
    T* data_;
    size_t size_;
    size_t expected_;
};

// ---------------------------------------------------------------------------
// Edge Property Maps: Checked and Unchecked (mutable and const)
// ---------------------------------------------------------------------------

template <typename T>
class UncheckedEdgePropertyMap {
public:
    UncheckedEdgePropertyMap(T* data, size_t size) : data_(data), size_(size) {}

    inline T& operator[](const Edge& e) noexcept { return data_[e.idx]; }
    inline const T& operator[](const Edge& e) const noexcept { return data_[e.idx]; }
    inline T& operator[](size_t idx) noexcept { return data_[idx]; }
    inline const T& operator[](size_t idx) const noexcept { return data_[idx]; }

    inline size_t size() const noexcept { return size_; }
    inline T* data() noexcept { return data_; }
    inline const T* data() const noexcept { return data_; }

private:
    T* data_;
    size_t size_;
};

template <typename T>
class EdgePropertyMap {
public:
    EdgePropertyMap(T* data, size_t size, size_t expected_edges)
        : owned_(nullptr), data_(data), size_(size), expected_(expected_edges) {}

    explicit EdgePropertyMap(size_t m, T init = T{})
        : owned_(std::make_shared<std::vector<T>>(m, init)),
          data_(owned_->data()), size_(m), expected_(m) {}

    inline T& operator[](const Edge& e) {
        return (*this)[e.idx];
    }

    inline const T& operator[](const Edge& e) const {
        return (*this)[e.idx];
    }

    inline T& operator[](size_t idx) {
        if (idx >= size_) {
            throw std::out_of_range("EdgePropertyMap index out of bounds: " +
                                    std::to_string(idx) + " >= " + std::to_string(size_));
        }
        return data_[idx];
    }

    inline const T& operator[](size_t idx) const {
        if (idx >= size_) {
            throw std::out_of_range("EdgePropertyMap index out of bounds: " +
                                    std::to_string(idx) + " >= " + std::to_string(size_));
        }
        return data_[idx];
    }

    inline T& at(const Edge& e) { return (*this)[e]; }
    inline const T& at(const Edge& e) const { return (*this)[e]; }

    UncheckedEdgePropertyMap<T> get_unchecked() {
        if (size_ < expected_) {
            throw std::out_of_range(
                "Cannot create UncheckedEdgePropertyMap: buffer length (" +
                std::to_string(size_) + ") is smaller than num_edges (" +
                std::to_string(expected_) + ")");
        }
        return UncheckedEdgePropertyMap<T>(data_, size_);
    }

    inline size_t size() const noexcept { return size_; }
    inline T* data() noexcept { return data_; }
    inline const T* data() const noexcept { return data_; }

private:
    std::shared_ptr<std::vector<T>> owned_;
    T* data_;
    size_t size_;
    size_t expected_;
};

template <typename T>
class UncheckedConstEdgePropertyMap {
public:
    UncheckedConstEdgePropertyMap(const T* data, size_t size) : data_(data), size_(size) {}

    inline const T& operator[](const Edge& e) const noexcept { return data_[e.idx]; }
    inline const T& operator[](size_t idx) const noexcept { return data_[idx]; }

    inline size_t size() const noexcept { return size_; }
    inline const T* data() const noexcept { return data_; }

private:
    const T* data_;
    size_t size_;
};

template <typename T>
class ConstEdgePropertyMap {
public:
    ConstEdgePropertyMap(const T* data, size_t size, size_t expected_edges)
        : data_(data), size_(size), expected_(expected_edges) {}

    inline const T& operator[](const Edge& e) const { return (*this)[e.idx]; }

    inline const T& operator[](size_t idx) const {
        if (idx >= size_) {
            throw std::out_of_range("ConstEdgePropertyMap index out of bounds: " +
                                    std::to_string(idx) + " >= " + std::to_string(size_));
        }
        return data_[idx];
    }

    inline const T& at(const Edge& e) const { return (*this)[e]; }

    UncheckedConstEdgePropertyMap<T> get_unchecked() const {
        if (size_ < expected_) {
            throw std::out_of_range(
                "Cannot create UncheckedConstEdgePropertyMap: buffer length (" +
                std::to_string(size_) + ") is smaller than num_edges (" +
                std::to_string(expected_) + ")");
        }
        return UncheckedConstEdgePropertyMap<T>(data_, size_);
    }

    inline size_t size() const noexcept { return size_; }
    inline const T* data() const noexcept { return data_; }

private:
    const T* data_;
    size_t size_;
    size_t expected_;
};

// ---------------------------------------------------------------------------
// Property Map factories
// ---------------------------------------------------------------------------

template <typename T>
inline VertexPropertyMap<T> make_vertex_property(const Graph& g, T init = T{}) {
    return VertexPropertyMap<T>(num_vertices(g), init);
}

template <typename T>
inline EdgePropertyMap<T> make_edge_property(const Graph& g, T init = T{}) {
    return EdgePropertyMap<T>(num_edges(g), init);
}

inline VertexPropertyMap<double>
vertex_property_view(Rcpp::NumericVector& vec, const Graph& g) {
    return VertexPropertyMap<double>(
        vec.size() > 0 ? REAL(vec) : nullptr,
        static_cast<size_t>(vec.size()),
        num_vertices(g)
    );
}

inline VertexPropertyMap<int>
vertex_property_view(Rcpp::IntegerVector& vec, const Graph& g) {
    return VertexPropertyMap<int>(
        vec.size() > 0 ? INTEGER(vec) : nullptr,
        static_cast<size_t>(vec.size()),
        num_vertices(g)
    );
}

inline EdgePropertyMap<double>
edge_property_view(Graph& g, const std::string& name) {
    auto it = g.edge_props.find(name);
    if (it == g.edge_props.end()) {
        throw std::out_of_range("Edge property '" + name + "' does not exist on graph");
    }
    auto& vec = it->second;
    return EdgePropertyMap<double>(vec.data(), vec.size(), num_edges(g));
}

inline ConstEdgePropertyMap<double>
edge_property_view(const Graph& g, const std::string& name) {
    auto it = g.edge_props.find(name);
    if (it == g.edge_props.end()) {
        throw std::out_of_range("Edge property '" + name + "' does not exist on graph");
    }
    const auto& vec = it->second;
    return ConstEdgePropertyMap<double>(vec.data(), vec.size(), num_edges(g));
}

} // namespace rgraphtool

#endif // RGRAPHTOOL_PUBLIC_H
