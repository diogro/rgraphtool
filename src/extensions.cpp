// rgraphtool: Built-in reference algorithms using the public <rgraphtool.h>
// C++ extension API (graph iteration, property maps, checked vs. unchecked).
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool k-core algorithm)
// Copyright (C) 2026 Diogo Melo (R/C++ port and extension API)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#include <rgraphtool.h>
#include <string>
#include <vector>

namespace {

template <typename CoreMap, typename DegMap, typename PosMap>
void kcore_impl(const rgraphtool::Graph& g,
                const std::string& deg_mode,
                CoreMap& core_map,
                DegMap& deg,
                PosMap& pos) {
    using rgraphtool::Vertex;
    std::vector<std::vector<Vertex>> bins;

    for (Vertex v : rgraphtool::vertices_range(g)) {
        size_t k = 0;
        if (deg_mode == "out" && rgraphtool::is_directed(g)) {
            k = rgraphtool::out_degree(v, g);
        } else if (deg_mode == "in" && rgraphtool::is_directed(g)) {
            k = rgraphtool::in_degree(v, g);
        } else {
            k = rgraphtool::total_degree(v, g);
        }
        deg[v] = k;
        if (k >= bins.size()) {
            bins.resize(k + 1);
        }
        pos[v] = bins[k].size();
        bins[k].push_back(v);
    }

    auto relax_neighbor = [&](Vertex v, Vertex u) {
        if (deg[u] > deg[v]) {
            size_t ku = deg[u];
            Vertex w = bins[ku].back();
            size_t pu = pos[u];
            bins[ku][pu] = w;
            pos[w] = pu;
            bins[ku].pop_back();
            pos[u] = bins[ku - 1].size();
            bins[ku - 1].push_back(u);
            deg[u] = ku - 1;
        }
    };

    for (size_t k = 0; k < bins.size(); ++k) {
        while (!bins[k].empty()) {
            Vertex v = bins[k].back();
            bins[k].pop_back();
            core_map[v] = static_cast<int>(k);
            if (deg_mode == "out" && rgraphtool::is_directed(g)) {
                for (Vertex u : rgraphtool::in_neighbors_range(v, g)) {
                    relax_neighbor(v, u);
                }
            } else if (deg_mode == "in" && rgraphtool::is_directed(g)) {
                for (Vertex u : rgraphtool::out_neighbors_range(v, g)) {
                    relax_neighbor(v, u);
                }
            } else {
                for (Vertex u : rgraphtool::all_neighbors_range(v, g)) {
                    relax_neighbor(v, u);
                }
            }
        }
    }
}

template <typename StrengthMap, typename WeightMap>
void strength_impl(const rgraphtool::Graph& g,
                   const std::string& mode,
                   StrengthMap& strength,
                   const WeightMap& weight) {
    using rgraphtool::Vertex;
    for (Vertex v : rgraphtool::vertices_range(g)) {
        double acc = 0.0;
        if (mode == "out" && rgraphtool::is_directed(g)) {
            for (const auto& e : rgraphtool::out_edges_range(v, g)) {
                acc += weight[e];
            }
        } else if (mode == "in" && rgraphtool::is_directed(g)) {
            for (const auto& e : rgraphtool::in_edges_range(v, g)) {
                acc += weight[e];
            }
        } else {
            for (const auto& e : rgraphtool::all_edges_range(v, g)) {
                acc += weight[e];
            }
        }
        strength[v] = acc;
    }
}

} // namespace

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_kcore_decomposition(SEXP graph_obj,
                                            std::string deg_mode = "total",
                                            bool checked = false) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph_obj);
    size_t n = rgraphtool::num_vertices(g);
    Rcpp::IntegerVector out(static_cast<R_xlen_t>(n));
    if (n == 0) return out;

    auto core_checked = rgraphtool::vertex_property_view(out, g);
    auto deg_checked = rgraphtool::make_vertex_property<size_t>(g, 0);
    auto pos_checked = rgraphtool::make_vertex_property<size_t>(g, 0);

    if (checked) {
        kcore_impl(g, deg_mode, core_checked, deg_checked, pos_checked);
    } else {
        auto core_u = core_checked.get_unchecked();
        auto deg_u = deg_checked.get_unchecked();
        auto pos_u = pos_checked.get_unchecked();
        kcore_impl(g, deg_mode, core_u, deg_u, pos_u);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_vertex_strength(SEXP graph_obj,
                                        std::string prop_name,
                                        std::string mode = "total",
                                        bool checked = false) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph_obj);
    size_t n = rgraphtool::num_vertices(g);
    Rcpp::NumericVector out(static_cast<R_xlen_t>(n), 0.0);
    if (n == 0) return out;

    auto strength_checked = rgraphtool::vertex_property_view(out, g);
    auto weight_checked = rgraphtool::edge_property_view(g, prop_name);

    if (checked) {
        strength_impl(g, mode, strength_checked, weight_checked);
    } else {
        auto strength_u = strength_checked.get_unchecked();
        auto weight_u = weight_checked.get_unchecked();
        strength_impl(g, mode, strength_u, weight_u);
    }
    return out;
}

// [[Rcpp::export]]
double cpp_extension_check_bounds(SEXP graph_obj,
                                  Rcpp::NumericVector vprop,
                                  int vertex_idx,
                                  bool use_unchecked = false) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph_obj);
    auto pm = rgraphtool::vertex_property_view(vprop, g);
    if (use_unchecked) {
        auto upm = pm.get_unchecked(); // validates vprop.size() >= num_vertices(g)
        rgraphtool::Vertex v = rgraphtool::checked_vertex(static_cast<size_t>(vertex_idx), g);
        return upm[v];
    }
    return pm[static_cast<size_t>(vertex_idx)];
}
