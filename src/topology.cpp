// rgraphtool: native C++23 centrality, shortest paths, connectivity,
// clustering, assortativity, classic random graph generators, and modularity.
//
// Copyright (C) 2026 Diogo Melo.
// Formulations and normalizations match graph-tool 3.9 (Tiago de Paula Peixoto).
// LGPL-3.0-or-later.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <stack>
#include <utility>
#include <vector>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include "graph_types.h"
#include "parallel.h"

namespace {

constexpr double TOP_INF = std::numeric_limits<double>::infinity();

struct TopRng {
    size_t index(size_t n) {
        return static_cast<size_t>(R_unif_index(static_cast<double>(n)));
    }
    double unif() { return unif_rand(); }
    template <class T>
    void shuffle(std::vector<T>& x) {
        for (size_t i = x.size(); i > 1; --i)
            std::swap(x[i - 1], x[index(i)]);
    }
};

struct AdjEdge {
    int to;
    int e_pos;    // 0-based position in boost::edges(g.adjacency)
    double w;
};

struct FlatAdj {
    int N = 0;
    int E = 0;
    bool directed = false;
    bool weighted = false;
    bool has_negative = false;
    std::vector<std::vector<AdjEdge>> out;
    std::vector<std::vector<AdjEdge>> in;
};

FlatAdj build_flat_adj(const Graph& g, Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
    FlatAdj fa;
    fa.N = static_cast<int>(boost::num_vertices(g.adjacency));
    fa.E = static_cast<int>(boost::num_edges(g.adjacency));
    fa.directed = g.directed;
    fa.out.resize(fa.N);
    if (fa.directed) fa.in.resize(fa.N);

    Rcpp::NumericVector wvec;
    if (weights.isNotNull()) {
        wvec = Rcpp::NumericVector(weights);
        if (wvec.size() != fa.E) {
            Rcpp::stop("weights length (%lld) must equal graph_num_edges(graph) (%d)",
                       static_cast<long long>(wvec.size()), fa.E);
        }
        fa.weighted = true;
    }

    auto erange = boost::edges(g.adjacency);
    int pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        double w = fa.weighted ? static_cast<double>(wvec[pos]) : 1.0;
        if (!std::isfinite(w)) Rcpp::stop("Edge weights must be finite numbers");
        if (w < 0.0) fa.has_negative = true;

        fa.out[u].push_back({v, pos, w});
        if (fa.directed) {
            fa.in[v].push_back({u, pos, w});
        } else if (u != v) {
            fa.out[v].push_back({u, pos, w});
        }
    }
    return fa;
}

// Single-source shortest distances (BFS if unweighted, Dijkstra if non-negative,
// Bellman-Ford if negative weights).
void single_source_distances(const FlatAdj& fa, int src, double max_dist,
                             std::vector<double>& dist) {
    dist.assign(fa.N, TOP_INF);
    if (src < 0 || src >= fa.N) return;
    dist[src] = 0.0;

    if (!fa.weighted) {
        std::vector<int> q;
        q.reserve(fa.N);
        q.push_back(src);
        size_t head = 0;
        while (head < q.size()) {
            int u = q[head++];
            double du = dist[u];
            if (du + 1.0 > max_dist) continue;
            for (const auto& e : fa.out[u]) {
                if (std::isinf(dist[e.to])) {
                    dist[e.to] = du + 1.0;
                    q.push_back(e.to);
                }
            }
        }
        return;
    }

    if (!fa.has_negative) {
        using P = std::pair<double, int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        pq.emplace(0.0, src);
        while (!pq.empty()) {
            auto [d, u] = pq.top();
            pq.pop();
            if (d > dist[u]) continue;
            if (d > max_dist) break;
            for (const auto& e : fa.out[u]) {
                double nd = d + e.w;
                if (nd <= max_dist && nd < dist[e.to] - 1e-15) {
                    dist[e.to] = nd;
                    pq.emplace(nd, e.to);
                }
            }
        }
        return;
    }

    // Bellman-Ford for graphs with negative edge weights
    for (int iter = 0; iter < fa.N - 1; ++iter) {
        bool changed = false;
        for (int u = 0; u < fa.N; ++u) {
            if (std::isinf(dist[u])) continue;
            for (const auto& e : fa.out[u]) {
                double nd = dist[u] + e.w;
                if (nd <= max_dist && nd < dist[e.to] - 1e-14) {
                    dist[e.to] = nd;
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    // Check negative cycle reachable from src
    for (int u = 0; u < fa.N; ++u) {
        if (std::isinf(dist[u])) continue;
        for (const auto& e : fa.out[u]) {
            if (dist[u] + e.w < dist[e.to] - 1e-12) {
                Rcpp::stop("Negative-weight cycle reachable from source vertex %d", src + 1);
            }
        }
    }
}

} // namespace

// ===========================================================================
// 1. Centrality Measures
// ===========================================================================

// [[Rcpp::export]]
Rcpp::List cpp_pagerank(SEXP pointer, double damping,
                        Rcpp::Nullable<Rcpp::NumericVector> pers,
                        Rcpp::Nullable<Rcpp::NumericVector> weights,
                        double tol, int max_iter) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    if (fa.has_negative) Rcpp::stop("PageRank requires non-negative edge weights");
    int N = fa.N;
    if (N == 0) {
        return Rcpp::List::create(Rcpp::Named("vector") = Rcpp::NumericVector(0),
                                  Rcpp::Named("iterations") = 0);
    }

    std::vector<double> p(N, 1.0 / N);
    if (pers.isNotNull()) {
        Rcpp::NumericVector pv(pers);
        if (pv.size() != N) Rcpp::stop("pers must have length graph_num_vertices(graph)");
        double psum = 0.0;
        for (int i = 0; i < N; ++i) {
            if (pv[i] < 0.0 || !std::isfinite(pv[i])) Rcpp::stop("pers entries must be finite and non-negative");
            psum += pv[i];
        }
        if (psum <= 0.0) Rcpp::stop("pers must have positive sum");
        for (int i = 0; i < N; ++i) p[i] = pv[i] / psum;
    }

    // Note: in undirected graphs with self-loops, graph-tool counts self-loop
    // weight twice in vertex strength for PageRank transition probabilities,
    // or once? Let's match boost/graph-tool: out_edges(u) for undirected visits
    // each loop once or twice? In boost::_all_edges_out, a self-loop appears TWICE!
    // Wait: in build_flat_adj above, for u == v we pushed once. Let's compute
    // exact out-strength accounting for undirected self-loops appearing twice in degree.
    std::vector<double> out_w(N, 0.0);
    for (int u = 0; u < N; ++u) {
        for (const auto& e : fa.out[u]) {
            out_w[u] += (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
        }
    }

    std::vector<double> r = p;
    std::vector<double> r_next(N, 0.0);
    int it = 0;
    for (it = 1; it <= max_iter; ++it) {
        double dangling = 0.0;
        for (int u = 0; u < N; ++u) {
            if (out_w[u] == 0.0) dangling += r[u];
        }
        double base_scale = (1.0 - damping) + damping * dangling;
        for (int v = 0; v < N; ++v) r_next[v] = base_scale * p[v];

        for (int u = 0; u < N; ++u) {
            if (out_w[u] <= 0.0) continue;
            double contrib = damping * r[u] / out_w[u];
            for (const auto& e : fa.out[u]) {
                double ew = (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
                r_next[e.to] += contrib * ew;
            }
        }

        double diff = 0.0;
        for (int v = 0; v < N; ++v) diff += std::abs(r_next[v] - r[v]);
        r.swap(r_next);
        if (diff < tol) break;
    }

    return Rcpp::List::create(
        Rcpp::Named("vector") = Rcpp::wrap(r),
        Rcpp::Named("iterations") = it
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_betweenness(SEXP pointer,
                           Rcpp::Nullable<Rcpp::NumericVector> weights,
                           bool normalized) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    if (fa.weighted) {
        for (int u = 0; u < fa.N; ++u) {
            for (const auto& e : fa.out[u]) {
                if (e.w <= 0.0) Rcpp::stop("Betweenness requires strictly positive edge weights");
            }
        }
    }
    int N = fa.N, E = fa.E;
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("vertex") = Rcpp::NumericVector(0),
            Rcpp::Named("edge") = Rcpp::NumericVector(0)
        );
    }

    int nw = rgt_parallel::resolve_threads(rgt_parallel::get_num_threads());
    std::vector<std::vector<double>> vb_tls(nw, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> eb_tls(nw, std::vector<double>(E, 0.0));

    rgt_parallel::parallel_for(static_cast<size_t>(N), nw, [&](int tid, size_t lo, size_t hi) {
        auto& vb = vb_tls[tid];
        auto& eb = eb_tls[tid];

        std::vector<int> stack_S;
        stack_S.reserve(N);
        std::vector<std::vector<std::pair<int, int>>> pred(N); // (pred_vertex, edge_pos)
        std::vector<double> sigma(N, 0.0);
        std::vector<double> dist(N, TOP_INF);
        std::vector<double> delta(N, 0.0);

        for (size_t s_sz = lo; s_sz < hi; ++s_sz) {
            int s = static_cast<int>(s_sz);
            stack_S.clear();
            for (int i = 0; i < N; ++i) pred[i].clear();
            std::fill(sigma.begin(), sigma.end(), 0.0);
            std::fill(dist.begin(), dist.end(), TOP_INF);
            sigma[s] = 1.0;
            dist[s] = 0.0;

            if (!fa.weighted) {
                std::vector<int> q;
                q.reserve(N);
                q.push_back(s);
                size_t head = 0;
                while (head < q.size()) {
                    int v = q[head++];
                    stack_S.push_back(v);
                    double dv = dist[v];
                    for (const auto& e : fa.out[v]) {
                        int w = e.to;
                        if (std::isinf(dist[w])) {
                            dist[w] = dv + 1.0;
                            q.push_back(w);
                        }
                        if (dist[w] == dv + 1.0) {
                            sigma[w] += sigma[v];
                            pred[w].emplace_back(v, e.e_pos);
                        }
                    }
                }
            } else {
                using P = std::pair<double, int>;
                std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
                pq.emplace(0.0, s);
                while (!pq.empty()) {
                    auto [dv, v] = pq.top();
                    pq.pop();
                    if (dv > dist[v] + 1e-15) continue;
                    stack_S.push_back(v);
                    for (const auto& e : fa.out[v]) {
                        int w = e.to;
                        double nd = dv + e.w;
                        double eps = 1e-13 * std::max(1.0, std::abs(nd));
                        if (nd < dist[w] - eps) {
                            dist[w] = nd;
                            pq.emplace(nd, w);
                            sigma[w] = sigma[v];
                            pred[w].clear();
                            pred[w].emplace_back(v, e.e_pos);
                        } else if (std::abs(nd - dist[w]) <= eps) {
                            sigma[w] += sigma[v];
                            pred[w].emplace_back(v, e.e_pos);
                        }
                    }
                }
            }

            std::fill(delta.begin(), delta.end(), 0.0);
            while (!stack_S.empty()) {
                int w = stack_S.back();
                stack_S.pop_back();
                for (const auto& pr : pred[w]) {
                    int v = pr.first;
                    int ep = pr.second;
                    double c = (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                    eb[ep] += c;
                    delta[v] += c;
                }
                if (w != s) vb[w] += delta[w];
            }
        }
    });

    Rcpp::NumericVector v_out(N, 0.0);
    Rcpp::NumericVector e_out(E, 0.0);
    for (int t = 0; t < nw; ++t) {
        for (int i = 0; i < N; ++i) v_out[i] += vb_tls[t][i];
        for (int e = 0; e < E; ++e) e_out[e] += eb_tls[t][e];
    }

    if (!fa.directed) {
        for (int i = 0; i < N; ++i) v_out[i] *= 0.5;
        for (int e = 0; e < E; ++e) e_out[e] *= 0.5;
    }

    if (normalized) {
        double v_norm = fa.directed
            ? (N > 2 ? static_cast<double>(N - 1) * (N - 2) : 1.0)
            : (N > 2 ? static_cast<double>(N - 1) * (N - 2) * 0.5 : 1.0);
        double e_norm = fa.directed
            ? (N > 1 ? static_cast<double>(N) * (N - 1) : 1.0)
            : (N > 1 ? static_cast<double>(N) * (N - 1) * 0.5 : 1.0);
        if (N > 2) for (int i = 0; i < N; ++i) v_out[i] /= v_norm;
        if (N > 1) for (int e = 0; e < E; ++e) e_out[e] /= e_norm;
    }

    return Rcpp::List::create(
        Rcpp::Named("vertex") = v_out,
        Rcpp::Named("edge") = e_out
    );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_closeness(SEXP pointer,
                                  Rcpp::Nullable<Rcpp::NumericVector> weights,
                                  bool normalized, bool harmonic) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    if (fa.has_negative) Rcpp::stop("Closeness requires non-negative edge weights");
    int N = fa.N;
    Rcpp::NumericVector out(N, 0.0);
    if (N <= 1) return out;

    std::vector<double> res(N, 0.0);
    int nw = rgt_parallel::resolve_threads(rgt_parallel::get_num_threads());
    rgt_parallel::parallel_for(static_cast<size_t>(N), nw, [&](int /*tid*/, size_t lo, size_t hi) {
        std::vector<double> dist;
        for (size_t s_sz = lo; s_sz < hi; ++s_sz) {
            int s = static_cast<int>(s_sz);
            single_source_distances(fa, s, TOP_INF, dist);
            if (harmonic) {
                double hsum = 0.0;
                for (int v = 0; v < N; ++v) {
                    if (v != s && std::isfinite(dist[v]) && dist[v] > 0.0) {
                        hsum += 1.0 / dist[v];
                    }
                }
                res[s] = normalized ? (hsum / static_cast<double>(N - 1)) : hsum;
            } else {
                double dsum = 0.0;
                int reached = 0;
                for (int v = 0; v < N; ++v) {
                    if (v != s && std::isfinite(dist[v])) {
                        dsum += dist[v];
                        reached++;
                    }
                }
                if (reached == 0 || dsum <= 0.0) {
                    res[s] = 0.0;
                } else {
                    res[s] = normalized ? (static_cast<double>(reached) / dsum) : (1.0 / dsum);
                }
            }
        }
    });

    for (int i = 0; i < N; ++i) out[i] = res[i];
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_eigenvector(SEXP pointer,
                           Rcpp::Nullable<Rcpp::NumericVector> weights,
                           double tol, int max_iter) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    if (N == 0) {
        return Rcpp::List::create(Rcpp::Named("eigenvalue") = NA_REAL,
                                  Rcpp::Named("vector") = Rcpp::NumericVector(0));
    }

    // Match graph-tool: x <- A * x (for directed, incoming adjacency or outgoing?
    // In graph-tool, eigenvector centrality uses incoming edges for directed graphs:
    // c_v = (1/lambda) sum_{u -> v} w_{uv} c_u).
    std::vector<double> x(N, 1.0 / std::sqrt(static_cast<double>(N)));
    std::vector<double> y(N, 0.0);
    double eig = 0.0;

    for (int it = 0; it < max_iter; ++it) {
        std::fill(y.begin(), y.end(), 0.0);
        for (int u = 0; u < N; ++u) {
            for (const auto& e : fa.out[u]) {
                double ew = (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
                if (fa.directed) {
                    y[e.to] += ew * x[u];
                } else {
                    y[u] += ew * x[e.to];
                }
            }
        }
        double norm2 = 0.0;
        for (int i = 0; i < N; ++i) norm2 += y[i] * y[i];
        eig = std::sqrt(norm2);
        if (eig <= 0.0) {
            std::fill(x.begin(), x.end(), 0.0);
            break;
        }
        double diff = 0.0;
        for (int i = 0; i < N; ++i) {
            double yn = y[i] / eig;
            diff += std::abs(yn - x[i]);
            x[i] = yn;
        }
        if (diff < tol) break;
    }

    return Rcpp::List::create(
        Rcpp::Named("eigenvalue") = eig,
        Rcpp::Named("vector") = Rcpp::wrap(x)
    );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_katz(SEXP pointer, double alpha,
                             Rcpp::Nullable<Rcpp::NumericVector> beta,
                             Rcpp::Nullable<Rcpp::NumericVector> weights,
                             bool normalized, double tol, int max_iter) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    if (N == 0) return Rcpp::NumericVector(0);

    std::vector<double> bvec(N, 1.0);
    if (beta.isNotNull()) {
        Rcpp::NumericVector bv(beta);
        if (bv.size() != N) Rcpp::stop("beta must have length graph_num_vertices(graph)");
        for (int i = 0; i < N; ++i) bvec[i] = bv[i];
    }

    std::vector<double> x = bvec;
    std::vector<double> y(N, 0.0);
    for (int it = 0; it < max_iter; ++it) {
        y = bvec;
        for (int u = 0; u < N; ++u) {
            for (const auto& e : fa.out[u]) {
                double ew = (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
                if (fa.directed) y[e.to] += alpha * ew * x[u];
                else y[u] += alpha * ew * x[e.to];
            }
        }
        double diff = 0.0;
        for (int i = 0; i < N; ++i) {
            diff += std::abs(y[i] - x[i]);
            if (!std::isfinite(y[i])) {
                Rcpp::stop("Katz centrality diverged; choose a smaller alpha (< 1 / lambda_max)");
            }
        }
        x.swap(y);
        if (diff < tol) break;
    }

    if (normalized) {
        double norm2 = 0.0;
        for (double v : x) norm2 += v * v;
        double nrm = std::sqrt(norm2);
        if (nrm > 0.0) for (double& v : x) v /= nrm;
    }
    return Rcpp::wrap(x);
}

// [[Rcpp::export]]
Rcpp::List cpp_hits(SEXP pointer,
                    Rcpp::Nullable<Rcpp::NumericVector> weights,
                    double tol, int max_iter) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("eigenvalue") = NA_REAL,
            Rcpp::Named("authority") = Rcpp::NumericVector(0),
            Rcpp::Named("hub") = Rcpp::NumericVector(0)
        );
    }

    std::vector<double> auth(N, 1.0 / std::sqrt(static_cast<double>(N)));
    std::vector<double> hub = auth;
    std::vector<double> auth_next(N, 0.0), hub_next(N, 0.0);
    double eig = 0.0;

    for (int it = 0; it < max_iter; ++it) {
        std::fill(auth_next.begin(), auth_next.end(), 0.0);
        for (int u = 0; u < N; ++u) {
            for (const auto& e : fa.out[u]) {
                double ew = (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
                auth_next[e.to] += ew * hub[u];
            }
        }
        double na = 0.0;
        for (double v : auth_next) na += v * v;
        na = std::sqrt(na);
        if (na > 0.0) for (double& v : auth_next) v /= na;

        std::fill(hub_next.begin(), hub_next.end(), 0.0);
        for (int u = 0; u < N; ++u) {
            for (const auto& e : fa.out[u]) {
                double ew = (!fa.directed && e.to == u) ? (2.0 * e.w) : e.w;
                hub_next[u] += ew * auth_next[e.to];
            }
        }
        double nh = 0.0;
        for (double v : hub_next) nh += v * v;
        nh = std::sqrt(nh);
        eig = nh;
        if (nh > 0.0) for (double& v : hub_next) v /= nh;

        double diff = 0.0;
        for (int i = 0; i < N; ++i) {
            diff += std::abs(auth_next[i] - auth[i]) + std::abs(hub_next[i] - hub[i]);
        }
        auth.swap(auth_next);
        hub.swap(hub_next);
        if (diff < tol) break;
    }

    return Rcpp::List::create(
        Rcpp::Named("eigenvalue") = eig,
        Rcpp::Named("authority") = Rcpp::wrap(auth),
        Rcpp::Named("hub") = Rcpp::wrap(hub)
    );
}

// ===========================================================================
// 2. Shortest Paths & Diameter
// ===========================================================================

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_shortest_distance(SEXP pointer,
                                          Rcpp::IntegerVector sources,
                                          Rcpp::Nullable<Rcpp::NumericVector> weights,
                                          double max_dist) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    int ns = sources.size();
    Rcpp::NumericMatrix out(ns, N);
    if (ns == 0 || N == 0) return out;

    for (int i = 0; i < ns; ++i) {
        if (sources[i] < 1 || sources[i] > N) {
            Rcpp::stop("source vertex IDs must be between 1 and %d", N);
        }
    }

    // Bellman-Ford throws Rcpp::stop on negative cycles, so run sequentially if negative weights
    int nw = fa.has_negative ? 1 : rgt_parallel::resolve_threads(rgt_parallel::get_num_threads());
    std::vector<std::vector<double>> rows(ns);
    rgt_parallel::parallel_for(static_cast<size_t>(ns), nw, [&](int /*tid*/, size_t lo, size_t hi) {
        for (size_t idx = lo; idx < hi; ++idx) {
            single_source_distances(fa, sources[static_cast<int>(idx)] - 1, max_dist, rows[idx]);
        }
    });

    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < N; ++j) out(i, j) = rows[i][j];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_shortest_path(SEXP pointer, int source, int target,
                             Rcpp::Nullable<Rcpp::NumericVector> weights) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    if (source < 1 || source > N || target < 1 || target > N) {
        Rcpp::stop("source and target must be valid 1-based vertex IDs in 1..%d", N);
    }
    int src = source - 1, tgt = target - 1;
    if (src == tgt) {
        return Rcpp::List::create(
            Rcpp::Named("vertices") = Rcpp::IntegerVector::create(source),
            Rcpp::Named("edges") = Rcpp::IntegerVector(0),
            Rcpp::Named("distance") = 0.0
        );
    }

    std::vector<double> dist(N, TOP_INF);
    std::vector<int> parent_v(N, -1);
    std::vector<int> parent_e(N, -1);
    dist[src] = 0.0;

    if (!fa.weighted) {
        std::vector<int> q;
        q.push_back(src);
        size_t head = 0;
        while (head < q.size()) {
            int u = q[head++];
            if (u == tgt) break;
            for (const auto& e : fa.out[u]) {
                if (std::isinf(dist[e.to])) {
                    dist[e.to] = dist[u] + 1.0;
                    parent_v[e.to] = u;
                    parent_e[e.to] = e.e_pos;
                    q.push_back(e.to);
                }
            }
        }
    } else if (!fa.has_negative) {
        using P = std::pair<double, int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        pq.emplace(0.0, src);
        while (!pq.empty()) {
            auto [d, u] = pq.top();
            pq.pop();
            if (d > dist[u]) continue;
            if (u == tgt) break;
            for (const auto& e : fa.out[u]) {
                double nd = d + e.w;
                if (nd < dist[e.to] - 1e-15) {
                    dist[e.to] = nd;
                    parent_v[e.to] = u;
                    parent_e[e.to] = e.e_pos;
                    pq.emplace(nd, e.to);
                }
            }
        }
    } else {
        for (int iter = 0; iter < N - 1; ++iter) {
            bool changed = false;
            for (int u = 0; u < N; ++u) {
                if (std::isinf(dist[u])) continue;
                for (const auto& e : fa.out[u]) {
                    double nd = dist[u] + e.w;
                    if (nd < dist[e.to] - 1e-14) {
                        dist[e.to] = nd;
                        parent_v[e.to] = u;
                        parent_e[e.to] = e.e_pos;
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
    }

    if (std::isinf(dist[tgt])) {
        return Rcpp::List::create(
            Rcpp::Named("vertices") = Rcpp::IntegerVector(0),
            Rcpp::Named("edges") = Rcpp::IntegerVector(0),
            Rcpp::Named("distance") = TOP_INF
        );
    }

    std::vector<int> vpath, epath;
    for (int cur = tgt; cur != src; cur = parent_v[cur]) {
        vpath.push_back(cur + 1);
        epath.push_back(parent_e[cur] + 1);
    }
    vpath.push_back(src + 1);
    std::reverse(vpath.begin(), vpath.end());
    std::reverse(epath.begin(), epath.end());

    return Rcpp::List::create(
        Rcpp::Named("vertices") = Rcpp::wrap(vpath),
        Rcpp::Named("edges") = Rcpp::wrap(epath),
        Rcpp::Named("distance") = dist[tgt]
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_pseudo_diameter(SEXP pointer, int source,
                               Rcpp::Nullable<Rcpp::NumericVector> weights) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g, weights);
    int N = fa.N;
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("diameter") = 0.0,
            Rcpp::Named("endpoints") = Rcpp::IntegerVector(0)
        );
    }
    if (source < 1 || source > N) Rcpp::stop("source must be in 1..%d", N);

    int cur = source - 1;
    double best_d = 0.0;
    int best_u = cur, best_v = cur;
    std::vector<double> dist;

    for (int iter = 0; iter < 64; ++iter) {
        single_source_distances(fa, cur, TOP_INF, dist);
        double max_d = -1.0;
        int far = cur;
        for (int v = 0; v < N; ++v) {
            if (std::isfinite(dist[v]) && dist[v] > max_d) {
                max_d = dist[v];
                far = v;
            }
        }
        if (max_d > best_d) {
            best_d = max_d;
            best_u = cur;
            best_v = far;
            cur = far;
        } else {
            break;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("diameter") = best_d,
        Rcpp::Named("endpoints") = Rcpp::IntegerVector::create(best_u + 1, best_v + 1)
    );
}

// ===========================================================================
// 3. Connectivity, Components, DAG, Clustering & Assortativity
// ===========================================================================

// [[Rcpp::export]]
Rcpp::List cpp_label_components(SEXP pointer) {
    const Graph& g = get_graph(pointer);
    FlatAdj fa = build_flat_adj(g);
    int N = fa.N;
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("membership") = Rcpp::IntegerVector(0),
            Rcpp::Named("sizes") = Rcpp::IntegerVector(0),
            Rcpp::Named("num_components") = 0,
            Rcpp::Named("attractors") = Rcpp::LogicalVector(0)
        );
    }

    std::vector<int> comp(N, -1);
    int n_comp = 0;

    if (!fa.directed) {
        std::vector<int> q;
        q.reserve(N);
        for (int s = 0; s < N; ++s) {
            if (comp[s] != -1) continue;
            comp[s] = n_comp;
            q.clear();
            q.push_back(s);
            size_t head = 0;
            while (head < q.size()) {
                int u = q[head++];
                for (const auto& e : fa.out[u]) {
                    if (comp[e.to] == -1) {
                        comp[e.to] = n_comp;
                        q.push_back(e.to);
                    }
                }
            }
            n_comp++;
        }
    } else {
        // Iterative Tarjan strongly connected components (avoids stack overflow on deep graphs)
        std::vector<int> disc(N, -1), low(N, 0), st;
        std::vector<uint8_t> on_stack(N, 0);
        int timer = 0;

        struct Frame { int u; size_t ei; };
        std::vector<Frame> call_stack;

        for (int start = 0; start < N; ++start) {
            if (disc[start] != -1) continue;
            disc[start] = low[start] = timer++;
            st.push_back(start);
            on_stack[start] = 1;
            call_stack.push_back({start, 0});

            while (!call_stack.empty()) {
                auto& fr = call_stack.back();
                int u = fr.u;
                if (fr.ei < fa.out[u].size()) {
                    int v = fa.out[u][fr.ei].to;
                    fr.ei++;
                    if (disc[v] == -1) {
                        disc[v] = low[v] = timer++;
                        st.push_back(v);
                        on_stack[v] = 1;
                        call_stack.push_back({v, 0});
                    } else if (on_stack[v]) {
                        low[u] = std::min(low[u], disc[v]);
                    }
                } else {
                    if (low[u] == disc[u]) {
                        while (true) {
                            int w = st.back();
                            st.pop_back();
                            on_stack[w] = 0;
                            comp[w] = n_comp;
                            if (w == u) break;
                        }
                        n_comp++;
                    }
                    call_stack.pop_back();
                    if (!call_stack.empty()) {
                        int parent = call_stack.back().u;
                        low[parent] = std::min(low[parent], low[u]);
                    }
                }
            }
        }
    }

    Rcpp::IntegerVector mem(N);
    Rcpp::IntegerVector sizes(n_comp, 0);
    for (int i = 0; i < N; ++i) {
        mem[i] = comp[i] + 1;
        sizes[comp[i]]++;
    }

    Rcpp::LogicalVector attractors(n_comp, TRUE);
    if (fa.directed) {
        for (int u = 0; u < N; ++u) {
            for (const auto& e : fa.out[u]) {
                if (comp[u] != comp[e.to]) {
                    attractors[comp[u]] = FALSE;
                }
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("membership") = mem,
        Rcpp::Named("sizes") = sizes,
        Rcpp::Named("num_components") = n_comp,
        Rcpp::Named("attractors") = attractors
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_biconnected_components(SEXP pointer) {
    const Graph& g = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(g.adjacency));
    int E = static_cast<int>(boost::num_edges(g.adjacency));

    // Treat graph as undirected for biconnectivity (standard Hopcroft-Tarjan)
    std::vector<std::vector<std::pair<int, int>>> adj(N);
    std::vector<std::pair<int, int>> edge_ends(E);
    auto erange = boost::edges(g.adjacency);
    int pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        edge_ends[pos] = {u, v};
        adj[u].emplace_back(v, pos);
        if (u != v) adj[v].emplace_back(u, pos);
    }

    std::vector<int> disc(N, -1), low(N, 0), parent(N, -1);
    std::vector<uint8_t> is_art(N, 0), edge_vis(E, 0);
    std::vector<int> edge_comp(E, 0);
    std::vector<int> edge_stack;
    int timer = 0, n_bicomp = 0;

    // Iterative Hopcroft-Tarjan DFS
    struct BFrame { int u; size_t ei; int children; };
    std::vector<BFrame> st;

    for (int root = 0; root < N; ++root) {
        if (disc[root] != -1) continue;
        disc[root] = low[root] = timer++;
        st.push_back({root, 0, 0});

        while (!st.empty()) {
            auto& fr = st.back();
            int u = fr.u;
            if (fr.ei < adj[u].size()) {
                auto [v, eid] = adj[u][fr.ei++];
                if (u == v) {
                    // Self-loop forms its own 1-edge component
                    if (!edge_vis[eid]) {
                        edge_vis[eid] = 1;
                        edge_comp[eid] = ++n_bicomp;
                    }
                    continue;
                }
                if (disc[v] == -1) {
                    parent[v] = u;
                    fr.children++;
                    edge_vis[eid] = 1;
                    edge_stack.push_back(eid);
                    disc[v] = low[v] = timer++;
                    st.push_back({v, 0, 0});
                } else if (!edge_vis[eid]) {
                    edge_vis[eid] = 1;
                    edge_stack.push_back(eid);
                    low[u] = std::min(low[u], disc[v]);
                }
            } else {
                int p = parent[u];
                int ch = fr.children;
                st.pop_back();
                if (p == -1) {
                    if (ch > 1) is_art[u] = 1;
                } else {
                    low[p] = std::min(low[p], low[u]);
                    if (parent[p] != -1 && low[u] >= disc[p]) {
                        is_art[p] = 1;
                    }
                    if (low[u] >= disc[p]) {
                        ++n_bicomp;
                        while (!edge_stack.empty()) {
                            int top_e = edge_stack.back();
                            edge_stack.pop_back();
                            edge_comp[top_e] = n_bicomp;
                            auto [eu, ev] = edge_ends[top_e];
                            if ((eu == p && ev == u) || (eu == u && ev == p)) break;
                        }
                    }
                }
            }
        }
        if (!edge_stack.empty()) {
            ++n_bicomp;
            while (!edge_stack.empty()) {
                edge_comp[edge_stack.back()] = n_bicomp;
                edge_stack.pop_back();
            }
        }
    }

    // Count edges per biconnected component to identify bridges (size-1 non-loop components)
    std::vector<int> comp_size(n_bicomp + 1, 0);
    for (int e = 0; e < E; ++e) {
        if (edge_comp[e] > 0) comp_size[edge_comp[e]]++;
    }
    std::vector<int> bridges;
    for (int e = 0; e < E; ++e) {
        int c = edge_comp[e];
        if (c > 0 && comp_size[c] == 1 && edge_ends[e].first != edge_ends[e].second) {
            bridges.push_back(e + 1);
        }
    }

    std::vector<int> art_pts;
    for (int v = 0; v < N; ++v) if (is_art[v]) art_pts.push_back(v + 1);

    return Rcpp::List::create(
        Rcpp::Named("edge_component") = Rcpp::wrap(edge_comp),
        Rcpp::Named("num_components") = n_bicomp,
        Rcpp::Named("articulation_points") = Rcpp::wrap(art_pts),
        Rcpp::Named("bridges") = Rcpp::wrap(bridges)
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_topological_sort(SEXP pointer) {
    const Graph& g = get_graph(pointer);
    if (!g.directed) Rcpp::stop("topological_sort requires a directed graph");
    FlatAdj fa = build_flat_adj(g);
    int N = fa.N;
    std::vector<int> indeg(N, 0);
    bool has_self_loop = false;
    for (int u = 0; u < N; ++u) {
        for (const auto& e : fa.out[u]) {
            if (u == e.to) has_self_loop = true;
            indeg[e.to]++;
        }
    }
    if (has_self_loop) {
        return Rcpp::List::create(Rcpp::Named("is_dag") = false,
                                  Rcpp::Named("order") = Rcpp::IntegerVector(0));
    }

    std::queue<int> q;
    for (int i = 0; i < N; ++i) if (indeg[i] == 0) q.push(i);

    std::vector<int> order;
    order.reserve(N);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        order.push_back(u + 1);
        for (const auto& e : fa.out[u]) {
            if (--indeg[e.to] == 0) q.push(e.to);
        }
    }

    bool is_dag = (static_cast<int>(order.size()) == N);
    Rcpp::IntegerVector ord_out = is_dag
        ? Rcpp::IntegerVector(order.begin(), order.end())
        : Rcpp::IntegerVector(0);
    return Rcpp::List::create(
        Rcpp::Named("is_dag") = is_dag,
        Rcpp::Named("order") = ord_out
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_clustering(SEXP pointer) {
    const Graph& g = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(g.adjacency));
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("local") = Rcpp::NumericVector(0),
            Rcpp::Named("global") = 0.0,
            Rcpp::Named("std_err") = 0.0
        );
    }

    // Match graph-tool: deduplicate parallel edges and ignore self-loops for simple
    // triangle/triple counting; for directed graphs, local clustering counts directed
    // edges between distinct neighbors of v divided by k_i * (k_i - 1), while global
    // clustering treats the underlying undirected simple graph.
    std::vector<boost::unordered_flat_set<int>> out_set(N), undir_set(N);
    auto erange = boost::edges(g.adjacency);
    for (auto it = erange.first; it != erange.second; ++it) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        if (u == v) continue;
        out_set[u].insert(v);
        if (!g.directed) out_set[v].insert(u);
        undir_set[u].insert(v);
        undir_set[v].insert(u);
    }

    Rcpp::NumericVector loc(N, 0.0);
    std::vector<double> tri_v(N, 0.0), triple_v(N, 0.0);
    double total_tri = 0.0, total_triples = 0.0;

    for (int v = 0; v < N; ++v) {
        // Local clustering (directed vs undirected matching graph-tool)
        if (!g.directed) {
            const auto& nbrs = undir_set[v];
            int k = static_cast<int>(nbrs.size());
            double tr = static_cast<double>(k) * (k - 1) * 0.5;
            triple_v[v] = tr;
            total_triples += tr;
            if (k >= 2) {
                std::vector<int> nv(nbrs.begin(), nbrs.end());
                int closed = 0;
                for (int i = 0; i < k; ++i) {
                    for (int j = i + 1; j < k; ++j) {
                        if (undir_set[nv[i]].contains(nv[j])) closed++;
                    }
                }
                tri_v[v] = static_cast<double>(closed);
                total_tri += static_cast<double>(closed);
                loc[v] = static_cast<double>(closed) / tr;
            }
        } else {
            // Directed local clustering: distinct neighbors (in U out), directed links between them
            const auto& nbrs = undir_set[v];
            int k = static_cast<int>(nbrs.size());
            if (k >= 2) {
                std::vector<int> nv(nbrs.begin(), nbrs.end());
                int links = 0;
                int undir_closed = 0;
                for (int i = 0; i < k; ++i) {
                    for (int j = i + 1; j < k; ++j) {
                        if (out_set[nv[i]].contains(nv[j])) links++;
                        if (out_set[nv[j]].contains(nv[i])) links++;
                        if (undir_set[nv[i]].contains(nv[j])) undir_closed++;
                    }
                }
                loc[v] = static_cast<double>(links) / (static_cast<double>(k) * (k - 1));
                double tr = static_cast<double>(k) * (k - 1) * 0.5;
                triple_v[v] = tr;
                tri_v[v] = static_cast<double>(undir_closed);
                total_triples += tr;
                total_tri += static_cast<double>(undir_closed);
            }
        }
    }

    double c_glob = (total_triples > 0.0) ? (total_tri / total_triples) : 0.0;
    // Jackknife standard error matching graph-tool's global_clustering()[1]
    double var_jk = 0.0;
    if (total_triples > 0.0 && N > 1) {
        // Removing vertex v removes 3 * tri_v[v] from total_tri (since each triangle has 3 vertices!)
        // and removes triple_v[v] + (triangles containing v as outer pair = 2 * tri_v[v])... wait:
        // Removing v removes triple_v[v] centered at v PLUS for each neighbor u of v, removes (deg(u)-1) triples!
        for (int v = 0; v < N; ++v) {
            double denom = total_triples - triple_v[v];
            double c_minus = (denom > 0.0) ? ((total_tri - tri_v[v]) / denom) : 0.0;
            double d = c_glob - c_minus;
            var_jk += d * d;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("local") = loc,
        Rcpp::Named("global") = c_glob,
        Rcpp::Named("std_err") = std::sqrt(var_jk)
    );
}

// [[Rcpp::export]]
double cpp_assortativity(SEXP pointer,
                         Rcpp::NumericVector x_src,
                         Rcpp::NumericVector x_tgt,
                         bool categorical,
                         Rcpp::Nullable<Rcpp::NumericVector> weights) {
    const Graph& g = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(g.adjacency));
    int E = static_cast<int>(boost::num_edges(g.adjacency));
    if (N == 0 || E == 0) return NA_REAL;

    Rcpp::NumericVector wvec;
    bool weighted = weights.isNotNull();
    if (weighted) {
        wvec = Rcpp::NumericVector(weights);
        if (wvec.size() != E) Rcpp::stop("weights must have length graph_num_edges(graph)");
    }

    if (categorical) {
        // Newman (2003) discrete assortativity r = (sum_i e_ii - sum_i a_i b_i) / (1 - sum_i a_i b_i)
        boost::unordered_flat_map<int64_t, double> a, b_sum;
        double e_diag = 0.0, W = 0.0;
        auto erange = boost::edges(g.adjacency);
        int pos = 0;
        for (auto it = erange.first; it != erange.second; ++it, ++pos) {
            int u = static_cast<int>(it->s);
            int v = static_cast<int>(it->t);
            double w = weighted ? static_cast<double>(wvec[pos]) : 1.0;
            int64_t cu = static_cast<int64_t>(x_src[u]);
            int64_t cv = static_cast<int64_t>(x_tgt[v]);
            if (g.directed) {
                a[cu] += w;
                b_sum[cv] += w;
                if (cu == cv) e_diag += w;
                W += w;
            } else {
                a[cu] += w;
                a[cv] += w;
                b_sum[cu] += w;
                b_sum[cv] += w;
                if (cu == cv) e_diag += 2.0 * w;
                W += 2.0 * w;
            }
        }
        if (W <= 0.0) return NA_REAL;
        double sum_ab = 0.0;
        for (const auto& kv : a) {
            auto it = b_sum.find(kv.first);
            if (it != b_sum.end()) sum_ab += (kv.second / W) * (it->second / W);
        }
        double denom = 1.0 - sum_ab;
        return (std::abs(denom) < 1e-15) ? NA_REAL : ((e_diag / W - sum_ab) / denom);
    }

    // Scalar Pearson correlation across oriented edges
    double W = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    auto erange = boost::edges(g.adjacency);
    int pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        double w = weighted ? static_cast<double>(wvec[pos]) : 1.0;
        auto add_pair = [&](double xv, double yv, double wt) {
            W += wt;
            sx += wt * xv;
            sy += wt * yv;
            sxx += wt * xv * xv;
            syy += wt * yv * yv;
            sxy += wt * xv * yv;
        };
        if (g.directed) {
            add_pair(x_src[u], x_tgt[v], w);
        } else {
            add_pair(x_src[u], x_tgt[v], w);
            add_pair(x_src[v], x_tgt[u], w);
        }
    }
    if (W <= 0.0) return NA_REAL;
    double mx = sx / W, my = sy / W;
    double vx = sxx / W - mx * mx;
    double vy = syy / W - my * my;
    if (vx <= 1e-15 || vy <= 1e-15) return NA_REAL;
    return (sxy / W - mx * my) / std::sqrt(vx * vy);
}

// ===========================================================================
// 4. Modularity & Multilevel Louvain Maximization
// ===========================================================================

// [[Rcpp::export]]
double cpp_modularity(SEXP pointer, Rcpp::IntegerVector blocks, double gamma,
                      Rcpp::Nullable<Rcpp::NumericVector> weights) {
    const Graph& g = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(g.adjacency));
    int E = static_cast<int>(boost::num_edges(g.adjacency));
    if (blocks.size() != N) Rcpp::stop("blocks must have length graph_num_vertices(graph)");
    if (N == 0 || E == 0) return 0.0;

    Rcpp::NumericVector wvec;
    bool weighted = weights.isNotNull();
    if (weighted) {
        wvec = Rcpp::NumericVector(weights);
        if (wvec.size() != E) Rcpp::stop("weights must have length graph_num_edges(graph)");
    }

    boost::unordered_flat_map<int, double> err, k_out, k_in;
    double W = 0.0;
    auto erange = boost::edges(g.adjacency);
    int pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        double w = weighted ? static_cast<double>(wvec[pos]) : 1.0;
        int ru = blocks[u], rv = blocks[v];
        if (g.directed) {
            W += w;
            k_out[ru] += w;
            k_in[rv] += w;
            if (ru == rv) err[ru] += w;
        } else {
            W += 2.0 * w; // W = 2m
            k_out[ru] += w;
            k_out[rv] += w;
            if (ru == rv) err[ru] += 2.0 * w;
        }
    }
    if (W == 0.0) return 0.0;

    double Q = 0.0;
    for (const auto& kv : k_out) {
        int r = kv.first;
        double e_in = 0.0;
        auto it_e = err.find(r);
        if (it_e != err.end()) e_in = it_e->second;
        if (g.directed) {
            double ki = 0.0;
            auto it_i = k_in.find(r);
            if (it_i != k_in.end()) ki = it_i->second;
            Q += (e_in - gamma * (kv.second * ki) / W) / W;
        } else {
            Q += (e_in - gamma * (kv.second * kv.second) / W) / W;
        }
    }
    if (g.directed) {
        // Include blocks that only have incoming edges (k_out == 0)
        for (const auto& kv : k_in) {
            if (k_out.find(kv.first) == k_out.end()) {
                auto it_e = err.find(kv.first);
                if (it_e != err.end()) Q += it_e->second / W;
            }
        }
    }
    return Q;
}

// [[Rcpp::export]]
Rcpp::List cpp_maximize_modularity(SEXP pointer, double gamma,
                                   Rcpp::Nullable<Rcpp::NumericVector> weights,
                                   int nruns, int max_levels) {
    const Graph& g = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(g.adjacency));
    int E = static_cast<int>(boost::num_edges(g.adjacency));
    if (N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("blocks") = Rcpp::IntegerVector(0),
            Rcpp::Named("modularity") = 0.0,
            Rcpp::Named("B") = 0
        );
    }
    if (E == 0) {
        Rcpp::IntegerVector b_triv(N, 1);
        return Rcpp::List::create(
            Rcpp::Named("blocks") = b_triv,
            Rcpp::Named("modularity") = 0.0,
            Rcpp::Named("B") = 1
        );
    }

    Rcpp::NumericVector wvec;
    bool weighted = weights.isNotNull();
    if (weighted) {
        wvec = Rcpp::NumericVector(weights);
        if (wvec.size() != E) Rcpp::stop("weights must have length graph_num_edges(graph)");
    }

    // Build base network representation for fast multilevel Louvain
    // Represent every graph in symmetric interaction form for modularity delta:
    // Notice for directed modularity: Q = (1/W) sum_r (e_rr - gamma * k_r^out * k_r^in / W).
    // Moving vertex u changes Q via (e_{u->s} + e_{s->u}) - gamma * (k_u^out * K_s^in + k_u^in * K_s^out) / W!
    // Undirected is the special case where W = 2m, k_u^out = k_u^in = k_u, and e_{u,s} = sum w_{u,v in s} * 2!
    struct LouvainGraph {
        int n = 0;
        bool directed = false;
        double W = 0.0;
        std::vector<double> self_w; // internal self-weight (2*w for undirected loop, w for directed loop)
        std::vector<double> k_out, k_in;
        // Undirected / combined neighbor interaction weight (already scaled: 2*w for undirected u!=v, w_uv + w_vu for directed u!=v)
        std::vector<std::vector<std::pair<int, double>>> nbrs;
        // Original directed edges for coarsening
        std::vector<std::tuple<int, int, double>> raw_edges;
    };

    LouvainGraph base;
    base.n = N;
    base.directed = g.directed;
    base.self_w.assign(N, 0.0);
    base.k_out.assign(N, 0.0);
    base.k_in.assign(N, 0.0);
    std::vector<boost::unordered_flat_map<int, double>> sym_adj(N);

    auto erange = boost::edges(g.adjacency);
    int pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        double w = weighted ? static_cast<double>(wvec[pos]) : 1.0;
        base.raw_edges.emplace_back(u, v, w);
        if (g.directed) {
            base.W += w;
            base.k_out[u] += w;
            base.k_in[v] += w;
            if (u == v) base.self_w[u] += w;
            else {
                sym_adj[u][v] += w;
                sym_adj[v][u] += w;
            }
        } else {
            base.W += 2.0 * w;
            base.k_out[u] += w;
            base.k_out[v] += w;
            base.k_in[u] += w;
            base.k_in[v] += w;
            if (u == v) base.self_w[u] += 2.0 * w;
            else {
                sym_adj[u][v] += 2.0 * w;
                sym_adj[v][u] += 2.0 * w;
            }
        }
    }
    base.nbrs.resize(N);
    for (int u = 0; u < N; ++u) {
        base.nbrs[u].assign(sym_adj[u].begin(), sym_adj[u].end());
    }

    Rcpp::RNGScope scope;
    TopRng rng;

    double best_Q = -TOP_INF;
    std::vector<int> best_membership(N, 1);

    for (int run = 0; run < nruns; ++run) {
        LouvainGraph cur = base;
        std::vector<int> global_b(N);
        std::iota(global_b.begin(), global_b.end(), 0);

        for (int lev = 0; lev < max_levels; ++lev) {
            int n = cur.n;
            std::vector<int> b(n);
            std::iota(b.begin(), b.end(), 0);
            std::vector<double> K_out = cur.k_out;
            std::vector<double> K_in = cur.k_in;

            std::vector<int> order(n);
            std::iota(order.begin(), order.end(), 0);

            bool any_level_move = false;
            for (int sweep = 0; sweep < 50; ++sweep) {
                rng.shuffle(order);
                bool moved_in_sweep = false;
                for (int u : order) {
                    int r = b[u];
                    double ku_out = cur.k_out[u];
                    double ku_in = cur.k_in[u];

                    // Temporarily remove u from r
                    K_out[r] -= ku_out;
                    K_in[r] -= ku_in;

                    boost::unordered_flat_map<int, double> e_to_c;
                    e_to_c[r] = 0.0;
                    for (const auto& [v, w_sym] : cur.nbrs[u]) {
                        e_to_c[b[v]] += w_sym;
                    }

                    // Gain of placing u into candidate community c (relative to isolated singleton)
                    auto gain = [&](int c, double e_uc) -> double {
                        if (cur.directed) {
                            return e_uc - gamma * (ku_out * K_in[c] + ku_in * K_out[c]) / cur.W;
                        } else {
                            // Note: for undirected, W = 2m, e_uc = 2 * sum_{v in c} w_uv,
                            // and delta(k_c^2 / 2m) = (2 * ku * K_c + ku^2) / 2m;
                            // subtracting singleton (ku^2 / 2m) leaves 2 * ku * K_c / W!
                            return e_uc - gamma * (2.0 * ku_out * K_out[c]) / cur.W;
                        }
                    };

                    double base_gain = gain(r, e_to_c[r]);
                    int best_c = r;
                    double best_diff = 1e-12;
                    for (const auto& [c, e_uc] : e_to_c) {
                        if (c == r) continue;
                        double diff = gain(c, e_uc) - base_gain;
                        if (diff > best_diff) {
                            best_diff = diff;
                            best_c = c;
                        }
                    }

                    b[u] = best_c;
                    K_out[best_c] += ku_out;
                    K_in[best_c] += ku_in;
                    if (best_c != r) {
                        moved_in_sweep = true;
                        any_level_move = true;
                    }
                }
                if (!moved_in_sweep) break;
            }

            if (!any_level_move) break;

            // Relabel communities 0..B_new-1
            std::vector<int> remap(n, -1);
            int B_new = 0;
            for (int u = 0; u < n; ++u) {
                if (remap[b[u]] == -1) remap[b[u]] = B_new++;
                b[u] = remap[b[u]];
            }
            for (int v = 0; v < N; ++v) global_b[v] = b[global_b[v]];
            if (B_new == n || B_new <= 1) break;

            // Coarsen cur -> next LouvainGraph
            LouvainGraph nxt;
            nxt.n = B_new;
            nxt.directed = cur.directed;
            nxt.W = cur.W;
            nxt.self_w.assign(B_new, 0.0);
            nxt.k_out.assign(B_new, 0.0);
            nxt.k_in.assign(B_new, 0.0);
            std::vector<boost::unordered_flat_map<int, double>> nxt_sym(B_new);
            boost::unordered_flat_map<uint64_t, double> nxt_raw;

            for (const auto& [u, v, w] : cur.raw_edges) {
                int ru = b[u], rv = b[v];
                uint64_t k = (static_cast<uint64_t>(static_cast<uint32_t>(ru)) << 32) |
                             static_cast<uint64_t>(static_cast<uint32_t>(rv));
                nxt_raw[k] += w;
            }
            for (const auto& [k, w] : nxt_raw) {
                int ru = static_cast<int>(k >> 32);
                int rv = static_cast<int>(k & 0xffffffffULL);
                nxt.raw_edges.emplace_back(ru, rv, w);
                if (nxt.directed) {
                    nxt.k_out[ru] += w;
                    nxt.k_in[rv] += w;
                    if (ru == rv) nxt.self_w[ru] += w;
                    else { nxt_sym[ru][rv] += w; nxt_sym[rv][ru] += w; }
                } else {
                    nxt.k_out[ru] += w;
                    nxt.k_out[rv] += w;
                    nxt.k_in[ru] += w;
                    nxt.k_in[rv] += w;
                    if (ru == rv) nxt.self_w[ru] += 2.0 * w;
                    else { nxt_sym[ru][rv] += 2.0 * w; nxt_sym[rv][ru] += 2.0 * w; }
                }
            }
            nxt.nbrs.resize(B_new);
            for (int r = 0; r < B_new; ++r) {
                nxt.nbrs[r].assign(nxt_sym[r].begin(), nxt_sym[r].end());
            }
            cur = std::move(nxt);
        }

        // Final fine-grain vertex refinement on original base graph
        {
            int B_cur = 0;
            for (int v = 0; v < N; ++v) B_cur = std::max(B_cur, global_b[v] + 1);
            std::vector<double> K_out(B_cur, 0.0), K_in(B_cur, 0.0);
            for (int v = 0; v < N; ++v) {
                K_out[global_b[v]] += base.k_out[v];
                K_in[global_b[v]] += base.k_in[v];
            }
            std::vector<int> order(N);
            std::iota(order.begin(), order.end(), 0);
            for (int sw = 0; sw < 15; ++sw) {
                rng.shuffle(order);
                bool moved = false;
                for (int u : order) {
                    int r = global_b[u];
                    double ku_out = base.k_out[u], ku_in = base.k_in[u];
                    K_out[r] -= ku_out;
                    K_in[r] -= ku_in;
                    boost::unordered_flat_map<int, double> e_to_c;
                    e_to_c[r] = 0.0;
                    for (const auto& [v, w_sym] : base.nbrs[u]) e_to_c[global_b[v]] += w_sym;
                    auto gain = [&](int c, double e_uc) {
                        return base.directed
                            ? (e_uc - gamma * (ku_out * K_in[c] + ku_in * K_out[c]) / base.W)
                            : (e_uc - gamma * (2.0 * ku_out * K_out[c]) / base.W);
                    };
                    double bg = gain(r, e_to_c[r]);
                    int best_c = r;
                    double best_d = 1e-12;
                    for (const auto& [c, e_uc] : e_to_c) {
                        if (c == r) continue;
                        double d = gain(c, e_uc) - bg;
                        if (d > best_d) { best_d = d; best_c = c; }
                    }
                    global_b[u] = best_c;
                    K_out[best_c] += ku_out;
                    K_in[best_c] += ku_in;
                    if (best_c != r) moved = true;
                }
                if (!moved) break;
            }
        }

        // Canonicalize 1..B and evaluate exact modularity
        std::vector<int> remap(N, -1);
        int B_fin = 0;
        Rcpp::IntegerVector b_cand(N);
        for (int v = 0; v < N; ++v) {
            if (remap[global_b[v]] == -1) remap[global_b[v]] = ++B_fin;
            b_cand[v] = remap[global_b[v]];
        }
        double Q = cpp_modularity(pointer, b_cand, gamma, weights);
        if (Q > best_Q) {
            best_Q = Q;
            best_membership.assign(b_cand.begin(), b_cand.end());
        }
    }

    int B_out = 0;
    for (int x : best_membership) B_out = std::max(B_out, x);
    return Rcpp::List::create(
        Rcpp::Named("blocks") = Rcpp::wrap(best_membership),
        Rcpp::Named("modularity") = best_Q,
        Rcpp::Named("B") = B_out
    );
}

// ===========================================================================
// 5. Fast Random Graph Generators
// ===========================================================================

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_erdos_renyi_gnp(int n, double p, bool directed, bool self_loops) {
    if (n <= 0 || p <= 0.0) return Rcpp::IntegerMatrix(0, 2);
    Rcpp::RNGScope scope;
    TopRng rng;

    std::vector<int> src, tgt;
    if (p >= 1.0) {
        for (int i = 0; i < n; ++i) {
            int j_start = directed ? 0 : (self_loops ? i : i + 1);
            for (int j = j_start; j < n; ++j) {
                if (!self_loops && i == j) continue;
                src.push_back(i + 1);
                tgt.push_back(j + 1);
            }
        }
    } else {
        // Batagelj-Brandes (2005) O(n + m) geometric skip algorithm
        const double lp = std::log1p(-p);
        int64_t total_pairs;
        if (directed) {
            total_pairs = self_loops ? (int64_t(n) * n) : (int64_t(n) * (n - 1));
        } else {
            total_pairs = self_loops ? (int64_t(n) * (n + 1) / 2) : (int64_t(n) * (n - 1) / 2);
        }

        int64_t idx = -1;
        while (true) {
            double u = std::max(1e-300, 1.0 - rng.unif());
            int64_t skip = static_cast<int64_t>(std::floor(std::log(u) / lp));
            idx += 1 + skip;
            if (idx >= total_pairs) break;

            int row, col;
            if (directed && self_loops) {
                row = static_cast<int>(idx / n);
                col = static_cast<int>(idx % n);
            } else if (directed && !self_loops) {
                row = static_cast<int>(idx / (n - 1));
                col = static_cast<int>(idx % (n - 1));
                if (col >= row) col++;
            } else if (!directed && !self_loops) {
                // Upper triangle i > j: idx = i*(i-1)/2 + j
                int i = static_cast<int>(std::floor((1.0 + std::sqrt(1.0 + 8.0 * static_cast<double>(idx))) * 0.5));
                while (int64_t(i) * (i - 1) / 2 > idx) i--;
                while (int64_t(i + 1) * i / 2 <= idx) i++;
                int j = static_cast<int>(idx - int64_t(i) * (i - 1) / 2);
                row = j;
                col = i;
            } else {
                // Upper triangle with diagonal i >= j: idx = i*(i+1)/2 + j
                int i = static_cast<int>(std::floor((-1.0 + std::sqrt(1.0 + 8.0 * static_cast<double>(idx))) * 0.5));
                while (int64_t(i) * (i + 1) / 2 > idx) i--;
                while (int64_t(i + 1) * (i + 2) / 2 <= idx) i++;
                int j = static_cast<int>(idx - int64_t(i) * (i + 1) / 2);
                row = j;
                col = i;
            }
            src.push_back(row + 1);
            tgt.push_back(col + 1);
        }
    }

    int m_out = static_cast<int>(src.size());
    Rcpp::IntegerMatrix out(m_out, 2);
    for (int k = 0; k < m_out; ++k) {
        out(k, 0) = src[k];
        out(k, 1) = tgt[k];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_erdos_renyi_gnm(int n, int m, bool directed, bool self_loops) {
    if (n <= 0 || m <= 0) return Rcpp::IntegerMatrix(0, 2);
    int64_t total_pairs;
    if (directed) {
        total_pairs = self_loops ? (int64_t(n) * n) : (int64_t(n) * (n - 1));
    } else {
        total_pairs = self_loops ? (int64_t(n) * (n + 1) / 2) : (int64_t(n) * (n - 1) / 2);
    }
    if (m > total_pairs) {
        Rcpp::stop("Requested edge count m (%d) exceeds maximum simple edges (%lld)",
                   m, static_cast<long long>(total_pairs));
    }

    Rcpp::RNGScope scope;
    TopRng rng;

    // Floyd's O(m) sample-without-replacement algorithm
    boost::unordered_flat_set<int64_t> chosen;
    chosen.reserve(static_cast<size_t>(m) * 2);
    for (int64_t j = total_pairs - m; j < total_pairs; ++j) {
        int64_t t = static_cast<int64_t>(rng.index(static_cast<size_t>(j + 1)));
        if (!chosen.insert(t).second) {
            chosen.insert(j);
        }
    }

    Rcpp::IntegerMatrix out(m, 2);
    int k = 0;
    for (int64_t idx : chosen) {
        int row, col;
        if (directed && self_loops) {
            row = static_cast<int>(idx / n);
            col = static_cast<int>(idx % n);
        } else if (directed && !self_loops) {
            row = static_cast<int>(idx / (n - 1));
            col = static_cast<int>(idx % (n - 1));
            if (col >= row) col++;
        } else if (!directed && !self_loops) {
            int i = static_cast<int>(std::floor((1.0 + std::sqrt(1.0 + 8.0 * static_cast<double>(idx))) * 0.5));
            while (int64_t(i) * (i - 1) / 2 > idx) i--;
            while (int64_t(i + 1) * i / 2 <= idx) i++;
            int j = static_cast<int>(idx - int64_t(i) * (i - 1) / 2);
            row = j; col = i;
        } else {
            int i = static_cast<int>(std::floor((-1.0 + std::sqrt(1.0 + 8.0 * static_cast<double>(idx))) * 0.5));
            while (int64_t(i) * (i + 1) / 2 > idx) i--;
            while (int64_t(i + 1) * (i + 2) / 2 <= idx) i++;
            int j = static_cast<int>(idx - int64_t(i) * (i + 1) / 2);
            row = j; col = i;
        }
        out(k, 0) = row + 1;
        out(k, 1) = col + 1;
        k++;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_barabasi_albert(int n, int m, double gamma_exp, bool directed) {
    if (n <= 1 || m <= 0) return Rcpp::IntegerMatrix(0, 2);
    int m0 = std::min(n, std::max(m + 1, 2));

    Rcpp::RNGScope scope;
    TopRng rng;

    std::vector<int> src, tgt;
    std::vector<int> deg(n, 0);

    // Seed clique / star on 0..m0-1
    for (int i = 1; i < m0; ++i) {
        for (int j = 0; j < i; ++j) {
            src.push_back(i + 1);
            tgt.push_back(j + 1);
            if (directed) deg[j]++;
            else { deg[i]++; deg[j]++; }
        }
    }

    // Add vertices m0..n-1 with attachment weight (deg[v]^gamma_exp + 1)
    for (int u = m0; u < n; ++u) {
        int k_pick = std::min(m, u);
        std::vector<double> w(u);
        double wsum = 0.0;
        for (int v = 0; v < u; ++v) {
            w[v] = (deg[v] > 0 ? std::pow(static_cast<double>(deg[v]), gamma_exp) : 0.0) + 1.0;
            wsum += w[v];
        }
        for (int step = 0; step < k_pick; ++step) {
            double r = rng.unif() * wsum;
            double acc = 0.0;
            int chosen = u - 1;
            for (int v = 0; v < u; ++v) {
                if (w[v] <= 0.0) continue;
                acc += w[v];
                if (r <= acc) { chosen = v; break; }
            }
            src.push_back(u + 1);
            tgt.push_back(chosen + 1);
            wsum -= w[chosen];
            w[chosen] = 0.0;
        }
        // Update degrees after selecting all k_pick targets
        for (size_t idx = src.size() - k_pick; idx < src.size(); ++idx) {
            int v = tgt[idx] - 1;
            if (directed) deg[v]++;
            else { deg[u]++; deg[v]++; }
        }
    }

    int me = static_cast<int>(src.size());
    Rcpp::IntegerMatrix out(me, 2);
    for (int i = 0; i < me; ++i) { out(i, 0) = src[i]; out(i, 1) = tgt[i]; }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_watts_strogatz(int n, int k, double p, bool circular) {
    if (n <= 1 || k <= 0) return Rcpp::IntegerMatrix(0, 2);
    int half_k = std::max(1, k / 2);
    if (half_k >= n) half_k = n - 1;

    Rcpp::RNGScope scope;
    TopRng rng;

    auto ekey = [](int a, int b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
               static_cast<uint64_t>(static_cast<uint32_t>(b));
    };

    std::vector<std::pair<int, int>> edges;
    boost::unordered_flat_set<uint64_t> eset;
    for (int u = 0; u < n; ++u) {
        for (int d = 1; d <= half_k; ++d) {
            if (!circular && u + d >= n) continue;
            int v = (u + d) % n;
            if (u == v) continue;
            uint64_t key = ekey(u, v);
            if (eset.insert(key).second) {
                edges.emplace_back(u, v);
            }
        }
    }

    if (p > 0.0) {
        int m = static_cast<int>(edges.size());
        for (int i = 0; i < m; ++i) {
            if (rng.unif() >= p) continue;
            int u = edges[i].first;
            int v_old = edges[i].second;
            for (int tries = 0; tries < 50; ++tries) {
                int w = static_cast<int>(rng.index(static_cast<size_t>(n)));
                if (w == u) continue;
                uint64_t nk = ekey(u, w);
                if (!eset.contains(nk)) {
                    eset.erase(ekey(u, v_old));
                    eset.insert(nk);
                    edges[i] = {u, w};
                    break;
                }
            }
        }
    }

    int me = static_cast<int>(edges.size());
    Rcpp::IntegerMatrix out(me, 2);
    for (int i = 0; i < me; ++i) {
        out(i, 0) = edges[i].first + 1;
        out(i, 1) = edges[i].second + 1;
    }
    return out;
}
