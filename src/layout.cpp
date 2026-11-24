// rgraphtool: native C++ SFDP force-directed layout, BFS traversal, and
// animation simulation engines (SIRS epidemics and spatial-segregation dynamic
// layouts), derived from graph-tool 3.9.
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
// Copyright (C) 2026 Diogo Melo (R port)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>
#include "graph_types.h"

namespace {

using Pos2 = std::array<double, 2>;

inline std::mt19937_64 make_rng_from_r() {
    // Seed from R's RNG so set.seed() gives deterministic results
    uint64_t hi = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
    uint64_t lo = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
    return std::mt19937_64((hi << 32) ^ lo ^ 0x9e3779b97f4a7c15ULL);
}

inline double dist2(const Pos2& a, const Pos2& b) {
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    return std::sqrt(dx * dx + dy * dy);
}

inline double get_diff2(const Pos2& p1, const Pos2& p2, Pos2& r) {
    r[0] = p1[0] - p2[0];
    r[1] = p1[1] - p2[1];
    double d2 = r[0] * r[0] + r[1] * r[1];
    if (d2 == 0.0) {
        r[0] = 0.0;
        r[1] = 0.0;
        return 0.0;
    }
    double d = std::sqrt(d2);
    r[0] /= d;
    r[1] /= d;
    return d;
}

inline double norm2_inplace(Pos2& x) {
    double d2 = x[0] * x[0] + x[1] * x[1];
    if (d2 <= 0.0) {
        x[0] = 0.0;
        x[1] = 0.0;
        return 0.0;
    }
    double d = std::sqrt(d2);
    x[0] /= d;
    x[1] /= d;
    return d;
}

inline double f_repulsive(double C, double K, double p, double d) {
    if (d <= 0.0) return 0.0;
    if (p == 2.0) return -C * K * K * K / (d * d);
    if (p == 1.0) return -C * K * K / d;
    return -C * std::pow(K, 1.0 + p) / std::pow(d, p);
}

inline double f_attractive(double K, double d) {
    return (d * d) / K;
}

class QuadTree2D {
public:
    QuadTree2D(const Pos2& ll, const Pos2& ur, int max_level)
        : ll_(ll), ur_(ur), cm_{0.0, 0.0}, count_(0.0), max_level_(max_level) {
        double dx = ur_[0] - ll_[0];
        double dy = ur_[1] - ll_[1];
        w_ = std::sqrt(dx * dx + dy * dy);
    }

    void put_pos(const Pos2& p, double w = 1.0) {
        count_ += w;
        cm_[0] += p[0] * w;
        cm_[1] += p[1] * w;
        if (max_level_ == 0) {
            dense_leafs_.push_back({p, w});
            return;
        }
        if (leafs_.empty()) {
            leafs_.reserve(4);
            double mx = ll_[0] + 0.5 * (ur_[0] - ll_[0]);
            double my = ll_[1] + 0.5 * (ur_[1] - ll_[1]);
            for (int q = 0; q < 4; ++q) {
                Pos2 lll = ll_, lur = ur_;
                if (q % 2) lll[0] = mx; else lur[0] = mx;
                if (q / 2) lll[1] = my; else lur[1] = my;
                leafs_.emplace_back(lll, lur, max_level_ - 1);
            }
        }
        int ix = (p[0] > ll_[0] + 0.5 * (ur_[0] - ll_[0])) ? 1 : 0;
        int iy = (p[1] > ll_[1] + 0.5 * (ur_[1] - ll_[1])) ? 1 : 0;
        leafs_[ix + 2 * iy].put_pos(p, w);
    }

    void accumulate_repulsion(const Pos2& pv, double C, double K, double p,
                              double theta, Pos2& ftot) const {
        if (count_ <= 0.0) return;
        Pos2 diff{0.0, 0.0};
        if (max_level_ == 0) {
            for (const auto& item : dense_leafs_) {
                double d = get_diff2(item.first, pv, diff);
                if (d == 0.0) continue;
                double f = f_repulsive(C, K, p, d) * item.second;
                ftot[0] += f * diff[0];
                ftot[1] += f * diff[1];
            }
            return;
        }
        Pos2 cm{cm_[0] / count_, cm_[1] / count_};
        double d = get_diff2(cm, pv, diff);
        if (w_ > theta * d && !leafs_.empty()) {
            for (const auto& leaf : leafs_) {
                if (leaf.count_ > 0.0)
                    leaf.accumulate_repulsion(pv, C, K, p, theta, ftot);
            }
        } else if (d > 0.0) {
            double f = f_repulsive(C, K, p, d) * count_;
            ftot[0] += f * diff[0];
            ftot[1] += f * diff[1];
        }
    }

private:
    Pos2 ll_, ur_;
    std::vector<QuadTree2D> leafs_;
    std::vector<std::pair<Pos2, double>> dense_leafs_;
    Pos2 cm_;
    double count_;
    int max_level_;
    double w_;
};

// Build symmetric adjacency list (for undirected spring attraction in SFDP)
std::vector<std::vector<size_t>> build_undirected_adj(const Graph& g) {
    size_t n = boost::num_vertices(g.adjacency);
    std::vector<std::vector<size_t>> adj(n);
    auto range = boost::edges(g.adjacency);
    for (auto it = range.first; it != range.second; ++it) {
        size_t u = boost::source(*it, g.adjacency);
        size_t v = boost::target(*it, g.adjacency);
        if (u == v) continue;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    return adj;
}

void sfdp_layout_core(const std::vector<std::vector<size_t>>& adj,
                      std::vector<Pos2>& pos,
                      double C, double K, double p, double theta,
                      double gamma, double mu, double mu_p,
                      double init_step, double cooling_step,
                      int max_level, double tol, int max_iter, bool simple,
                      const std::vector<int>& group,
                      const std::vector<uint8_t>& pinned,
                      std::mt19937_64& rng) {
    size_t n = pos.size();
    if (n <= 1) return;

    int num_groups = 1;
    for (int g_id : group) {
        if (g_id + 1 > num_groups) num_groups = g_id + 1;
    }

    std::vector<Pos2> group_cm(static_cast<size_t>(num_groups), Pos2{0.0, 0.0});
    std::vector<double> group_size(static_cast<size_t>(num_groups), 0.0);
    std::vector<size_t> movable;
    movable.reserve(n);

    Pos2 ll{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Pos2 ur{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};

    for (size_t v = 0; v < n; ++v) {
        if (!pinned[v]) movable.push_back(v);
        ll[0] = std::min(ll[0], pos[v][0]);
        ll[1] = std::min(ll[1], pos[v][1]);
        ur[0] = std::max(ur[0], pos[v][0]);
        ur[1] = std::max(ur[1], pos[v][1]);
        size_t s = static_cast<size_t>(group[v]);
        group_size[s] += 1.0;
        group_cm[s][0] += pos[v][0];
        group_cm[s][1] += pos[v][1];
    }
    if (movable.empty()) return;

    for (size_t s = 0; s < group_size.size(); ++s) {
        if (group_size[s] > 0.0) {
            group_cm[s][0] /= group_size[s];
            group_cm[s][1] /= group_size[s];
        }
    }

    double delta = tol * K + 1.0;
    double E = 0.0;
    double E0 = std::numeric_limits<double>::max();
    int n_iter = 0;
    double step = init_step;
    int progress = 0;

    while (delta > tol * K && (max_iter == 0 || n_iter < max_iter)) {
        delta = 0.0;
        E0 = E;
        E = 0.0;

        Pos2 qll = ll, qur = ur;
        if (qur[0] - qll[0] < 1e-6) { qll[0] -= 0.5e-6; qur[0] += 0.5e-6; }
        if (qur[1] - qll[1] < 1e-6) { qll[1] -= 0.5e-6; qur[1] += 0.5e-6; }

        bool use_qt = (n > 128 && max_level > 0);
        QuadTree2D qt(qll, qur, use_qt ? max_level : 0);
        if (use_qt) {
            for (size_t v = 0; v < n; ++v)
                qt.put_pos(pos[v], 1.0);
        }

        std::shuffle(movable.begin(), movable.end(), rng);

        Pos2 nll{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        Pos2 nur{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
        for (size_t v = 0; v < n; ++v) {
            if (pinned[v]) {
                nll[0] = std::min(nll[0], pos[v][0]);
                nll[1] = std::min(nll[1], pos[v][1]);
                nur[0] = std::max(nur[0], pos[v][0]);
                nur[1] = std::max(nur[1], pos[v][1]);
            }
        }

        for (size_t v : movable) {
            Pos2 ftot{0.0, 0.0};
            Pos2 diff{0.0, 0.0};

            // Global repulsive forces
            if (use_qt) {
                qt.accumulate_repulsion(pos[v], C, K, p, theta, ftot);
            } else {
                for (size_t u = 0; u < n; ++u) {
                    if (u == v) continue;
                    double d = get_diff2(pos[u], pos[v], diff);
                    if (d == 0.0) {
                        // Tiny deterministic nudge to separate coincident vertices
                        double angle = 2.0 * M_PI * R::unif_rand();
                        diff[0] = std::cos(angle);
                        diff[1] = std::sin(angle);
                        double f = f_repulsive(C, K, p, 1e-3 * K);
                        ftot[0] += f * diff[0];
                        ftot[1] += f * diff[1];
                        continue;
                    }
                    double f = f_repulsive(C, K, p, d);
                    ftot[0] += f * diff[0];
                    ftot[1] += f * diff[1];
                }
            }

            // Local attractive forces along edges
            for (size_t u : adj[v]) {
                if (u == v) continue;
                double d = get_diff2(pos[u], pos[v], diff);
                if (d == 0.0) continue;
                double f = f_attractive(K, d);
                ftot[0] += f * diff[0];
                ftot[1] += f * diff[1];
            }

            // Inter-group forces
            size_t gv = static_cast<size_t>(group[v]);
            if (gamma > 0.0) {
                double Kp = K * static_cast<double>(n) * static_cast<double>(n);
                for (size_t s = 0; s < group_size.size(); ++s) {
                    if (group_size[s] <= 0.0 || s == gv) continue;
                    double d = get_diff2(group_cm[s], pos[v], diff);
                    if (d == 0.0) continue;
                    double f = f_attractive(Kp, d) * gamma * group_size[s];
                    ftot[0] += f * diff[0];
                    ftot[1] += f * diff[1];
                }
            } else if (gamma < 0.0) {
                for (size_t s = 0; s < group_size.size(); ++s) {
                    if (group_size[s] <= 0.0 || s == gv) continue;
                    double d = get_diff2(group_cm[s], pos[v], diff);
                    if (d == 0.0) continue;
                    double f = f_repulsive(C, K, p, d) * group_size[s] * std::abs(gamma);
                    ftot[0] += f * diff[0];
                    ftot[1] += f * diff[1];
                }
            }

            // Intra-group attractive forces
            if (mu > 0.0 && group_size[gv] > 1.0) {
                double d = get_diff2(group_cm[gv], pos[v], diff);
                if (d > 0.0) {
                    double Kp = K * std::pow(group_size[gv], mu_p);
                    double f = f_attractive(Kp, d) * mu * group_size[gv];
                    ftot[0] += f * diff[0];
                    ftot[1] += f * diff[1];
                }
            }

            double fmag = norm2_inplace(ftot);
            E += fmag * fmag;

            if (fmag > 0.0) {
                double gs = group_size[gv];
                group_cm[gv][0] = group_cm[gv][0] * gs - pos[v][0];
                group_cm[gv][1] = group_cm[gv][1] * gs - pos[v][1];

                double dx = ftot[0] * step;
                double dy = ftot[1] * step;
                pos[v][0] += dx;
                pos[v][1] += dy;

                group_cm[gv][0] = (group_cm[gv][0] + pos[v][0]) / gs;
                group_cm[gv][1] = (group_cm[gv][1] + pos[v][1]) / gs;

                delta += std::sqrt(dx * dx + dy * dy);
            }

            nll[0] = std::min(nll[0], pos[v][0]);
            nll[1] = std::min(nll[1], pos[v][1]);
            nur[0] = std::max(nur[0], pos[v][0]);
            nur[1] = std::max(nur[1], pos[v][1]);
        }

        n_iter++;
        ll = nll;
        ur = nur;
        delta /= static_cast<double>(movable.size());

        if (simple) {
            step *= cooling_step;
        } else {
            if (E < E0) {
                ++progress;
                if (progress >= 5) {
                    progress = 0;
                    step /= cooling_step;
                }
            } else {
                progress = 0;
                step *= cooling_step;
            }
        }
    }
}

inline uint64_t pack_edge(size_t u, size_t v) {
    if (u > v) std::swap(u, v);
    return (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_sfdp_layout(
    SEXP pointer,
    Rcpp::Nullable<Rcpp::NumericMatrix> init_pos = R_NilValue,
    double C = 0.2,
    double K = 1.0,
    double p = 2.0,
    double theta = 0.6,
    double gamma = 1.0,
    double mu = 0.0,
    double mu_p = 1.0,
    double init_step = -1.0,
    double cooling_step = 0.95,
    int max_level = 10,
    double tol = 1e-4,
    int max_iter = 200,
    bool simple = false,
    Rcpp::Nullable<Rcpp::IntegerVector> groups = R_NilValue,
    Rcpp::Nullable<Rcpp::LogicalVector> pin = R_NilValue) {

    const Graph& g = get_graph(pointer);
    size_t n = boost::num_vertices(g.adjacency);
    Rcpp::NumericMatrix out(static_cast<R_xlen_t>(n), 2);
    Rcpp::colnames(out) = Rcpp::CharacterVector::create("x", "y");
    if (n == 0) return out;
    if (n == 1) {
        if (init_pos.isNotNull()) {
            Rcpp::NumericMatrix ip(init_pos);
            out(0, 0) = ip(0, 0);
            out(0, 1) = ip(0, 1);
        } else {
            out(0, 0) = 0.0;
            out(0, 1) = 0.0;
        }
        return out;
    }

    std::mt19937_64 rng = make_rng_from_r();
    std::vector<Pos2> pos(n);
    if (init_pos.isNotNull()) {
        Rcpp::NumericMatrix ip(init_pos);
        if (static_cast<size_t>(ip.nrow()) != n || ip.ncol() != 2)
            Rcpp::stop("pos must be an N x 2 numeric matrix");
        for (size_t i = 0; i < n; ++i) {
            pos[i][0] = ip(static_cast<R_xlen_t>(i), 0);
            pos[i][1] = ip(static_cast<R_xlen_t>(i), 1);
        }
    } else {
        double scale = K * std::sqrt(static_cast<double>(n));
        for (size_t i = 0; i < n; ++i) {
            pos[i][0] = (R::unif_rand() - 0.5) * scale;
            pos[i][1] = (R::unif_rand() - 0.5) * scale;
        }
    }

    std::vector<int> group(n, 0);
    if (groups.isNotNull()) {
        Rcpp::IntegerVector gv(groups);
        if (static_cast<size_t>(gv.size()) != n)
            Rcpp::stop("groups must have length equal to graph_num_vertices(graph)");
        // Map arbitrary integer labels into 0..(G-1)
        std::map<int, int> gmap;
        int next_id = 0;
        for (size_t i = 0; i < n; ++i) {
            int val = gv[static_cast<R_xlen_t>(i)];
            if (val == NA_INTEGER) Rcpp::stop("groups must not contain NA");
            auto it = gmap.find(val);
            if (it == gmap.end()) {
                gmap[val] = next_id;
                group[i] = next_id++;
            } else {
                group[i] = it->second;
            }
        }
    }

    std::vector<uint8_t> pinned(n, 0);
    if (pin.isNotNull()) {
        Rcpp::LogicalVector pv(pin);
        if (static_cast<size_t>(pv.size()) != n)
            Rcpp::stop("pin must have length equal to graph_num_vertices(graph)");
        for (size_t i = 0; i < n; ++i)
            pinned[i] = (pv[static_cast<R_xlen_t>(i)] == TRUE) ? 1 : 0;
    }

    double step = (init_step > 0.0) ? init_step : std::max(K, 0.1 * K * std::sqrt(static_cast<double>(n)));
    auto adj = build_undirected_adj(g);
    sfdp_layout_core(adj, pos, C, K, p, theta, gamma, mu, mu_p,
                     step, cooling_step, max_level, tol, max_iter, simple,
                     group, pinned, rng);

    for (size_t i = 0; i < n; ++i) {
        out(static_cast<R_xlen_t>(i), 0) = pos[i][0];
        out(static_cast<R_xlen_t>(i), 1) = pos[i][1];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_graph_bfs(SEXP pointer, int source, int max_dist = -1, bool directed = false) {
    const Graph& g = get_graph(pointer);
    size_t n = boost::num_vertices(g.adjacency);
    if (n == 0) Rcpp::stop("Graph has no vertices");
    if (source < 1 || static_cast<size_t>(source) > n)
        Rcpp::stop("source vertex must be between 1 and graph_num_vertices(graph)");

    bool use_directed = directed && g.directed;
    std::vector<std::vector<size_t>> adj(n);
    auto erange = boost::edges(g.adjacency);
    std::vector<std::pair<size_t, size_t>> edge_list;
    edge_list.reserve(boost::num_edges(g.adjacency));
    for (auto it = erange.first; it != erange.second; ++it) {
        size_t u = boost::source(*it, g.adjacency);
        size_t v = boost::target(*it, g.adjacency);
        edge_list.push_back({u, v});
        if (u == v) continue;
        adj[u].push_back(v);
        if (!use_directed) adj[v].push_back(u);
    }

    size_t s = static_cast<size_t>(source - 1);
    std::vector<int> dist(n, -1);
    std::vector<int> parent(n, NA_INTEGER);
    std::vector<int> order;
    order.reserve(n);

    std::queue<size_t> q;
    dist[s] = 0;
    q.push(s);
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        order.push_back(static_cast<int>(u + 1));
        if (max_dist >= 0 && dist[u] >= max_dist) continue;
        for (size_t v : adj[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                parent[v] = static_cast<int>(u + 1);
                q.push(v);
            }
        }
    }

    Rcpp::NumericVector dist_out(static_cast<R_xlen_t>(n));
    Rcpp::IntegerVector parent_out(static_cast<R_xlen_t>(n));
    for (size_t i = 0; i < n; ++i) {
        dist_out[static_cast<R_xlen_t>(i)] = (dist[i] >= 0) ? static_cast<double>(dist[i]) : R_PosInf;
        parent_out[static_cast<R_xlen_t>(i)] = parent[i];
    }

    // Identify tree edges aligned with graph_edges(graph)
    size_t m = edge_list.size();
    Rcpp::LogicalVector tree_edges(static_cast<R_xlen_t>(m), false);
    std::vector<uint8_t> matched_child(n, 0);
    for (size_t idx = 0; idx < m; ++idx) {
        size_t u = edge_list[idx].first;
        size_t v = edge_list[idx].second;
        if (u == v) continue;
        if (parent[v] != NA_INTEGER && static_cast<size_t>(parent[v] - 1) == u && !matched_child[v]) {
            tree_edges[static_cast<R_xlen_t>(idx)] = true;
            matched_child[v] = 1;
        } else if (!use_directed && parent[u] != NA_INTEGER &&
                   static_cast<size_t>(parent[u] - 1) == v && !matched_child[u]) {
            tree_edges[static_cast<R_xlen_t>(idx)] = true;
            matched_child[u] = 1;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("source") = source,
        Rcpp::Named("max_dist") = (max_dist >= 0) ? Rcpp::wrap(max_dist) : Rcpp::wrap(R_PosInf),
        Rcpp::Named("directed") = use_directed,
        Rcpp::Named("distance") = dist_out,
        Rcpp::Named("parent") = parent_out,
        Rcpp::Named("order") = Rcpp::wrap(order),
        Rcpp::Named("tree_edges") = tree_edges
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_sirs_animation_sim(SEXP pointer,
                                  int n_steps,
                                  double x_spont,
                                  double r_rec,
                                  double s_resusc,
                                  Rcpp::IntegerVector init_state,
                                  bool directed = false) {
    const Graph& g = get_graph(pointer);
    size_t n = boost::num_vertices(g.adjacency);
    if (static_cast<size_t>(init_state.size()) != n)
        Rcpp::stop("init_state must have length equal to graph_num_vertices(graph)");
    if (n_steps < 0)
        Rcpp::stop("n_steps must be non-negative");

    // Neighbor list for infection transmission (in-neighbors if directed, all neighbors if undirected)
    bool use_directed = directed && g.directed;
    std::vector<std::vector<size_t>> nbrs(n);
    auto erange = boost::edges(g.adjacency);
    for (auto it = erange.first; it != erange.second; ++it) {
        size_t u = boost::source(*it, g.adjacency);
        size_t v = boost::target(*it, g.adjacency);
        if (u == v) continue;
        if (use_directed) {
            nbrs[v].push_back(u); // u -> v infects v
        } else {
            nbrs[u].push_back(v);
            nbrs[v].push_back(u);
        }
    }

    std::vector<int> state(n);
    for (size_t i = 0; i < n; ++i) {
        int val = init_state[static_cast<R_xlen_t>(i)];
        if (val < 0 || val > 2)
            Rcpp::stop("init_state values must be 0 (S), 1 (I), or 2 (R)");
        state[i] = val;
    }

    int n_frames = n_steps + 1;
    Rcpp::IntegerMatrix states_mat(static_cast<R_xlen_t>(n), n_frames);
    Rcpp::LogicalMatrix newly_mat(static_cast<R_xlen_t>(n), n_frames);
    Rcpp::LogicalMatrix fear_mat(static_cast<R_xlen_t>(n), n_frames);
    Rcpp::IntegerMatrix counts_mat(n_frames, 3);
    Rcpp::colnames(counts_mat) = Rcpp::CharacterVector::create("S", "I", "R");

    auto record_frame = [&](int t, const std::vector<uint8_t>& newly) {
        int cS = 0, cI = 0, cR = 0;
        for (size_t v = 0; v < n; ++v) {
            states_mat(static_cast<R_xlen_t>(v), t) = state[v];
            newly_mat(static_cast<R_xlen_t>(v), t) = (newly[v] != 0);
            if (state[v] == 0) {
                ++cS;
                bool has_inf = false;
                for (size_t w : nbrs[v]) {
                    if (state[w] == 1) { has_inf = true; break; }
                }
                fear_mat(static_cast<R_xlen_t>(v), t) = has_inf;
            } else if (state[v] == 1) {
                ++cI;
                fear_mat(static_cast<R_xlen_t>(v), t) = false;
            } else {
                ++cR;
                fear_mat(static_cast<R_xlen_t>(v), t) = false;
            }
        }
        counts_mat(t, 0) = cS;
        counts_mat(t, 1) = cI;
        counts_mat(t, 2) = cR;
    };

    std::vector<uint8_t> newly(n, 0);
    record_frame(0, newly);

    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::mt19937_64 rng = make_rng_from_r();

    for (int t = 1; t <= n_steps; ++t) {
        std::fill(newly.begin(), newly.end(), 0);
        std::shuffle(order.begin(), order.end(), rng);
        for (size_t v : order) {
            if (state[v] == 1) {
                if (R::unif_rand() < r_rec) {
                    state[v] = 2;
                }
            } else if (state[v] == 0) {
                if (R::unif_rand() < x_spont) {
                    state[v] = 1;
                    newly[v] = 1;
                } else if (!nbrs[v].empty()) {
                    size_t idx = static_cast<size_t>(std::floor(R::unif_rand() * static_cast<double>(nbrs[v].size())));
                    if (idx >= nbrs[v].size()) idx = nbrs[v].size() - 1;
                    size_t w = nbrs[v][idx];
                    if (state[w] == 1) {
                        state[v] = 1;
                        newly[v] = 1;
                    }
                }
            } else {
                if (R::unif_rand() < s_resusc) {
                    state[v] = 0;
                }
            }
        }
        record_frame(t, newly);
    }

    return Rcpp::List::create(
        Rcpp::Named("states") = states_mat,
        Rcpp::Named("newly_infected") = newly_mat,
        Rcpp::Named("fear") = fear_mat,
        Rcpp::Named("counts") = counts_mat
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_dynamic_layout_sim(SEXP pointer,
                                  Rcpp::NumericMatrix init_pos,
                                  int n_steps = 60,
                                  double K = 0.5,
                                  double C = 0.2,
                                  double p = 2.0,
                                  double step = 0.005,
                                  int n_rewire = 10,
                                  double theta = 0.6,
                                  int max_level = 10) {
    const Graph& g = get_graph(pointer);
    size_t n = boost::num_vertices(g.adjacency);
    if (static_cast<size_t>(init_pos.nrow()) != n || init_pos.ncol() != 2)
        Rcpp::stop("init_pos must be an N x 2 numeric matrix");
    if (n_steps < 0)
        Rcpp::stop("n_steps must be non-negative");

    std::vector<Pos2> pos(n);
    for (size_t i = 0; i < n; ++i) {
        pos[i][0] = init_pos(static_cast<R_xlen_t>(i), 0);
        pos[i][1] = init_pos(static_cast<R_xlen_t>(i), 1);
    }

    // Build simple undirected edge list + fast lookup set + degree vector
    std::vector<std::pair<size_t, size_t>> edges;
    std::unordered_set<uint64_t> edge_set;
    std::vector<size_t> deg(n, 0);

    auto erange = boost::edges(g.adjacency);
    for (auto it = erange.first; it != erange.second; ++it) {
        size_t u = boost::source(*it, g.adjacency);
        size_t v = boost::target(*it, g.adjacency);
        if (u == v) continue;
        uint64_t key = pack_edge(u, v);
        if (edge_set.insert(key).second) {
            edges.push_back({std::min(u, v), std::max(u, v)});
            deg[u]++;
            deg[v]++;
        }
    }

    int n_frames = n_steps + 1;
    Rcpp::List pos_list(n_frames);
    Rcpp::List edge_list(n_frames);
    Rcpp::NumericVector mean_len(n_frames);
    Rcpp::IntegerVector rewired_vec(n_frames, 0);

    auto record_frame = [&](int t, int n_acc) {
        Rcpp::NumericMatrix pm(static_cast<R_xlen_t>(n), 2);
        Rcpp::colnames(pm) = Rcpp::CharacterVector::create("x", "y");
        for (size_t i = 0; i < n; ++i) {
            pm(static_cast<R_xlen_t>(i), 0) = pos[i][0];
            pm(static_cast<R_xlen_t>(i), 1) = pos[i][1];
        }
        Rcpp::IntegerMatrix em(static_cast<R_xlen_t>(edges.size()), 2);
        Rcpp::colnames(em) = Rcpp::CharacterVector::create("from", "to");
        double total_len = 0.0;
        for (size_t idx = 0; idx < edges.size(); ++idx) {
            size_t u = edges[idx].first;
            size_t v = edges[idx].second;
            em(static_cast<R_xlen_t>(idx), 0) = static_cast<int>(u + 1);
            em(static_cast<R_xlen_t>(idx), 1) = static_cast<int>(v + 1);
            total_len += dist2(pos[u], pos[v]);
        }
        pos_list[t] = pm;
        edge_list[t] = em;
        mean_len[t] = edges.empty() ? 0.0 : (total_len / static_cast<double>(edges.size()));
        rewired_vec[t] = n_acc;
    };

    record_frame(0, 0);

    std::vector<int> group(n, 0);
    std::vector<uint8_t> pinned(n, 0);
    std::mt19937_64 rng = make_rng_from_r();

    for (int t = 1; t <= n_steps; ++t) {
        // 1. One incremental SFDP step
        std::vector<std::vector<size_t>> adj(n);
        for (const auto& e : edges) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
        sfdp_layout_core(adj, pos, C, K, p, theta, 0.0, 0.0, 1.0,
                         step, 1.0, max_level, 0.0, 1, true,
                         group, pinned, rng);

        // 2. Spatial-segregation distance-shortening rewiring (animation_dancing.py)
        int accepted = 0;
        if (n >= 3 && !edges.empty() && n_rewire > 0) {
            for (int r = 0; r < n_rewire; ++r) {
                size_t eidx = static_cast<size_t>(std::floor(R::unif_rand() * static_cast<double>(edges.size())));
                if (eidx >= edges.size()) eidx = edges.size() - 1;
                size_t s1 = edges[eidx].first;
                size_t t1 = edges[eidx].second;
                if (R::unif_rand() < 0.5) std::swap(s1, t1);
                if (deg[t1] <= 1) continue;

                size_t t2 = static_cast<size_t>(std::floor(R::unif_rand() * static_cast<double>(n)));
                if (t2 >= n) t2 = n - 1;
                if (t2 == s1 || t2 == t1) continue;

                uint64_t new_key = pack_edge(s1, t2);
                if (edge_set.find(new_key) != edge_set.end()) continue;

                double d_old = dist2(pos[s1], pos[t1]);
                double d_new = dist2(pos[s1], pos[t2]);
                if (d_new <= d_old) {
                    uint64_t old_key = pack_edge(s1, t1);
                    edge_set.erase(old_key);
                    edge_set.insert(new_key);
                    edges[eidx] = {std::min(s1, t2), std::max(s1, t2)};
                    deg[t1]--;
                    deg[t2]++;
                    accepted++;
                }
            }
        }

        // Re-center positions to keep viewport stable
        if (n > 0) {
            double cx = 0.0, cy = 0.0;
            for (size_t i = 0; i < n; ++i) {
                cx += pos[i][0];
                cy += pos[i][1];
            }
            cx /= static_cast<double>(n);
            cy /= static_cast<double>(n);
            for (size_t i = 0; i < n; ++i) {
                pos[i][0] -= cx;
                pos[i][1] -= cy;
            }
        }

        record_frame(t, accepted);
    }

    return Rcpp::List::create(
        Rcpp::Named("positions") = pos_list,
        Rcpp::Named("edges") = edge_list,
        Rcpp::Named("mean_edge_length") = mean_len,
        Rcpp::Named("rewired") = rewired_vec
    );
}
