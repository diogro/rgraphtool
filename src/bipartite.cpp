// rgraphtool: native C++23 flat, nested, and weighted bipartite stochastic
// block models (Yen & Larremore 2020, Phys. Rev. E 102, 032309; Peixoto 2018;
// graph-tool 3.9 pclabel parity).
//
// Copyright (C) 2026 Diogo Melo.
// Algorithms follow graph-tool by Tiago de Paula Peixoto and bipartiteSBM by
// Tzu-Chi Yen & Daniel B. Larremore. LGPL-3.0-or-later.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include <boost/unordered/unordered_flat_map.hpp>
#include "graph_types.h"
#include "parallel.h"
#include "logq.h"

namespace {

template <class K, class V, class H = boost::hash<K>>
using bp_hash_map = boost::unordered_flat_map<K, V, H>;

constexpr double BP_INF = std::numeric_limits<double>::infinity();

struct BpRng {
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

struct BpMath {
    std::vector<double> lg, lg1, ln;
    explicit BpMath(size_t max_n = 262144) { ensure(max_n); }
    void ensure(size_t n) {
        if (n < lg.size()) return;
        size_t sz = std::max(n + 1, lg.empty() ? size_t(4096) : lg.size() * 2);
        size_t old = lg.size();
        lg.resize(sz);
        lg1.resize(sz);
        ln.resize(sz);
        for (size_t i = old; i < sz; ++i) {
            lg[i] = (i == 0) ? 0.0 : std::lgamma(static_cast<double>(i));
            lg1[i] = std::lgamma(static_cast<double>(i) + 1.0);
            ln[i] = (i == 0) ? 0.0 : std::log(static_cast<double>(i));
        }
    }
    inline double lgam1(int64_t n) const {
        if (n < 0) return 0.0;
        if (static_cast<size_t>(n) < lg1.size()) return lg1[n];
        return std::lgamma(static_cast<double>(n) + 1.0);
    }
    inline double logi(int64_t n) const {
        if (n <= 1) return 0.0;
        if (static_cast<size_t>(n) < ln.size()) return ln[n];
        return std::log(static_cast<double>(n));
    }
    inline double lbinom(int64_t n, int64_t k) const {
        if (k <= 0 || k >= n) return 0.0;
        if (n < k) return 0.0;
        return lgam1(n) - lgam1(k) - lgam1(n - k);
    }
};

BpMath& bp_math() {
    static BpMath m;
    return m;
}

inline uint64_t bp_pair_key(int r, int s) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(r)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(s));
}

inline double bp_safe_lgamma(double x) {
    return (x <= 0.0) ? 0.0 : std::lgamma(x);
}

inline double bp_lbinom_d(double n, double k) {
    if (k <= 0.0 || k >= n || n < k) return 0.0;
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

inline double bp_lbeta(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

// Shared log q(n, k) table (src/logq.h).
using BpLogQ = rgt_logq::LogQ;
inline BpLogQ& bp_log_q() { return rgt_logq::shared(); }

// ---------------------------------------------------------------------------
// Weighted covariate specifications (exact Peixoto 2018 / graph-tool 3.9)

enum class BpCovType : int {
    REAL_EXPONENTIAL = 0,
    REAL_NORMAL = 1,
    DISCRETE_GEOMETRIC = 2,
    DISCRETE_POISSON = 3,
    DISCRETE_BINOMIAL = 4
};

struct BpCovSpec {
    BpCovType type;
    bool micro = true;
    double p0 = NAN, p1 = NAN, p2 = NAN, p3 = NAN;
    double eps = 1e-8;
    double const_term = 0.0;
};

inline double bp_cov_neg_log_P(const BpCovSpec& c, double m, double x, double x2) {
    if (m <= 0.0) return 0.0;
    switch (c.type) {
    case BpCovType::REAL_EXPONENTIAL: {
        if (std::isnan(c.p0) && std::isnan(c.p1)) {
            if (x < c.eps || m <= 1.0) return 0.0;
            return -(bp_safe_lgamma(m) - (m - 1.0) * std::log(x));
        }
        return -(bp_safe_lgamma(m + c.p0) - bp_safe_lgamma(c.p0) +
                 c.p0 * std::log(c.p1) - (c.p0 + m) * std::log(c.p1 + x));
    }
    case BpCovType::REAL_NORMAL: {
        if (std::isnan(c.p0) && std::isnan(c.p1)) {
            double smu1 = x * (x / m);
            if (m < 2.0 || smu1 >= x2 || (x2 - smu1) < c.eps * c.eps) return 0.0;
            return -(bp_safe_lgamma((m - 1.0) / 2.0) + 0.5 * std::log(m) -
                     ((static_cast<int>(m) - 3) / 2.0) * std::log(x2 - smu1) -
                     ((m - 1.0) / 2.0) * std::log(M_PI));
        }
        double v = std::max(0.0, x2 - x * (x / m));
        double k_n = c.p1 + m;
        double nu_n = c.p3 + m;
        double diff = c.p0 - x / m;
        double nu_v_n = c.p2 * c.p3 + v + ((m * c.p1) / k_n) * diff * diff;
        return -(bp_safe_lgamma(nu_n / 2.0) - bp_safe_lgamma(c.p3 / 2.0) +
                 0.5 * (std::log(c.p1) - std::log(k_n)) +
                 0.5 * c.p3 * std::log(c.p3 * c.p2) - 0.5 * nu_n * std::log(nu_v_n) -
                 0.5 * m * std::log(M_PI));
    }
    case BpCovType::DISCRETE_GEOMETRIC: {
        if (std::isnan(c.p0) && std::isnan(c.p1))
            return bp_lbinom_d((m - 1.0) + x, x);
        return -(bp_lbeta(m + c.p0, x + c.p1) - bp_lbeta(c.p0, c.p1));
    }
    case BpCovType::DISCRETE_POISSON: {
        if (std::isnan(c.p0) && std::isnan(c.p1))
            return -(bp_safe_lgamma(x + 1.0) - x * std::log(m));
        return -(bp_safe_lgamma(x + c.p0) - (x + c.p0) * std::log(m + c.p1) -
                 bp_safe_lgamma(c.p0) + c.p0 * std::log(c.p1));
    }
    case BpCovType::DISCRETE_BINOMIAL: {
        if (std::isnan(c.p1) && std::isnan(c.p2))
            return bp_lbinom_d(m * c.p0, x);
        return -(bp_lbeta(x + c.p1, m * c.p0 - x + c.p2) - bp_lbeta(c.p1, c.p2));
    }
    }
    return 0.0;
}

std::vector<BpCovSpec> parse_cov_specs(const Rcpp::CharacterVector& types,
                                       const Rcpp::List& params,
                                       const Rcpp::NumericMatrix& W) {
    int M = types.size();
    std::vector<BpCovSpec> specs(M);
    int E = W.nrow();
    for (int m = 0; m < M; ++m) {
        std::string tp = Rcpp::as<std::string>(types[m]);
        BpCovSpec s;
        if (tp == "real-exponential") s.type = BpCovType::REAL_EXPONENTIAL;
        else if (tp == "real-normal") s.type = BpCovType::REAL_NORMAL;
        else if (tp == "discrete-geometric") s.type = BpCovType::DISCRETE_GEOMETRIC;
        else if (tp == "discrete-poisson") s.type = BpCovType::DISCRETE_POISSON;
        else if (tp == "discrete-binomial") s.type = BpCovType::DISCRETE_BINOMIAL;
        else Rcpp::stop("Unknown covariate type: %s", tp.c_str());

        Rcpp::NumericVector p = params[m];
        if (p.size() >= 1) s.p0 = p[0];
        if (p.size() >= 2) s.p1 = p[1];
        if (p.size() >= 3) s.p2 = p[2];
        if (p.size() >= 4) s.p3 = p[3];

        double min_nz = BP_INF;
        double cterm = 0.0;
        for (int i = 0; i < E; ++i) {
            double w = W(i, m);
            if (w != 0.0 && std::abs(w) < min_nz) min_nz = std::abs(w);
            if (s.type == BpCovType::DISCRETE_POISSON) {
                cterm += bp_safe_lgamma(w + 1.0);
            } else if (s.type == BpCovType::DISCRETE_BINOMIAL) {
                cterm -= bp_lbinom_d(s.p0, w);
            }
        }
        s.eps = std::isfinite(min_nz) ? (0.1 * min_nz) : 1e-8;
        s.const_term = cterm;
        specs[m] = s;
    }
    return specs;
}

// ---------------------------------------------------------------------------
// Compact bipartite graph representation (supports both base graphs and
// upper-level block multigraphs with edge multiplicities and optional weights)

struct BpEdgeEntry {
    int nbr;
    int64_t mult;
    int w_offset; // index into edge weight arrays (only for M > 0)
};

struct BpGraphData {
    int N = 0;
    int N_type[2] = {0, 0};
    int64_t E = 0;
    bool directed = false;
    std::vector<int> vtype;          // length N, values in {0, 1}
    std::vector<int64_t> deg_out;    // length N (total degree if undirected)
    std::vector<int64_t> deg_in;     // length N (unused if undirected)
    std::vector<std::vector<BpEdgeEntry>> out_adj; // for undirected: all incident edges
    std::vector<std::vector<BpEdgeEntry>> in_adj;  // only used if directed
    double edge_factor = 0.0;        // sum lgamma(A_ij + 1)
    double degree_factor = 0.0;      // sum lgamma(k_i + 1) (or in + out)

    // Optional covariates (M covariates per edge)
    int M = 0;
    std::vector<BpCovSpec> cov_specs;
    std::vector<double> edge_sum_w;   // flat [edge_id * M + m]
    std::vector<double> edge_sum_w2;  // flat [edge_id * M + m]
};

BpGraphData build_bp_graph(const Graph& g, const Rcpp::IntegerVector& type_vec,
                           const Rcpp::NumericMatrix& W,
                           const Rcpp::CharacterVector& cov_types,
                           const Rcpp::List& cov_params) {
    BpGraphData d;
    d.N = static_cast<int>(boost::num_vertices(g.adjacency));
    d.directed = g.directed;
    if (type_vec.size() != d.N)
        Rcpp::stop("type vector length must equal number of vertices");
    d.vtype.resize(d.N);
    for (int i = 0; i < d.N; ++i) {
        int t = type_vec[i];
        if (t != 0 && t != 1) Rcpp::stop("type values must be 0 or 1");
        d.vtype[i] = t;
        d.N_type[t]++;
    }

    d.M = cov_types.size();
    if (d.M > 0) {
        d.cov_specs = parse_cov_specs(cov_types, cov_params, W);
    }

    d.deg_out.assign(d.N, 0);
    d.deg_in.assign(d.N, 0);
    d.out_adj.resize(d.N);
    if (d.directed) d.in_adj.resize(d.N);

    // Aggregate parallel edges between same (u, v) for faster neighbor loops
    bp_hash_map<uint64_t, int> pair_to_id;
    std::vector<std::pair<int, int>> uniq_edges;
    std::vector<int64_t> uniq_mult;

    auto erange = boost::edges(g.adjacency);
    int e_pos = 0;
    for (auto it = erange.first; it != erange.second; ++it, ++e_pos) {
        int u = static_cast<int>(it->s);
        int v = static_cast<int>(it->t);
        if (d.vtype[u] == d.vtype[v]) {
            Rcpp::stop("Graph is not bipartite: edge (%d, %d) connects same-type vertices", u + 1, v + 1);
        }
        int cu = u, cv = v;
        if (!d.directed && d.vtype[cu] == 1) std::swap(cu, cv);
        uint64_t key = bp_pair_key(cu, cv);
        auto found = pair_to_id.find(key);
        int id;
        if (found == pair_to_id.end()) {
            id = static_cast<int>(uniq_edges.size());
            pair_to_id.emplace(key, id);
            uniq_edges.emplace_back(cu, cv);
            uniq_mult.push_back(1);
            if (d.M > 0) {
                for (int m = 0; m < d.M; ++m) {
                    double w = W(e_pos, m);
                    d.edge_sum_w.push_back(w);
                    d.edge_sum_w2.push_back(w * w);
                }
            }
        } else {
            id = found->second;
            uniq_mult[id]++;
            if (d.M > 0) {
                for (int m = 0; m < d.M; ++m) {
                    double w = W(e_pos, m);
                    d.edge_sum_w[id * d.M + m] += w;
                    d.edge_sum_w2[id * d.M + m] += w * w;
                }
            }
        }
        d.E++;
        if (d.directed) {
            d.deg_out[u]++;
            d.deg_in[v]++;
        } else {
            d.deg_out[u]++;
            d.deg_out[v]++;
        }
    }

    bp_math().ensure(static_cast<size_t>(std::max<int64_t>(d.N + 2, d.E + 2)));

    for (size_t id = 0; id < uniq_edges.size(); ++id) {
        int u = uniq_edges[id].first;
        int v = uniq_edges[id].second;
        int64_t mult = uniq_mult[id];
        d.edge_factor += bp_math().lgam1(mult);
        int woff = (d.M > 0) ? static_cast<int>(id) * d.M : -1;
        d.out_adj[u].push_back({v, mult, woff});
        if (d.directed) {
            d.in_adj[v].push_back({u, mult, woff});
        } else {
            d.out_adj[v].push_back({u, mult, woff});
        }
    }

    for (int i = 0; i < d.N; ++i) {
        if (d.directed) {
            d.degree_factor += bp_math().lgam1(d.deg_out[i]) + bp_math().lgam1(d.deg_in[i]);
        } else {
            d.degree_factor += bp_math().lgam1(d.deg_out[i]);
        }
    }
    return d;
}

// ---------------------------------------------------------------------------
// Pair statistics for block-pair (r, s)

struct PairStat {
    int64_t e = 0;
    std::vector<double> sw;
    std::vector<double> sw2;
    PairStat() = default;
    explicit PairStat(int M) : e(0), sw(M, 0.0), sw2(M, 0.0) {}
};

// ---------------------------------------------------------------------------
// Flat Bipartite Block State (with exact O(deg(v)) incremental move_delta)

class BipartiteBlockState {
public:
    const BpGraphData* g;
    bool dc;
    bool dist;
    bool bp_edge_prior;
    bool include_edges_dl;
    bool multiset_adj;

    int B_cap = 0;
    std::vector<int> b;          // vertex -> block id in [0, B_cap)
    std::vector<int> btype;      // block -> type in {0, 1} (-1 if free)
    std::vector<int64_t> wr;     // block vertex count
    std::vector<int64_t> mr_out; // block out-degree sum (total if undirected)
    std::vector<int64_t> mr_in;  // block in-degree sum (directed only)
    int B_type[2] = {0, 0};      // active blocks per type
    std::vector<int> free_blocks;
    std::vector<int> active_by_type[2];
    std::vector<int> active_pos; // block -> index in active_by_type[btype[r]]

    // Block-pair edge counts & covariate sums.
    // Undirected: key is (r_type0, s_type1).
    // Directed:   key is (r_src, s_tgt) with btype[r_src] != btype[s_tgt].
    bp_hash_map<uint64_t, PairStat> mrs;

    // Per-block adjacency lists of neighbor blocks (for fast merge & proposal)
    std::vector<bp_hash_map<int, int64_t>> blk_nbrs;

    // Degree-distribution histograms per block (for DC distributed prior)
    // Undirected: deg_hist_out stores total degree counts.
    // Directed:   deg_hist_out stores out-degree counts, deg_hist_in stores in-degree counts.
    std::vector<bp_hash_map<int64_t, int64_t>> deg_hist_out;
    std::vector<bp_hash_map<int64_t, int64_t>> deg_hist_in;

    BipartiteBlockState(const BpGraphData* g_, const std::vector<int>& init_b,
                        bool dc_, bool dist_, bool bp_edge_prior_,
                        bool include_edges_dl_ = true,
                        bool multiset_adj_ = false)
        : g(g_), dc(dc_), dist(dist_), bp_edge_prior(bp_edge_prior_),
          include_edges_dl(include_edges_dl_), multiset_adj(multiset_adj_) {
        init_from_blocks(init_b);
    }

    inline uint64_t pair_key_for(int r, int s) const {
        if (!g->directed && btype[r] == 1) return bp_pair_key(s, r);
        return bp_pair_key(r, s);
    }

    void init_from_blocks(const std::vector<int>& init_b) {
        int N = g->N;
        B_cap = std::max(N + 2, 4);
        int max_b = 0;
        for (int v = 0; v < N; ++v) max_b = std::max(max_b, init_b[v]);
        if (max_b + 2 > B_cap) B_cap = max_b + 2;

        b.assign(N, 0);
        btype.assign(B_cap, -1);
        wr.assign(B_cap, 0);
        mr_out.assign(B_cap, 0);
        mr_in.assign(B_cap, 0);
        active_pos.assign(B_cap, -1);
        active_by_type[0].clear();
        active_by_type[1].clear();
        B_type[0] = B_type[1] = 0;
        mrs.clear();
        blk_nbrs.assign(B_cap, {});
        deg_hist_out.assign(B_cap, {});
        if (g->directed) deg_hist_in.assign(B_cap, {});
        free_blocks.clear();

        // Validate and assign vertex blocks
        for (int v = 0; v < N; ++v) {
            int r = init_b[v];
            int t = g->vtype[v];
            if (btype[r] == -1) {
                btype[r] = t;
                activate_block(r, t);
            } else if (btype[r] != t) {
                Rcpp::stop("Block %d contains vertices of both bipartite types", r + 1);
            }
            b[v] = r;
            wr[r]++;
            mr_out[r] += g->deg_out[v];
            if (g->directed) mr_in[r] += g->deg_in[v];
            if (dc && dist) {
                deg_hist_out[r][g->deg_out[v]]++;
                if (g->directed) deg_hist_in[r][g->deg_in[v]]++;
            }
        }

        for (int r = B_cap - 1; r >= 0; --r) {
            if (wr[r] == 0) {
                btype[r] = -1;
                free_blocks.push_back(r);
            }
        }

        // Populate mrs and blk_nbrs
        for (int u = 0; u < N; ++u) {
            if (!g->directed && g->vtype[u] == 1) continue; // visit each undirected edge once
            int ru = b[u];
            for (const auto& e : g->out_adj[u]) {
                int v = e.nbr;
                int rv = b[v];
                uint64_t key = pair_key_for(ru, rv);
                auto it = mrs.find(key);
                if (it == mrs.end()) {
                    it = mrs.emplace(key, PairStat(g->M)).first;
                }
                it->second.e += e.mult;
                for (int m = 0; m < g->M; ++m) {
                    it->second.sw[m] += g->edge_sum_w[e.w_offset + m];
                    it->second.sw2[m] += g->edge_sum_w2[e.w_offset + m];
                }
                blk_nbrs[ru][rv] += e.mult;
                blk_nbrs[rv][ru] += e.mult;
            }
        }
    }

    void activate_block(int r, int t) {
        btype[r] = t;
        active_pos[r] = static_cast<int>(active_by_type[t].size());
        active_by_type[t].push_back(r);
        B_type[t]++;
    }

    void deactivate_block(int r) {
        int t = btype[r];
        int pos = active_pos[r];
        int last = active_by_type[t].back();
        active_by_type[t][pos] = last;
        active_pos[last] = pos;
        active_by_type[t].pop_back();
        active_pos[r] = -1;
        btype[r] = -1;
        B_type[t]--;
        free_blocks.push_back(r);
    }

    int pop_free_block() {
        if (free_blocks.empty()) {
            int old_cap = B_cap;
            B_cap *= 2;
            btype.resize(B_cap, -1);
            wr.resize(B_cap, 0);
            mr_out.resize(B_cap, 0);
            mr_in.resize(B_cap, 0);
            active_pos.resize(B_cap, -1);
            blk_nbrs.resize(B_cap);
            deg_hist_out.resize(B_cap);
            if (g->directed) deg_hist_in.resize(B_cap);
            for (int r = B_cap - 1; r >= old_cap; --r) free_blocks.push_back(r);
        }
        int r = free_blocks.back();
        free_blocks.pop_back();
        return r;
    }

    inline double edges_dl(int B0, int B1) const {
        if (!include_edges_dl || g->E == 0) return 0.0;
        if (bp_edge_prior) {
            int64_t slots = g->directed ? (2LL * B0 * B1) : (int64_t(B0) * B1);
            if (slots <= 0) return 0.0;
            return bp_lbinom_d(static_cast<double>(slots + g->E - 1),
                               static_cast<double>(g->E));
        } else {
            int64_t B = int64_t(B0) + B1;
            int64_t slots = g->directed ? (B * B) : (B * (B + 1) / 2);
            if (slots <= 0) return 0.0;
            return bp_lbinom_d(static_cast<double>(slots + g->E - 1),
                               static_cast<double>(g->E));
        }
    }

    inline double part_dl_type(int c, int Bc) const {
        int Nc = g->N_type[c];
        if (Nc == 0 || Bc == 0) return 0.0;
        return bp_math().logi(Nc) + bp_math().lbinom(Nc - 1, Bc - 1) + bp_math().lgam1(Nc);
    }

    inline double block_deg_dl(int r) const {
        if (!dc || wr[r] == 0) return 0.0;
        if (dist) {
            double s = bp_log_q()(mr_out[r], wr[r]) + bp_math().lgam1(wr[r]);
            for (const auto& kv : deg_hist_out[r]) {
                if (kv.second > 1) s -= bp_math().lgam1(kv.second);
            }
            if (g->directed) {
                s += bp_log_q()(mr_in[r], wr[r]) + bp_math().lgam1(wr[r]);
                for (const auto& kv : deg_hist_in[r]) {
                    if (kv.second > 1) s -= bp_math().lgam1(kv.second);
                }
            }
            return s;
        } else {
            if (g->directed) {
                return bp_math().lbinom(wr[r] + mr_out[r] - 1, mr_out[r]) +
                       bp_math().lbinom(wr[r] + mr_in[r] - 1, mr_in[r]);
            }
            return bp_math().lbinom(wr[r] + mr_out[r] - 1, mr_out[r]);
        }
    }

    inline double block_lik_term(int64_t nr, int64_t er_out, int64_t er_in) const {
        if (nr == 0) return 0.0;
        if (dc) {
            if (g->directed) return bp_math().lgam1(er_out) + bp_math().lgam1(er_in);
            return bp_math().lgam1(er_out);
        } else {
            int64_t etot = g->directed ? (er_out + er_in) : er_out;
            return (etot > 0) ? (static_cast<double>(etot) * bp_math().logi(nr)) : 0.0;
        }
    }

    inline double pair_entropy(const PairStat& ps) const {
        if (ps.e <= 0) return 0.0;
        double s = -bp_math().lgam1(ps.e);
        for (int m = 0; m < g->M; ++m) {
            s += bp_cov_neg_log_P(g->cov_specs[m], static_cast<double>(ps.e),
                                  ps.sw[m], ps.sw2[m]);
        }
        return s;
    }

    inline double multiset_pair_dl(int64_t na, int64_t nb, int64_t e) const {
        if (e <= 0 || na <= 0 || nb <= 0) return 0.0;
        return bp_lbinom_d(static_cast<double>(na * nb + e - 1), static_cast<double>(e));
    }

    double entropy() const {
        if (multiset_adj) {
            double S = part_dl_type(0, B_type[0]) + part_dl_type(1, B_type[1]);
            S += edges_dl(B_type[0], B_type[1]);
            for (const auto& kv : mrs) {
                if (kv.second.e <= 0) continue;
                int a = static_cast<int>(kv.first >> 32);
                int b_id = static_cast<int>(kv.first & 0xffffffffULL);
                S += multiset_pair_dl(wr[a], wr[b_id], kv.second.e);
            }
            for (int c = 0; c < 2; ++c) {
                for (int r : active_by_type[c]) {
                    S -= bp_math().lgam1(wr[r]);
                }
            }
            return S;
        }

        double S = g->edge_factor - (dc ? g->degree_factor : 0.0);
        for (int m = 0; m < g->M; ++m) S += g->cov_specs[m].const_term;

        // Pairwise edge counts & covariate likelihoods
        for (const auto& kv : mrs) {
            if (kv.second.e > 0) S += pair_entropy(kv.second);
        }

        // Per-block likelihood, partition, and degree terms
        S += part_dl_type(0, B_type[0]) + part_dl_type(1, B_type[1]);
        S += edges_dl(B_type[0], B_type[1]);

        for (int c = 0; c < 2; ++c) {
            for (int r : active_by_type[c]) {
                S += block_lik_term(wr[r], mr_out[r], mr_in[r]);
                S -= bp_math().lgam1(wr[r]); // -log(n_r!) from partition prior
                S += block_deg_dl(r);
            }
        }
        return S;
    }

    // Exact O(deg(v)) incremental entropy change for moving vertex v from r=b[v]
    // to target block s (s == -1 means moving to a new empty block of type vtype[v]).
    double move_delta(int v, int s) const {
        int r = b[v];
        int t = g->vtype[v];
        if (s == r) return 0.0;
        bool is_new = (s < 0 || wr[s] == 0);
        if (!is_new && btype[s] != t) return BP_INF;
        if (is_new && wr[r] == 1) return 0.0; // moving singleton to empty block is a no-op

        int64_t k_out = g->deg_out[v];
        int64_t k_in = g->directed ? g->deg_in[v] : 0;

        int64_t nr = wr[r], nr_new = nr - 1;
        int64_t ns = is_new ? 0 : wr[s], ns_new = ns + 1;

        if (multiset_adj) {
            double dS = (-bp_math().lgam1(nr_new) + bp_math().lgam1(nr))
                      + (-bp_math().lgam1(ns_new) + bp_math().lgam1(ns));
            int dB = (is_new ? 1 : 0) - (nr_new == 0 ? 1 : 0);
            if (dB != 0) {
                int B0 = B_type[0], B1 = B_type[1];
                int B0_new = B0 + (t == 0 ? dB : 0);
                int B1_new = B1 + (t == 1 ? dB : 0);
                dS += part_dl_type(t, B_type[t] + dB) - part_dl_type(t, B_type[t]);
                dS += edges_dl(B0_new, B1_new) - edges_dl(B0, B1);
            }
            // Collect delta edges from v per opposite-type block q (and direction)
            auto multiset_dir = [&](const std::vector<BpEdgeEntry>& adj, bool v_is_src) {
                bp_hash_map<int, int64_t> dq;
                for (const auto& e : adj) dq[b[e.nbr]] += e.mult;
                bp_hash_map<int, bool> touched;
                for (const auto& qkv : blk_nbrs[r]) touched[qkv.first] = true;
                if (!is_new) for (const auto& qkv : blk_nbrs[s]) touched[qkv.first] = true;
                for (const auto& qkv : dq) touched[qkv.first] = true;

                for (const auto& tkv : touched) {
                    int q = tkv.first;
                    int64_t wq = wr[q];
                    int64_t dv = 0;
                    auto it_d = dq.find(q);
                    if (it_d != dq.end()) dv = it_d->second;

                    uint64_t krq = (!g->directed)
                        ? (t == 0 ? bp_pair_key(r, q) : bp_pair_key(q, r))
                        : (v_is_src ? bp_pair_key(r, q) : bp_pair_key(q, r));
                    auto it_rq = mrs.find(krq);
                    int64_t erq = (it_rq != mrs.end()) ? it_rq->second.e : 0;
                    if (erq > 0 || dv > 0) {
                        dS -= multiset_pair_dl(nr, wq, erq);
                        dS += multiset_pair_dl(nr_new, wq, erq - dv);
                    }

                    int64_t esq = 0;
                    if (!is_new) {
                        uint64_t ksq = (!g->directed)
                            ? (t == 0 ? bp_pair_key(s, q) : bp_pair_key(q, s))
                            : (v_is_src ? bp_pair_key(s, q) : bp_pair_key(q, s));
                        auto it_sq = mrs.find(ksq);
                        if (it_sq != mrs.end()) esq = it_sq->second.e;
                    }
                    if (esq > 0 || dv > 0) {
                        dS -= multiset_pair_dl(ns, wq, esq);
                        dS += multiset_pair_dl(ns_new, wq, esq + dv);
                    }
                }
            };
            multiset_dir(g->out_adj[v], true);
            if (g->directed) multiset_dir(g->in_adj[v], false);
            return dS;
        }

        int64_t er_out = mr_out[r], er_out_new = er_out - k_out;
        int64_t er_in = mr_in[r], er_in_new = er_in - k_in;
        int64_t es_out = is_new ? 0 : mr_out[s], es_out_new = es_out + k_out;
        int64_t es_in = is_new ? 0 : mr_in[s], es_in_new = es_in + k_in;

        double dS = 0.0;

        // 1. Block-level likelihood + partition (-lgam1(nr))
        dS += block_lik_term(nr_new, er_out_new, er_in_new) - block_lik_term(nr, er_out, er_in);
        dS += block_lik_term(ns_new, es_out_new, es_in_new) - block_lik_term(ns, es_out, es_in);
        dS += (-bp_math().lgam1(nr_new) + bp_math().lgam1(nr));
        dS += (-bp_math().lgam1(ns_new) + bp_math().lgam1(ns));

        // 2. Partition prior & edges_dl change if B_type[t] changes
        int dB = (is_new ? 1 : 0) - (nr_new == 0 ? 1 : 0);
        if (dB != 0) {
            int B0 = B_type[0], B1 = B_type[1];
            int B0_new = B0 + (t == 0 ? dB : 0);
            int B1_new = B1 + (t == 1 ? dB : 0);
            dS += part_dl_type(t, B_type[t] + dB) - part_dl_type(t, B_type[t]);
            dS += edges_dl(B0_new, B1_new) - edges_dl(B0, B1);
        }

        // 3. Degree prior change (DC only)
        if (dc) {
            if (dist) {
                auto one_hist_delta = [&](const std::vector<bp_hash_map<int64_t, int64_t>>& hist,
                                          int64_t kv_val, int64_t er_old, int64_t er_new_v,
                                          int64_t es_old, int64_t es_new_v) -> double {
                    double d = 0.0;
                    d -= bp_log_q()(er_old, nr);
                    if (nr_new > 0) d += bp_log_q()(er_new_v, nr_new);
                    d += bp_math().lgam1(nr_new) - bp_math().lgam1(nr);
                    auto it_r = hist[r].find(kv_val);
                    int64_t c_r = (it_r != hist[r].end()) ? it_r->second : 0;
                    d += bp_math().lgam1(c_r) - bp_math().lgam1(c_r - 1);

                    if (!is_new) d -= bp_log_q()(es_old, ns);
                    d += bp_log_q()(es_new_v, ns_new);
                    d += bp_math().lgam1(ns_new) - bp_math().lgam1(ns);
                    int64_t c_s = 0;
                    if (!is_new) {
                        auto it_s = hist[s].find(kv_val);
                        if (it_s != hist[s].end()) c_s = it_s->second;
                    }
                    d += bp_math().lgam1(c_s) - bp_math().lgam1(c_s + 1);
                    return d;
                };
                dS += one_hist_delta(deg_hist_out, k_out, er_out, er_out_new, es_out, es_out_new);
                if (g->directed) {
                    dS += one_hist_delta(deg_hist_in, k_in, er_in, er_in_new, es_in, es_in_new);
                }
            } else {
                if (g->directed) {
                    dS -= bp_math().lbinom(nr + er_out - 1, er_out) + bp_math().lbinom(nr + er_in - 1, er_in);
                    if (nr_new > 0) dS += bp_math().lbinom(nr_new + er_out_new - 1, er_out_new) +
                                          bp_math().lbinom(nr_new + er_in_new - 1, er_in_new);
                    if (!is_new) dS -= bp_math().lbinom(ns + es_out - 1, es_out) + bp_math().lbinom(ns + es_in - 1, es_in);
                    dS += bp_math().lbinom(ns_new + es_out_new - 1, es_out_new) +
                          bp_math().lbinom(ns_new + es_in_new - 1, es_in_new);
                } else {
                    dS -= bp_math().lbinom(nr + er_out - 1, er_out);
                    if (nr_new > 0) dS += bp_math().lbinom(nr_new + er_out_new - 1, er_out_new);
                    if (!is_new) dS -= bp_math().lbinom(ns + es_out - 1, es_out);
                    dS += bp_math().lbinom(ns_new + es_out_new - 1, es_out_new);
                }
            }
        }

        // 4. Pairwise edge-count & covariate delta (neighbors of v all have type 1 - t!)
        // Accumulate deltas per neighbor block q (and direction if directed)
        auto accum_dir = [&](const std::vector<BpEdgeEntry>& adj, bool v_is_src) {
            bp_hash_map<int, PairStat> delta_q;
            for (const auto& e : adj) {
                int q = b[e.nbr];
                auto it = delta_q.find(q);
                if (it == delta_q.end()) it = delta_q.emplace(q, PairStat(g->M)).first;
                it->second.e += e.mult;
                for (int m = 0; m < g->M; ++m) {
                    it->second.sw[m] += g->edge_sum_w[e.w_offset + m];
                    it->second.sw2[m] += g->edge_sum_w2[e.w_offset + m];
                }
            }
            PairStat tmp(g->M);
            for (const auto& kv : delta_q) {
                int q = kv.first;
                const PairStat& dp = kv.second;
                // Key for (r, q)
                uint64_t krq = (!g->directed)
                    ? (t == 0 ? bp_pair_key(r, q) : bp_pair_key(q, r))
                    : (v_is_src ? bp_pair_key(r, q) : bp_pair_key(q, r));
                auto it_rq = mrs.find(krq);
                if (it_rq != mrs.end()) {
                    const PairStat& old_rq = it_rq->second;
                    dS -= pair_entropy(old_rq);
                    tmp.e = old_rq.e - dp.e;
                    for (int m = 0; m < g->M; ++m) {
                        tmp.sw[m] = old_rq.sw[m] - dp.sw[m];
                        tmp.sw2[m] = old_rq.sw2[m] - dp.sw2[m];
                    }
                    dS += pair_entropy(tmp);
                }
                // Key for (s, q)
                if (is_new) {
                    dS += pair_entropy(dp);
                } else {
                    uint64_t ksq = (!g->directed)
                        ? (t == 0 ? bp_pair_key(s, q) : bp_pair_key(q, s))
                        : (v_is_src ? bp_pair_key(s, q) : bp_pair_key(q, s));
                    auto it_sq = mrs.find(ksq);
                    if (it_sq != mrs.end()) {
                        const PairStat& old_sq = it_sq->second;
                        dS -= pair_entropy(old_sq);
                        tmp.e = old_sq.e + dp.e;
                        for (int m = 0; m < g->M; ++m) {
                            tmp.sw[m] = old_sq.sw[m] + dp.sw[m];
                            tmp.sw2[m] = old_sq.sw2[m] + dp.sw2[m];
                        }
                        dS += pair_entropy(tmp);
                    } else {
                        dS += pair_entropy(dp);
                    }
                }
            }
        };

        accum_dir(g->out_adj[v], true);
        if (g->directed) accum_dir(g->in_adj[v], false);

        return dS;
    }

    // Apply move of vertex v to block s (s == -1 allocates a free block).
    // Returns the actual block id s assigned.
    int move_vertex(int v, int s) {
        int r = b[v];
        int t = g->vtype[v];
        if (s == r) return r;
        if (s < 0) {
            if (wr[r] == 1) return r;
            s = pop_free_block();
            activate_block(s, t);
        } else if (wr[s] == 0) {
            // s was free: remove from free_blocks and activate
            auto fit = std::find(free_blocks.begin(), free_blocks.end(), s);
            if (fit != free_blocks.end()) free_blocks.erase(fit);
            activate_block(s, t);
        }

        int64_t k_out = g->deg_out[v];
        int64_t k_in = g->directed ? g->deg_in[v] : 0;

        auto apply_dir = [&](const std::vector<BpEdgeEntry>& adj, bool v_is_src) {
            for (const auto& e : adj) {
                int q = b[e.nbr];
                uint64_t krq = (!g->directed)
                    ? (t == 0 ? bp_pair_key(r, q) : bp_pair_key(q, r))
                    : (v_is_src ? bp_pair_key(r, q) : bp_pair_key(q, r));
                uint64_t ksq = (!g->directed)
                    ? (t == 0 ? bp_pair_key(s, q) : bp_pair_key(q, s))
                    : (v_is_src ? bp_pair_key(s, q) : bp_pair_key(q, s));

                auto it_rq = mrs.find(krq);
                if (it_rq != mrs.end()) {
                    it_rq->second.e -= e.mult;
                    for (int m = 0; m < g->M; ++m) {
                        it_rq->second.sw[m] -= g->edge_sum_w[e.w_offset + m];
                        it_rq->second.sw2[m] -= g->edge_sum_w2[e.w_offset + m];
                    }
                    if (it_rq->second.e <= 0) mrs.erase(it_rq);
                }
                auto it_sq = mrs.find(ksq);
                if (it_sq == mrs.end()) it_sq = mrs.emplace(ksq, PairStat(g->M)).first;
                it_sq->second.e += e.mult;
                for (int m = 0; m < g->M; ++m) {
                    it_sq->second.sw[m] += g->edge_sum_w[e.w_offset + m];
                    it_sq->second.sw2[m] += g->edge_sum_w2[e.w_offset + m];
                }

                // Update blk_nbrs
                auto bn_rq = blk_nbrs[r].find(q);
                if (bn_rq != blk_nbrs[r].end()) {
                    bn_rq->second -= e.mult;
                    if (bn_rq->second <= 0) blk_nbrs[r].erase(bn_rq);
                }
                auto bn_qr = blk_nbrs[q].find(r);
                if (bn_qr != blk_nbrs[q].end()) {
                    bn_qr->second -= e.mult;
                    if (bn_qr->second <= 0) blk_nbrs[q].erase(bn_qr);
                }
                blk_nbrs[s][q] += e.mult;
                blk_nbrs[q][s] += e.mult;
            }
        };

        apply_dir(g->out_adj[v], true);
        if (g->directed) apply_dir(g->in_adj[v], false);

        if (dc && dist) {
            auto it_r = deg_hist_out[r].find(k_out);
            if (it_r != deg_hist_out[r].end()) {
                it_r->second--;
                if (it_r->second <= 0) deg_hist_out[r].erase(it_r);
            }
            deg_hist_out[s][k_out]++;
            if (g->directed) {
                auto it_ri = deg_hist_in[r].find(k_in);
                if (it_ri != deg_hist_in[r].end()) {
                    it_ri->second--;
                    if (it_ri->second <= 0) deg_hist_in[r].erase(it_ri);
                }
                deg_hist_in[s][k_in]++;
            }
        }

        b[v] = s;
        wr[r]--;
        wr[s]++;
        mr_out[r] -= k_out;
        mr_out[s] += k_out;
        if (g->directed) {
            mr_in[r] -= k_in;
            mr_in[s] += k_in;
        }
        if (wr[r] == 0) deactivate_block(r);
        return s;
    }

    // Exact merge delta for merging entire block r into block s (same type t)
    double merge_delta(int r, int s) const {
        if (r == s || wr[r] == 0 || wr[s] == 0 || btype[r] != btype[s]) return BP_INF;
        int t = btype[r];
        int64_t nr = wr[r], ns = wr[s], n_new = nr + ns;
        if (multiset_adj) {
            double dS = -bp_math().lgam1(n_new) + bp_math().lgam1(nr) + bp_math().lgam1(ns);
            int B0_new = B_type[0] - (t == 0 ? 1 : 0);
            int B1_new = B_type[1] - (t == 1 ? 1 : 0);
            dS += part_dl_type(t, B_type[t] - 1) - part_dl_type(t, B_type[t]);
            dS += edges_dl(B0_new, B1_new) - edges_dl(B_type[0], B_type[1]);

            bp_hash_map<int, bool> touched;
            for (const auto& qkv : blk_nbrs[r]) touched[qkv.first] = true;
            for (const auto& qkv : blk_nbrs[s]) touched[qkv.first] = true;

            auto do_multiset = [&](int q, uint64_t krq, uint64_t ksq) {
                int64_t wq = wr[q];
                auto it_rq = mrs.find(krq);
                auto it_sq = mrs.find(ksq);
                int64_t erq = (it_rq != mrs.end()) ? it_rq->second.e : 0;
                int64_t esq = (it_sq != mrs.end()) ? it_sq->second.e : 0;
                if (erq > 0 || esq > 0) {
                    dS -= multiset_pair_dl(nr, wq, erq) + multiset_pair_dl(ns, wq, esq);
                    dS += multiset_pair_dl(n_new, wq, erq + esq);
                }
            };
            for (const auto& tkv : touched) {
                int q = tkv.first;
                if (!g->directed) {
                    do_multiset(q, t == 0 ? bp_pair_key(r, q) : bp_pair_key(q, r),
                                   t == 0 ? bp_pair_key(s, q) : bp_pair_key(q, s));
                } else {
                    do_multiset(q, bp_pair_key(r, q), bp_pair_key(s, q));
                    do_multiset(q, bp_pair_key(q, r), bp_pair_key(q, s));
                }
            }
            return dS;
        }

        int64_t er_out = mr_out[r], es_out = mr_out[s], e_out_new = er_out + es_out;
        int64_t er_in = mr_in[r], es_in = mr_in[s], e_in_new = er_in + es_in;

        double dS = 0.0;
        dS += block_lik_term(n_new, e_out_new, e_in_new)
              - block_lik_term(nr, er_out, er_in)
              - block_lik_term(ns, es_out, es_in);
        dS += -bp_math().lgam1(n_new) + bp_math().lgam1(nr) + bp_math().lgam1(ns);

        int B0_new = B_type[0] - (t == 0 ? 1 : 0);
        int B1_new = B_type[1] - (t == 1 ? 1 : 0);
        dS += part_dl_type(t, B_type[t] - 1) - part_dl_type(t, B_type[t]);
        dS += edges_dl(B0_new, B1_new) - edges_dl(B_type[0], B_type[1]);

        if (dc) {
            if (dist) {
                dS -= block_deg_dl(r) + block_deg_dl(s);
                auto merge_hist_dl = [&](const std::vector<bp_hash_map<int64_t, int64_t>>& hist,
                                         int64_t e_new_v) -> double {
                    double sn = bp_log_q()(e_new_v, n_new) + bp_math().lgam1(n_new);
                    const auto& h_small = (hist[r].size() <= hist[s].size()) ? hist[r] : hist[s];
                    const auto& h_large = (hist[r].size() <= hist[s].size()) ? hist[s] : hist[r];
                    for (const auto& kv : h_large) {
                        auto it = h_small.find(kv.first);
                        int64_t cnt = kv.second + (it != h_small.end() ? it->second : 0);
                        if (cnt > 1) sn -= bp_math().lgam1(cnt);
                    }
                    for (const auto& kv : h_small) {
                        if (h_large.find(kv.first) == h_large.end() && kv.second > 1) {
                            sn -= bp_math().lgam1(kv.second);
                        }
                    }
                    return sn;
                };
                dS += merge_hist_dl(deg_hist_out, e_out_new);
                if (g->directed) dS += merge_hist_dl(deg_hist_in, e_in_new);
            } else {
                dS -= block_deg_dl(r) + block_deg_dl(s);
                if (g->directed) {
                    dS += bp_math().lbinom(n_new + e_out_new - 1, e_out_new) +
                          bp_math().lbinom(n_new + e_in_new - 1, e_in_new);
                } else {
                    dS += bp_math().lbinom(n_new + e_out_new - 1, e_out_new);
                }
            }
        }

        // Pairwise terms: only neighbor blocks q of r are affected!
        PairStat tmp(g->M);
        for (const auto& qkv : blk_nbrs[r]) {
            int q = qkv.first;
            auto do_pair = [&](uint64_t krq, uint64_t ksq) {
                auto it_rq = mrs.find(krq);
                if (it_rq == mrs.end()) return;
                dS -= pair_entropy(it_rq->second);
                auto it_sq = mrs.find(ksq);
                if (it_sq != mrs.end()) {
                    dS -= pair_entropy(it_sq->second);
                    tmp.e = it_rq->second.e + it_sq->second.e;
                    for (int m = 0; m < g->M; ++m) {
                        tmp.sw[m] = it_rq->second.sw[m] + it_sq->second.sw[m];
                        tmp.sw2[m] = it_rq->second.sw2[m] + it_sq->second.sw2[m];
                    }
                    dS += pair_entropy(tmp);
                } else {
                    dS += pair_entropy(it_rq->second);
                }
            };
            if (!g->directed) {
                uint64_t krq = (t == 0) ? bp_pair_key(r, q) : bp_pair_key(q, r);
                uint64_t ksq = (t == 0) ? bp_pair_key(s, q) : bp_pair_key(q, s);
                do_pair(krq, ksq);
            } else {
                do_pair(bp_pair_key(r, q), bp_pair_key(s, q));
                do_pair(bp_pair_key(q, r), bp_pair_key(q, s));
            }
        }
        return dS;
    }

    // Merge all vertices of block r into block s (same type).
    void merge_blocks(int r, int s, const std::vector<std::vector<int>>& blk_verts) {
        for (int v : blk_verts[r]) {
            move_vertex(v, s);
        }
    }

    // Fast 2-hop candidate block proposal for vertex v (same type t = vtype[v]).
    // Returns candidate block s in active_by_type[t], or -1 for new empty block.
    int propose_block(int v, double c_prop, double d_empty, BpRng& rng,
                      int B_max_t) const {
        int t = g->vtype[v];
        const auto& act = active_by_type[t];
        int Bt = static_cast<int>(act.size());
        if (Bt == 0) return -1;

        if (d_empty > 0.0 && Bt < B_max_t && rng.unif() < d_empty) {
            return -1;
        }

        // 2-hop proposal via a random neighbor u -> neighbor's block q -> neighbor block s of q
        int64_t deg_v = g->directed ? (g->deg_out[v] + g->deg_in[v]) : g->deg_out[v];
        if (deg_v > 0 && (c_prop == 0.0 || rng.unif() > c_prop / (c_prop + static_cast<double>(deg_v)))) {
            // Pick a random incident edge of v
            int64_t pick = static_cast<int64_t>(rng.index(static_cast<size_t>(deg_v)));
            int u = -1;
            int64_t acc = 0;
            for (const auto& e : g->out_adj[v]) {
                acc += e.mult;
                if (pick < acc) { u = e.nbr; break; }
            }
            if (u < 0 && g->directed) {
                for (const auto& e : g->in_adj[v]) {
                    acc += e.mult;
                    if (pick < acc) { u = e.nbr; break; }
                }
            }
            if (u >= 0) {
                int q = b[u];
                int64_t eq = g->directed ? (mr_out[q] + mr_in[q]) : mr_out[q];
                if (eq > 0 && rng.unif() > c_prop / (c_prop + static_cast<double>(eq))) {
                    int64_t epick = static_cast<int64_t>(rng.index(static_cast<size_t>(eq)));
                    int64_t eacc = 0;
                    for (const auto& kv : blk_nbrs[q]) {
                        eacc += kv.second;
                        if (epick < eacc) return kv.first;
                    }
                }
            }
        }
        return act[rng.index(static_cast<size_t>(Bt))];
    }

    // Exact log proposal probability log P(s | v, state) for uniform-in-type + empty
    // proposal (used for exact Metropolis-Hastings detailed balance in finite-T MCMC).
    double uniform_proposal_lprob(int v, int s, double d_empty,
                                  int B_min_t, int B_max_t) const {
        int t = g->vtype[v];
        int Bt = B_type[t];
        bool is_new = (s < 0 || wr[s] == 0);
        bool can_new = (d_empty > 0.0 && Bt < B_max_t && wr[b[v]] > 1);
        int n_exist = Bt - 1; // excluding current block b[v]
        if (is_new) {
            if (!can_new) return -BP_INF;
            return std::log(d_empty);
        }
        if (n_exist <= 0) return -BP_INF;
        double p_exist = can_new ? (1.0 - d_empty) : 1.0;
        return std::log(p_exist) - std::log(static_cast<double>(n_exist));
    }

    // Relabel active blocks consecutively: type 0 -> 1..B0, type 1 -> (B0+1)..(B0+B1)
    Rcpp::IntegerVector canonical_blocks() const {
        std::vector<int> map_id(B_cap, 0);
        int next_id = 1;
        for (int c = 0; c < 2; ++c) {
            for (int r : active_by_type[c]) {
                map_id[r] = next_id++;
            }
        }
        Rcpp::IntegerVector out(g->N);
        for (int v = 0; v < g->N; ++v) out[v] = map_id[b[v]];
        return out;
    }

    // Greedy zero-temperature single-node refinement sweeps
    int greedy_sweep(int niter, double c_prop, BpRng& rng,
                     int B0_min, int B0_max, int B1_min, int B1_max) {
        int total_moved = 0;
        std::vector<int> verts(g->N);
        std::iota(verts.begin(), verts.end(), 0);
        for (int it = 0; it < niter; ++it) {
            rng.shuffle(verts);
            int moved = 0;
            for (int v : verts) {
                int r = b[v];
                int t = g->vtype[v];
                int Bmin_t = (t == 0) ? B0_min : B1_min;
                int Bmax_t = (t == 0) ? B0_max : B1_max;
                if (wr[r] == 1 && B_type[t] <= Bmin_t) continue;
                int s = propose_block(v, c_prop, 0.01, rng, Bmax_t);
                if (s == r) continue;
                if (s < 0 && B_type[t] >= Bmax_t) continue;
                double dS = move_delta(v, s);
                if (dS < -1e-12) {
                    move_vertex(v, s);
                    moved++;
                }
            }
            total_moved += moved;
            if (moved == 0) break;
        }
        return total_moved;
    }

    // Parallel agglomerative merge sweep within type c (reduces B_type[c] toward target_Bc)
    bool merge_sweep_type(int c, int target_Bc, int n_cand, BpRng& rng) {
        const auto& act = active_by_type[c];
        int Bc = static_cast<int>(act.size());
        if (Bc <= target_Bc || Bc <= 1) return false;

        int n_merges_goal = std::max(1, Bc - target_Bc);

        // Build per-block vertex lists once
        std::vector<std::vector<int>> blk_verts(B_cap);
        for (int v = 0; v < g->N; ++v) {
            if (g->vtype[v] == c) blk_verts[b[v]].push_back(v);
        }

        // Pre-generate candidate merge targets per block in act (deterministic via RRng)
        std::vector<std::vector<int>> cands(Bc);
        for (int i = 0; i < Bc; ++i) {
            int r = act[i];
            if (Bc <= n_cand + 1) {
                for (int j = 0; j < Bc; ++j) {
                    if (act[j] != r) cands[i].push_back(act[j]);
                }
            } else {
                // 2-hop neighbor blocks of r + random blocks in act
                bp_hash_map<int, bool> seen;
                seen[r] = true;
                for (const auto& qkv : blk_nbrs[r]) {
                    int q = qkv.first;
                    for (const auto& skv : blk_nbrs[q]) {
                        int s = skv.first;
                        if (seen.find(s) == seen.end()) {
                            seen[s] = true;
                            cands[i].push_back(s);
                            if (static_cast<int>(cands[i].size()) >= n_cand) break;
                        }
                    }
                    if (static_cast<int>(cands[i].size()) >= n_cand) break;
                }
                int tries = 0;
                while (static_cast<int>(cands[i].size()) < n_cand && tries < n_cand * 3) {
                    int s = act[rng.index(static_cast<size_t>(Bc))];
                    if (seen.find(s) == seen.end()) {
                        seen[s] = true;
                        cands[i].push_back(s);
                    }
                    tries++;
                }
            }
        }

        // Parallel read-only evaluation of best merge partner for each block r in act
        struct MergeProposal {
            double dS = BP_INF;
            int r = -1;
            int s = -1;
        };
        std::vector<MergeProposal> best(Bc);
        rgt_parallel::parallel_for(static_cast<size_t>(Bc), rgt_parallel::get_num_threads(),
                                   [&](int /*tid*/, size_t lo, size_t hi) {
            for (size_t i = lo; i < hi; ++i) {
                int r = act[i];
                double best_dS = BP_INF;
                int best_s = -1;
                for (int s : cands[i]) {
                    double dS = merge_delta(r, s);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best_s = s;
                    }
                }
                best[i] = {best_dS, r, best_s};
            }
        });

        // Sort proposals by dS ascending and commit greedily (re-checking validity)
        std::sort(best.begin(), best.end(), [](const MergeProposal& a, const MergeProposal& b) {
            return a.dS < b.dS;
        });

        int merged = 0;
        std::vector<int> redirect(B_cap, -1);
        auto resolve = [&](int x) {
            while (redirect[x] != -1) x = redirect[x];
            return x;
        };

        for (const auto& prop : best) {
            if (merged >= n_merges_goal || B_type[c] <= target_Bc) break;
            if (prop.r < 0 || prop.s < 0 || !std::isfinite(prop.dS)) continue;
            int r = resolve(prop.r);
            int s = resolve(prop.s);
            if (r == s || wr[r] == 0 || wr[s] == 0) continue;
            // Merge smaller into larger to minimize vertex moves
            if (blk_verts[r].size() > blk_verts[s].size()) std::swap(r, s);
            merge_blocks(r, s, blk_verts);
            blk_verts[s].insert(blk_verts[s].end(), blk_verts[r].begin(), blk_verts[r].end());
            blk_verts[r].clear();
            redirect[r] = s;
            merged++;
        }
        return merged > 0;
    }

    // Finite-temperature Metropolis-Hastings sweep (single-node + optional merge-split)
    void mcmc_sweep(double beta, double d_empty, bool multiflip,
                    int B0_min, int B0_max, int B1_min, int B1_max,
                    BpRng& rng, int& n_attempts, int& n_accepted) {
        std::vector<int> verts(g->N);
        std::iota(verts.begin(), verts.end(), 0);
        rng.shuffle(verts);

        for (int v : verts) {
            int r = b[v];
            int t = g->vtype[v];
            int Bmin_t = (t == 0) ? B0_min : B1_min;
            int Bmax_t = (t == 0) ? B0_max : B1_max;
            int Bt = B_type[t];

            bool can_new = (d_empty > 0.0 && Bt < Bmax_t && wr[r] > 1);
            int n_exist = Bt - 1;
            if (!can_new && n_exist <= 0) continue;

            int s;
            if (can_new && rng.unif() < d_empty) {
                s = -1;
            } else {
                if (n_exist <= 0) continue;
                int idx = static_cast<int>(rng.index(static_cast<size_t>(n_exist)));
                const auto& act = active_by_type[t];
                s = act[idx];
                if (s == r) s = act[n_exist];
            }

            if (wr[r] == 1 && Bt <= Bmin_t) continue;

            n_attempts++;
            double dS = move_delta(v, s);
            double lp_fwd = uniform_proposal_lprob(v, s, d_empty, Bmin_t, Bmax_t);

            // Reverse proposal probability after hypothetical move
            int Bt_after = Bt + (s < 0 ? 1 : 0) - (wr[r] == 1 ? 1 : 0);
            int wr_s_after = (s < 0 ? 1 : wr[s] + 1);
            bool rev_is_new = (wr[r] == 1);
            bool rev_can_new = (d_empty > 0.0 && Bt_after < Bmax_t && wr_s_after > 1);
            int rev_n_exist = Bt_after - 1;
            double lp_rev;
            if (rev_is_new) {
                lp_rev = rev_can_new ? std::log(d_empty) : -BP_INF;
            } else {
                double p_exist = rev_can_new ? (1.0 - d_empty) : 1.0;
                lp_rev = (rev_n_exist > 0) ? (std::log(p_exist) - std::log(static_cast<double>(rev_n_exist))) : -BP_INF;
            }

            double log_acc;
            if (std::isinf(beta)) {
                log_acc = (dS < -1e-12) ? 1.0 : -1.0;
            } else {
                log_acc = -beta * dS + (lp_rev - lp_fwd);
            }

            if (log_acc >= 0.0 || (std::isfinite(log_acc) && std::log(rng.unif()) < log_acc)) {
                move_vertex(v, s);
                n_accepted++;
            }
        }

        // Optional within-type merge / split proposal for faster mixing
        if (multiflip) {
            for (int c = 0; c < 2; ++c) {
                int Bc = B_type[c];
                int Bmin_c = (c == 0) ? B0_min : B1_min;
                int Bmax_c = (c == 0) ? B0_max : B1_max;
                if (Bc >= 2 && Bc > Bmin_c && rng.unif() < 0.5) {
                    // Propose random merge of two blocks in type c
                    const auto& act = active_by_type[c];
                    int i = static_cast<int>(rng.index(static_cast<size_t>(Bc)));
                    int j = static_cast<int>(rng.index(static_cast<size_t>(Bc - 1)));
                    if (j >= i) j++;
                    int r = act[i], s = act[j];
                    n_attempts++;
                    double dS = merge_delta(r, s);
                    // Balance with random split reverse probability:
                    // Merge picks unordered pair {r,s} with prob 1 / C(Bc, 2).
                    // Reverse split picks merged block (1/(Bc-1)) and random non-empty
                    // bipartition of its (nr+ns) vertices (1 / (2^(nr+ns-1) - 1)).
                    int64_t n_tot = wr[r] + wr[s];
                    double lp_fwd = -bp_math().lbinom(Bc, 2);
                    double lp_rev = -bp_math().logi(Bc - 1) -
                                    std::log(std::ldexp(1.0, static_cast<int>(std::min<int64_t>(n_tot - 1, 50))) - 1.0);
                    if (n_tot - 1 > 50) {
                        lp_rev = -bp_math().logi(Bc - 1) - static_cast<double>(n_tot - 1) * std::log(2.0);
                    }
                    double log_acc = std::isinf(beta) ? ((dS < -1e-12) ? 1.0 : -1.0)
                                                      : (-beta * dS + (lp_rev - lp_fwd));
                    if (log_acc >= 0.0 || (std::isfinite(log_acc) && std::log(rng.unif()) < log_acc)) {
                        std::vector<std::vector<int>> bv(B_cap);
                        for (int v = 0; v < g->N; ++v) if (b[v] == r) bv[r].push_back(v);
                        merge_blocks(r, s, bv);
                        n_accepted++;
                    }
                } else if (Bc < Bmax_c) {
                    // Propose random split of a block with size >= 2 in type c
                    std::vector<int> splittable;
                    for (int r : active_by_type[c]) if (wr[r] >= 2) splittable.push_back(r);
                    if (splittable.empty()) continue;
                    int r = splittable[rng.index(splittable.size())];
                    std::vector<int> verts_r;
                    for (int v = 0; v < g->N; ++v) if (b[v] == r) verts_r.push_back(v);
                    // Random non-empty proper subset
                    std::vector<int> move_set;
                    do {
                        move_set.clear();
                        for (int v : verts_r) if (rng.unif() < 0.5) move_set.push_back(v);
                    } while (move_set.empty() || move_set.size() == verts_r.size());

                    n_attempts++;
                    double S_before = entropy();
                    int s_new = move_vertex(move_set[0], -1);
                    for (size_t k = 1; k < move_set.size(); ++k) move_vertex(move_set[k], s_new);
                    double S_after = entropy();
                    double dS = S_after - S_before;

                    int64_t n_tot = static_cast<int64_t>(verts_r.size());
                    double lp_fwd = -bp_math().logi(static_cast<int64_t>(splittable.size())) -
                                    ((n_tot - 1 <= 50)
                                     ? std::log(std::ldexp(1.0, static_cast<int>(n_tot - 1)) - 1.0)
                                     : static_cast<double>(n_tot - 1) * std::log(2.0));
                    double lp_rev = -bp_math().lbinom(Bc + 1, 2);
                    double log_acc = std::isinf(beta) ? ((dS < -1e-12) ? 1.0 : -1.0)
                                                      : (-beta * dS + (lp_rev - lp_fwd));
                    if (log_acc >= 0.0 || (std::isfinite(log_acc) && std::log(rng.unif()) < log_acc)) {
                        n_accepted++;
                    } else {
                        // Revert split
                        for (int v : move_set) move_vertex(v, r);
                    }
                }
            }
        }
    }
};

// Build upper-level bipartite block multigraph from a lower-level state
BpGraphData build_upper_bp_graph(const BipartiteBlockState& lower) {
    BpGraphData u;
    int B0 = lower.B_type[0];
    int B1 = lower.B_type[1];
    u.N = B0 + B1;
    u.N_type[0] = B0;
    u.N_type[1] = B1;
    u.directed = lower.g->directed;
    u.vtype.resize(u.N);
    std::vector<int> old_to_new(lower.B_cap, -1);
    int idx = 0;
    for (int r : lower.active_by_type[0]) {
        u.vtype[idx] = 0;
        old_to_new[r] = idx++;
    }
    for (int r : lower.active_by_type[1]) {
        u.vtype[idx] = 1;
        old_to_new[r] = idx++;
    }

    u.deg_out.assign(u.N, 0);
    u.deg_in.assign(u.N, 0);
    u.out_adj.resize(u.N);
    if (u.directed) u.in_adj.resize(u.N);

    for (const auto& kv : lower.mrs) {
        if (kv.second.e <= 0) continue;
        int r_old = static_cast<int>(kv.first >> 32);
        int s_old = static_cast<int>(kv.first & 0xffffffffULL);
        int r = old_to_new[r_old];
        int s = old_to_new[s_old];
        int64_t mult = kv.second.e;
        u.E += mult;
        u.edge_factor += bp_math().lgam1(mult);
        u.out_adj[r].push_back({s, mult, -1});
        if (u.directed) {
            u.in_adj[s].push_back({r, mult, -1});
            u.deg_out[r] += mult;
            u.deg_in[s] += mult;
        } else {
            u.out_adj[s].push_back({r, mult, -1});
            u.deg_out[r] += mult;
            u.deg_out[s] += mult;
        }
    }
    return u;
}

} // namespace

// ===========================================================================
// Exported Rcpp entry points
// ===========================================================================

// [[Rcpp::export]]
double cpp_bipartite_entropy(SEXP pointer,
                             Rcpp::IntegerVector type_vec,
                             Rcpp::IntegerVector blocks,
                             bool degree_corrected,
                             bool distributed,
                             bool bp_edge_prior,
                             Rcpp::NumericMatrix W,
                             Rcpp::CharacterVector cov_types,
                             Rcpp::List cov_params,
                             double q_max) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    BpGraphData gd = build_bp_graph(g, type_vec, W, cov_types, cov_params);
    std::vector<int> b0(gd.N);
    for (int i = 0; i < gd.N; ++i) b0[i] = blocks[i] - 1;
    BipartiteBlockState st(&gd, b0, degree_corrected, distributed, bp_edge_prior, true);
    return st.entropy();
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_bipartite_check_moves(SEXP pointer,
                                              Rcpp::IntegerVector type_vec,
                                              Rcpp::IntegerVector blocks,
                                              bool degree_corrected,
                                              bool distributed,
                                              bool bp_edge_prior,
                                              Rcpp::NumericMatrix W,
                                              Rcpp::CharacterVector cov_types,
                                              Rcpp::List cov_params,
                                              int nmoves,
                                              double q_max) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    BpGraphData gd = build_bp_graph(g, type_vec, W, cov_types, cov_params);
    std::vector<int> b0(gd.N);
    for (int i = 0; i < gd.N; ++i) b0[i] = blocks[i] - 1;
    BipartiteBlockState st(&gd, b0, degree_corrected, distributed, bp_edge_prior, true);

    Rcpp::RNGScope scope;
    BpRng rng;
    Rcpp::NumericVector diffs(nmoves);
    for (int m = 0; m < nmoves; ++m) {
        if (gd.N == 0) { diffs[m] = 0.0; continue; }
        int v = static_cast<int>(rng.index(static_cast<size_t>(gd.N)));
        int t = gd.vtype[v];
        int s;
        if (rng.unif() < 0.2) {
            s = -1;
        } else {
            const auto& act = st.active_by_type[t];
            s = act[rng.index(act.size())];
        }
        double S0 = st.entropy();
        double dS = st.move_delta(v, s);
        st.move_vertex(v, s);
        double S1 = st.entropy();
        diffs[m] = (S1 - S0) - dS;
    }
    return diffs;
}

// [[Rcpp::export]]
Rcpp::List cpp_bipartite_minimize(SEXP pointer,
                                  Rcpp::IntegerVector type_vec,
                                  bool degree_corrected,
                                  bool distributed,
                                  bool bp_edge_prior,
                                  Rcpp::NumericMatrix W,
                                  Rcpp::CharacterVector cov_types,
                                  Rcpp::List cov_params,
                                  int B0_min, int B0_max,
                                  int B1_min, int B1_max,
                                  double r_shrink,
                                  int n_cand,
                                  bool refine,
                                  double q_max,
                                  Rcpp::Nullable<Rcpp::IntegerVector> init_blocks = R_NilValue) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    BpGraphData gd = build_bp_graph(g, type_vec, W, cov_types, cov_params);

    Rcpp::RNGScope scope;
    BpRng rng;

    if (gd.N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("blocks") = Rcpp::IntegerVector(0),
            Rcpp::Named("description_length") = 0.0,
            Rcpp::Named("B0") = 0,
            Rcpp::Named("B1") = 0
        );
    }

    int N0 = gd.N_type[0], N1 = gd.N_type[1];
    B0_min = std::max(N0 > 0 ? 1 : 0, std::min(B0_min, N0));
    B0_max = std::max(B0_min, std::min(B0_max, N0));
    B1_min = std::max(N1 > 0 ? 1 : 0, std::min(B1_min, N1));
    B1_max = std::max(B1_min, std::min(B1_max, N1));

    // Initialize: either user-supplied init_blocks or singleton blocks
    std::vector<int> b0(gd.N);
    if (init_blocks.isNotNull()) {
        Rcpp::IntegerVector ib(init_blocks);
        for (int i = 0; i < gd.N; ++i) b0[i] = ib[i] - 1;
    } else {
        std::iota(b0.begin(), b0.end(), 0);
    }

    BipartiteBlockState st(&gd, b0, degree_corrected, distributed, bp_edge_prior, true);

    double best_S = BP_INF;
    Rcpp::IntegerVector best_b = st.canonical_blocks();
    int best_B0 = st.B_type[0], best_B1 = st.B_type[1];

    auto check_best = [&]() {
        if (st.B_type[0] >= B0_min && st.B_type[0] <= B0_max &&
            st.B_type[1] >= B1_min && st.B_type[1] <= B1_max) {
            double S = st.entropy();
            if (S < best_S) {
                best_S = S;
                best_b = st.canonical_blocks();
                best_B0 = st.B_type[0];
                best_B1 = st.B_type[1];
            }
        }
    };

    check_best();

    // Agglomerative multilevel reduction alternating/jointly shrinking B0 and B1
    while (st.B_type[0] > B0_min || st.B_type[1] > B1_min) {
        Rcpp::checkUserInterrupt();
        bool progressed = false;

        // Shrink whichever side has more remaining headroom (or both proportionally)
        for (int c = 0; c < 2; ++c) {
            int Bc = st.B_type[c];
            int Bmin_c = (c == 0) ? B0_min : B1_min;
            if (Bc <= Bmin_c) continue;
            int target = std::max(Bmin_c, static_cast<int>(std::floor(Bc / std::max(1.05, r_shrink))));
            if (target >= Bc) target = Bc - 1;
            if (st.merge_sweep_type(c, target, n_cand, rng)) {
                progressed = true;
            }
        }

        // Greedy single-node refinement sweeps at this scale
        st.greedy_sweep(5, 0.5, rng, B0_min, B0_max, B1_min, B1_max);
        check_best();

        if (!progressed) break;
    }

    // Final refinement from best found state
    if (refine) {
        std::vector<int> b_best(gd.N);
        for (int i = 0; i < gd.N; ++i) b_best[i] = best_b[i] - 1;
        BipartiteBlockState st_ref(&gd, b_best, degree_corrected, distributed, bp_edge_prior, true);
        st_ref.greedy_sweep(20, 0.5, rng, B0_min, B0_max, B1_min, B1_max);
        // Also try strictly improving merges on each side
        bool improved = true;
        while (improved) {
            improved = false;
            for (int c = 0; c < 2; ++c) {
                int Bmin_c = (c == 0) ? B0_min : B1_min;
                if (st_ref.B_type[c] <= Bmin_c) continue;
                const auto& act = st_ref.active_by_type[c];
                int Bc = static_cast<int>(act.size());
                double min_dS = -1e-12;
                int br = -1, bs = -1;
                if (Bc <= 64) {
                    for (int i = 0; i < Bc; ++i) {
                        for (int j = i + 1; j < Bc; ++j) {
                            double dS = st_ref.merge_delta(act[i], act[j]);
                            if (dS < min_dS) { min_dS = dS; br = act[i]; bs = act[j]; }
                        }
                    }
                }
                if (br >= 0 && bs >= 0) {
                    std::vector<std::vector<int>> bv(st_ref.B_cap);
                    for (int v = 0; v < gd.N; ++v) if (st_ref.b[v] == br) bv[br].push_back(v);
                    st_ref.merge_blocks(br, bs, bv);
                    st_ref.greedy_sweep(10, 0.5, rng, B0_min, B0_max, B1_min, B1_max);
                    improved = true;
                }
            }
        }
        double S_ref = st_ref.entropy();
        if (S_ref < best_S &&
            st_ref.B_type[0] >= B0_min && st_ref.B_type[0] <= B0_max &&
            st_ref.B_type[1] >= B1_min && st_ref.B_type[1] <= B1_max) {
            best_S = S_ref;
            best_b = st_ref.canonical_blocks();
            best_B0 = st_ref.B_type[0];
            best_B1 = st_ref.B_type[1];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("blocks") = best_b,
        Rcpp::Named("description_length") = best_S,
        Rcpp::Named("B0") = best_B0,
        Rcpp::Named("B1") = best_B1
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_bipartite_mcmc(SEXP pointer,
                              Rcpp::IntegerVector type_vec,
                              Rcpp::IntegerVector blocks,
                              bool degree_corrected,
                              bool distributed,
                              bool bp_edge_prior,
                              Rcpp::NumericMatrix W,
                              Rcpp::CharacterVector cov_types,
                              Rcpp::List cov_params,
                              Rcpp::List mcmc_opts,
                              Rcpp::List eq_opts,
                              double q_max) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    BpGraphData gd = build_bp_graph(g, type_vec, W, cov_types, cov_params);
    std::vector<int> b0(gd.N);
    for (int i = 0; i < gd.N; ++i) b0[i] = blocks[i] - 1;
    BipartiteBlockState st(&gd, b0, degree_corrected, distributed, bp_edge_prior, true);

    Rcpp::RNGScope scope;
    BpRng rng;

    double beta = Rcpp::as<double>(mcmc_opts["beta"]);
    double d_empty = Rcpp::as<double>(mcmc_opts["d"]);
    int niter = static_cast<int>(Rcpp::as<double>(mcmc_opts["niter"]));
    bool multiflip = Rcpp::as<bool>(mcmc_opts["multiflip"]);
    int B0_min = gd.N_type[0] > 0 ? 1 : 0;
    int B0_max = gd.N_type[0];
    int B1_min = gd.N_type[1] > 0 ? 1 : 0;
    int B1_max = gd.N_type[1];

    std::string mode = eq_opts.containsElementNamed("mode")
        ? Rcpp::as<std::string>(eq_opts["mode"]) : "sweep";
    bool do_eq = (mode == "equilibrate" || mode == "anneal");
    double wait_d = eq_opts.containsElementNamed("wait") ? Rcpp::as<double>(eq_opts["wait"]) : 0.0;
    int wait = std::isfinite(wait_d) ? static_cast<int>(wait_d) : 0;
    double nbreaks_d = eq_opts.containsElementNamed("nbreaks") ? Rcpp::as<double>(eq_opts["nbreaks"]) : 1.0;
    int nbreaks = std::isfinite(nbreaks_d) ? static_cast<int>(nbreaks_d) : 1;
    double force_d = eq_opts.containsElementNamed("force_niter") ? Rcpp::as<double>(eq_opts["force_niter"]) : -1.0;
    double max_d = eq_opts.containsElementNamed("max_niter") ? Rcpp::as<double>(eq_opts["max_niter"]) : static_cast<double>(niter);
    bool collect = eq_opts.containsElementNamed("collect") && Rcpp::as<bool>(eq_opts["collect"]);
    double thin_d = eq_opts.containsElementNamed("thin") ? Rcpp::as<double>(eq_opts["thin"]) : 1.0;
    int thin = std::max(1, std::isfinite(thin_d) ? static_cast<int>(thin_d) : 1);
    bool record_history = eq_opts.containsElementNamed("history") && Rcpp::as<bool>(eq_opts["history"]);

    int total_sweeps;
    bool fixed_niter = false;
    if (do_eq && force_d > 0.0 && std::isfinite(force_d)) {
        total_sweeps = static_cast<int>(force_d);
        fixed_niter = true;
    } else if (do_eq) {
        total_sweeps = std::isfinite(max_d) ? static_cast<int>(max_d) : 1000000;
    } else {
        total_sweeps = niter;
        fixed_niter = true;
    }

    int n_attempts = 0, n_accepted = 0;
    double S_init = st.entropy();
    double S_cur = S_init;
    double S_min = S_cur, S_max = S_cur;
    int break_count = 0, no_break = 0;

    std::vector<double> hist_entropy;
    std::vector<int> hist_iter;
    std::vector<double> sample_entropies;
    Rcpp::List samples;
    int actual_sweeps = 0;

    for (int sw = 1; sw <= total_sweeps; ++sw) {
        Rcpp::checkUserInterrupt();
        st.mcmc_sweep(beta, d_empty, multiflip, B0_min, B0_max, B1_min, B1_max,
                      rng, n_attempts, n_accepted);
        S_cur = st.entropy();
        actual_sweeps = sw;
        if (record_history) {
            hist_iter.push_back(sw);
            hist_entropy.push_back(S_cur);
        }
        if (collect && (sw % thin == 0)) {
            samples.push_back(st.canonical_blocks());
            sample_entropies.push_back(S_cur);
        }
        if (do_eq && !fixed_niter && wait > 0) {
            if (S_cur < S_min - 1e-9 || S_cur > S_max + 1e-9) {
                S_min = std::min(S_min, S_cur);
                S_max = std::max(S_max, S_cur);
                break_count++;
                no_break = 0;
            } else {
                no_break++;
                if (no_break >= wait && break_count >= nbreaks) break;
            }
        }
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("blocks") = st.canonical_blocks(),
        Rcpp::Named("entropy") = S_cur,
        Rcpp::Named("dS") = S_cur - S_init,
        Rcpp::Named("B0") = st.B_type[0],
        Rcpp::Named("B1") = st.B_type[1],
        Rcpp::Named("nattempts") = static_cast<double>(n_attempts),
        Rcpp::Named("nmoves") = static_cast<double>(n_accepted),
        Rcpp::Named("niter") = static_cast<double>(actual_sweeps)
    );
    if (record_history) {
        out["history"] = Rcpp::DataFrame::create(
            Rcpp::Named("iter") = hist_iter,
            Rcpp::Named("entropy") = hist_entropy
        );
    }
    if (collect) {
        out["samples"] = samples;
        out["sample_entropies"] = sample_entropies;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_bipartite_entropy(SEXP pointer,
                                                 Rcpp::IntegerVector type_vec,
                                                 Rcpp::List levels,
                                                 bool degree_corrected,
                                                 bool distributed,
                                                 bool bp_edge_prior,
                                                 Rcpp::NumericMatrix W,
                                                 Rcpp::CharacterVector cov_types,
                                                 Rcpp::List cov_params,
                                                 double q_max) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    int L = levels.size();
    Rcpp::NumericVector lev_S(L);

    BpGraphData cur_g = build_bp_graph(g, type_vec, W, cov_types, cov_params);
    for (int l = 0; l < L; ++l) {
        Rcpp::IntegerVector bl = levels[l];
        if (bl.size() != cur_g.N) {
            Rcpp::stop("Level %d block vector length (%d) does not match level vertex count (%d)",
                       l + 1, static_cast<int>(bl.size()), cur_g.N);
        }
        std::vector<int> b0(cur_g.N);
        for (int i = 0; i < cur_g.N; ++i) b0[i] = bl[i] - 1;
        bool is_top = (l == L - 1);
        BipartiteBlockState st(&cur_g, b0, (l == 0) ? degree_corrected : false,
                               distributed, bp_edge_prior, is_top, l > 0);
        lev_S[l] = st.entropy();
        if (!is_top) {
            cur_g = build_upper_bp_graph(st);
        }
    }
    return lev_S;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_bipartite_minimize(SEXP pointer,
                                         Rcpp::IntegerVector type_vec,
                                         bool degree_corrected,
                                         bool distributed,
                                         bool bp_edge_prior,
                                         Rcpp::NumericMatrix W,
                                         Rcpp::CharacterVector cov_types,
                                         Rcpp::List cov_params,
                                         double r_shrink,
                                         int n_cand,
                                         double q_max) {
    bp_log_q().init((!std::isfinite(q_max) || q_max > 10000) ? 10000 : static_cast<size_t>(std::max(100.0, q_max)));
    const Graph& g = get_graph(pointer);
    BpGraphData cur_g = build_bp_graph(g, type_vec, W, cov_types, cov_params);

    Rcpp::RNGScope scope;
    BpRng rng;

    Rcpp::List levels;
    std::vector<double> level_entropies;

    if (cur_g.N == 0) {
        return Rcpp::List::create(
            Rcpp::Named("levels") = levels,
            Rcpp::Named("level_entropies") = Rcpp::NumericVector(0),
            Rcpp::Named("description_length") = 0.0
        );
    }

    int max_levels = 20;
    for (int l = 0; l < max_levels; ++l) {
        Rcpp::checkUserInterrupt();
        int N0 = cur_g.N_type[0], N1 = cur_g.N_type[1];
        int B0_min = N0 > 0 ? 1 : 0;
        int B1_min = N1 > 0 ? 1 : 0;

        // Fit level l with flat bipartite minimization
        std::vector<int> b0(cur_g.N);
        std::iota(b0.begin(), b0.end(), 0);
        BipartiteBlockState st(&cur_g, b0, (l == 0) ? degree_corrected : false,
                               distributed, bp_edge_prior, true, l > 0);

        double best_S = st.entropy();
        Rcpp::IntegerVector best_b = st.canonical_blocks();
        int best_B0 = st.B_type[0], best_B1 = st.B_type[1];

        while (st.B_type[0] > B0_min || st.B_type[1] > B1_min) {
            bool progressed = false;
            for (int c = 0; c < 2; ++c) {
                int Bc = st.B_type[c];
                int Bmin_c = (c == 0) ? B0_min : B1_min;
                if (Bc <= Bmin_c) continue;
                int target = std::max(Bmin_c, static_cast<int>(std::floor(Bc / std::max(1.05, r_shrink))));
                if (target >= Bc) target = Bc - 1;
                if (st.merge_sweep_type(c, target, n_cand, rng)) progressed = true;
            }
            st.greedy_sweep(5, 0.5, rng, B0_min, N0, B1_min, N1);
            double S = st.entropy();
            if (S < best_S) {
                best_S = S;
                best_b = st.canonical_blocks();
                best_B0 = st.B_type[0];
                best_B1 = st.B_type[1];
            }
            if (!progressed) break;
        }

        // Reconstruct best state at level l and run greedy + merge refinement
        std::vector<int> b_best(cur_g.N);
        for (int i = 0; i < cur_g.N; ++i) b_best[i] = best_b[i] - 1;
        BipartiteBlockState st_best(&cur_g, b_best, (l == 0) ? degree_corrected : false,
                                    distributed, bp_edge_prior, true, l > 0);
        st_best.greedy_sweep(20, 0.5, rng, B0_min, N0, B1_min, N1);
        bool improved = true;
        while (improved) {
            improved = false;
            for (int c = 0; c < 2; ++c) {
                int Bmin_c = (c == 0) ? B0_min : B1_min;
                if (st_best.B_type[c] <= Bmin_c) continue;
                const auto& act = st_best.active_by_type[c];
                int Bc = static_cast<int>(act.size());
                double min_dS = -1e-12;
                int br = -1, bs = -1;
                if (Bc <= 64) {
                    for (int i = 0; i < Bc; ++i) {
                        for (int j = i + 1; j < Bc; ++j) {
                            double dS = st_best.merge_delta(act[i], act[j]);
                            if (dS < min_dS) { min_dS = dS; br = act[i]; bs = act[j]; }
                        }
                    }
                }
                if (br >= 0 && bs >= 0) {
                    std::vector<std::vector<int>> bv(st_best.B_cap);
                    for (int v = 0; v < cur_g.N; ++v) if (st_best.b[v] == br) bv[br].push_back(v);
                    st_best.merge_blocks(br, bs, bv);
                    st_best.greedy_sweep(10, 0.5, rng, B0_min, N0, B1_min, N1);
                    improved = true;
                }
            }
        }
        best_b = st_best.canonical_blocks();
        best_B0 = st_best.B_type[0];
        best_B1 = st_best.B_type[1];

        levels.push_back(best_b);

        // Stop if we reached (1, 1) or no further compression occurred
        if (best_B0 <= B0_min && best_B1 <= B1_min) {
            break;
        }
        if (l > 0 && best_B0 == N0 && best_B1 == N1) {
            // Force a final top level merging each type to 1 superblock
            Rcpp::IntegerVector top_b(best_B0 + best_B1);
            for (int i = 0; i < best_B0; ++i) top_b[i] = 1;
            for (int i = 0; i < best_B1; ++i) top_b[best_B0 + i] = (best_B0 > 0) ? 2 : 1;
            levels.push_back(top_b);
            break;
        }

        cur_g = build_upper_bp_graph(st_best);
    }

    // Compute exact nested level entropies (only top level includes edges_dl)
    Rcpp::NumericVector lev_S = cpp_nested_bipartite_entropy(
        pointer, type_vec, levels, degree_corrected, distributed, bp_edge_prior,
        W, cov_types, cov_params, q_max
    );
    double total_S = 0.0;
    for (R_xlen_t i = 0; i < lev_S.size(); ++i) total_S += lev_S[i];

    return Rcpp::List::create(
        Rcpp::Named("levels") = levels,
        Rcpp::Named("level_entropies") = lev_S,
        Rcpp::Named("description_length") = total_S
    );
}
