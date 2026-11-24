// graph-tool-R: native flat stochastic block model state and multilevel
// agglomerative description-length minimization.
//
// The entropy terms follow graph-tool 3.9's microcanonical BlockState with
// dense=False, multigraph=True, exact=True and the default partition, edge and
// degree priors. The search follows the structure of graph-tool's
// multilevel_mcmc_sweep (loops/multilevel.hh, support/bisection_search.hh) at
// beta = infinity with global moves, using sequential sweeps and merges.
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
// Copyright (C) 2026 Diogo Melo (R port)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <vector>
#include "graph_types.h"
#include "parallel.h"
#include <rgraphtool/vendor/hash_map_wrap.hh>

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double LOG2 = 0.693147180559945309417232121458176568;

// ---------------------------------------------------------------------------
// Random numbers come from R's generator so that set.seed() controls the search.

struct RRng {
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

// ---------------------------------------------------------------------------
// Fast integer special-function tables (exact IEEE-754 match to std::lgamma /
// std::log for cached integer arguments).

struct MathTables {
    static constexpr size_t MAX_TABLE = 2097152;
    std::vector<double> lgamma_tab;
    std::vector<double> log_tab;

    void ensure(size_t n) {
        size_t target = std::min(n + 1, MAX_TABLE);
        if (target <= lgamma_tab.size())
            return;
        size_t old = lgamma_tab.size();
        size_t cap = std::max<size_t>(4096, old);
        while (cap < target)
            cap = std::min(cap * 2, MAX_TABLE);
        lgamma_tab.resize(cap);
        log_tab.resize(cap);
        for (size_t i = old; i < cap; ++i) {
            lgamma_tab[i] = (i == 0) ? inf : std::lgamma(static_cast<double>(i));
            log_tab[i] = (i == 0) ? -inf : std::log(static_cast<double>(i));
        }
    }
};

MathTables& math_tables() {
    static MathTables t;
    return t;
}

inline double lgam_int(int64_t x) {
    const auto& tab = math_tables().lgamma_tab;
    if (x >= 1 && static_cast<size_t>(x) < tab.size())
        return tab[static_cast<size_t>(x)];
    return std::lgamma(static_cast<double>(x));
}

inline double log_int(int64_t x) {
    const auto& tab = math_tables().log_tab;
    if (x >= 1 && static_cast<size_t>(x) < tab.size())
        return tab[static_cast<size_t>(x)];
    return std::log(static_cast<double>(x));
}

inline double lbinom_int(int64_t n, int64_t k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return lgam_int(n + 1) - lgam_int(k + 1) - lgam_int(n - k + 1);
}

inline double lgam(double x) { return std::lgamma(x); }

inline double lbinom(double n, double k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return lgam(n + 1) - lgam(k + 1) - lgam(n - k + 1);
}

// Fast integer pair key mixer for flat hash maps.
inline uint64_t pair_key(uint64_t a, uint64_t b) {
    return (a << 32) | (b & 0xffffffffULL);
}

// log of the number of partitions of n into at most k parts. Exact values are
// tabulated (as in graph-tool's q cache) and an asymptotic approximation is used
// beyond the table (memoized in a bounded flat hash table).
class LogQ {
public:
    void ensure(size_t n_max) {
        if (n_max <= _n_max && !_table.empty())
            return;
        _n_max = n_max;
        _table.assign(n_max * (n_max + 1) / 2 + 1, -inf);
        for (size_t n = 1; n <= n_max; ++n) {
            for (size_t k = 1; k <= n; ++k) {
                double a = (k == 1) ? -inf : get(n, k - 1);
                double b = get_raw(n - k, k);
                double hi = std::max(a, b), lo = std::min(a, b);
                at(n, k) = (hi == -inf) ? -inf : hi + std::log1p(std::exp(lo - hi));
            }
        }
    }
    size_t n_max() const { return _n_max; }
    double operator()(long n, long k) const {
        if (k > n)
            k = n;
        if (n == 0 && k == 0)
            return 0;
        if (n < 0 || k <= 0)
            return -inf;
        if (static_cast<size_t>(n) <= _n_max)
            return get(n, k);
        thread_local gt_hash_map<uint64_t, double> approx_cache;
        uint64_t key = pair_key(static_cast<uint64_t>(n), static_cast<uint64_t>(k));
        auto it = approx_cache.find(key);
        if (it != approx_cache.end())
            return it->second;
        double val = approx(static_cast<double>(n), static_cast<double>(k));
        if (approx_cache.size() >= 65536)
            approx_cache.clear();
        approx_cache.emplace(key, val);
        return val;
    }
    // graph-tool's log_q_approx: the Szekeres saddle-point formula, with the
    // small-k limit for k < n^(1/4).
    static double approx(double n, double k) {
        if (k < std::pow(n, 0.25))
            return lbinom(n - 1, k - 1) - lgam(k + 1);
        double u = k / std::sqrt(n);
        double v = u, delta = 1;
        while (delta > 1e-8) {
            double n_v = u * std::sqrt(spence(std::exp(-v)));
            delta = std::abs(n_v - v);
            v = n_v;
        }
        double lf = std::log(v) - std::log1p(-std::exp(-v) * (1 + u * u / 2)) / 2 -
                    LOG2 * 1.5 - std::log(u) - std::log(M_PI);
        double g = 2 * v / u - u * std::log1p(-std::exp(-v));
        return lf - std::log(n) + std::sqrt(n) * g;
    }
    // spence(x) = Li2(1 - x) for 0 <= x <= 1 (scipy/cephes convention)
    static double spence(double x) {
        if (x <= 0)
            return M_PI * M_PI / 6;
        if (x >= 1)
            return 0;
        if (x < 0.5)
            return M_PI * M_PI / 6 - std::log1p(-x) * std::log(x) - dilog(x);
        return dilog(1 - x);
    }
    static double dilog(double z) {
        if (z <= 0)
            return 0;
        if (z >= 1)
            return M_PI * M_PI / 6;
        if (z > 0.5)
            return M_PI * M_PI / 6 - std::log(z) * std::log1p(-z) - dilog(1 - z);
        double sum = 0, term = z;
        for (int k = 1; k < 200; ++k) {
            double add = term / (double(k) * k);
            sum += add;
            if (add < 1e-18 * sum)
                break;
            term *= z;
        }
        return sum;
    }
private:
    size_t _n_max = 0;
    std::vector<double> _table;
    double& at(size_t n, size_t k) { return _table[n * (n - 1) / 2 + (k - 1)]; }
    double get(size_t n, size_t k) const { return _table[n * (n - 1) / 2 + (k - 1)]; }
    double get_raw(size_t n, size_t k) const {
        if (k > n)
            k = n;
        if (n == 0)
            return 0;
        return get(n, k);
    }
};

LogQ& log_q_cache() {
    static LogQ cache;
    return cache;
}

// ---------------------------------------------------------------------------
// Flat block state

struct Delta {
    int x, y;
    long d;
};

// Accumulates group-pair count changes, normalising undirected keys to x <= y.
class DeltaAcc {
public:
    void reset(bool directed) {
        _directed = directed;
        clear();
        _idx.reserve(64);
    }
    void clear() {
        _v.clear();
        _idx.clear();
    }
    void add(int x, int y, long d) {
        if (!_directed && x > y)
            std::swap(x, y);
        uint64_t key = pair_key(static_cast<uint64_t>(x), static_cast<uint64_t>(y));
        auto it = _idx.find(key);
        if (it == _idx.end()) {
            _idx.emplace(key, _v.size());
            _v.push_back({x, y, d});
        } else {
            _v[it->second].d += d;
        }
    }
    const std::vector<Delta>& items() const { return _v; }
private:
    bool _directed = false;
    std::vector<Delta> _v;
    gt_hash_map<uint64_t, size_t> _idx;
};

// Changes passed from one level to the level above: group-pair count changes
// (node-pair changes for the upper level) and nodes whose weight changes.
struct ChangeSet {
    std::vector<Delta> edges;
    std::vector<std::pair<int, int>> weights;
    void clear() {
        edges.clear();
        weights.clear();
    }
    bool empty() const { return edges.empty() && weights.empty(); }
};

// Per-thread scratch buffers for lock-free read-only virtual_move / virtual_merge.
struct MergeScratch {
    std::vector<Delta> deltas;
    std::vector<ChangeSet> cs;
    DeltaAcc acc, pairs;
    std::vector<std::pair<int, int>> wd;
    gt_hash_map<int, long> cnt_out, cnt_in;

    void init(bool directed, int levels) {
        acc.reset(directed);
        pairs.reset(directed);
        cs.resize(static_cast<size_t>(levels) + 1);
    }
};

class SBM {
public:
    SBM(const Graph& graph, bool dc, bool dist, size_t q_cache_max, bool edges_dl = true)
        : _directed(graph.directed), _dc(dc), _dist(dist), _edges_dl(edges_dl) {
        const Adjacency& g = graph.adjacency;
        _N = static_cast<int>(boost::num_vertices(g));
        _E = static_cast<long>(boost::num_edges(g));
        math_tables().ensure(static_cast<size_t>(std::max<int64_t>(2 * _E + _N + 16, 4096)));

        _out.resize(_N);
        _in.resize(_N);
        _all.resize(_N);
        _loops.assign(_N, 0);
        _kout.assign(_N, 0);
        _kin.assign(_N, 0);
        gt_hash_map<uint64_t, long> mult;
        mult.reserve(static_cast<size_t>(_E));
        auto range = boost::edges(g);
        for (auto it = range.first; it != range.second; ++it) {
            int s = static_cast<int>(it->s), t = static_cast<int>(it->t);
            if (s == t) {
                _loops[s]++;
                _all[s].push_back(s);
                _all[s].push_back(s);
            } else {
                _out[s].push_back(t);
                _all[s].push_back(t);
                _all[t].push_back(s);
                if (_directed)
                    _in[t].push_back(s);
                else
                    _out[t].push_back(s);
            }
            if (_directed) {
                _kout[s]++;
                _kin[t]++;
            } else {
                _kout[s]++;
                _kout[t]++;
            }
            uint64_t a = static_cast<uint64_t>(s), c = static_cast<uint64_t>(t);
            if (!_directed && a > c)
                std::swap(a, c);
            mult[pair_key(a, c)]++;
        }
        _const = 0;
        for (const auto& [key, m] : mult) {
            _const += lgam_int(m + 1);
            if (!_directed && (key >> 32) == (key & 0xffffffffULL))
                _const += m * LOG2;
        }
        if (_dc) {
            for (int v = 0; v < _N; ++v)
                _const -= lgam_int(_kout[v] + 1) + (_directed ? lgam_int(_kin[v] + 1) : 0);
        }
        if (_dc && _dist) {
            size_t n_max = std::min<size_t>(std::max<long>(2 * _E, 100), 10000);
            n_max = std::min(n_max, q_cache_max);
            log_q_cache().ensure(n_max);
        }

        // Pre-allocate state storage once; set_partition() resets only occupied groups.
        _b.assign(_N, 0);
        _wr.assign(_N, 0);
        _ep.assign(_N, 0);
        _em.assign(_N, 0);
        _m_out.resize(_N);
        if (_directed)
            _m_in.resize(_N);
        if (_dc && _dist) {
            _h_out.resize(_N);
            if (_directed)
                _h_in.resize(_N);
        }
        _members.resize(_N);
        _pos.assign(_N, 0);
        _occ_pos.assign(_N, -1);
        _empty_pos.assign(_N, -1);
        _empty.reserve(_N);
        for (int r = _N - 1; r >= 0; --r)
            push_empty(r);
        _cnt_out.assign(_N, 0);
        if (_directed)
            _cnt_in.assign(_N, 0);
    }

    int N() const { return _N; }
    int B() const { return static_cast<int>(_occupied.size()); }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    bool directed() const { return _directed; }
    long E() const { return _E; }
    // Group-level statistics, used as the adjacency of the level above.
    const std::vector<gt_hash_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<gt_hash_map<int, long>>& m_in() const { return _m_in; }
    const std::vector<long>& ep() const { return _ep; }
    const std::vector<long>& em() const { return _em; }
    const std::vector<int>& wr() const { return _wr; }
    // Group-pair count changes of the last virtual_move(), move(), virtual_merge(), or merge().
    const std::vector<Delta>& deltas() const { return _deltas; }

    void set_partition(const std::vector<int>& b) {
        invalidate_cache();
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _m_out[r].clear();
            if (_directed)
                _m_in[r].clear();
            if (_dc && _dist) {
                _h_out[r].clear();
                if (_directed)
                    _h_in[r].clear();
            }
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        _b = b;
        for (int v = 0; v < _N; ++v) {
            int r = _b[v];
            if (r < 0 || r >= _N)
                Rcpp::stop("Internal error: block label out of range");
            if (_wr[r] == 0)
                take_empty(r);
            add_member(v, r);
            _ep[r] += _kout[v];
            _em[r] += _kin[v];
            if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                if (_directed)
                    _h_in[r][_kin[v]]++;
            }
        }
        for (int v = 0; v < _N; ++v) {
            int r = _b[v];
            auto& row = _m_out[r];
            for (int u : _out[v]) {
                int bu = _b[u];
                row[bu]++;
                if (_directed)
                    _m_in[bu][r]++;
            }
            if (_loops[v] > 0) {
                if (_directed) {
                    row[r] += _loops[v];
                    _m_in[r][r] += _loops[v];
                } else {
                    row[r] += 2 * _loops[v];
                }
            }
        }
    }

    double global_term(long B) const {
        if (_N == 0)
            return 0;
        int64_t pairs = _directed ? int64_t(B) * B : int64_t(B) * (B + 1) / 2;
        return lbinom_int(_N - 1, B - 1) + lgam_int(_N + 1) + log_int(_N) +
               (_edges_dl ? lbinom_int(pairs + _E - 1, _E) : 0);
    }

    double group_term(long n, long ep, long em) const {
        if (n == 0)
            return 0;
        double S = -lgam_int(n + 1);
        if (_dc) {
            S += lgam_int(ep + 1) + (_directed ? lgam_int(em + 1) : 0);
            if (_dist) {
                const LogQ& q = log_q_cache();
                S += q(ep, n) + (_directed ? q(em, n) + 2 * lgam_int(n + 1) : lgam_int(n + 1));
            } else {
                S += lbinom_int(n + ep - 1, ep) + (_directed ? lbinom_int(n + em - 1, em) : 0);
            }
        } else {
            S += (ep + (_directed ? em : 0)) * log_int(n);
        }
        return S;
    }

    double edge_term(int x, int y, long m) const {
        if (m <= 0)
            return 0;
        if (!_directed && x == y) {
            long half = m >> 1;
            return lgam_int(half + 1) + half * LOG2;
        }
        return lgam_int(m + 1);
    }

    long get_m(int x, int y) const {
        const auto& row = _m_out[x];
        auto it = row.find(y);
        return it == row.end() ? 0 : it->second;
    }

    double entropy() const {
        if (_N == 0)
            return 0;
        double S = _const + global_term(B());
        for (int r : _occupied) {
            S += group_term(_wr[r], _ep[r], _em[r]);
            if (_dc && _dist) {
                for (const auto& kc : _h_out[r])
                    S -= lgam_int(kc.second + 1);
                if (_directed)
                    for (const auto& kc : _h_in[r])
                        S -= lgam_int(kc.second + 1);
            }
            for (const auto& [t, m] : _m_out[r])
                if (_directed || t >= r)
                    S -= edge_term(r, t, m);
        }
        return S;
    }

    double virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s)
            return 0;
        double dS = 0;
        int dB = (_wr[s] == 0 ? 1 : 0) - (_wr[r] == 1 ? 1 : 0);
        if (dB != 0)
            dS += global_term(B() + dB) - global_term(B());
        long ko = _kout[v], ki = _kin[v];
        dS += group_term(_wr[r] - 1, _ep[r] - ko, _em[r] - ki) -
              group_term(_wr[r], _ep[r], _em[r]);
        dS += group_term(_wr[s] + 1, _ep[s] + ko, _em[s] + ki) -
              group_term(_wr[s], _ep[s], _em[s]);
        if (_dc && _dist) {
            dS += hist_delta(_h_out, r, s, ko);
            if (_directed)
                dS += hist_delta(_h_in, r, s, ki);
        }
        collect_deltas(v, r, s);
        _last_v = v;
        _last_r = r;
        _last_s = s;
        _last_mr = -1;
        for (const Delta& e : _deltas) {
            long m = get_m(e.x, e.y);
            dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
        }
        return dS;
    }

    void move(int v, int s) {
        int r = _b[v];
        if (r == s)
            return;
        if (_last_v != v || _last_r != r || _last_s != s)
            collect_deltas(v, r, s);
        invalidate_cache();
        for (const Delta& e : _deltas) {
            add_m(_m_out[e.x], e.y, e.d);
            if (_directed)
                add_m(_m_in[e.y], e.x, e.d);
            else if (e.x != e.y)
                add_m(_m_out[e.y], e.x, e.d);
        }
        long ko = _kout[v], ki = _kin[v];
        _ep[r] -= ko;
        _em[r] -= ki;
        _ep[s] += ko;
        _em[s] += ki;
        if (_dc && _dist) {
            hist_move(_h_out, r, s, ko);
            if (_directed)
                hist_move(_h_in, r, s, ki);
        }
        remove_member(v, r);
        if (_wr[s] == 0)
            take_empty(s);
        add_member(v, s);
        _b[v] = s;
        if (_wr[r] == 0)
            push_empty(r);
    }

    // Direct O(deg_group(r) + |hist(r)|) evaluation of merging group r into s.
    double virtual_merge(int r, int s) {
        if (r == s)
            return 0;
        double dS = global_term(B() - 1) - global_term(B());
        dS += group_term(_wr[s] + _wr[r], _ep[s] + _ep[r], _em[s] + _em[r]) -
              group_term(_wr[s], _ep[s], _em[s]) -
              group_term(_wr[r], _ep[r], _em[r]);
        if (_dc && _dist) {
            for (const auto& [k, cr] : _h_out[r]) {
                auto it = _h_out[s].find(k);
                if (it != _h_out[s].end()) {
                    long cs = it->second;
                    dS += lgam_int(cr + 1) + lgam_int(cs + 1) - lgam_int(cr + cs + 1);
                }
            }
            if (_directed) {
                for (const auto& [k, cr] : _h_in[r]) {
                    auto it = _h_in[s].find(k);
                    if (it != _h_in[s].end()) {
                        long cs = it->second;
                        dS += lgam_int(cr + 1) + lgam_int(cs + 1) - lgam_int(cr + cs + 1);
                    }
                }
            }
        }
        collect_merge_deltas(r, s);
        _last_v = -1;
        _last_mr = r;
        _last_ms = s;
        for (const Delta& e : _deltas) {
            if (e.d < 0) {
                // Every negative merge delta zeroes out an edge incident to r.
                dS += edge_term(e.x, e.y, -e.d);
            } else {
                long m = get_m(e.x, e.y);
                dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
            }
        }
        return dS;
    }

    void merge(int r, int s) {
        if (r == s)
            return;
        if (_last_mr != r || _last_ms != s)
            collect_merge_deltas(r, s);
        invalidate_cache();
        for (const Delta& e : _deltas) {
            add_m(_m_out[e.x], e.y, e.d);
            if (_directed)
                add_m(_m_in[e.y], e.x, e.d);
            else if (e.x != e.y)
                add_m(_m_out[e.y], e.x, e.d);
        }
        _ep[s] += _ep[r];
        _em[s] += _em[r];
        _ep[r] = 0;
        _em[r] = 0;
        if (_dc && _dist) {
            for (const auto& [k, cr] : _h_out[r])
                _h_out[s][k] += cr;
            _h_out[r].clear();
            if (_directed) {
                for (const auto& [k, cr] : _h_in[r])
                    _h_in[s][k] += cr;
                _h_in[r].clear();
            }
        }
        auto& vr = _members[r];
        auto& vs = _members[s];
        vs.reserve(vs.size() + vr.size());
        for (int v : vr) {
            _pos[v] = static_cast<int>(vs.size());
            vs.push_back(v);
            _b[v] = s;
        }
        _wr[s] += _wr[r];
        vr.clear();
        _wr[r] = 0;
        remove_occupied(r);
        push_empty(r);
    }

    // Thread-safe read-only virtual_move using caller-supplied scratch (zero O(N) allocs).
    double virtual_move_const(int v, int s, MergeScratch& sc) const {
        int r = _b[v];
        if (r == s)
            return 0;
        double dS = 0;
        int dB = (_wr[s] == 0 ? 1 : 0) - (_wr[r] == 1 ? 1 : 0);
        if (dB != 0)
            dS += global_term(B() + dB) - global_term(B());
        long ko = _kout[v], ki = _kin[v];
        dS += group_term(_wr[r] - 1, _ep[r] - ko, _em[r] - ki) -
              group_term(_wr[r], _ep[r], _em[r]);
        dS += group_term(_wr[s] + 1, _ep[s] + ko, _em[s] + ki) -
              group_term(_wr[s], _ep[s], _em[s]);
        if (_dc && _dist) {
            dS += hist_delta(_h_out, r, s, ko);
            if (_directed)
                dS += hist_delta(_h_in, r, s, ki);
        }
        collect_deltas_const(v, r, s, sc);
        for (const Delta& e : sc.deltas) {
            long m = get_m(e.x, e.y);
            dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
        }
        return dS;
    }

    // Thread-safe read-only virtual_merge using caller-supplied scratch.
    double virtual_merge_const(int r, int s, MergeScratch& sc) const {
        if (r == s)
            return 0;
        double dS = global_term(B() - 1) - global_term(B());
        dS += group_term(_wr[s] + _wr[r], _ep[s] + _ep[r], _em[s] + _em[r]) -
              group_term(_wr[s], _ep[s], _em[s]) -
              group_term(_wr[r], _ep[r], _em[r]);
        if (_dc && _dist) {
            for (const auto& [k, cr] : _h_out[r]) {
                auto it = _h_out[s].find(k);
                if (it != _h_out[s].end()) {
                    long cs = it->second;
                    dS += lgam_int(cr + 1) + lgam_int(cs + 1) - lgam_int(cr + cs + 1);
                }
            }
            if (_directed) {
                for (const auto& [k, cr] : _h_in[r]) {
                    auto it = _h_in[s].find(k);
                    if (it != _h_in[s].end()) {
                        long cs = it->second;
                        dS += lgam_int(cr + 1) + lgam_int(cs + 1) - lgam_int(cr + cs + 1);
                    }
                }
            }
        }
        collect_merge_deltas_into(r, s, sc.deltas);
        for (const Delta& e : sc.deltas) {
            if (e.d < 0) {
                dS += edge_term(e.x, e.y, -e.d);
            } else {
                long m = get_m(e.x, e.y);
                dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
            }
        }
        return dS;
    }

    template <class Rng>
    int random_occupied(Rng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    long m_sym(int x, int y) const {
        return get_m(x, y) + (_directed ? get_m(y, x) : 0);
    }

    long e_tot(int r) const { return _ep[r] + (_directed ? _em[r] : 0); }

    int sample_group_neighbor(int t, RRng& rng) const {
        double total = static_cast<double>(e_tot(t));
        if (total <= 0)
            return random_occupied(rng);
        double u = rng.unif() * total;
        for (const auto& [s, m] : _m_out[t]) {
            u -= m;
            if (u < 0)
                return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in[t]) {
                u -= m;
                if (u < 0)
                    return s;
            }
        }
        return t;
    }

    // Exact neighbor-based proposal for finite-temperature MCMC (Peixoto 2014).
    int sample_exact_group(int v, double c, RRng& rng) const {
        if (std::isinf(c) || _all[v].empty())
            return random_occupied(rng);
        int u = _all[v][rng.index(_all[v].size())];
        int t = _b[u];
        double et = static_cast<double>(e_tot(t));
        double B_occ = static_cast<double>(_occupied.size());
        double p_rand = (c > 0) ? (c * B_occ / (et + c * B_occ)) : 0.0;
        if (et <= 0 || (c > 0 && rng.unif() < p_rand))
            return random_occupied(rng);
        return sample_group_neighbor(t, rng);
    }

    // Exact probability q(r -> s | v) (reverse=false) or post-move q'(s -> r | v)
    // (reverse=true, assuming r stays occupied after moving v: r -> s).
    double move_prob(int v, int r, int s, double c, bool reverse) {
        double kv = static_cast<double>(_all[v].size());
        if (!reverse) {
            double B_occ = static_cast<double>(_occupied.size());
            if (kv <= 0 || std::isinf(c))
                return 1.0 / B_occ;
            _touched_out.clear();
            auto add_w = [&](int t, long w) {
                if (_cnt_out[t] == 0)
                    _touched_out.push_back(t);
                _cnt_out[t] += w;
            };
            for (int u : _out[v])
                add_w(_b[u], 1);
            if (_directed) {
                for (int u : _in[v])
                    add_w(_b[u], 1);
            }
            if (_loops[v] > 0)
                add_w(r, 2 * _loops[v]);
            double prob = 0;
            for (int t : _touched_out) {
                double wt = static_cast<double>(_cnt_out[t]);
                _cnt_out[t] = 0;
                double et = static_cast<double>(e_tot(t));
                double mts = static_cast<double>(m_sym(t, s));
                prob += (wt / kv) * ((mts + c) / (et + c * B_occ));
            }
            return prob;
        }
        double B_new = static_cast<double>(_occupied.size() + (_wr[s] == 0 ? 1 : 0));
        if (kv <= 0 || std::isinf(c))
            return 1.0 / B_new;
        _touched_out.clear();
        auto add_w = [&](int t, long w) {
            if (_cnt_out[t] == 0)
                _touched_out.push_back(t);
            _cnt_out[t] += w;
        };
        for (int u : _out[v])
            add_w(_b[u], 1);
        if (_directed) {
            for (int u : _in[v])
                add_w(_b[u], 1);
        }
        long lp_w = 2 * _loops[v];
        long wnl_r = _cnt_out[r];
        long wnl_s = _cnt_out[s];
        if (lp_w > 0 && _cnt_out[s] == 0)
            _touched_out.push_back(s);
        double prob = 0;
        for (int t : _touched_out) {
            long wnl = _cnt_out[t];
            _cnt_out[t] = 0;
            double wt = static_cast<double>(wnl + (t == s ? lp_w : 0));
            if (wt <= 0)
                continue;
            double et_new, mtr_new;
            if (t == r) {
                et_new = static_cast<double>(e_tot(r) - static_cast<long>(kv));
                mtr_new = static_cast<double>(m_sym(r, r) - 2 * wnl_r - lp_w);
            } else if (t == s) {
                et_new = static_cast<double>(e_tot(s) + static_cast<long>(kv));
                mtr_new = static_cast<double>(m_sym(s, r) + wnl_r - wnl_s);
            } else {
                et_new = static_cast<double>(e_tot(t));
                mtr_new = static_cast<double>(m_sym(t, r) - wnl);
            }
            prob += (wt / kv) * ((mtr_new + c) / (et_new + c * B_new));
        }
        return prob;
    }

    // Group-level merge target sampling and exact unordered pair probability.
    int sample_merge_target(int r, double c, RRng& rng) const {
        size_t B_occ = _occupied.size();
        if (B_occ <= 1)
            return r;
        auto rand_other = [&]() {
            int s = _occupied[rng.index(B_occ - 1)];
            return (s == r) ? _occupied[B_occ - 1] : s;
        };
        double ext = static_cast<double>(e_tot(r) - m_sym(r, r));
        if (ext <= 0 || std::isinf(c))
            return rand_other();
        double B_other = static_cast<double>(B_occ - 1);
        double p_rand = (c > 0) ? (c * B_other / (ext + c * B_other)) : 0.0;
        if (c > 0 && rng.unif() < p_rand)
            return rand_other();
        double u = rng.unif() * ext;
        for (const auto& [s, m] : _m_out[r]) {
            if (s == r)
                continue;
            u -= m;
            if (u < 0)
                return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in[r]) {
                if (s == r)
                    continue;
                u -= m;
                if (u < 0)
                    return s;
            }
        }
        return rand_other();
    }

    double merge_prob(int r, int s, double c) const {
        size_t B_occ = _occupied.size();
        if (B_occ <= 1 || r == s)
            return 0;
        double B_other = static_cast<double>(B_occ - 1);
        double mrs = static_cast<double>(m_sym(r, s));
        auto q_dir = [&](int a) {
            double ext = static_cast<double>(e_tot(a) - m_sym(a, a));
            if (ext <= 0 || std::isinf(c))
                return 1.0 / B_other;
            return (mrs + c) / (ext + c * B_other);
        };
        return (q_dir(r) + q_dir(s)) / static_cast<double>(B_occ);
    }

    // Neighbour-based proposal of graph-tool's BlockState::sample_group with
    // d = 0 and the fast (non-edge-group) proposal.
    template <class Rng>
    int sample_group(int v, double c, Rng& rng) const {
        if (std::isinf(c) || _all[v].empty())
            return random_occupied(rng);
        int u = _all[v][rng.index(_all[v].size())];
        int t = _b[u];
        double p_rand = 0;
        if (c > 0) {
            double B = _occupied.size();
            double et = _ep[t] + (_directed ? _em[t] : 0);
            p_rand = c * B / (et + c * B);
        }
        if (c == 0 || rng.unif() >= p_rand) {
            const auto& vs = _members[t];
            int w = vs[rng.index(vs.size())];
            if (!_all[w].empty())
                w = _all[w][rng.index(_all[w].size())];
            return _b[w];
        }
        return random_occupied(rng);
    }

private:
    bool _directed, _dc, _dist, _edges_dl;
    int _N;
    long _E;
    double _const;
    std::vector<std::vector<int>> _out, _in, _all;
    std::vector<long> _loops, _kout, _kin;

    std::vector<int> _b, _wr, _pos, _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<long> _ep, _em;
    std::vector<gt_hash_map<int, long>> _m_out, _m_in;
    std::vector<gt_hash_map<long, long>> _h_out, _h_in;
    std::vector<std::vector<int>> _members;

    std::vector<Delta> _deltas;
    std::vector<long> _cnt_out, _cnt_in;
    std::vector<int> _touched_out, _touched_in;
    int _last_v = -1, _last_r = -1, _last_s = -1;
    int _last_mr = -1, _last_ms = -1;

    void invalidate_cache() {
        _last_v = -1;
        _last_mr = -1;
    }

    static void add_m(gt_hash_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) {
            row.emplace(y, d);
        } else {
            it->second += d;
            if (it->second == 0)
                row.erase(it);
        }
    }

    double hist_delta(const std::vector<gt_hash_map<long, long>>& h,
                      int r, int s, long k) const {
        auto itr = h[r].find(k);
        long cr = itr == h[r].end() ? 0L : itr->second;
        auto its = h[s].find(k);
        long cs = its == h[s].end() ? 0L : its->second;
        // -(lgam(cr) - lgam(cr + 1)) - (lgam(cs + 2) - lgam(cs + 1)) = log(cr) - log(cs + 1)
        return log_int(cr) - log_int(cs + 1);
    }

    static void hist_move(std::vector<gt_hash_map<long, long>>& h,
                          int r, int s, long k) {
        auto it = h[r].find(k);
        if (--(it->second) == 0)
            h[r].erase(it);
        h[s][k]++;
    }

    void add_member(int v, int r) {
        if (_wr[r] == 0) {
            _occ_pos[r] = static_cast<int>(_occupied.size());
            _occupied.push_back(r);
        }
        _pos[v] = static_cast<int>(_members[r].size());
        _members[r].push_back(v);
        _wr[r]++;
    }

    void remove_occupied(int r) {
        int i = _occ_pos[r];
        int other = _occupied.back();
        _occupied[i] = other;
        _occ_pos[other] = i;
        _occupied.pop_back();
        _occ_pos[r] = -1;
    }

    void remove_member(int v, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[v]] = last;
        _pos[last] = _pos[v];
        vs.pop_back();
        if (--_wr[r] == 0)
            remove_occupied(r);
    }

    void push_empty(int r) {
        _empty_pos[r] = static_cast<int>(_empty.size());
        _empty.push_back(r);
    }

    void take_empty(int s) {
        int i = _empty_pos[s];
        if (i >= 0) {
            int last = _empty.back();
            _empty[i] = last;
            _empty_pos[last] = i;
            _empty.pop_back();
            _empty_pos[s] = -1;
        }
    }

    void collect_deltas(int v, int r, int s) {
        _deltas.clear();
        _touched_out.clear();
        if (_directed) {
            _touched_in.clear();
            for (int u : _out[v]) {
                int t = _b[u];
                if (_cnt_out[t] == 0)
                    _touched_out.push_back(t);
                _cnt_out[t]++;
            }
            for (int u : _in[v]) {
                int t = _b[u];
                if (_cnt_in[t] == 0)
                    _touched_in.push_back(t);
                _cnt_in[t]++;
            }
            long co_r = _cnt_out[r], co_s = _cnt_out[s];
            long ci_r = _cnt_in[r], ci_s = _cnt_in[s];
            long lp = _loops[v];
            for (int t : _touched_out) {
                long c = _cnt_out[t];
                _cnt_out[t] = 0;
                if (t == r || t == s)
                    continue;
                _deltas.push_back({r, t, -c});
                _deltas.push_back({s, t, c});
            }
            for (int t : _touched_in) {
                long c = _cnt_in[t];
                _cnt_in[t] = 0;
                if (t == r || t == s)
                    continue;
                _deltas.push_back({t, r, -c});
                _deltas.push_back({t, s, c});
            }
            long d_rr = -co_r - ci_r - lp;
            long d_rs = -co_s + ci_r;
            long d_sr = co_r - ci_s;
            long d_ss = co_s + ci_s + lp;
            if (d_rr != 0) _deltas.push_back({r, r, d_rr});
            if (d_rs != 0) _deltas.push_back({r, s, d_rs});
            if (d_sr != 0) _deltas.push_back({s, r, d_sr});
            if (d_ss != 0) _deltas.push_back({s, s, d_ss});
        } else {
            for (int u : _out[v]) {
                int t = _b[u];
                if (_cnt_out[t] == 0)
                    _touched_out.push_back(t);
                _cnt_out[t]++;
            }
            long cr = _cnt_out[r];
            long cs = _cnt_out[s];
            long lp = _loops[v];
            for (int t : _touched_out) {
                long c = _cnt_out[t];
                _cnt_out[t] = 0;
                if (t == r || t == s)
                    continue;
                int x1 = r < t ? r : t, y1 = r < t ? t : r;
                int x2 = s < t ? s : t, y2 = s < t ? t : s;
                _deltas.push_back({x1, y1, -c});
                _deltas.push_back({x2, y2, c});
            }
            long d_rr = -2 * cr - 2 * lp;
            long d_rs = cr - cs;
            long d_ss = 2 * cs + 2 * lp;
            if (d_rr != 0) _deltas.push_back({r, r, d_rr});
            if (d_rs != 0) {
                int x = r < s ? r : s, y = r < s ? s : r;
                _deltas.push_back({x, y, d_rs});
            }
            if (d_ss != 0) _deltas.push_back({s, s, d_ss});
        }
    }

    void collect_deltas_const(int v, int r, int s, MergeScratch& sc) const {
        auto& deltas = sc.deltas;
        deltas.clear();
        if (_directed) {
            auto& co = sc.cnt_out;
            auto& ci = sc.cnt_in;
            co.clear();
            ci.clear();
            for (int u : _out[v]) co[_b[u]]++;
            for (int u : _in[v]) ci[_b[u]]++;
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k);
                return it == mp.end() ? 0L : it->second;
            };
            long co_r = get_c(co, r), co_s = get_c(co, s);
            long ci_r = get_c(ci, r), ci_s = get_c(ci, s);
            long lp = _loops[v];
            for (const auto& [t, c] : co) {
                if (t == r || t == s) continue;
                deltas.push_back({r, t, -c});
                deltas.push_back({s, t, c});
            }
            for (const auto& [t, c] : ci) {
                if (t == r || t == s) continue;
                deltas.push_back({t, r, -c});
                deltas.push_back({t, s, c});
            }
            long d_rr = -co_r - ci_r - lp;
            long d_rs = -co_s + ci_r;
            long d_sr = co_r - ci_s;
            long d_ss = co_s + ci_s + lp;
            if (d_rr != 0) deltas.push_back({r, r, d_rr});
            if (d_rs != 0) deltas.push_back({r, s, d_rs});
            if (d_sr != 0) deltas.push_back({s, r, d_sr});
            if (d_ss != 0) deltas.push_back({s, s, d_ss});
        } else {
            auto& co = sc.cnt_out;
            co.clear();
            for (int u : _out[v]) co[_b[u]]++;
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k);
                return it == mp.end() ? 0L : it->second;
            };
            long cr = get_c(co, r), cs = get_c(co, s);
            long lp = _loops[v];
            for (const auto& [t, c] : co) {
                if (t == r || t == s) continue;
                int x1 = r < t ? r : t, y1 = r < t ? t : r;
                int x2 = s < t ? s : t, y2 = s < t ? t : s;
                deltas.push_back({x1, y1, -c});
                deltas.push_back({x2, y2, c});
            }
            long d_rr = -2 * cr - 2 * lp;
            long d_rs = cr - cs;
            long d_ss = 2 * cs + 2 * lp;
            if (d_rr != 0) deltas.push_back({r, r, d_rr});
            if (d_rs != 0) {
                int x = r < s ? r : s, y = r < s ? s : r;
                deltas.push_back({x, y, d_rs});
            }
            if (d_ss != 0) deltas.push_back({s, s, d_ss});
        }
    }

    void collect_merge_deltas_into(int r, int s, std::vector<Delta>& deltas) const {
        deltas.clear();
        if (_directed) {
            long m_rr = 0, m_rs = 0, m_sr = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) {
                    m_rr = m;
                } else if (t == s) {
                    m_rs = m;
                } else {
                    deltas.push_back({r, t, -m});
                    deltas.push_back({s, t, m});
                }
            }
            for (const auto& [t, m] : _m_in[r]) {
                if (t == r) {
                    continue;
                } else if (t == s) {
                    m_sr = m;
                } else {
                    deltas.push_back({t, r, -m});
                    deltas.push_back({t, s, m});
                }
            }
            if (m_rr != 0) deltas.push_back({r, r, -m_rr});
            if (m_rs != 0) deltas.push_back({r, s, -m_rs});
            if (m_sr != 0) deltas.push_back({s, r, -m_sr});
            long d_ss = m_rr + m_rs + m_sr;
            if (d_ss != 0) deltas.push_back({s, s, d_ss});
        } else {
            long m_rr = 0, m_rs = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) {
                    m_rr = m;
                } else if (t == s) {
                    m_rs = m;
                } else {
                    int x1 = r < t ? r : t, y1 = r < t ? t : r;
                    int x2 = s < t ? s : t, y2 = s < t ? t : s;
                    deltas.push_back({x1, y1, -m});
                    deltas.push_back({x2, y2, m});
                }
            }
            if (m_rr != 0) deltas.push_back({r, r, -m_rr});
            if (m_rs != 0) {
                int x = r < s ? r : s, y = r < s ? s : r;
                deltas.push_back({x, y, -m_rs});
            }
            long d_ss = m_rr + 2 * m_rs;
            if (d_ss != 0) deltas.push_back({s, s, d_ss});
        }
    }

    void collect_merge_deltas(int r, int s) {
        collect_merge_deltas_into(r, s, _deltas);
    }
};

// ---------------------------------------------------------------------------
// Nested (hierarchical) block model
//
// Level 0 is the SBM above without the edge-count prior. Level l >= 1 is a
// non-degree-corrected, dense, multigraph block state whose nodes are the
// level-(l-1) groups, with unit weight for nonempty groups and edge
// multiplicities equal to the level-(l-1) group-pair counts. Only the top level
// carries the edge-count prior. This is graph-tool's NestedBlockState with the
// default BlockState arguments.

// The statistics of the level below, which form this level's weighted graph.
struct LowerView {
    const std::vector<gt_hash_map<int, long>>* m_out;
    const std::vector<gt_hash_map<int, long>>* m_in;
    const std::vector<long>* ep;
    const std::vector<long>* em;
    const std::vector<int>* wr;
    const std::vector<int>* occupied;
};

class BlockLevel {
public:
    BlockLevel(int cap, bool directed, long E, bool top)
        : _cap(cap), _directed(directed), _top(top), _E(E),
          _b(cap, 0), _wr(cap, 0), _pos(cap, 0), _occ_pos(cap, -1),
          _empty_pos(cap, -1), _ep(cap, 0), _em(cap, 0),
          _m_out(cap), _m_in(directed ? cap : 0), _members(cap),
          _cnt_out(cap, 0), _cnt_in(directed ? cap : 0) {
        _acc.reset(directed);
        _pairs.reset(directed);
        _empty.reserve(cap);
        for (int r = _cap - 1; r >= 0; --r)
            push_empty(r);
    }

    int B() const { return static_cast<int>(_occupied.size()); }
    int N() const { return _N; }
    int group(int x) const { return _b[x]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    void ensure_cap(int need_cap) {
        if (need_cap <= _cap)
            return;
        int new_cap = std::max(need_cap, _cap * 2);
        int old_cap = _cap;
        _cap = new_cap;
        _b.resize(new_cap, 0);
        _wr.resize(new_cap, 0);
        _pos.resize(new_cap, 0);
        _occ_pos.resize(new_cap, -1);
        _empty_pos.resize(new_cap, -1);
        _ep.resize(new_cap, 0);
        _em.resize(new_cap, 0);
        _m_out.resize(new_cap);
        if (_directed) _m_in.resize(new_cap);
        _members.resize(new_cap);
        _cnt_out.resize(new_cap, 0);
        if (_directed) _cnt_in.resize(new_cap, 0);
        _empty.reserve(new_cap);
        for (int r = new_cap - 1; r >= old_cap; --r)
            push_empty(r);
    }

    template <class Rng>
    int random_occupied(Rng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    LowerView view() const { return {&_m_out, &_m_in, &_ep, &_em, &_wr, &_occupied}; }

    void build(const std::vector<int>& b, const LowerView& lower) {
        int max_id = 0;
        for (int x : *lower.occupied) {
            if (x > max_id) max_id = x;
            int g = (x < static_cast<int>(b.size())) ? b[x] : 0;
            if (g > max_id) max_id = g;
        }
        ensure_cap(max_id + 1);
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _m_out[r].clear();
            if (_directed)
                _m_in[r].clear();
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        for (int x : *lower.occupied) {
            int g = (x < static_cast<int>(b.size())) ? b[x] : 0;
            _b[x] = g;
        }
        _N = 0;
        for (int x : *lower.occupied) {
            int g = _b[x];
            if (_wr[g] == 0)
                take_empty(g);
            add_member(x, g);
            ++_N;
        }
        for (int x : *lower.occupied) {
            int g = _b[x];
            auto& row = _m_out[g];
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                int h = _b[y];
                row[h] += val;
                _ep[g] += val;
                if (_directed) {
                    _m_in[h][g] += val;
                    _em[h] += val;
                }
            }
        }
    }

    void build_sparse(const std::vector<std::pair<int, int>>& entries, const LowerView& lower) {
        int max_id = 0;
        for (const auto& [x, g] : entries) {
            if (x > max_id) max_id = x;
            if (g > max_id) max_id = g;
        }
        for (int x : *lower.occupied) {
            if (x > max_id) max_id = x;
        }
        ensure_cap(max_id + 1);
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _m_out[r].clear();
            if (_directed)
                _m_in[r].clear();
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        for (const auto& [x, g] : entries)
            _b[x] = g;
        _N = 0;
        for (int x : *lower.occupied) {
            int g = _b[x];
            if (_wr[g] == 0)
                take_empty(g);
            add_member(x, g);
            ++_N;
        }
        for (int x : *lower.occupied) {
            int g = _b[x];
            auto& row = _m_out[g];
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                int h = _b[y];
                row[h] += val;
                _ep[g] += val;
                if (_directed) {
                    _m_in[h][g] += val;
                    _em[h] += val;
                }
            }
        }
    }

    double partition_term(long N, long B) const {
        if (N == 0)
            return 0;
        return lbinom_int(N - 1, B - 1) + lgam_int(N + 1) + log_int(N);
    }

    double edges_term(long B) const {
        if (!_top)
            return 0;
        int64_t pairs = _directed ? int64_t(B) * B : int64_t(B) * (B + 1) / 2;
        return lbinom_int(pairs + _E - 1, _E);
    }

    double dense_term(int x, int y, long m, long wx, long wy) const {
        if (m <= 0)
            return 0;
        int64_t nn, e = m;
        if (_directed || x != y) {
            nn = int64_t(wx) * wy;
        } else {
            nn = int64_t(wx) * (wx + 1) / 2;
            e = m >> 1;
        }
        return lbinom_int(nn + e - 1, e);
    }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0 : it->second;
    }

    double entropy() const {
        double S = partition_term(_N, B()) + edges_term(B());
        for (int r : _occupied) {
            S -= lgam_int(_wr[r] + 1);
            for (const auto& [t, m] : _m_out[r])
                if (_directed || t >= r)
                    S += dense_term(r, t, m, _wr[r], _wr[t]);
        }
        return S;
    }

    // Changes arriving from the level below.
    double lower_changes(const ChangeSet& in, ChangeSet& out, bool apply) {
        _acc.clear();
        for (const Delta& e : in.edges) {
            int gx = _b[e.x], gy = _b[e.y];
            if (_directed)
                _acc.add(gx, gy, e.d);
            else if (e.x == e.y)
                _acc.add(gx, gx, e.d);
            else if (gx == gy)
                _acc.add(gx, gx, 2 * e.d);
            else
                _acc.add(gx, gy, e.d);
        }
        _wd.clear();
        int dN = 0;
        for (const auto& [x, dw] : in.weights) {
            add_weight(_b[x], dw);
            dN += dw;
        }
        out.clear();
        for (const Delta& e : _acc.items())
            if (e.d != 0)
                out.edges.push_back(e);
        double dS = process(dN, out, apply);
        if (apply) {
            for (const auto& [x, dw] : in.weights) {
                if (dw > 0)
                    push_node(x, _b[x]);
                else
                    pop_node(x, _b[x]);
            }
        }
        return dS;
    }

    // Move node x (a nonempty group of the level below) to group h.
    double node_move(int x, int h, const LowerView& lower, ChangeSet& out, bool apply) {
        int g = _b[x];
        out.clear();
        _touched_out.clear();
        if (_directed) {
            _touched_in.clear();
            long lp = 0;
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) {
                    lp = val;
                } else {
                    int t = _b[y];
                    if (_cnt_out[t] == 0)
                        _touched_out.push_back(t);
                    _cnt_out[t] += val;
                }
            }
            for (const auto& [y, val] : (*lower.m_in)[x]) {
                if (y == x)
                    continue;
                int t = _b[y];
                if (_cnt_in[t] == 0)
                    _touched_in.push_back(t);
                _cnt_in[t] += val;
            }
            long co_g = _cnt_out[g], co_h = _cnt_out[h];
            long ci_g = _cnt_in[g], ci_h = _cnt_in[h];
            for (int t : _touched_out) {
                long c = _cnt_out[t];
                _cnt_out[t] = 0;
                if (t == g || t == h)
                    continue;
                out.edges.push_back({g, t, -c});
                out.edges.push_back({h, t, c});
            }
            for (int t : _touched_in) {
                long c = _cnt_in[t];
                _cnt_in[t] = 0;
                if (t == g || t == h)
                    continue;
                out.edges.push_back({t, g, -c});
                out.edges.push_back({t, h, c});
            }
            long d_gg = -co_g - ci_g - lp;
            long d_gh = -co_h + ci_g;
            long d_hg = co_g - ci_h;
            long d_hh = co_h + ci_h + lp;
            if (d_gg != 0) out.edges.push_back({g, g, d_gg});
            if (d_gh != 0) out.edges.push_back({g, h, d_gh});
            if (d_hg != 0) out.edges.push_back({h, g, d_hg});
            if (d_hh != 0) out.edges.push_back({h, h, d_hh});
        } else {
            long lp = 0;
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) {
                    lp = val;
                } else {
                    int t = _b[y];
                    if (_cnt_out[t] == 0)
                        _touched_out.push_back(t);
                    _cnt_out[t] += val;
                }
            }
            long cg = _cnt_out[g];
            long ch = _cnt_out[h];
            for (int t : _touched_out) {
                long c = _cnt_out[t];
                _cnt_out[t] = 0;
                if (t == g || t == h)
                    continue;
                int x1 = g < t ? g : t, y1 = g < t ? t : g;
                int x2 = h < t ? h : t, y2 = h < t ? t : h;
                out.edges.push_back({x1, y1, -c});
                out.edges.push_back({x2, y2, c});
            }
            long d_gg = -2 * cg - lp;
            long d_gh = cg - ch;
            long d_hh = 2 * ch + lp;
            if (d_gg != 0) out.edges.push_back({g, g, d_gg});
            if (d_gh != 0) {
                int a = g < h ? g : h, b = g < h ? h : g;
                out.edges.push_back({a, b, d_gh});
            }
            if (d_hh != 0) out.edges.push_back({h, h, d_hh});
        }
        _wd.clear();
        _wd.emplace_back(g, -1);
        _wd.emplace_back(h, +1);
        double dS = process(0, out, apply);
        if (apply) {
            pop_node(x, g);
            _b[x] = h;
            push_node(x, h);
        }
        return dS;
    }

    // Direct group-level merge of r into s at this level.
    double group_merge(int r, int s, ChangeSet& out, bool apply) {
        if (r == s)
            return 0;
        out.clear();
        if (_directed) {
            long m_rr = 0, m_rs = 0, m_sr = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) {
                    m_rr = m;
                } else if (t == s) {
                    m_rs = m;
                } else {
                    out.edges.push_back({r, t, -m});
                    out.edges.push_back({s, t, m});
                }
            }
            for (const auto& [t, m] : _m_in[r]) {
                if (t == r) {
                    continue;
                } else if (t == s) {
                    m_sr = m;
                } else {
                    out.edges.push_back({t, r, -m});
                    out.edges.push_back({t, s, m});
                }
            }
            if (m_rr != 0) out.edges.push_back({r, r, -m_rr});
            if (m_rs != 0) out.edges.push_back({r, s, -m_rs});
            if (m_sr != 0) out.edges.push_back({s, r, -m_sr});
            long d_ss = m_rr + m_rs + m_sr;
            if (d_ss != 0) out.edges.push_back({s, s, d_ss});
        } else {
            long m_rr = 0, m_rs = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) {
                    m_rr = m;
                } else if (t == s) {
                    m_rs = m;
                } else {
                    int x1 = r < t ? r : t, y1 = r < t ? t : r;
                    int x2 = s < t ? s : t, y2 = s < t ? t : s;
                    out.edges.push_back({x1, y1, -m});
                    out.edges.push_back({x2, y2, m});
                }
            }
            if (m_rr != 0) out.edges.push_back({r, r, -m_rr});
            if (m_rs != 0) {
                int x = r < s ? r : s, y = r < s ? s : r;
                out.edges.push_back({x, y, -m_rs});
            }
            long d_ss = m_rr + 2 * m_rs;
            if (d_ss != 0) out.edges.push_back({s, s, d_ss});
        }
        _wd.clear();
        _wd.emplace_back(r, -_wr[r]);
        _wd.emplace_back(s, _wr[r]);
        double dS = process(0, out, apply);
        if (apply) {
            auto& vr = _members[r];
            auto& vs = _members[s];
            vs.reserve(vs.size() + vr.size());
            for (int x : vr) {
                _pos[x] = static_cast<int>(vs.size());
                vs.push_back(x);
                _b[x] = s;
            }
            _wr[s] += _wr[r];
            vr.clear();
            _wr[r] = 0;
            remove_occupied(r);
            push_empty(r);
        }
        return dS;
    }

    // Assign an empty (weight-zero) node to a group; this changes no statistics.
    void set_empty_label(int x, int g) {
        ensure_cap(std::max(x, g) + 1);
        _b[x] = g;
    }

    // Thread-safe read-only evaluation of lower_changes, node_move, and group_merge.
    double lower_changes_eval(const ChangeSet& in, ChangeSet& out, MergeScratch& sc) const {
        sc.acc.clear();
        for (const Delta& e : in.edges) {
            int gx = _b[e.x], gy = _b[e.y];
            if (_directed)
                sc.acc.add(gx, gy, e.d);
            else if (e.x == e.y)
                sc.acc.add(gx, gx, e.d);
            else if (gx == gy)
                sc.acc.add(gx, gx, 2 * e.d);
            else
                sc.acc.add(gx, gy, e.d);
        }
        sc.wd.clear();
        int dN = 0;
        for (const auto& [x, dw] : in.weights) {
            int g = _b[x];
            bool found = false;
            for (auto& p : sc.wd) {
                if (p.first == g) { p.second += dw; found = true; break; }
            }
            if (!found) sc.wd.emplace_back(g, dw);
            dN += dw;
        }
        out.clear();
        for (const Delta& e : sc.acc.items())
            if (e.d != 0)
                out.edges.push_back(e);
        return process_eval(dN, out, sc);
    }

    double node_move_eval(int x, int h, const LowerView& lower, ChangeSet& out,
                          MergeScratch& sc) const {
        int g = _b[x];
        out.clear();
        if (_directed) {
            auto& co = sc.cnt_out;
            auto& ci = sc.cnt_in;
            co.clear();
            ci.clear();
            long lp = 0;
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) lp = val;
                else co[_b[y]] += val;
            }
            for (const auto& [y, val] : (*lower.m_in)[x]) {
                if (y == x) continue;
                ci[_b[y]] += val;
            }
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k);
                return it == mp.end() ? 0L : it->second;
            };
            long co_g = get_c(co, g), co_h = get_c(co, h);
            long ci_g = get_c(ci, g), ci_h = get_c(ci, h);
            for (const auto& [t, c] : co) {
                if (t == g || t == h) continue;
                out.edges.push_back({g, t, -c});
                out.edges.push_back({h, t, c});
            }
            for (const auto& [t, c] : ci) {
                if (t == g || t == h) continue;
                out.edges.push_back({t, g, -c});
                out.edges.push_back({t, h, c});
            }
            long d_gg = -co_g - ci_g - lp;
            long d_gh = -co_h + ci_g;
            long d_hg = co_g - ci_h;
            long d_hh = co_h + ci_h + lp;
            if (d_gg != 0) out.edges.push_back({g, g, d_gg});
            if (d_gh != 0) out.edges.push_back({g, h, d_gh});
            if (d_hg != 0) out.edges.push_back({h, g, d_hg});
            if (d_hh != 0) out.edges.push_back({h, h, d_hh});
        } else {
            auto& co = sc.cnt_out;
            co.clear();
            long lp = 0;
            for (const auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) lp = val;
                else co[_b[y]] += val;
            }
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k);
                return it == mp.end() ? 0L : it->second;
            };
            long cg = get_c(co, g), ch = get_c(co, h);
            for (const auto& [t, c] : co) {
                if (t == g || t == h) continue;
                int x1 = g < t ? g : t, y1 = g < t ? t : g;
                int x2 = h < t ? h : t, y2 = h < t ? t : h;
                out.edges.push_back({x1, y1, -c});
                out.edges.push_back({x2, y2, c});
            }
            long d_gg = -2 * cg - lp;
            long d_gh = cg - ch;
            long d_hh = 2 * ch + lp;
            if (d_gg != 0) out.edges.push_back({g, g, d_gg});
            if (d_gh != 0) {
                int a = g < h ? g : h, b = g < h ? h : g;
                out.edges.push_back({a, b, d_gh});
            }
            if (d_hh != 0) out.edges.push_back({h, h, d_hh});
        }
        sc.wd.clear();
        sc.wd.emplace_back(g, -1);
        sc.wd.emplace_back(h, +1);
        return process_eval(0, out, sc);
    }

    double group_merge_eval(int r, int s, ChangeSet& out, MergeScratch& sc) const {
        if (r == s) return 0;
        out.clear();
        if (_directed) {
            long m_rr = 0, m_rs = 0, m_sr = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) m_rr = m;
                else if (t == s) m_rs = m;
                else {
                    out.edges.push_back({r, t, -m});
                    out.edges.push_back({s, t, m});
                }
            }
            for (const auto& [t, m] : _m_in[r]) {
                if (t == r) continue;
                else if (t == s) m_sr = m;
                else {
                    out.edges.push_back({t, r, -m});
                    out.edges.push_back({t, s, m});
                }
            }
            if (m_rr != 0) out.edges.push_back({r, r, -m_rr});
            if (m_rs != 0) out.edges.push_back({r, s, -m_rs});
            if (m_sr != 0) out.edges.push_back({s, r, -m_sr});
            long d_ss = m_rr + m_rs + m_sr;
            if (d_ss != 0) out.edges.push_back({s, s, d_ss});
        } else {
            long m_rr = 0, m_rs = 0;
            for (const auto& [t, m] : _m_out[r]) {
                if (t == r) m_rr = m;
                else if (t == s) m_rs = m;
                else {
                    int x1 = r < t ? r : t, y1 = r < t ? t : r;
                    int x2 = s < t ? s : t, y2 = s < t ? t : s;
                    out.edges.push_back({x1, y1, -m});
                    out.edges.push_back({x2, y2, m});
                }
            }
            if (m_rr != 0) out.edges.push_back({r, r, -m_rr});
            if (m_rs != 0) {
                int x = r < s ? r : s, y = r < s ? s : r;
                out.edges.push_back({x, y, -m_rs});
            }
            long d_ss = m_rr + 2 * m_rs;
            if (d_ss != 0) out.edges.push_back({s, s, d_ss});
        }
        sc.wd.clear();
        sc.wd.emplace_back(r, -_wr[r]);
        sc.wd.emplace_back(s, _wr[r]);
        return process_eval(0, out, sc);
    }

    template <class Rng>
    int random_neighbor(int x, const LowerView& lower, Rng& rng) const {
        double total = (*lower.ep)[x] + (_directed ? (*lower.em)[x] : 0);
        if (total <= 0)
            return -1;
        double u = rng.unif() * total;
        for (const auto& [y, val] : (*lower.m_out)[x]) {
            u -= val;
            if (u < 0)
                return y;
        }
        if (_directed) {
            for (const auto& [y, val] : (*lower.m_in)[x]) {
                u -= val;
                if (u < 0)
                    return y;
            }
        }
        return x;
    }

    long m_sym(int x, int y) const {
        return get_m(x, y) + (_directed ? get_m(y, x) : 0);
    }

    long e_tot(int r) const { return _ep[r] + (_directed ? _em[r] : 0); }

    int sample_group_neighbor(int t, RRng& rng) const {
        double total = static_cast<double>(e_tot(t));
        if (total <= 0)
            return random_occupied(rng);
        double u = rng.unif() * total;
        for (const auto& [s, m] : _m_out[t]) {
            u -= m;
            if (u < 0)
                return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in[t]) {
                u -= m;
                if (u < 0)
                    return s;
            }
        }
        return t;
    }

    int sample_exact_group(int x, double c, const LowerView& lower, RRng& rng) const {
        double kx = static_cast<double>((*lower.ep)[x] + (_directed ? (*lower.em)[x] : 0));
        if (std::isinf(c) || kx <= 0)
            return random_occupied(rng);
        int u = random_neighbor(x, lower, rng);
        if (u < 0)
            return random_occupied(rng);
        int t = _b[u];
        double et = static_cast<double>(e_tot(t));
        double B_occ = static_cast<double>(_occupied.size());
        double p_rand = (c > 0) ? (c * B_occ / (et + c * B_occ)) : 0.0;
        if (et <= 0 || (c > 0 && rng.unif() < p_rand))
            return random_occupied(rng);
        return sample_group_neighbor(t, rng);
    }

    double move_prob(int x, int r, int s, double c, const LowerView& lower, bool reverse) {
        double kx = static_cast<double>((*lower.ep)[x] + (_directed ? (*lower.em)[x] : 0));
        if (!reverse) {
            double B_occ = static_cast<double>(_occupied.size());
            if (kx <= 0 || std::isinf(c))
                return 1.0 / B_occ;
            _touched_out.clear();
            auto add_w = [&](int t, long w) {
                if (_cnt_out[t] == 0)
                    _touched_out.push_back(t);
                _cnt_out[t] += w;
            };
            for (const auto& [y, val] : (*lower.m_out)[x])
                add_w(y == x ? r : _b[y], val);
            if (_directed) {
                for (const auto& [y, val] : (*lower.m_in)[x])
                    add_w(y == x ? r : _b[y], val);
            }
            double prob = 0;
            for (int t : _touched_out) {
                double wt = static_cast<double>(_cnt_out[t]);
                _cnt_out[t] = 0;
                double et = static_cast<double>(e_tot(t));
                double mts = static_cast<double>(m_sym(t, s));
                prob += (wt / kx) * ((mts + c) / (et + c * B_occ));
            }
            return prob;
        }
        double B_new = static_cast<double>(_occupied.size() + (_wr[s] == 0 ? 1 : 0));
        if (kx <= 0 || std::isinf(c))
            return 1.0 / B_new;
        _touched_out.clear();
        auto add_w = [&](int t, long w) {
            if (_cnt_out[t] == 0)
                _touched_out.push_back(t);
            _cnt_out[t] += w;
        };
        long lp_w = 0;
        for (const auto& [y, val] : (*lower.m_out)[x]) {
            if (y == x)
                lp_w += val;
            else
                add_w(_b[y], val);
        }
        if (_directed) {
            for (const auto& [y, val] : (*lower.m_in)[x]) {
                if (y == x)
                    lp_w += val;
                else
                    add_w(_b[y], val);
            }
        }
        long wnl_r = _cnt_out[r];
        long wnl_s = _cnt_out[s];
        if (lp_w > 0 && _cnt_out[s] == 0)
            _touched_out.push_back(s);
        double prob = 0;
        for (int t : _touched_out) {
            long wnl = _cnt_out[t];
            _cnt_out[t] = 0;
            double wt = static_cast<double>(wnl + (t == s ? lp_w : 0));
            if (wt <= 0)
                continue;
            double et_new, mtr_new;
            if (t == r) {
                et_new = static_cast<double>(e_tot(r) - static_cast<long>(kx));
                mtr_new = static_cast<double>(m_sym(r, r) - 2 * wnl_r - lp_w);
            } else if (t == s) {
                et_new = static_cast<double>(e_tot(s) + static_cast<long>(kx));
                mtr_new = static_cast<double>(m_sym(s, r) + wnl_r - wnl_s);
            } else {
                et_new = static_cast<double>(e_tot(t));
                mtr_new = static_cast<double>(m_sym(t, r) - wnl);
            }
            prob += (wt / kx) * ((mtr_new + c) / (et_new + c * B_new));
        }
        return prob;
    }

    int sample_merge_target(int r, double c, RRng& rng) const {
        size_t B_occ = _occupied.size();
        if (B_occ <= 1)
            return r;
        auto rand_other = [&]() {
            int s = _occupied[rng.index(B_occ - 1)];
            return (s == r) ? _occupied[B_occ - 1] : s;
        };
        double ext = static_cast<double>(e_tot(r) - m_sym(r, r));
        if (ext <= 0 || std::isinf(c))
            return rand_other();
        double B_other = static_cast<double>(B_occ - 1);
        double p_rand = (c > 0) ? (c * B_other / (ext + c * B_other)) : 0.0;
        if (c > 0 && rng.unif() < p_rand)
            return rand_other();
        double u = rng.unif() * ext;
        for (const auto& [s, m] : _m_out[r]) {
            if (s == r)
                continue;
            u -= m;
            if (u < 0)
                return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in[r]) {
                if (s == r)
                    continue;
                u -= m;
                if (u < 0)
                    return s;
            }
        }
        return rand_other();
    }

    double merge_prob(int r, int s, double c) const {
        size_t B_occ = _occupied.size();
        if (B_occ <= 1 || r == s)
            return 0;
        double B_other = static_cast<double>(B_occ - 1);
        double mrs = static_cast<double>(m_sym(r, s));
        auto q_dir = [&](int a) {
            double ext = static_cast<double>(e_tot(a) - m_sym(a, a));
            if (ext <= 0 || std::isinf(c))
                return 1.0 / B_other;
            return (mrs + c) / (ext + c * B_other);
        };
        return (q_dir(r) + q_dir(s)) / static_cast<double>(B_occ);
    }

    template <class Rng>
    int sample_group(int x, double c, const LowerView& lower, Rng& rng) const {
        int u = std::isinf(c) ? -1 : random_neighbor(x, lower, rng);
        if (u < 0)
            return random_occupied(rng);
        int t = _b[u];
        double p_rand = 0;
        if (c > 0) {
            double B = _occupied.size();
            double et = _ep[t] + (_directed ? _em[t] : 0);
            p_rand = c * B / (et + c * B);
        }
        if (c == 0 || rng.unif() >= p_rand) {
            const auto& vs = _members[t];
            int w = vs[rng.index(vs.size())];
            int z = random_neighbor(w, lower, rng);
            return _b[z < 0 ? w : z];
        }
        return random_occupied(rng);
    }

private:
    int _cap;
    bool _directed, _top;
    long _E;
    int _N = 0;
    std::vector<int> _b, _wr, _pos, _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<long> _ep, _em;
    std::vector<gt_hash_map<int, long>> _m_out, _m_in;
    std::vector<std::vector<int>> _members;
    DeltaAcc _acc, _pairs;
    std::vector<std::pair<int, int>> _wd;
    std::vector<long> _cnt_out, _cnt_in;
    std::vector<int> _touched_out, _touched_in;

    void add_weight(int g, int dw) {
        for (auto& p : _wd) {
            if (p.first == g) {
                p.second += dw;
                return;
            }
        }
        _wd.emplace_back(g, dw);
    }

    int new_weight(int g) const {
        int w = _wr[g];
        for (const auto& p : _wd)
            if (p.first == g)
                w += p.second;
        return w;
    }

    double process_eval(int dN, ChangeSet& out, MergeScratch& sc) const {
        bool any_wd = false;
        long B0 = B(), B1 = B0;
        for (const auto& [g, dw] : sc.wd) {
            if (dw == 0)
                continue;
            any_wd = true;
            int w0 = _wr[g], w1 = w0 + dw;
            if (w0 == 0 && w1 > 0) {
                ++B1;
                out.weights.emplace_back(g, +1);
            } else if (w0 > 0 && w1 == 0) {
                --B1;
                out.weights.emplace_back(g, -1);
            }
        }
        if (dN == 0 && !any_wd && out.edges.empty())
            return 0;

        double dS = 0;
        if (dN != 0 || B1 != B0) {
            dS += partition_term(_N + dN, B1) - partition_term(_N, B0);
            dS += edges_term(B1) - edges_term(B0);
        }
        for (const auto& [g, dw] : sc.wd) {
            if (dw != 0)
                dS += lgam_int(_wr[g] + 1) - lgam_int(_wr[g] + dw + 1);
        }
        sc.pairs.clear();
        for (const Delta& e : out.edges)
            sc.pairs.add(e.x, e.y, e.d);
        for (const auto& [g, dw] : sc.wd) {
            if (dw == 0)
                continue;
            for (const auto& [t, m] : _m_out[g])
                sc.pairs.add(g, t, 0);
            if (_directed)
                for (const auto& [t, m] : _m_in[g])
                    sc.pairs.add(t, g, 0);
        }
        auto new_w = [&](int g) {
            int w = _wr[g];
            for (const auto& p : sc.wd)
                if (p.first == g) w += p.second;
            return w;
        };
        for (const Delta& p : sc.pairs.items()) {
            long m0 = get_m(p.x, p.y);
            dS += dense_term(p.x, p.y, m0 + p.d, new_w(p.x), new_w(p.y)) -
                  dense_term(p.x, p.y, m0, _wr[p.x], _wr[p.y]);
        }
        return dS;
    }

    double process(int dN, ChangeSet& out, bool apply) {
        bool any_wd = false;
        long B0 = B(), B1 = B0;
        for (const auto& [g, dw] : _wd) {
            if (dw == 0)
                continue;
            any_wd = true;
            int w0 = _wr[g], w1 = w0 + dw;
            if (w0 == 0 && w1 > 0) {
                ++B1;
                out.weights.emplace_back(g, +1);
            } else if (w0 > 0 && w1 == 0) {
                --B1;
                out.weights.emplace_back(g, -1);
            }
        }
        if (dN == 0 && !any_wd && out.edges.empty())
            return 0;

        double dS = 0;
        if (!apply) {
            if (dN != 0 || B1 != B0) {
                dS += partition_term(_N + dN, B1) - partition_term(_N, B0);
                dS += edges_term(B1) - edges_term(B0);
            }
            for (const auto& [g, dw] : _wd) {
                if (dw != 0)
                    dS += lgam_int(_wr[g] + 1) - lgam_int(_wr[g] + dw + 1);
            }
            _pairs.clear();
            for (const Delta& e : out.edges)
                _pairs.add(e.x, e.y, e.d);
            for (const auto& [g, dw] : _wd) {
                if (dw == 0)
                    continue;
                for (const auto& [t, m] : _m_out[g])
                    _pairs.add(g, t, 0);
                if (_directed)
                    for (const auto& [t, m] : _m_in[g])
                        _pairs.add(t, g, 0);
            }
            for (const Delta& p : _pairs.items()) {
                long m0 = get_m(p.x, p.y);
                dS += dense_term(p.x, p.y, m0 + p.d, new_weight(p.x), new_weight(p.y)) -
                      dense_term(p.x, p.y, m0, _wr[p.x], _wr[p.y]);
            }
            return dS;
        }

        for (const Delta& e : out.edges) {
            add_m(_m_out[e.x], e.y, e.d);
            if (_directed) {
                add_m(_m_in[e.y], e.x, e.d);
                _ep[e.x] += e.d;
                _em[e.y] += e.d;
            } else {
                if (e.x != e.y)
                    add_m(_m_out[e.y], e.x, e.d);
                _ep[e.x] += e.d;
                if (e.x != e.y)
                    _ep[e.y] += e.d;
            }
        }
        _N += dN;
        return 0;
    }

    static void add_m(gt_hash_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) {
            row.emplace(y, d);
        } else {
            it->second += d;
            if (it->second == 0)
                row.erase(it);
        }
    }

    void add_member(int x, int r) {
        if (_wr[r] == 0) {
            _occ_pos[r] = static_cast<int>(_occupied.size());
            _occupied.push_back(r);
        }
        _pos[x] = static_cast<int>(_members[r].size());
        _members[r].push_back(x);
        _wr[r]++;
    }

    void remove_occupied(int r) {
        int i = _occ_pos[r];
        int other = _occupied.back();
        _occupied[i] = other;
        _occ_pos[other] = i;
        _occupied.pop_back();
        _occ_pos[r] = -1;
    }

    void push_empty(int r) {
        _empty_pos[r] = static_cast<int>(_empty.size());
        _empty.push_back(r);
    }

    void take_empty(int s) {
        int i = _empty_pos[s];
        if (i >= 0) {
            int last = _empty.back();
            _empty[i] = last;
            _empty_pos[last] = i;
            _empty.pop_back();
            _empty_pos[s] = -1;
        }
    }

    void push_node(int x, int r) {
        if (_wr[r] == 0)
            take_empty(r);
        add_member(x, r);
    }

    void pop_node(int x, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[x]] = last;
        _pos[last] = _pos[x];
        vs.pop_back();
        if (--_wr[r] == 0) {
            remove_occupied(r);
            push_empty(r);
        }
    }
};

class Nested {
public:
    typedef std::vector<std::vector<int>> Hierarchy;

    // Compact snapshot: level 0 has length N; upper levels store only (node, group)
    // pairs for occupied lower-level groups, reducing cache memory by ~15-20x on
    // large networks (N >= 10^4).
    struct CompactHierarchy {
        std::vector<int> b0;
        std::vector<std::vector<std::pair<int, int>>> upper;
    };

    Nested(const Graph& graph, bool dc, bool dist, size_t q_cache_max, int levels)
        : _base(graph, dc, dist, q_cache_max, levels <= 1) {
        int max_n = std::max(_base.N(), 1);
        _up.reserve(levels > 1 ? levels - 1 : 0);
        for (int l = 1; l < levels; ++l) {
            // Start compact and grow dynamically on demand via ensure_cap().
            int init_cap = (l == 1) ? std::min(max_n, 1024) : std::min(max_n, 64);
            _up.emplace_back(std::max(init_cap, 1), _base.directed(), _base.E(), l == levels - 1);
        }
        _cs.resize(levels + 1);
    }

    int levels() const { return static_cast<int>(_up.size()) + 1; }
    bool directed() const { return _base.directed(); }
    SBM& base() { return _base; }
    const SBM& base() const { return _base; }
    BlockLevel& upper(int l) { return _up[l - 1]; }
    const BlockLevel& upper(int l) const { return _up[l - 1]; }
    LowerView lower(int l) const {
        if (l == 1)
            return {&_base.m_out(), &_base.m_in(), &_base.ep(), &_base.em(), &_base.wr(),
                    &_base.occupied()};
        return _up[l - 2].view();
    }

    void build(const Hierarchy& bs) {
        _base.set_partition(bs[0]);
        for (size_t l = 1; l < bs.size(); ++l)
            _up[l - 1].build(bs[l], lower(l));
    }

    CompactHierarchy save_compact() const {
        CompactHierarchy snap;
        snap.b0 = _base.blocks();
        snap.upper.resize(_up.size());
        for (size_t l = 0; l < _up.size(); ++l) {
            const auto& locc = (l == 0) ? _base.occupied() : _up[l - 1].occupied();
            auto& dst = snap.upper[l];
            dst.reserve(locc.size());
            const auto& bl = _up[l].blocks();
            for (int x : locc)
                dst.emplace_back(x, bl[x]);
        }
        return snap;
    }

    void restore_compact(const CompactHierarchy& snap) {
        _base.set_partition(snap.b0);
        for (size_t l = 0; l < _up.size(); ++l)
            _up[l].build_sparse(snap.upper[l], lower(static_cast<int>(l) + 1));
    }

    void rebuild_from(int l, const std::vector<int>& bl,
                      const std::vector<int>* b_next = nullptr) {
        if (l == 0)
            _base.set_partition(bl);
        else
            _up[l - 1].build(bl, lower(l));
        for (int k = l + 1; k < levels(); ++k) {
            if (k == l + 1 && b_next != nullptr)
                _up[k - 1].build(*b_next, lower(k));
            else
                _up[k - 1].build(_up[k - 1].blocks(), lower(k));
        }
    }

    Hierarchy hierarchy() const {
        Hierarchy bs;
        bs.reserve(levels());
        bs.push_back(_base.blocks());
        for (const auto& u : _up)
            bs.push_back(u.blocks());
        return bs;
    }

    double entropy() const {
        double S = _base.entropy();
        for (const auto& u : _up)
            S += u.entropy();
        return S;
    }

    std::vector<double> level_entropies() const {
        std::vector<double> S{_base.entropy()};
        for (const auto& u : _up)
            S.push_back(u.entropy());
        return S;
    }

    int group(int l, int v) const { return l == 0 ? _base.group(v) : _up[l - 1].group(v); }
    int size(int l, int r) const { return l == 0 ? _base.size(r) : _up[l - 1].size(r); }

    double virtual_move(int l, int v, int s) { return change(l, v, s, false); }
    void move(int l, int v, int s) { change(l, v, s, true); }

    double virtual_merge(int l, int r, int s) { return merge_groups(l, r, s, false); }
    void merge(int l, int r, int s) { merge_groups(l, r, s, true); }

    // Thread-safe read-only virtual_move and virtual_merge using per-thread scratch.
    double virtual_move_const(int l, int v, int s, MergeScratch& sc) const {
        int r = group(l, v);
        if (r == s) return 0;
        double dS = 0;
        ChangeSet& cs = sc.cs[l];
        if (l == 0) {
            bool empties = _base.size(r) == 1, fills = _base.size(s) == 0;
            dS += _base.virtual_move_const(v, s, sc);
            if (levels() == 1) return dS;
            if (!empties && !fills && _up[0].group(r) == _up[0].group(s))
                return dS;
            cs.clear();
            cs.edges = sc.deltas;
            if (empties) cs.weights.emplace_back(r, -1);
            if (fills) cs.weights.emplace_back(s, +1);
        } else {
            bool empties = _up[l - 1].size(r) == 1, fills = _up[l - 1].size(s) == 0;
            dS += _up[l - 1].node_move_eval(v, s, lower(l), cs, sc);
            if (l + 1 < levels() && !empties && !fills &&
                _up[l].group(r) == _up[l].group(s))
                return dS;
        }
        for (int k = l + 1; k < levels(); ++k) {
            if (sc.cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes_eval(sc.cs[k - 1], sc.cs[k], sc);
        }
        return dS;
    }

    double virtual_merge_const(int l, int r, int s, MergeScratch& sc) const {
        if (r == s) return 0;
        double dS = 0;
        ChangeSet& cs = sc.cs[l];
        if (l == 0) {
            cs.clear();
            dS += _base.virtual_merge_const(r, s, sc);
            if (levels() == 1) return dS;
            cs.edges = sc.deltas;
            cs.weights.emplace_back(r, -1);
        } else {
            dS += _up[l - 1].group_merge_eval(r, s, cs, sc);
        }
        for (int k = l + 1; k < levels(); ++k) {
            if (sc.cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes_eval(sc.cs[k - 1], sc.cs[k], sc);
        }
        return dS;
    }

    // A new group t at level l takes the upper-level label of group r.
    void inherit(int l, int t, int r) {
        if (l + 1 < levels())
            _up[l].set_empty_label(t, _up[l].group(r));
    }

    bool same_parent(int l, int r, int s) const {
        if (l + 1 >= levels())
            return true;
        return _up[l].group(r) == _up[l].group(s);
    }

    int empty_group(int l) const {
        return l == 0 ? _base.empty_group() : _up[l - 1].empty_group();
    }

private:
    SBM _base;
    std::vector<BlockLevel> _up;
    std::vector<ChangeSet> _cs;

    double change(int l, int v, int s, bool apply) {
        int r = group(l, v);
        if (r == s)
            return 0;
        double dS = 0;
        ChangeSet& cs = _cs[l];
        if (l == 0) {
            bool empties = _base.size(r) == 1, fills = _base.size(s) == 0;
            if (apply)
                _base.move(v, s);
            else
                dS += _base.virtual_move(v, s);
            if (levels() == 1)
                return dS;
            // Fast path: if r and s share the same upper-level group and neither
            // empties nor fills, all upper-level edge deltas cancel identically.
            if (!empties && !fills && _up[0].group(r) == _up[0].group(s))
                return dS;
            cs.clear();
            cs.edges = _base.deltas();
            if (empties)
                cs.weights.emplace_back(r, -1);
            if (fills)
                cs.weights.emplace_back(s, +1);
        } else {
            bool empties = _up[l - 1].size(r) == 1, fills = _up[l - 1].size(s) == 0;
            dS += _up[l - 1].node_move(v, s, lower(l), cs, apply);
            if (l + 1 < levels() && !empties && !fills &&
                _up[l].group(r) == _up[l].group(s))
                return dS;
        }
        for (int k = l + 1; k < levels(); ++k) {
            if (_cs[k - 1].empty())
                break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k], apply);
        }
        return dS;
    }

    double merge_groups(int l, int r, int s, bool apply) {
        if (r == s)
            return 0;
        double dS = 0;
        ChangeSet& cs = _cs[l];
        if (l == 0) {
            cs.clear();
            if (apply)
                _base.merge(r, s);
            else
                dS += _base.virtual_merge(r, s);
            if (levels() == 1)
                return dS;
            cs.edges = _base.deltas();
            cs.weights.emplace_back(r, -1);
        } else {
            dS += _up[l - 1].group_merge(r, s, cs, apply);
        }
        for (int k = l + 1; k < levels(); ++k) {
            if (_cs[k - 1].empty())
                break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k], apply);
        }
        return dS;
    }
};

class LevelView {
public:
    typedef Nested::CompactHierarchy Snapshot;
    LevelView(Nested& n, int l) : _n(n), _l(l) {
        if (l == 0) {
            _nodes.resize(n.base().N());
            std::iota(_nodes.begin(), _nodes.end(), 0);
        } else {
            _nodes = *n.lower(l).occupied;
        }
    }
    bool directed() const { return _n.directed(); }
    int num_levels() const { return _n.levels(); }
    size_t N() const { return _nodes.size(); }
    const std::vector<int>& nodes() const { return _nodes; }
    int B() const { return _l == 0 ? _n.base().B() : _n.upper(_l).B(); }
    int group(int v) const { return _n.group(_l, v); }
    int size(int r) const { return _n.size(_l, r); }
    const std::vector<int>& occupied() const {
        return _l == 0 ? _n.base().occupied() : _n.upper(_l).occupied();
    }
    const std::vector<int>& members(int r) const {
        return _l == 0 ? _n.base().members(r) : _n.upper(_l).members(r);
    }
    double virtual_move(int v, int s) { return _n.virtual_move(_l, v, s); }
    double virtual_move_const(int v, int s, MergeScratch& sc) const {
        return _n.virtual_move_const(_l, v, s, sc);
    }
    void move(int v, int s) { _n.move(_l, v, s); }
    double virtual_merge(int r, int s) { return _n.virtual_merge(_l, r, s); }
    double virtual_merge_const(int r, int s, MergeScratch& sc) const {
        return _n.virtual_merge_const(_l, r, s, sc);
    }
    void merge(int r, int s) { _n.merge(_l, r, s); }
    template <class Rng>
    int sample_group(int v, double c, Rng& rng) const {
        if (_l == 0)
            return _n.base().sample_group(v, c, rng);
        return _n.upper(_l).sample_group(v, c, _n.lower(_l), rng);
    }
    int empty_group() const { return _n.empty_group(_l); }
    void inherit(int t, int r) { _n.inherit(_l, t, r); }
    bool same_parent(int r, int s) const { return _n.same_parent(_l, r, s); }
    int sample_exact_group(int v, double c, RRng& rng) const {
        if (_l == 0)
            return _n.base().sample_exact_group(v, c, rng);
        return _n.upper(_l).sample_exact_group(v, c, _n.lower(_l), rng);
    }
    double move_prob(int v, int r, int s, double c, bool reverse) {
        if (_l == 0)
            return _n.base().move_prob(v, r, s, c, reverse);
        return _n.upper(_l).move_prob(v, r, s, c, _n.lower(_l), reverse);
    }
    int sample_merge_target(int r, double c, RRng& rng) const {
        if (_l == 0)
            return _n.base().sample_merge_target(r, c, rng);
        return _n.upper(_l).sample_merge_target(r, c, rng);
    }
    double merge_prob(int r, int s, double c) const {
        if (_l == 0)
            return _n.base().merge_prob(r, s, c);
        return _n.upper(_l).merge_prob(r, s, c);
    }
    template <class Rng>
    int random_occupied(Rng& rng) const {
        return _l == 0 ? _n.base().random_occupied(rng) : _n.upper(_l).random_occupied(rng);
    }
    double entropy() const { return _n.entropy(); }
    Snapshot save() const { return _n.save_compact(); }
    void restore(const Snapshot& bs) { _n.restore_compact(bs); }
    void spread() {
        int max_v = 0;
        for (int v : _nodes) if (v > max_v) max_v = v;
        int sz = (_l == 0) ? _n.base().N() : (max_v + 1);
        std::vector<int> bl(std::max(sz, 1), 0);
        std::vector<int> b_next;
        bool has_upper = (_l + 1 < _n.levels());
        if (has_upper) {
            b_next = _n.upper(_l + 1).blocks();
            const std::vector<int> old_next = b_next;
            if (static_cast<int>(b_next.size()) <= max_v)
                b_next.resize(max_v + 1, 0);
            for (int v : _nodes) {
                int r = group(v);
                bl[v] = v;
                b_next[v] = (r < static_cast<int>(old_next.size())) ? old_next[r] : 0;
            }
            _n.rebuild_from(_l, bl, &b_next);
        } else {
            for (int v : _nodes)
                bl[v] = v;
            _n.rebuild_from(_l, bl, nullptr);
        }
    }
    void collapse(int t) {
        int max_v = 0;
        for (int v : _nodes) if (v > max_v) max_v = v;
        int sz = (_l == 0) ? _n.base().N() : (max_v + 1);
        std::vector<int> bl(std::max(sz, 1), t);
        _n.rebuild_from(_l, bl, nullptr);
    }
private:
    Nested& _n;
    int _l;
    std::vector<int> _nodes;
};

// ---------------------------------------------------------------------------
// Views used by the multilevel search. A view exposes the nodes of one level,
// their groups, moves, and whole-model snapshots.

class FlatView {
public:
    typedef std::vector<int> Snapshot;
    explicit FlatView(SBM& s) : _s(s), _nodes(s.N()) {
        std::iota(_nodes.begin(), _nodes.end(), 0);
    }
    bool directed() const { return _s.directed(); }
    int num_levels() const { return 1; }
    size_t N() const { return _nodes.size(); }
    const std::vector<int>& nodes() const { return _nodes; }
    int B() const { return _s.B(); }
    int group(int v) const { return _s.group(v); }
    int size(int r) const { return _s.size(r); }
    const std::vector<int>& occupied() const { return _s.occupied(); }
    const std::vector<int>& members(int r) const { return _s.members(r); }
    double virtual_move(int v, int s) { return _s.virtual_move(v, s); }
    double virtual_move_const(int v, int s, MergeScratch& sc) const {
        return _s.virtual_move_const(v, s, sc);
    }
    void move(int v, int s) { _s.move(v, s); }
    double virtual_merge(int r, int s) { return _s.virtual_merge(r, s); }
    double virtual_merge_const(int r, int s, MergeScratch& sc) const {
        return _s.virtual_merge_const(r, s, sc);
    }
    void merge(int r, int s) { _s.merge(r, s); }
    template <class Rng>
    int sample_group(int v, double c, Rng& rng) const { return _s.sample_group(v, c, rng); }
    int empty_group() const { return _s.empty_group(); }
    void inherit(int, int) {}
    bool same_parent(int, int) const { return true; }
    int sample_exact_group(int v, double c, RRng& rng) const {
        return _s.sample_exact_group(v, c, rng);
    }
    double move_prob(int v, int r, int s, double c, bool reverse) {
        return _s.move_prob(v, r, s, c, reverse);
    }
    int sample_merge_target(int r, double c, RRng& rng) const {
        return _s.sample_merge_target(r, c, rng);
    }
    double merge_prob(int r, int s, double c) const { return _s.merge_prob(r, s, c); }
    template <class Rng>
    int random_occupied(Rng& rng) const { return _s.random_occupied(rng); }
    double entropy() const { return _s.entropy(); }
    Snapshot save() const { return _s.blocks(); }
    void restore(const Snapshot& b) { _s.set_partition(b); }
    void spread() {
        std::vector<int> b(_nodes.size());
        std::iota(b.begin(), b.end(), 0);
        _s.set_partition(b);
    }
    void collapse(int t) { _s.set_partition(std::vector<int>(_nodes.size(), t)); }
private:
    SBM& _s;
    std::vector<int> _nodes;
};

// ---------------------------------------------------------------------------
// Multilevel agglomerative search

struct Params {
    double r = 0.9;
    size_t merge_sweeps = 10;
    size_t mh_sweeps = 10;
    double init_r = 0.99;
    size_t init_min_iter = 8;
    size_t init_min_iter_beta = 4;
    double init_beta = 1.0;
    double c = 0.5;
    size_t B_min = 1;
    size_t B_max = std::numeric_limits<size_t>::max();
    double epsilon = 1e-3;
    size_t max_refine = 1000;
    bool bisection = true;
    int nthreads = 1;
    bool verbose = false;
};

template <class State>
class Multilevel {
public:
    Multilevel(State& state, const Params& p) : _s(state), _p(p) {}

    bool metropolis(double dS, double beta) {
        if (std::isinf(beta))
            return dS < 0;
        double a = -beta * dS;
        return a > 0 || _rng.unif() < std::exp(a);
    }

    double mh_sweep(double beta, size_t B_min, size_t B_max) {
        size_t N = _s.N();
        size_t B = _s.B();
        if (B == 1 || (B == N && B_min == B))
            return 0;
        std::vector<int> vs = _s.nodes();
        _rng.shuffle(vs);
        double S = 0;

        int nw = rgt_parallel::resolve_threads(_p.nthreads);
        if (std::isinf(beta) && nw > 1 && N >= 2048) {
            // Parallel batch proposal filtering for zero-temperature sweeps:
            // read-only workers screen out non-improving moves across all cores
            // with chunk-indexed deterministic PRNGs, then the main thread
            // re-verifies and commits exact dS < 0 moves in deterministic order.
            size_t batch = std::min(N, std::max<size_t>(2048, static_cast<size_t>(nw) * 512));
            std::vector<MergeScratch> scratches(static_cast<size_t>(nw));
            for (int t = 0; t < nw; ++t)
                scratches[t].init(_s.directed(), _s.num_levels());
            std::vector<int> cand_s;

            for (size_t base = 0; base < N; base += batch) {
                size_t len = std::min(batch, N - base);
                cand_s.assign(len, -1);
                uint64_t batch_seed = static_cast<uint64_t>(_rng.index(1ULL << 31)) ^
                                      (static_cast<uint64_t>(_rng.index(1ULL << 31)) << 31) ^
                                      0x9e3779b97f4a7c15ULL;
                rgt_parallel::parallel_for(len, nw, [&](int tid, size_t lo, size_t hi) {
                    rgt_parallel::ThreadRng rng(batch_seed ^
                                                (static_cast<uint64_t>(lo + 1) * 0xbf58476d1ce4e5b9ULL));
                    auto& sc = scratches[tid];
                    for (size_t idx = lo; idx < hi; ++idx) {
                        int v = vs[base + idx];
                        int r = _s.group(v);
                        int s = _s.sample_group(v, 0, rng);
                        if (s == r) continue;
                        double dS = _s.virtual_move_const(v, s, sc);
                        if (dS < 0.0) cand_s[idx] = s;
                    }
                }, 64);

                for (size_t idx = 0; idx < len; ++idx) {
                    int s = cand_s[idx];
                    if (s < 0) continue;
                    int v = vs[base + idx];
                    int r = _s.group(v);
                    if (s == r) continue;
                    int dB = (_s.size(s) == 0 ? 1 : 0) - (_s.size(r) == 1 ? 1 : 0);
                    size_t Bc = _s.B();
                    if ((dB < 0 && Bc + dB < B_min) || (dB > 0 && Bc + dB > B_max))
                        continue;
                    double dS = _s.virtual_move(v, s);
                    if (dS < 0.0) {
                        _s.move(v, s);
                        S += dS;
                    }
                }
            }
            return S;
        }

        for (int v : vs) {
            int r = _s.group(v);
            int s = _s.sample_group(v, 0, _rng);
            if (s == r)
                continue;
            int dB = (_s.size(s) == 0 ? 1 : 0) - (_s.size(r) == 1 ? 1 : 0);
            size_t Bc = _s.B();
            if ((dB < 0 && Bc + dB < B_min) || (dB > 0 && Bc + dB > B_max))
                continue;
            double dS = _s.virtual_move(v, s);
            if (metropolis(dS, beta)) {
                _s.move(v, s);
                S += dS;
            }
        }
        return S;
    }

    double merge_sweep(size_t B_target, bool strict, size_t B_max) {
        const auto& occ = _s.occupied();
        size_t B_occ = occ.size();
        if (B_occ <= 1)
            return 0;
        std::vector<int> rlist = occ;
        std::vector<std::pair<int, int>> pairs;
        std::vector<double> dSs;
        pairs.reserve(B_occ);
        dSs.reserve(B_occ);

        int nw = rgt_parallel::resolve_threads(_p.nthreads);
        if (nw > 1 && B_occ >= 32) {
            // Lock-free parallel candidate merge evaluation with chunk-indexed
            // deterministic PRNGs (bit-for-bit reproducible for fixed nthreads).
            std::vector<int> best_tgt(B_occ, -1);
            std::vector<double> best_val(B_occ, inf);
            uint64_t sweep_seed = static_cast<uint64_t>(_rng.index(1ULL << 31)) ^
                                  (static_cast<uint64_t>(_rng.index(1ULL << 31)) << 31) ^
                                  0xbf58476d1ce4e5b9ULL;
            std::vector<MergeScratch> scratches(static_cast<size_t>(nw));
            std::vector<gt_hash_set<int>> pasts(static_cast<size_t>(nw));
            for (int t = 0; t < nw; ++t) {
                scratches[t].init(_s.directed(), _s.num_levels());
                pasts[t].reserve(_p.merge_sweeps * 2);
            }

            rgt_parallel::parallel_for(B_occ, nw, [&](int tid, size_t lo, size_t hi) {
                rgt_parallel::ThreadRng rng(sweep_seed ^
                                            (static_cast<uint64_t>(lo + 1) * 0x94d049bb133111ebULL));
                auto& sc = scratches[tid];
                auto& past = pasts[tid];
                for (size_t idx = lo; idx < hi; ++idx) {
                    int r = rlist[idx];
                    int best = r;
                    double best_dS = inf;
                    past.clear();
                    auto try_target = [&](int s) {
                        if (s != r && past.find(s) == past.end()) {
                            double dS = _s.virtual_merge_const(r, s, sc);
                            if (dS < best_dS) {
                                best = s;
                                best_dS = dS;
                            }
                            past.insert(s);
                        }
                    };
                    auto find_candidates = [&](double c_val) {
                        const auto& vs = _s.members(r);
                        for (size_t i = 0; i < _p.merge_sweeps; ++i) {
                            int v = vs[rng.index(vs.size())];
                            try_target(_s.sample_group(v, c_val, rng));
                        }
                    };
                    find_candidates(0);
                    if (best == r)
                        find_candidates(_p.c);
                    if (best == r) {
                        size_t tries = std::min(_p.merge_sweeps, B_occ - 1);
                        for (size_t i = 0; i < tries; ++i) {
                            int s = occ[rng.index(B_occ - 1)];
                            if (s == r)
                                s = occ[B_occ - 1];
                            try_target(s);
                        }
                    }
                    if (best != r) {
                        best_tgt[idx] = best;
                        best_val[idx] = best_dS;
                    }
                }
            }, 8);

            for (size_t idx = 0; idx < B_occ; ++idx) {
                if (best_tgt[idx] >= 0) {
                    pairs.emplace_back(rlist[idx], best_tgt[idx]);
                    dSs.push_back(best_val[idx]);
                }
            }
        } else {
            gt_hash_set<int> past;
            past.reserve(_p.merge_sweeps * 2);
            for (int r : rlist) {
                int best = r;
                double best_dS = inf;
                past.clear();
                auto try_target = [&](int s) {
                    if (s != r && past.find(s) == past.end()) {
                        double dS = _s.virtual_merge(r, s);
                        if (dS < best_dS) {
                            best = s;
                            best_dS = dS;
                        }
                        past.insert(s);
                    }
                };
                auto find_candidates = [&](double c_val) {
                    const auto& vs = _s.members(r);
                    for (size_t i = 0; i < _p.merge_sweeps; ++i) {
                        int v = vs[_rng.index(vs.size())];
                        try_target(_s.sample_group(v, c_val, _rng));
                    }
                };
                find_candidates(0);
                if (best == r)
                    find_candidates(_p.c);
                if (best == r) {
                    // Guaranteed fallback to uniform random other occupied group (c = inf).
                    size_t tries = std::min(_p.merge_sweeps, B_occ - 1);
                    for (size_t i = 0; i < tries; ++i) {
                        int s = occ[_rng.index(B_occ - 1)];
                        if (s == r)
                            s = occ[B_occ - 1];
                        try_target(s);
                    }
                }
                if (best != r) {
                    pairs.emplace_back(r, best);
                    dSs.push_back(best_dS);
                }
            }
        }

        std::vector<size_t> pis(pairs.size());
        std::iota(pis.begin(), pis.end(), 0);
        _rng.shuffle(pis);
        auto cmp = [&](size_t i, size_t j) { return dSs[i] > dSs[j]; };
        std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)>
            queue(cmp, std::move(pis));

        gt_hash_map<int, int> root;
        root.reserve(B_occ);
        auto get_root = [&](int r) {
            int x = r;
            while (true) {
                auto it = root.find(x);
                if (it == root.end() || it->second == x)
                    break;
                x = it->second;
            }
            root[r] = x;
            return x;
        };

        double S = 0;
        while (static_cast<size_t>(_s.B()) > B_target && !queue.empty()) {
            size_t i = queue.top();
            queue.pop();
            int s = pairs[i].first;
            int t = get_root(pairs[i].second);
            if (s == t || _s.size(s) == 0 || _s.size(t) == 0)
                continue;
            double ndS = _s.virtual_merge(s, t);
            if (!queue.empty() && ndS > dSs[queue.top()]) {
                dSs[i] = ndS;
                queue.push(i);
                continue;
            }
            if (ndS > 0 && strict && static_cast<size_t>(_s.B()) < B_max)
                break;
            root[s] = t;
            _s.merge(s, t);
            S += ndS;
        }
        return S;
    }

    double merge_and_sweep(bool strict, size_t B_min, size_t B_max) {
        size_t B = _s.B();
        size_t Bnext = std::min(std::max(std::min(B - 1, size_t(std::round(B * _p.r))), B_min),
                                B_max);
        double S = 0;
        size_t stalled = 0;
        while (static_cast<size_t>(_s.B()) > Bnext) {
            size_t before = _s.B();
            S += merge_sweep(Bnext, strict, B_max);
            if (strict && static_cast<size_t>(_s.B()) <= B_max)
                break;
            if (static_cast<size_t>(_s.B()) == before && ++stalled > 100)
                break;
            Rcpp::checkUserInterrupt();
        }
        for (size_t i = 0; i < _p.mh_sweeps; ++i) {
            double dS = mh_sweep(inf, Bnext, B_max);
            S += dS;
            if (std::abs(dS) < 1e-8)
                break;
        }
        return S;
    }

    double single_move_init(size_t B_min, size_t B_max) {
        size_t i = 0, Bprev;
        double S = 0;
        do {
            double beta = (i < _p.init_min_iter_beta) ? _p.init_beta : inf;
            Bprev = _s.B();
            S += mh_sweep(beta, B_min, B_max);
            ++i;
            Rcpp::checkUserInterrupt();
        } while ((i < _p.init_min_iter || _s.B() / double(Bprev) < _p.init_r) &&
                 static_cast<size_t>(_s.B()) > B_min);
        return S;
    }

    void put_cache(double S) {
        _cache[_s.B()] = {S, _s.save()};
        _trace_B.push_back(_s.B());
        _trace_S.push_back(S);
    }

    double get_cache(size_t B) {
        auto& c = _cache.at(B);
        _s.restore(c.second);
        return c.first;
    }

    double get_S(size_t B, size_t B_max) {
        auto it = _cache.lower_bound(B);
        if (it == _cache.end())
            Rcpp::stop("Internal error: no cached partition at or above B = %d", int(B));
        if (it->first == B)
            return it->second.first;
        get_cache(it->first);
        while (static_cast<size_t>(_s.B()) > B) {
            size_t before = _s.B();
            merge_and_sweep(false, B, B_max);
            if (static_cast<size_t>(_s.B()) == before)
                Rcpp::stop("Internal error: multilevel search could not merge groups");
            put_cache(_s.entropy());
            if (_p.verbose)
                Rprintf("    %d -> %d: %.4f\n", int(before), _s.B(), _cache[_s.B()].first);
        }
        return _cache[B].first;
    }

    size_t random_mid(size_t a, size_t b) {
        if (a == b || b - a < 2)
            return a;
        return a + 1 + _rng.index(b - a - 1);
    }

    // BisectionSearch::search with random midpoints (random_bisection = True)
    void bisection(size_t x_min, size_t x_max, size_t B_max) {
        auto f = [&](size_t B) { return get_S(B, B_max); };
        size_t x_mid = random_mid(x_min, x_max);
        double f_max = f(x_max);
        double f_min = f(x_min);
        double f_mid = f(x_mid);
        while (f_mid > f_min || f_mid > f_max) {
            if (f_min < f_max) {
                x_max = x_mid;
                f_max = f_mid;
                x_mid = random_mid(x_min, x_mid);
            } else {
                x_min = x_mid;
                f_min = f_mid;
                x_mid = random_mid(x_mid, x_max);
            }
            f_mid = f(x_mid);
            if (x_min == x_mid && (x_max - x_mid) <= 1)
                break;
        }
        while (x_max - x_mid > 1) {
            bool upper = x_max - x_mid > x_mid - x_min;
            size_t x = upper ? random_mid(x_mid, x_max) : random_mid(x_min, x_mid);
            double f_x = f(x);
            if (f_x < f_mid) {
                if (upper) {
                    x_min = x_mid;
                    f_min = f_mid;
                } else {
                    x_max = x_mid;
                    f_max = f_mid;
                }
                x_mid = x;
                f_mid = f_x;
            } else {
                if (upper) {
                    x_max = x;
                    f_max = f_x;
                } else {
                    x_min = x;
                    f_min = f_x;
                }
            }
        }
    }

    // Strictly agglomerative geometric coarsening without bisection backtracking.
    // Retains only a single best snapshot in memory (O(N) peak snapshot memory)
    // and shrinks B -> round(B * r) -> ... -> B_min in O(log(N) / log(1/r)) steps.
    double stage_agglomerative(size_t B_min, size_t B_max, bool strict_target = false) {
        size_t N = _s.N();
        double S_init = _s.entropy();
        if (N <= 1)
            return S_init;
        B_max = std::min(N, B_max);
        B_min = std::min(std::max<size_t>(1, B_min), B_max);
        _cache.clear();

        double best_S = inf;
        typename State::Snapshot best_snap;
        auto consider = [&](double S_cur) {
            size_t Bc = static_cast<size_t>(_s.B());
            _trace_B.push_back(static_cast<double>(Bc));
            _trace_S.push_back(S_cur);
            if (Bc >= B_min && Bc <= B_max && S_cur < best_S) {
                best_S = S_cur;
                best_snap = _s.save();
            }
        };

        if (B_min == 1 && !strict_target) {
            auto initial = _s.save();
            _s.collapse(_s.group(_s.nodes()[_rng.index(N)]));
            consider(_s.entropy());
            _s.restore(initial);
        }

        _s.spread();
        if (!strict_target)
            consider(_s.entropy());
        single_move_init(B_min, B_max);
        double S_cur = _s.entropy();
        if (!strict_target)
            consider(S_cur);
        if (_p.verbose)
            Rprintf("init (agglom): B = %d, S = %.4f\n", _s.B(), S_cur);

        while (static_cast<size_t>(_s.B()) > B_min) {
            size_t before = _s.B();
            merge_and_sweep(false, B_min, B_max);
            if (static_cast<size_t>(_s.B()) == before)
                break;
            S_cur = _s.entropy();
            if (!strict_target)
                consider(S_cur);
            else {
                _trace_B.push_back(static_cast<double>(_s.B()));
                _trace_S.push_back(S_cur);
            }
            if (_p.verbose)
                Rprintf("    %d -> %d: %.4f\n", int(before), _s.B(), S_cur);
        }

        if (strict_target)
            return _s.entropy();
        if (std::isfinite(best_S)) {
            _s.restore(best_snap);
            return best_S;
        }
        return _s.entropy();
    }

    // stage_multilevel() over all vertices at beta = infinity
    double stage_multilevel() {
        if (!_p.bisection)
            return stage_agglomerative(_p.B_min, _p.B_max, false);
        size_t N = _s.N();
        double S_init = _s.entropy();
        if (N <= 1)
            return S_init;
        size_t B_max = std::min(N, _p.B_max);
        size_t B_min = std::min(std::max<size_t>(1, _p.B_min), B_max);
        _cache.clear();
        auto initial = _s.save();
        if (B_min == 1) {
            _s.collapse(_s.group(_s.nodes()[_rng.index(N)]));
            put_cache(_s.entropy());
            _s.restore(initial);
        }
        _s.spread();
        put_cache(_s.entropy());
        single_move_init(B_min, B_max);
        B_max = std::min<size_t>(_s.B(), B_max);
        if (_cache.find(_s.B()) == _cache.end())
            put_cache(_s.entropy());
        if (_p.verbose)
            Rprintf("init: B = %d, S = %.4f\n", _s.B(), _s.entropy());
        size_t B_max_search = std::max(B_max, B_min);
        bisection(B_min, B_max_search, B_max_search);

        double best = inf;
        size_t best_B = 0;
        for (const auto& [B, c] : _cache) {
            if (B < B_min || B > B_max_search)
                continue;
            if (c.first < best) {
                best = c.first;
                best_B = B;
            }
        }
        return get_cache(best_B);
    }

    double refine(double S) {
        size_t iter = 0;
        while (iter++ < _p.max_refine) {
            double dS = merge_and_sweep(true, _p.B_min, _p.B_max);
            S = _s.entropy();
            if (_p.verbose)
                Rprintf("refine: B = %d, dS = %.6f, S = %.4f\n", _s.B(), dS, S);
            if (S == 0 || std::abs(dS / S) <= _p.epsilon)
                break;
        }
        return S;
    }

    std::vector<double> _trace_B, _trace_S;

private:
    State& _s;
    Params _p;
    RRng _rng;
    std::map<size_t, std::pair<double, typename State::Snapshot>> _cache;
};

Params read_params(const Rcpp::List& control) {
    Params p;
    auto num = [&](const char* name, double def) {
        return control.containsElementNamed(name) ? Rcpp::as<double>(control[name]) : def;
    };
    p.r = num("r", p.r);
    p.merge_sweeps = num("merge_sweeps", p.merge_sweeps);
    p.mh_sweeps = num("mh_sweeps", p.mh_sweeps);
    p.init_r = num("init_r", p.init_r);
    p.init_min_iter = num("init_min_iter", p.init_min_iter);
    p.init_min_iter_beta = num("init_min_iter_beta", p.init_min_iter_beta);
    p.init_beta = num("init_beta", p.init_beta);
    p.c = num("c", p.c);
    p.epsilon = num("epsilon", p.epsilon);
    p.max_refine = num("max_refine", p.max_refine);
    p.B_min = num("B_min", 1);
    double bmax = num("B_max", inf);
    p.B_max = std::isinf(bmax) ? std::numeric_limits<size_t>::max() : size_t(bmax);
    p.bisection = !control.containsElementNamed("bisection") || Rcpp::as<bool>(control["bisection"]);
    int nthr = control.containsElementNamed("nthreads") ? Rcpp::as<int>(control["nthreads"]) : 0;
    p.nthreads = rgt_parallel::resolve_threads(nthr);
    p.verbose = control.containsElementNamed("verbose") && Rcpp::as<bool>(control["verbose"]);
    return p;
}

std::vector<int> zero_based(const Graph& graph, const Rcpp::IntegerVector& blocks) {
    size_t N = boost::num_vertices(graph.adjacency);
    if (static_cast<size_t>(blocks.size()) != N)
        Rcpp::stop("blocks must have one entry per vertex");
    // relabel to 0..B-1 preserving first-appearance order
    gt_hash_map<int, int> map;
    std::vector<int> b(N);
    for (size_t v = 0; v < N; ++v) {
        auto it = map.find(blocks[v]);
        if (it == map.end())
            it = map.emplace(blocks[v], static_cast<int>(map.size())).first;
        b[v] = it->second;
    }
    return b;
}

Rcpp::IntegerVector one_based(const std::vector<int>& b) {
    gt_hash_map<int, int> map;
    Rcpp::IntegerVector out(b.size());
    for (size_t v = 0; v < b.size(); ++v) {
        auto it = map.find(b[v]);
        if (it == map.end())
            it = map.emplace(b[v], static_cast<int>(map.size()) + 1).first;
        out[v] = it->second;
    }
    return out;
}

size_t q_cache_max(double x) {
    return std::isinf(x) ? 10000 : static_cast<size_t>(x);
}

}  // namespace

// [[Rcpp::export]]
double cpp_sbm_entropy(SEXP pointer, Rcpp::IntegerVector blocks, bool degree_corrected,
                       bool distributed, double q_max) {
    Graph& graph = get_graph(pointer);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    state.set_partition(zero_based(graph, blocks));
    return state.entropy();
}

// Applies random single-vertex moves and group merges and returns, for each
// step, the difference between the incremental and the recomputed entropy change.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_sbm_check_moves(SEXP pointer, Rcpp::IntegerVector blocks,
                                        bool degree_corrected, bool distributed,
                                        int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    state.set_partition(zero_based(graph, blocks));
    RRng rng;
    Rcpp::NumericVector err(nmoves);
    if (state.N() == 0)
        return err;
    double S = state.entropy();
    for (int i = 0; i < nmoves; ++i) {
        if (rng.unif() < 0.15 && state.B() > 1) {
            const auto& occ = state.occupied();
            int r = occ[rng.index(occ.size())];
            int s = occ[rng.index(occ.size() - 1)];
            if (s == r)
                s = occ.back();
            double dS = state.virtual_merge(r, s);
            state.merge(r, s);
            double S_new = state.entropy();
            err[i] = dS - (S_new - S);
            S = S_new;
        } else {
            int v = rng.index(state.N());
            int s = (rng.unif() < 0.25 && state.empty_group() >= 0)
                        ? state.empty_group() : state.random_occupied(rng);
            double dS = state.virtual_move(v, s);
            state.move(v, s);
            double S_new = state.entropy();
            err[i] = dS - (S_new - S);
            S = S_new;
        }
    }
    return err;
}

// [[Rcpp::export]]
int cpp_max_threads() {
    return rgt_parallel::hardware_threads();
}

// [[Rcpp::export]]
int cpp_get_num_threads() {
    return rgt_parallel::get_num_threads();
}

// [[Rcpp::export]]
int cpp_set_num_threads(int nthreads) {
    return rgt_parallel::set_num_threads(nthreads);
}

// [[Rcpp::export]]
Rcpp::List cpp_sbm_minimize(SEXP pointer, bool degree_corrected, bool distributed,
                            Rcpp::List control, bool refine, double q_max,
                            Rcpp::Nullable<Rcpp::IntegerVector> init_blocks = R_NilValue) {
    Graph& graph = get_graph(pointer);
    Params p = read_params(control);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    FlatView view(state);
    Multilevel<FlatView> ml(view, p);
    double S = 0.0;
    double S_stage = 0.0;
    if (init_blocks.isNotNull()) {
        Rcpp::IntegerVector b0(init_blocks);
        state.set_partition(zero_based(graph, b0));
        if (static_cast<size_t>(state.B()) > p.B_max) {
            while (static_cast<size_t>(state.B()) > p.B_max) {
                size_t before = state.B();
                ml.merge_and_sweep(false, p.B_max, p.B_max);
                if (static_cast<size_t>(state.B()) == before) break;
            }
        }
        S = state.entropy();
        S_stage = S;
        ml._trace_B.push_back(static_cast<double>(state.B()));
        ml._trace_S.push_back(S);
    } else {
        std::vector<int> b(state.N(), 0);
        state.set_partition(b);
        S = ml.stage_multilevel();
        S_stage = S;
    }
    if (refine && state.N() > 1)
        S = ml.refine(S);
    return Rcpp::List::create(Rcpp::Named("blocks") = one_based(state.blocks()),
                              Rcpp::Named("entropy") = S,
                              Rcpp::Named("entropy_multilevel") = S_stage,
                              Rcpp::Named("trace") = Rcpp::DataFrame::create(
                                  Rcpp::Named("B") = ml._trace_B,
                                  Rcpp::Named("entropy") = ml._trace_S));
}

// [[Rcpp::export]]
Rcpp::List cpp_sbm_refine(SEXP pointer, Rcpp::IntegerVector blocks, bool degree_corrected,
                          bool distributed, Rcpp::List control, double q_max) {
    Graph& graph = get_graph(pointer);
    Params p = read_params(control);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    state.set_partition(zero_based(graph, blocks));
    FlatView view(state);
    Multilevel<FlatView> ml(view, p);
    double S = state.entropy();
    if (state.N() > 1)
        S = ml.refine(S);
    return Rcpp::List::create(Rcpp::Named("blocks") = one_based(state.blocks()),
                              Rcpp::Named("entropy") = S);
}

// [[Rcpp::export]]
double cpp_log_q(double n, double k, double q_max) {
    LogQ& q = log_q_cache();
    q.ensure(std::min<size_t>(q_cache_max(q_max), std::max<double>(n, 1)));
    return q(static_cast<long>(n), static_cast<long>(k));
}

// ---------------------------------------------------------------------------
// Nested model interface. Hierarchies cross the R boundary as lists of
// 1-based contiguous labels: element 1 labels vertices, element l labels the
// groups of level l - 1.

namespace {

Nested::Hierarchy internal_hierarchy(int N, const Rcpp::List& bs) {
    Nested::Hierarchy out;
    out.reserve(bs.size());
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        std::vector<int> b(l == 0 ? N : std::max<int>(static_cast<int>(x.size()), 1), 0);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            b[i] = x[i] - 1;
        out.push_back(std::move(b));
    }
    return out;
}

Rcpp::List external_hierarchy(const Nested::Hierarchy& bs) {
    Rcpp::List out(bs.size());
    std::vector<int> lower;  // internal labels of the previous level, in order
    for (size_t l = 0; l < bs.size(); ++l) {
        const std::vector<int>& nodes = (l == 0) ? bs[0] : bs[l];
        std::vector<int> items;
        if (l == 0) {
            items = bs[0];
        } else {
            for (int g : lower)
                items.push_back(nodes[g]);
        }
        gt_hash_map<int, int> map;
        std::vector<int> order;
        Rcpp::IntegerVector x(items.size());
        for (size_t i = 0; i < items.size(); ++i) {
            auto it = map.find(items[i]);
            if (it == map.end()) {
                it = map.emplace(items[i], static_cast<int>(map.size())).first;
                order.push_back(items[i]);
            }
            x[i] = it->second + 1;
        }
        out[l] = x;
        lower = order;
    }
    return out;
}

void nested_refine_loop(Nested& n, const Params& p, double c, std::vector<double>& trace) {
    int L = n.levels();
    auto params = [&](int l) {
        Params q = p;
        q.c = c * std::pow(2.0, l);
        if (l > 0) {
            q.B_min = 1;
            q.B_max = std::numeric_limits<size_t>::max();
        }
        return q;
    };
    for (size_t iter = 0; iter < p.max_refine; ++iter) {
        double dS = 0;
        for (int k = 0; k < L; ++k) {
            LevelView view(n, k);
            if (view.N() <= 1)
                continue;
            Params q = params(k);
            Multilevel<LevelView> ml(view, q);
            dS += ml.merge_and_sweep(true, q.B_min, q.B_max);
        }
        double S = n.entropy();
        trace.push_back(S);
        if (p.verbose)
            Rprintf("refine: dS = %.6f, S = %.4f\n", dS, S);
        if (S == 0 || std::abs(dS / S) <= p.epsilon)
            break;
        Rcpp::checkUserInterrupt();
    }
}

void nested_minimize_simple_init(Nested& n, Params p, double c, bool refine,
                                 bool simple_init_reverse, double base_factor,
                                 int B_min_base, double top_factor, bool strict,
                                 std::vector<double>& trace) {
    int L = n.levels();
    auto params = [&](int l) {
        Params q = p;
        q.c = c * std::pow(2.0, l);
        q.bisection = false;
        if (l > 0) {
            q.B_min = 1;
            q.B_max = std::numeric_limits<size_t>::max();
        }
        return q;
    };

    std::vector<int> level_order(L);
    std::iota(level_order.begin(), level_order.end(), 0);
    if (simple_init_reverse) {
        // Seed initial geometric hierarchy top-down before reverse sweep
        int N0 = n.base().N();
        Nested::Hierarchy bs(L);
        size_t cur_n = static_cast<size_t>(N0);
        for (int l = 0; l < L; ++l) {
            size_t B_tgt;
            if (l == 0) {
                double raw = std::max(cur_n / base_factor, static_cast<double>(B_min_base));
                B_tgt = static_cast<size_t>(std::max(1.0, std::floor(raw)));
                B_tgt = std::max(p.B_min, std::min(std::min(p.B_max, cur_n), B_tgt));
            } else {
                B_tgt = std::max<size_t>(1, static_cast<size_t>(std::floor(cur_n / top_factor)));
            }
            bs[l].resize(std::max<size_t>(cur_n, 1), 0);
            for (size_t i = 0; i < cur_n; ++i)
                bs[l][i] = static_cast<int>(i % B_tgt);
            cur_n = B_tgt;
        }
        n.build(bs);
        std::reverse(level_order.begin(), level_order.end());
    }

    for (int l : level_order) {
        LevelView view(n, l);
        if (view.N() <= 1)
            continue;
        size_t B_tgt;
        if (l == 0) {
            double raw = std::max(view.N() / base_factor, static_cast<double>(B_min_base));
            B_tgt = static_cast<size_t>(std::max(1.0, std::floor(raw)));
            B_tgt = std::max(p.B_min, std::min(std::min(p.B_max, view.N()), B_tgt));
        } else {
            B_tgt = std::max<size_t>(1, static_cast<size_t>(std::floor(view.N() / top_factor)));
        }
        Params q = params(l);
        Multilevel<LevelView> ml(view, q);
        ml.stage_agglomerative(B_tgt, strict ? B_tgt : view.N(), strict);
        double S = n.entropy();
        trace.push_back(S);
        if (p.verbose)
            Rprintf("simple_init level %d: B = %d (target %d), S = %.4f\n",
                    l, view.B(), int(B_tgt), S);
        Rcpp::checkUserInterrupt();
    }

    if (refine)
        nested_refine_loop(n, p, c, trace);
}

void nested_minimize(Nested& n, Params p, double c, bool refine, std::vector<double>& trace) {
    int L = n.levels();
    auto params = [&](int l) {
        Params q = p;
        q.c = c * std::pow(2.0, l);
        if (l > 0) {
            q.B_min = 1;
            q.B_max = std::numeric_limits<size_t>::max();
        }
        return q;
    };
    int l = 0;
    bool forced_pending;
    {
        size_t B0 = n.base().B();
        forced_pending = B0 < p.B_min || B0 > p.B_max;
    }
    while (l >= 0) {
        LevelView view(n, l);
        auto before = view.save();
        double S_old = n.entropy();
        Multilevel<LevelView> ml(view, params(l));
        double S_new = ml.stage_multilevel();
        double dS = S_new - S_old;
        // Unlike graph-tool, accept the first level-0 result when the starting
        // partition violates the requested group-count bounds.
        size_t B_now = view.B();
        bool forced = l == 0 && forced_pending;
        forced_pending = false;
        if (forced && (B_now < p.B_min || B_now > p.B_max))
            forced = false;
        if (!(dS < 0) && !forced) {
            view.restore(before);
            dS = 0;
        }
        if (p.verbose)
            Rprintf("level %d: dS = %.6f, S = %.4f\n", l, dS, n.entropy());
        trace.push_back(n.entropy());
        if (std::abs(dS) < 1e-8 && !forced)
            --l;
        else
            l = std::min(l + 1, L - 1);
        Rcpp::checkUserInterrupt();
    }
    if (refine)
        nested_refine_loop(n, p, c, trace);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_entropy(SEXP pointer, Rcpp::List bs, bool degree_corrected,
                                       bool distributed, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = boost::num_vertices(graph.adjacency);
    Nested n(graph, degree_corrected, distributed, q_cache_max(q_max), bs.size());
    n.build(internal_hierarchy(N, bs));
    std::vector<double> S = n.level_entropies();
    return Rcpp::NumericVector(S.begin(), S.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_check_moves(SEXP pointer, Rcpp::List bs, bool degree_corrected,
                                           bool distributed, int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = boost::num_vertices(graph.adjacency);
    Nested n(graph, degree_corrected, distributed, q_cache_max(q_max), bs.size());
    n.build(internal_hierarchy(N, bs));
    RRng rng;
    Rcpp::NumericVector err(nmoves);
    double S = n.entropy();
    for (int i = 0; i < nmoves; ++i) {
        int l = rng.index(n.levels());
        LevelView view(n, l);
        if (view.N() == 0)
            continue;
        if (rng.unif() < 0.15 && view.B() > 1) {
            const auto& occ = view.occupied();
            int r = occ[rng.index(occ.size())];
            int s = occ[rng.index(occ.size() - 1)];
            if (s == r)
                s = occ.back();
            double dS = view.virtual_merge(r, s);
            view.merge(r, s);
            double S_new = n.entropy();
            err[i] = dS - (S_new - S);
            S = S_new;
        } else {
            int v = view.nodes()[rng.index(view.N())];
            int s;
            int empty = l == 0 ? n.base().empty_group() : n.upper(l).empty_group();
            if (rng.unif() < 0.25 && empty >= 0) {
                s = empty;
                n.inherit(l, s, view.group(v));
            } else {
                const auto& occ = view.occupied();
                s = occ[rng.index(occ.size())];
            }
            double dS = view.virtual_move(v, s);
            view.move(v, s);
            double S_new = n.entropy();
            err[i] = dS - (S_new - S);
            S = S_new;
        }
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_minimize(SEXP pointer, bool degree_corrected, bool distributed,
                               Rcpp::List control, bool refine, double q_max,
                               bool simple_init = false, bool simple_init_reverse = false,
                               double base_factor = 10.0, int B_min_base = 1000,
                               double top_factor = 2.0, bool strict = false,
                               Rcpp::Nullable<Rcpp::List> init_levels = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = boost::num_vertices(graph.adjacency);
    Params p = read_params(control);
    double c = control.containsElementNamed("c") ? Rcpp::as<double>(control["c"]) : 1.0;

    int L;
    if (init_levels.isNotNull()) {
        Rcpp::List il(init_levels);
        L = std::max(static_cast<int>(il.size()), 1);
    } else if (simple_init && N > 1) {
        // Depth matching the requested geometric hierarchy shape
        L = 1;
        size_t cur_b = static_cast<size_t>(std::max(1.0, std::floor(
            std::max(N / base_factor, static_cast<double>(std::min(B_min_base, N))))));
        cur_b = std::max(p.B_min, std::min(std::min(p.B_max, static_cast<size_t>(N)), cur_b));
        while (cur_b > 1) {
            ++L;
            cur_b = static_cast<size_t>(std::max(1.0, std::floor(cur_b / top_factor)));
        }
        L = std::max(L, 2);
    } else {
        L = N > 1 ? static_cast<int>(std::ceil(std::log2(double(N)))) + 1 : 1;
    }

    Nested n(graph, degree_corrected, distributed, q_cache_max(q_max), L);
    std::vector<double> trace;
    if (init_levels.isNotNull()) {
        Rcpp::List il(init_levels);
        n.build(internal_hierarchy(N, il));
        trace.push_back(n.entropy());
        if (N > 1 && refine)
            nested_refine_loop(n, p, c, trace);
    } else {
        Nested::Hierarchy bs(L);
        bs[0].assign(N, 0);
        for (int l = 1; l < L; ++l)
            bs[l].assign(1, 0);
        n.build(bs);
        if (N > 1) {
            if (simple_init) {
                nested_minimize_simple_init(n, p, c, refine, simple_init_reverse,
                                            base_factor, B_min_base, top_factor,
                                            strict, trace);
            } else {
                nested_minimize(n, p, c, refine, trace);
            }
        }
    }
    std::vector<double> S = n.level_entropies();
    return Rcpp::List::create(Rcpp::Named("bs") = external_hierarchy(n.hierarchy()),
                              Rcpp::Named("level_entropy") = S,
                              Rcpp::Named("trace") = trace);
}


// ---------------------------------------------------------------------------
// Exact finite-temperature single-node and merge-split MCMC (Peixoto 2014,
// 2020), equilibration, and simulated annealing.

namespace {

inline double log_sum_exp(double a, double b) {
    if (std::isinf(a) && a < 0)
        return b;
    if (std::isinf(b) && b < 0)
        return a;
    double hi = std::max(a, b);
    double lo = std::min(a, b);
    return hi + std::log1p(std::exp(lo - hi));
}

struct MCMCParams {
    double beta = 1.0;
    double c = 0.5;
    double d = 0.01;
    size_t niter = 1;
    size_t B_min = 1;
    size_t B_max = std::numeric_limits<size_t>::max();
    bool multiflip = true;
    double psingle = -1.0;  // < 0 => defaults to N
    double psplit = 1.0;
    double pmerge = 1.0;
    double pmergesplit = 1.0;
    size_t gibbs_sweeps = 10;
    bool ls_shuffle = true;
};

struct EquilParams {
    std::string mode = "sweep";  // "sweep", "equilibrate", "anneal"
    size_t wait = 1000;
    size_t nbreaks = 2;
    size_t max_niter = std::numeric_limits<size_t>::max();
    long force_niter = -1;
    double epsilon = 0.0;
    bool history = false;
    bool collect = false;
    size_t thin = 1;
    bool verbose = false;
    // Annealing parameters
    double beta_start = 1.0;
    double beta_end = 10.0;
    size_t anneal_niter = 100;
};

struct SweepResult {
    double dS = 0;
    size_t nattempts = 0;
    size_t nmoves = 0;
};

MCMCParams read_mcmc_params(const Rcpp::List& opts, double default_c) {
    MCMCParams p;
    p.c = default_c;
    auto num = [&](const char* name, double def) {
        return opts.containsElementNamed(name) ? Rcpp::as<double>(opts[name]) : def;
    };
    auto flag = [&](const char* name, bool def) {
        return opts.containsElementNamed(name) ? Rcpp::as<bool>(opts[name]) : def;
    };
    p.beta = num("beta", p.beta);
    p.c = num("c", p.c);
    p.d = num("d", p.d);
    p.niter = static_cast<size_t>(std::max(0.0, num("niter", static_cast<double>(p.niter))));
    p.B_min = static_cast<size_t>(std::max(1.0, num("B_min", 1.0)));
    double bmax = num("B_max", inf);
    p.B_max = std::isinf(bmax) ? std::numeric_limits<size_t>::max() : static_cast<size_t>(bmax);
    p.multiflip = flag("multiflip", p.multiflip);
    p.psingle = num("psingle", p.psingle);
    p.psplit = num("psplit", p.psplit);
    p.pmerge = num("pmerge", p.pmerge);
    p.pmergesplit = num("pmergesplit", p.pmergesplit);
    p.gibbs_sweeps = static_cast<size_t>(std::max(1.0, num("gibbs_sweeps", 10.0)));
    p.ls_shuffle = flag("ls_shuffle", p.ls_shuffle);
    return p;
}

EquilParams read_equil_params(const Rcpp::List& opts) {
    EquilParams e;
    auto num = [&](const char* name, double def) {
        return opts.containsElementNamed(name) ? Rcpp::as<double>(opts[name]) : def;
    };
    auto flag = [&](const char* name, bool def) {
        return opts.containsElementNamed(name) ? Rcpp::as<bool>(opts[name]) : def;
    };
    if (opts.containsElementNamed("mode"))
        e.mode = Rcpp::as<std::string>(opts["mode"]);
    e.wait = static_cast<size_t>(std::max(1.0, num("wait", 1000.0)));
    e.nbreaks = static_cast<size_t>(std::max(1.0, num("nbreaks", 2.0)));
    double mx = num("max_niter", inf);
    e.max_niter = std::isinf(mx) ? std::numeric_limits<size_t>::max() : static_cast<size_t>(mx);
    e.force_niter = static_cast<long>(num("force_niter", -1.0));
    e.epsilon = num("epsilon", 0.0);
    e.history = flag("history", false);
    e.collect = flag("collect", false);
    e.thin = static_cast<size_t>(std::max(1.0, num("thin", 1.0)));
    e.verbose = flag("verbose", false);
    e.beta_start = num("beta_start", 1.0);
    e.beta_end = num("beta_end", 10.0);
    e.anneal_niter = static_cast<size_t>(std::max(1.0, num("anneal_niter", 100.0)));
    return e;
}

template <class State>
class MCMC {
public:
    MCMC(State& state, RRng& rng) : _s(state), _rng(rng) {}

    std::pair<double, bool> single_move(int v, const MCMCParams& p) {
        size_t N = _s.N();
        if (N <= 1)
            return {0.0, false};
        int r = _s.group(v);
        size_t B = static_cast<size_t>(_s.B());
        bool can_create = (B < N && B < p.B_max && p.d > 0);
        double p_new = can_create ? p.d : 0.0;

        if (can_create && _rng.unif() < p_new) {
            if (_s.size(r) == 1)
                return {0.0, false};
            int s = _s.empty_group();
            if (s < 0)
                return {0.0, false};
            _s.inherit(s, r);
            double dS = _s.virtual_move(v, s);
            bool accept = false;
            if (std::isinf(p.beta)) {
                accept = (dS < 0);
            } else {
                bool rev_can_create = (B + 1 < N && B + 1 < p.B_max && p.d > 0);
                double p_new_rev = rev_can_create ? p.d : 0.0;
                double prev = (1.0 - p_new_rev) * _s.move_prob(v, r, s, p.c, true);
                double pfwd = p_new;
                if (prev > 0 && pfwd > 0) {
                    double log_a = -p.beta * dS + std::log(prev) - std::log(pfwd);
                    accept = (log_a >= 0 || _rng.unif() < std::exp(log_a));
                }
            }
            if (accept) {
                _s.move(v, s);
                return {dS, true};
            }
            return {0.0, false};
        }

        int s = _s.sample_exact_group(v, p.c, _rng);
        if (s == r)
            return {0.0, false};
        bool empties = (_s.size(r) == 1);
        if (empties) {
            if (B <= p.B_min || !_s.same_parent(r, s))
                return {0.0, false};
        }
        double dS = _s.virtual_move(v, s);
        bool accept = false;
        if (std::isinf(p.beta)) {
            accept = (dS < 0);
        } else {
            double pfwd = (1.0 - p_new) * _s.move_prob(v, r, s, p.c, false);
            double prev = 0.0;
            if (empties) {
                // Reverse move creates an empty group from state with B - 1 groups.
                prev = (p.d > 0) ? p.d : 0.0;
            } else {
                prev = (1.0 - p_new) * _s.move_prob(v, r, s, p.c, true);
            }
            if (prev > 0 && pfwd > 0) {
                double log_a = -p.beta * dS + std::log(prev) - std::log(pfwd);
                accept = (log_a >= 0 || _rng.unif() < std::exp(log_a));
            }
        }
        if (accept) {
            _s.move(v, s);
            return {dS, true};
        }
        return {0.0, false};
    }

    SweepResult mcmc_sweep(const MCMCParams& p) {
        SweepResult res;
        size_t N = _s.N();
        if (N <= 1)
            return res;
        std::vector<int> vs = _s.nodes();
        _rng.shuffle(vs);
        for (int v : vs) {
            ++res.nattempts;
            auto [dS, moved] = single_move(v, p);
            if (moved) {
                res.dS += dS;
                ++res.nmoves;
            }
        }
        return res;
    }

    SweepResult multiflip_sweep(const MCMCParams& p) {
        SweepResult res;
        size_t N = _s.N();
        if (N <= 1)
            return res;
        const auto& nodes = _s.nodes();
        double ps = (p.psingle < 0) ? static_cast<double>(N) : p.psingle;
        double psp = std::max(0.0, p.psplit);
        double pm = std::max(0.0, p.pmerge);
        double pms = std::max(0.0, p.pmergesplit);
        double ptot = ps + psp + pm + pms;
        if (ptot <= 0)
            return res;

        size_t nsteps = std::max<size_t>(1, static_cast<size_t>(std::llround(ptot)));
        double c1 = ps / ptot;
        double c2 = (ps + psp) / ptot;
        double c3 = (ps + psp + pm) / ptot;

        for (size_t step = 0; step < nsteps; ++step) {
            double u = _rng.unif();
            if (u < c1) {
                int v = nodes[_rng.index(N)];
                ++res.nattempts;
                auto [dS, moved] = single_move(v, p);
                if (moved) {
                    res.dS += dS;
                    ++res.nmoves;
                }
            } else if (u < c2) {
                ++res.nattempts;
                auto [dS, nm] = split_move(p, psp, pm);
                res.dS += dS;
                res.nmoves += nm;
            } else if (u < c3) {
                ++res.nattempts;
                auto [dS, nm] = merge_move(p, psp, pm);
                res.dS += dS;
                res.nmoves += nm;
            } else {
                ++res.nattempts;
                auto [dS, nm] = mergesplit_move(p);
                res.dS += dS;
                res.nmoves += nm;
            }
        }
        return res;
    }

private:
    State& _s;
    RRng& _rng;

    double gibbs_step(int v, int r1, int r2, double beta, bool sample, int target) {
        int cur = _s.group(v);
        int other = (cur == r1) ? r2 : r1;
        double dS_other = _s.virtual_move(v, other);
        if (std::isinf(beta)) {
            if (sample) {
                if (dS_other < 0 || (dS_other == 0 && _rng.unif() < 0.5))
                    _s.move(v, other);
            } else if (target == other) {
                _s.move(v, other);
            }
            return 0.0;
        }
        double x = -beta * dS_other;
        double log_p_other, log_p_cur;
        if (x > 0) {
            double e = std::exp(-x);
            log_p_other = -std::log1p(e);
            log_p_cur = -x - std::log1p(e);
        } else {
            double e = std::exp(x);
            log_p_other = x - std::log1p(e);
            log_p_cur = -std::log1p(e);
        }
        if (sample) {
            if (_rng.unif() < std::exp(log_p_other)) {
                _s.move(v, other);
                return log_p_other;
            }
            return log_p_cur;
        }
        if (target == other) {
            _s.move(v, other);
            return log_p_other;
        }
        return log_p_cur;
    }

    template <class Item, class GetV>
    void make_aux_split(std::vector<Item>& items, GetV get_v, int r1, int r2,
                        double beta, size_t gibbs_sweeps, std::vector<int>& aux) {
        _rng.shuffle(items);
        for (const auto& it : items) {
            int v = get_v(it);
            int dst = (_rng.unif() < 0.5) ? r1 : r2;
            if (_s.group(v) != dst)
                _s.move(v, dst);
        }
        for (size_t sw = 1; sw < gibbs_sweeps; ++sw) {
            _rng.shuffle(items);
            for (const auto& it : items)
                gibbs_step(get_v(it), r1, r2, beta, true, 0);
        }
        _rng.shuffle(items);
        aux.resize(items.size());
        for (size_t i = 0; i < items.size(); ++i)
            aux[i] = _s.group(get_v(items[i]));
    }

    template <class Item, class GetV>
    void restore_assign(const std::vector<Item>& items, GetV get_v,
                        const std::vector<int>& assign) {
        for (size_t i = 0; i < items.size(); ++i) {
            int v = get_v(items[i]);
            if (_s.group(v) != assign[i])
                _s.move(v, assign[i]);
        }
    }

    std::pair<double, size_t> split_move(const MCMCParams& p, double psplit, double pmerge) {
        size_t B_before = static_cast<size_t>(_s.B());
        if (B_before >= _s.N() || B_before >= p.B_max || psplit <= 0)
            return {0.0, 0};
        if (!std::isinf(p.beta) && pmerge <= 0)
            return {0.0, 0};
        int r1 = _s.random_occupied(_rng);
        if (_s.size(r1) < 2)
            return {0.0, 0};
        int r2 = _s.empty_group();
        if (r2 < 0)
            return {0.0, 0};
        _s.inherit(r2, r1);

        std::vector<int> vs = _s.members(r1);
        auto get_v = [](int v) { return v; };
        std::vector<int> aux;
        make_aux_split(vs, get_v, r1, r2, p.beta, p.gibbs_sweeps, aux);

        double log_p_z = 0.0;
        std::vector<int> prop(vs.size());
        for (size_t i = 0; i < vs.size(); ++i) {
            log_p_z += gibbs_step(vs[i], r1, r2, p.beta, true, 0);
            prop[i] = _s.group(vs[i]);
        }

        if (_s.size(r1) == 0 || _s.size(r2) == 0) {
            for (int v : vs)
                if (_s.group(v) != r1)
                    _s.move(v, r1);
            return {0.0, 0};
        }

        double dS = -_s.virtual_merge(r2, r1);
        if (std::isinf(p.beta)) {
            if (dS < 0)
                return {dS, static_cast<size_t>(_s.size(r2))};
            _s.merge(r2, r1);
            return {0.0, 0};
        }

        double P_merge = _s.merge_prob(r1, r2, p.c);
        restore_assign(vs, get_v, aux);
        double log_p_zbar = 0.0;
        for (size_t i = 0; i < vs.size(); ++i) {
            int tgt = (prop[i] == r1) ? r2 : r1;
            log_p_zbar += gibbs_step(vs[i], r1, r2, p.beta, false, tgt);
        }

        double log_pfwd = -std::log(static_cast<double>(B_before)) +
                          log_sum_exp(log_p_z, log_p_zbar);
        double log_prev = std::log(P_merge);
        double log_a = -p.beta * dS + std::log(pmerge / psplit) + log_prev - log_pfwd;
        if (P_merge > 0 && (log_a >= 0 || _rng.unif() < std::exp(log_a)))
            return {dS, static_cast<size_t>(_s.size(r2))};
        _s.merge(r2, r1);
        return {0.0, 0};
    }

    std::pair<double, size_t> merge_move(const MCMCParams& p, double psplit, double pmerge) {
        size_t B_before = static_cast<size_t>(_s.B());
        if (B_before <= 1 || B_before <= p.B_min || pmerge <= 0)
            return {0.0, 0};
        if (!std::isinf(p.beta) && psplit <= 0)
            return {0.0, 0};
        int r1 = _s.random_occupied(_rng);
        int r2 = _s.sample_merge_target(r1, p.c, _rng);
        if (r2 == r1 || !_s.same_parent(r1, r2))
            return {0.0, 0};

        double dS = _s.virtual_merge(r2, r1);
        size_t nm = static_cast<size_t>(_s.size(r2));
        if (std::isinf(p.beta)) {
            if (dS < 0) {
                _s.merge(r2, r1);
                return {dS, nm};
            }
            return {0.0, 0};
        }

        double P_merge = _s.merge_prob(r1, r2, p.c);
        std::vector<std::pair<int, int>> vo;
        vo.reserve(_s.size(r1) + _s.size(r2));
        for (int v : _s.members(r1))
            vo.emplace_back(v, r1);
        for (int v : _s.members(r2))
            vo.emplace_back(v, r2);
        auto get_v = [](const std::pair<int, int>& x) { return x.first; };

        std::vector<int> aux;
        make_aux_split(vo, get_v, r1, r2, p.beta, p.gibbs_sweeps, aux);

        double log_p_z = 0.0;
        for (size_t i = 0; i < vo.size(); ++i)
            log_p_z += gibbs_step(vo[i].first, r1, r2, p.beta, false, vo[i].second);

        restore_assign(vo, get_v, aux);
        double log_p_zbar = 0.0;
        for (size_t i = 0; i < vo.size(); ++i) {
            int tgt = (vo[i].second == r1) ? r2 : r1;
            log_p_zbar += gibbs_step(vo[i].first, r1, r2, p.beta, false, tgt);
        }

        double log_prev = -std::log(static_cast<double>(B_before - 1)) +
                          log_sum_exp(log_p_z, log_p_zbar);
        double log_pfwd = std::log(P_merge);
        double log_a = -p.beta * dS + std::log(psplit / pmerge) + log_prev - log_pfwd;
        if (P_merge > 0 && (log_a >= 0 || _rng.unif() < std::exp(log_a))) {
            _s.merge(r2, r1);
            return {dS, nm};
        }
        // Restore exact original assignment.
        for (const auto& [v, g] : vo)
            if (_s.group(v) != g)
                _s.move(v, g);
        return {0.0, 0};
    }

    std::pair<double, size_t> mergesplit_move(const MCMCParams& p) {
        size_t B_before = static_cast<size_t>(_s.B());
        if (B_before <= 1)
            return {0.0, 0};
        int r1 = _s.random_occupied(_rng);
        int r2 = _s.sample_merge_target(r1, p.c, _rng);
        if (r2 == r1 || !_s.same_parent(r1, r2))
            return {0.0, 0};

        double P_before = _s.merge_prob(r1, r2, p.c);
        double dS_to_merged = _s.virtual_merge(r2, r1);

        std::vector<std::pair<int, int>> vo;
        vo.reserve(_s.size(r1) + _s.size(r2));
        for (int v : _s.members(r1))
            vo.emplace_back(v, r1);
        for (int v : _s.members(r2))
            vo.emplace_back(v, r2);
        auto get_v = [](const std::pair<int, int>& x) { return x.first; };

        std::vector<int> aux;
        make_aux_split(vo, get_v, r1, r2, p.beta, p.gibbs_sweeps, aux);

        double log_p_prop1 = 0.0;
        std::vector<int> prop(vo.size());
        for (size_t i = 0; i < vo.size(); ++i) {
            log_p_prop1 += gibbs_step(vo[i].first, r1, r2, p.beta, true, 0);
            prop[i] = _s.group(vo[i].first);
        }

        if (_s.size(r1) == 0 || _s.size(r2) == 0) {
            for (const auto& [v, g] : vo)
                if (_s.group(v) != g)
                    _s.move(v, g);
            return {0.0, 0};
        }

        double dS_prop_to_merged = _s.virtual_merge(r2, r1);
        double dS = dS_to_merged - dS_prop_to_merged;
        size_t nm = 0;
        for (size_t i = 0; i < vo.size(); ++i)
            if (prop[i] != vo[i].second)
                ++nm;

        if (std::isinf(p.beta)) {
            if (dS < 0)
                return {dS, nm};
            for (const auto& [v, g] : vo)
                if (_s.group(v) != g)
                    _s.move(v, g);
            return {0.0, 0};
        }

        double P_after = _s.merge_prob(r1, r2, p.c);

        restore_assign(vo, get_v, aux);
        double log_p_prop2 = 0.0;
        for (size_t i = 0; i < vo.size(); ++i) {
            int tgt = (prop[i] == r1) ? r2 : r1;
            log_p_prop2 += gibbs_step(vo[i].first, r1, r2, p.beta, false, tgt);
        }

        restore_assign(vo, get_v, aux);
        double log_p_orig1 = 0.0;
        for (size_t i = 0; i < vo.size(); ++i)
            log_p_orig1 += gibbs_step(vo[i].first, r1, r2, p.beta, false, vo[i].second);

        restore_assign(vo, get_v, aux);
        double log_p_orig2 = 0.0;
        for (size_t i = 0; i < vo.size(); ++i) {
            int tgt = (vo[i].second == r1) ? r2 : r1;
            log_p_orig2 += gibbs_step(vo[i].first, r1, r2, p.beta, false, tgt);
        }

        double log_p_prop = log_sum_exp(log_p_prop1, log_p_prop2);
        double log_p_orig = log_sum_exp(log_p_orig1, log_p_orig2);
        double log_a = -p.beta * dS + std::log(P_after) + log_p_orig -
                       std::log(P_before) - log_p_prop;
        if (P_before > 0 && P_after > 0 && (log_a >= 0 || _rng.unif() < std::exp(log_a))) {
            restore_assign(vo, get_v, prop);
            return {dS, nm};
        }
        for (const auto& [v, g] : vo)
            if (_s.group(v) != g)
                _s.move(v, g);
        return {0.0, 0};
    }
};

struct RunHistory {
    std::vector<double> nattempts;
    std::vector<double> nmoves;
    std::vector<double> entropy;
    std::vector<double> beta;
};

template <class StepFn, class EntropyFn, class CollectFn>
SweepResult run_equilibrate(StepFn step_fn, EntropyFn entropy_fn, CollectFn collect_fn,
                            const EquilParams& eq, double beta_val,
                            RunHistory& hist, size_t& niter_out) {
    SweepResult total;
    size_t count = 0;
    size_t break_count = 0;
    size_t niter = 0;
    double S = entropy_fn();
    double min_S = S, max_S = S;
    constexpr double m_eps = 1e-6;

    while (count < eq.wait) {
        SweepResult r = step_fn();
        S += r.dS;
        ++niter;
        total.dS += r.dS;
        total.nattempts += r.nattempts;
        total.nmoves += r.nmoves;

        if (eq.force_niter >= 0) {
            if (S > max_S) max_S = S;
            if (S < min_S) min_S = S;
        } else {
            if (std::abs(r.dS) >= std::abs(S - r.dS) * eq.epsilon) {
                if (S > max_S + m_eps) {
                    max_S = S;
                    count = 0;
                } else if (S < min_S - m_eps) {
                    min_S = S;
                    count = 0;
                } else {
                    ++count;
                }
            } else {
                ++count;
            }
            if (count >= eq.wait) {
                ++break_count;
                if (break_count < eq.nbreaks) {
                    count = 0;
                    min_S = max_S = S;
                }
            }
        }

        if (eq.collect && (niter % eq.thin == 0))
            collect_fn(S);

        if (eq.history) {
            hist.nattempts.push_back(static_cast<double>(r.nattempts));
            hist.nmoves.push_back(static_cast<double>(r.nmoves));
            hist.entropy.push_back(S);
            hist.beta.push_back(beta_val);
        }

        if (eq.verbose) {
            Rprintf("niter: %5zu  count: %4zu  breaks: %2zu  min_S: %.6f  max_S: %.6f  S: %.6f  dS: %.6f  moves: %zu\n",
                    niter, count, break_count, min_S, max_S, S, r.dS, r.nmoves);
        }

        Rcpp::checkUserInterrupt();

        if (eq.force_niter >= 0 && niter >= static_cast<size_t>(eq.force_niter))
            break;
        if (niter >= eq.max_niter)
            break;
    }
    niter_out += niter;
    return total;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_sbm_mcmc(SEXP pointer, Rcpp::IntegerVector blocks, bool degree_corrected,
                        bool distributed, Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                        double q_max) {
    Graph& graph = get_graph(pointer);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    state.set_partition(zero_based(graph, blocks));
    FlatView view(state);
    RRng rng;
    MCMC<FlatView> mcmc(view, rng);
    MCMCParams mp = read_mcmc_params(mcmc_opts, 0.5);
    EquilParams eq = read_equil_params(eq_opts);

    SweepResult total;
    RunHistory hist;
    std::vector<Rcpp::IntegerVector> samples;
    std::vector<double> sample_entropies;
    size_t niter_done = 0;

    auto do_step = [&](const MCMCParams& cur_p) {
        SweepResult r;
        for (size_t i = 0; i < cur_p.niter; ++i) {
            SweepResult s = cur_p.multiflip ? mcmc.multiflip_sweep(cur_p)
                                            : mcmc.mcmc_sweep(cur_p);
            r.dS += s.dS;
            r.nattempts += s.nattempts;
            r.nmoves += s.nmoves;
        }
        return r;
    };

    if (eq.mode == "sweep") {
        total = do_step(mp);
        niter_done = mp.niter;
    } else if (eq.mode == "equilibrate") {
        total = run_equilibrate(
            [&]() { return do_step(mp); },
            [&]() { return view.entropy(); },
            [&](double cur_S) {
                samples.push_back(one_based(state.blocks()));
                sample_entropies.push_back(cur_S);
            },
            eq, mp.beta, hist, niter_done);
    } else if (eq.mode == "anneal") {
        double beta = eq.beta_start;
        double speed = std::exp((std::log(eq.beta_end) - std::log(eq.beta_start)) /
                                static_cast<double>(eq.anneal_niter));
        for (size_t step = 0; step <= eq.anneal_niter; ++step) {
            MCMCParams cur_p = mp;
            cur_p.beta = (step == eq.anneal_niter) ? eq.beta_end : beta;
            SweepResult r = run_equilibrate(
                [&]() { return do_step(cur_p); },
                [&]() { return view.entropy(); },
                [&](double cur_S) {
                    samples.push_back(one_based(state.blocks()));
                    sample_entropies.push_back(cur_S);
                },
                eq, cur_p.beta, hist, niter_done);
            total.dS += r.dS;
            total.nattempts += r.nattempts;
            total.nmoves += r.nmoves;
            beta *= speed;
        }
    } else {
        Rcpp::stop("Unknown MCMC mode: %s", eq.mode.c_str());
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("blocks") = one_based(state.blocks()),
        Rcpp::Named("entropy") = state.entropy(),
        Rcpp::Named("dS") = total.dS,
        Rcpp::Named("nattempts") = static_cast<double>(total.nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(total.nmoves),
        Rcpp::Named("niter") = static_cast<double>(niter_done));

    if (eq.history) {
        out["history"] = Rcpp::DataFrame::create(
            Rcpp::Named("nattempts") = hist.nattempts,
            Rcpp::Named("nmoves") = hist.nmoves,
            Rcpp::Named("entropy") = hist.entropy,
            Rcpp::Named("beta") = hist.beta);
    }
    if (eq.collect) {
        Rcpp::List slist(samples.size());
        for (size_t i = 0; i < samples.size(); ++i)
            slist[i] = samples[i];
        out["samples"] = slist;
        out["sample_entropies"] = sample_entropies;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_mcmc(SEXP pointer, Rcpp::List bs, bool degree_corrected,
                           bool distributed, Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                           double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    Nested n(graph, degree_corrected, distributed, q_cache_max(q_max), bs.size());
    n.build(internal_hierarchy(N, bs));
    RRng rng;
    MCMCParams mp = read_mcmc_params(mcmc_opts, 1.0);
    EquilParams eq = read_equil_params(eq_opts);

    int L = n.levels();
    std::vector<double> c_levels(L);
    if (mcmc_opts.containsElementNamed("c_levels")) {
        Rcpp::NumericVector cv = mcmc_opts["c_levels"];
        for (int l = 0; l < L; ++l)
            c_levels[l] = cv[std::min<R_xlen_t>(l, cv.size() - 1)];
    } else {
        for (int l = 0; l < L; ++l)
            c_levels[l] = mp.c * std::pow(2.0, l);
    }

    std::vector<int> lrange(L);
    std::iota(lrange.begin(), lrange.end(), 0);

    auto do_step = [&](const MCMCParams& cur_p) {
        SweepResult r;
        if (cur_p.ls_shuffle)
            rng.shuffle(lrange);
        for (int l : lrange) {
            LevelView view(n, l);
            if (view.N() <= 1)
                continue;
            MCMCParams lp = cur_p;
            lp.c = c_levels[l];
            if (l > 0) {
                lp.B_min = 1;
                lp.B_max = std::numeric_limits<size_t>::max();
            }
            if (lp.psingle < 0)
                lp.psingle = static_cast<double>(view.N());
            MCMC<LevelView> mcmc(view, rng);
            for (size_t i = 0; i < lp.niter; ++i) {
                SweepResult s = lp.multiflip ? mcmc.multiflip_sweep(lp)
                                             : mcmc.mcmc_sweep(lp);
                r.dS += s.dS;
                r.nattempts += s.nattempts;
                r.nmoves += s.nmoves;
            }
        }
        return r;
    };

    SweepResult total;
    RunHistory hist;
    std::vector<Rcpp::List> samples;
    std::vector<double> sample_entropies;
    size_t niter_done = 0;

    if (eq.mode == "sweep") {
        total = do_step(mp);
        niter_done = mp.niter;
    } else if (eq.mode == "equilibrate") {
        total = run_equilibrate(
            [&]() { return do_step(mp); },
            [&]() { return n.entropy(); },
            [&](double cur_S) {
                samples.push_back(external_hierarchy(n.hierarchy()));
                sample_entropies.push_back(cur_S);
            },
            eq, mp.beta, hist, niter_done);
    } else if (eq.mode == "anneal") {
        double beta = eq.beta_start;
        double speed = std::exp((std::log(eq.beta_end) - std::log(eq.beta_start)) /
                                static_cast<double>(eq.anneal_niter));
        for (size_t step = 0; step <= eq.anneal_niter; ++step) {
            MCMCParams cur_p = mp;
            cur_p.beta = (step == eq.anneal_niter) ? eq.beta_end : beta;
            SweepResult r = run_equilibrate(
                [&]() { return do_step(cur_p); },
                [&]() { return n.entropy(); },
                [&](double cur_S) {
                    samples.push_back(external_hierarchy(n.hierarchy()));
                    sample_entropies.push_back(cur_S);
                },
                eq, cur_p.beta, hist, niter_done);
            total.dS += r.dS;
            total.nattempts += r.nattempts;
            total.nmoves += r.nmoves;
            beta *= speed;
        }
    } else {
        Rcpp::stop("Unknown MCMC mode: %s", eq.mode.c_str());
    }

    std::vector<double> S = n.level_entropies();
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("bs") = external_hierarchy(n.hierarchy()),
        Rcpp::Named("level_entropy") = S,
        Rcpp::Named("entropy") = std::accumulate(S.begin(), S.end(), 0.0),
        Rcpp::Named("dS") = total.dS,
        Rcpp::Named("nattempts") = static_cast<double>(total.nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(total.nmoves),
        Rcpp::Named("niter") = static_cast<double>(niter_done));

    if (eq.history) {
        out["history"] = Rcpp::DataFrame::create(
            Rcpp::Named("nattempts") = hist.nattempts,
            Rcpp::Named("nmoves") = hist.nmoves,
            Rcpp::Named("entropy") = hist.entropy,
            Rcpp::Named("beta") = hist.beta);
    }
    if (eq.collect) {
        Rcpp::List slist(samples.size());
        for (size_t i = 0; i < samples.size(); ++i)
            slist[i] = samples[i];
        out["samples"] = slist;
        out["sample_entropies"] = sample_entropies;
    }
    return out;
}

// Verifies exact proposal normalization sum_s q(r->s|v) == 1 and exact
// reverse probability move_prob(v, r, s, reverse=true) == actual post-move
// move_prob(v, s, r, reverse=false).
// [[Rcpp::export]]
Rcpp::NumericVector cpp_sbm_check_move_prob(SEXP pointer, Rcpp::IntegerVector blocks,
                                            double c, int nmoves) {
    Graph& graph = get_graph(pointer);
    SBM state(graph, false, false, 100);
    state.set_partition(zero_based(graph, blocks));
    RRng rng;
    double max_sum_err = 0.0;
    double max_rev_err = 0.0;
    double max_merge_sum_err = 0.0;
    for (int i = 0; i < nmoves; ++i) {
        int v = rng.index(state.N());
        int r = state.group(v);
        double sum_q = 0.0;
        for (int s : state.occupied())
            sum_q += state.move_prob(v, r, s, c, false);
        max_sum_err = std::max(max_sum_err, std::abs(sum_q - 1.0));

        if (state.B() > 1) {
            const auto& occ = state.occupied();
            double sum_m = 0.0;
            for (size_t a = 0; a < occ.size(); ++a)
                for (size_t b = a + 1; b < occ.size(); ++b)
                    sum_m += state.merge_prob(occ[a], occ[b], c);
            max_merge_sum_err = std::max(max_merge_sum_err, std::abs(sum_m - 1.0));
        }

        int s = (rng.unif() < 0.2 && state.empty_group() >= 0)
                    ? state.empty_group()
                    : state.random_occupied(rng);
        if (s == r)
            continue;
        bool stays = state.size(r) > 1;
        double pred_rev = stays ? state.move_prob(v, r, s, c, true) : 0.0;
        state.move(v, s);
        if (stays) {
            double actual_rev = state.move_prob(v, s, r, c, false);
            max_rev_err = std::max(max_rev_err, std::abs(pred_rev - actual_rev));
        }
    }
    return Rcpp::NumericVector::create(
        Rcpp::Named("sum_err") = max_sum_err,
        Rcpp::Named("rev_err") = max_rev_err,
        Rcpp::Named("merge_sum_err") = max_merge_sum_err);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_check_move_prob(SEXP pointer, Rcpp::List bs,
                                               double c, int nmoves) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    Nested n(graph, false, false, 100, bs.size());
    n.build(internal_hierarchy(N, bs));
    RRng rng;
    double max_sum_err = 0.0;
    double max_rev_err = 0.0;
    double max_merge_sum_err = 0.0;
    for (int i = 0; i < nmoves; ++i) {
        int l = rng.index(n.levels());
        LevelView view(n, l);
        if (view.N() == 0)
            continue;
        int v = view.nodes()[rng.index(view.N())];
        int r = view.group(v);
        double sum_q = 0.0;
        for (int s : view.occupied())
            sum_q += view.move_prob(v, r, s, c, false);
        max_sum_err = std::max(max_sum_err, std::abs(sum_q - 1.0));

        if (view.B() > 1) {
            const auto& occ = view.occupied();
            double sum_m = 0.0;
            for (size_t a = 0; a < occ.size(); ++a)
                for (size_t b = a + 1; b < occ.size(); ++b)
                    sum_m += view.merge_prob(occ[a], occ[b], c);
            max_merge_sum_err = std::max(max_merge_sum_err, std::abs(sum_m - 1.0));
        }

        int s = (rng.unif() < 0.2 && view.empty_group() >= 0)
                    ? view.empty_group()
                    : view.random_occupied(rng);
        if (s == r)
            continue;
        bool stays = view.size(r) > 1;
        if (view.size(s) == 0)
            view.inherit(s, r);
        double pred_rev = stays ? view.move_prob(v, r, s, c, true) : 0.0;
        view.move(v, s);
        if (stays) {
            double actual_rev = view.move_prob(v, s, r, c, false);
            max_rev_err = std::max(max_rev_err, std::abs(pred_rev - actual_rev));
        }
    }
    return Rcpp::NumericVector::create(
        Rcpp::Named("sum_err") = max_sum_err,
        Rcpp::Named("rev_err") = max_rev_err,
        Rcpp::Named("merge_sum_err") = max_merge_sum_err);
}
