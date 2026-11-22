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
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "graph_types.h"

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

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
// Special functions

inline double lgam(double x) { return std::lgamma(x); }

inline double lbinom(double n, double k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return lgam(n + 1) - lgam(k + 1) - lgam(n - k + 1);
}

// log of the number of partitions of n into at most k parts. Exact values are
// tabulated (as in graph-tool's q cache) and an asymptotic approximation is used
// beyond the table.
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
        return approx(n, k);
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
                    std::log(2.0) * 3 / 2. - std::log(u) - std::log(M_PI);
        double g = 2 * v / u - u * std::log1p(-std::exp(-v));
        return lf - std::log(n) + std::sqrt(n) * g;
    }
    // spence(x) = Li2(1 - x) for 0 < x <= 1 (scipy/cephes convention)
    static double spence(double x) {
        return dilog(1 - x);
    }
    static double dilog(double z) {
        if (z == 0)
            return 0;
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

class SBM {
public:
    SBM(const Graph& graph, bool dc, bool dist, size_t q_cache_max, bool edges_dl = true)
        : _directed(graph.directed), _dc(dc), _dist(dist), _edges_dl(edges_dl) {
        const Adjacency& g = graph.adjacency;
        _N = static_cast<int>(boost::num_vertices(g));
        _E = static_cast<long>(boost::num_edges(g));
        _out.resize(_N);
        _in.resize(_N);
        _all.resize(_N);
        _loops.assign(_N, 0);
        _kout.assign(_N, 0);
        _kin.assign(_N, 0);
        std::map<std::pair<size_t, size_t>, long> mult;
        auto range = boost::edges(g);
        for (auto it = range.first; it != range.second; ++it) {
            int s = static_cast<int>(it->s), t = static_cast<int>(it->t);
            if (s == t) {
                _loops[s]++;
                _all[s].push_back(s);
                if (!_directed)
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
            size_t a = s, c = t;
            if (!_directed && a > c)
                std::swap(a, c);
            mult[{a, c}]++;
        }
        _const = 0;
        for (auto& [key, m] : mult) {
            _const += lgam(m + 1);
            if (!_directed && key.first == key.second)
                _const += m * std::log(2.0);
        }
        if (_dc) {
            for (int v = 0; v < _N; ++v)
                _const -= lgam(_kout[v] + 1) + (_directed ? lgam(_kin[v] + 1) : 0);
        }
        if (_dc && _dist) {
            size_t n_max = std::min<size_t>(std::max<long>(2 * _E, 100), 10000);
            n_max = std::min(n_max, q_cache_max);
            log_q_cache().ensure(n_max);
        }
        _key.reserve(64);
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
    const std::vector<std::unordered_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<std::unordered_map<int, long>>& m_in() const { return _m_in; }
    const std::vector<long>& ep() const { return _ep; }
    const std::vector<long>& em() const { return _em; }
    const std::vector<int>& wr() const { return _wr; }
    // Group-pair count changes of the last virtual_move() or move().
    const std::vector<Delta>& deltas() const { return _deltas; }

    void set_partition(const std::vector<int>& b) {
        _b = b;
        _wr.assign(_N, 0);
        _ep.assign(_N, 0);
        _em.assign(_N, 0);
        _m_out.assign(_N, {});
        if (_directed)
            _m_in.assign(_N, {});
        _h_out.assign(_N, {});
        if (_directed)
            _h_in.assign(_N, {});
        _members.assign(_N, {});
        _pos.assign(_N, 0);
        _occ_pos.assign(_N, -1);
        _occupied.clear();
        for (int v = 0; v < _N; ++v) {
            int r = _b[v];
            if (r < 0 || r >= _N)
                Rcpp::stop("Internal error: block label out of range");
            add_member(v, r);
            _ep[r] += _kout[v];
            _em[r] += _kin[v];
            _h_out[r][_kout[v]]++;
            if (_directed)
                _h_in[r][_kin[v]]++;
        }
        for (int v = 0; v < _N; ++v) {
            int r = _b[v];
            for (int u : _out[v]) {
                if (_directed) {
                    _m_out[r][_b[u]]++;
                    _m_in[_b[u]][r]++;
                } else {
                    _m_out[r][_b[u]]++;  // each undirected edge is seen from both ends
                }
            }
            if (_loops[v] > 0) {
                if (_directed) {
                    _m_out[r][r] += _loops[v];
                    _m_in[r][r] += _loops[v];
                } else {
                    _m_out[r][r] += 2 * _loops[v];
                }
            }
        }
        _empty.clear();
        for (int r = _N - 1; r >= 0; --r)
            if (_wr[r] == 0)
                _empty.push_back(r);
    }

    double global_term(long B) const {
        if (_N == 0)
            return 0;
        double pairs = _directed ? double(B) * B : double(B) * (B + 1) / 2;
        return lbinom(_N - 1, B - 1) + lgam(_N + 1) + std::log(double(_N)) +
               (_edges_dl ? lbinom(pairs + _E - 1, _E) : 0);
    }

    double group_term(long n, long ep, long em) const {
        if (n == 0)
            return 0;
        double S = -lgam(n + 1);
        if (_dc) {
            S += lgam(ep + 1) + (_directed ? lgam(em + 1) : 0);
            if (_dist) {
                const LogQ& q = log_q_cache();
                S += q(ep, n) + (_directed ? q(em, n) + 2 * lgam(n + 1) : lgam(n + 1));
            } else {
                S += lbinom(n + ep - 1, ep) + (_directed ? lbinom(n + em - 1, em) : 0);
            }
        } else {
            S += (ep + (_directed ? em : 0)) * std::log(double(n));
        }
        return S;
    }

    double edge_term(int x, int y, long m) const {
        if (m <= 0)
            return 0;
        if (!_directed && x == y)
            return lgam(m / 2.0 + 1) + (m / 2.0) * std::log(2.0);
        return lgam(m + 1);
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
                for (auto& kc : _h_out[r])
                    S -= lgam(kc.second + 1);
                if (_directed)
                    for (auto& kc : _h_in[r])
                        S -= lgam(kc.second + 1);
            }
            for (auto& [t, m] : _m_out[r])
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
        for (const Delta& e : _deltas) {
            if (e.d == 0)
                continue;
            long m = get_m(e.x, e.y);
            dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
        }
        return dS;
    }

    void move(int v, int s) {
        int r = _b[v];
        if (r == s)
            return;
        collect_deltas(v, r, s);
        for (const Delta& e : _deltas) {
            if (e.d == 0)
                continue;
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
        hist_move(_h_out, r, s, ko);
        if (_directed)
            hist_move(_h_in, r, s, ki);
        remove_member(v, r);
        if (_wr[s] == 0)
            take_empty(s);
        add_member(v, s);
        _b[v] = s;
        if (_wr[r] == 0)
            _empty.push_back(r);
    }

    int random_occupied(RRng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    // Neighbour-based proposal of graph-tool's BlockState::sample_group with
    // d = 0 and the fast (non-edge-group) proposal.
    int sample_group(int v, double c, RRng& rng) {
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

    std::vector<int> _b, _wr, _pos, _occ_pos, _occupied, _empty;
    std::vector<long> _ep, _em;
    std::vector<std::unordered_map<int, long>> _m_out, _m_in;
    std::vector<std::unordered_map<long, long>> _h_out, _h_in;
    std::vector<std::vector<int>> _members;

    std::vector<Delta> _deltas;
    std::unordered_map<long long, size_t> _key;

    static void add_m(std::unordered_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) {
            row.emplace(y, d);
        } else {
            it->second += d;
            if (it->second == 0)
                row.erase(it);
        }
    }

    double hist_delta(const std::vector<std::unordered_map<long, long>>& h,
                      int r, int s, long k) const {
        auto count = [&](int x) {
            auto it = h[x].find(k);
            return it == h[x].end() ? 0L : it->second;
        };
        long cr = count(r), cs = count(s);
        // the histogram term enters with a negative sign
        return -(lgam(cr) - lgam(cr + 1)) - (lgam(cs + 2) - lgam(cs + 1));
    }

    static void hist_move(std::vector<std::unordered_map<long, long>>& h,
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

    void remove_member(int v, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[v]] = last;
        _pos[last] = _pos[v];
        vs.pop_back();
        if (--_wr[r] == 0) {
            int i = _occ_pos[r];
            int other = _occupied.back();
            _occupied[i] = other;
            _occ_pos[other] = i;
            _occupied.pop_back();
            _occ_pos[r] = -1;
        }
    }

    void take_empty(int s) {
        auto it = std::find(_empty.rbegin(), _empty.rend(), s);
        if (it != _empty.rend())
            _empty.erase(std::next(it).base());
    }

    void add_delta(int x, int y, long d) {
        if (!_directed && x > y)
            std::swap(x, y);
        long long key = static_cast<long long>(x) * _N + y;
        auto it = _key.find(key);
        if (it == _key.end()) {
            _key.emplace(key, _deltas.size());
            _deltas.push_back({x, y, d});
        } else {
            _deltas[it->second].d += d;
        }
    }

    void collect_deltas(int v, int r, int s) {
        _deltas.clear();
        _key.clear();
        if (_directed) {
            for (int u : _out[v]) {
                int t = _b[u];
                add_delta(r, t, -1);
                add_delta(s, t, +1);
            }
            for (int u : _in[v]) {
                int t = _b[u];
                add_delta(t, r, -1);
                add_delta(t, s, +1);
            }
            if (_loops[v] > 0) {
                add_delta(r, r, -_loops[v]);
                add_delta(s, s, +_loops[v]);
            }
        } else {
            for (int u : _out[v]) {
                int t = _b[u];
                add_delta(r, t, t == r ? -2 : -1);
                add_delta(s, t, t == s ? +2 : +1);
            }
            if (_loops[v] > 0) {
                add_delta(r, r, -2 * _loops[v]);
                add_delta(s, s, +2 * _loops[v]);
            }
        }
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

// Accumulates group-pair count changes, normalising undirected keys to x <= y.
class DeltaAcc {
public:
    void reset(int n, bool directed) {
        _n = n;
        _directed = directed;
        clear();
    }
    void clear() {
        _v.clear();
        _idx.clear();
    }
    void add(int x, int y, long d) {
        if (!_directed && x > y)
            std::swap(x, y);
        long long key = static_cast<long long>(x) * _n + y;
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
    int _n = 0;
    bool _directed = false;
    std::vector<Delta> _v;
    std::unordered_map<long long, size_t> _idx;
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

// The statistics of the level below, which form this level's weighted graph.
struct LowerView {
    const std::vector<std::unordered_map<int, long>>* m_out;
    const std::vector<std::unordered_map<int, long>>* m_in;
    const std::vector<long>* ep;
    const std::vector<long>* em;
    const std::vector<int>* wr;
};

class BlockLevel {
public:
    BlockLevel(int cap, bool directed, long E, bool top)
        : _cap(cap), _directed(directed), _top(top), _E(E) {
        _acc.reset(cap, directed);
        _pairs.reset(cap, directed);
    }

    int B() const { return static_cast<int>(_occupied.size()); }
    int N() const { return _N; }
    int group(int x) const { return _b[x]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }
    int random_occupied(RRng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    LowerView view() const { return {&_m_out, &_m_in, &_ep, &_em, &_wr}; }

    void build(const std::vector<int>& b, const LowerView& lower) {
        _b = b;
        _wr.assign(_cap, 0);
        _ep.assign(_cap, 0);
        _em.assign(_cap, 0);
        _m_out.assign(_cap, {});
        if (_directed)
            _m_in.assign(_cap, {});
        _members.assign(_cap, {});
        _pos.assign(_cap, 0);
        _occ_pos.assign(_cap, -1);
        _occupied.clear();
        _N = 0;
        for (int x = 0; x < _cap; ++x) {
            if ((*lower.wr)[x] > 0) {
                add_member(x, _b[x]);
                ++_N;
            }
        }
        for (int x = 0; x < _cap; ++x) {
            if ((*lower.wr)[x] == 0)
                continue;
            int g = _b[x];
            for (auto& [y, val] : (*lower.m_out)[x]) {
                int h = _b[y];
                _m_out[g][h] += val;
                _ep[g] += val;
                if (_directed) {
                    _m_in[h][g] += val;
                    _em[h] += val;
                }
            }
        }
        _empty.clear();
        for (int r = _cap - 1; r >= 0; --r)
            if (_wr[r] == 0)
                _empty.push_back(r);
    }

    double partition_term(long N, long B) const {
        if (N == 0)
            return 0;
        return lbinom(N - 1, B - 1) + lgam(N + 1) + std::log(double(N));
    }

    double edges_term(long B) const {
        if (!_top)
            return 0;
        double pairs = _directed ? double(B) * B : double(B) * (B + 1) / 2;
        return lbinom(pairs + _E - 1, _E);
    }

    double dense_term(int x, int y, long m, long wx, long wy) const {
        if (m <= 0)
            return 0;
        double nn, e = m;
        if (_directed || x != y) {
            nn = double(wx) * wy;
        } else {
            nn = double(wx) * (wx + 1) / 2;
            e = m / 2.0;
        }
        return lbinom(nn + e - 1, e);
    }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0 : it->second;
    }

    double entropy() const {
        double S = partition_term(_N, B()) + edges_term(B());
        for (int r : _occupied) {
            S -= lgam(_wr[r] + 1);
            for (auto& [t, m] : _m_out[r])
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
        for (auto& [x, dw] : in.weights) {
            add_weight(_b[x], dw);
            dN += dw;
        }
        double dS = process(dN, out, apply);
        if (apply) {
            for (auto& [x, dw] : in.weights) {
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
        _acc.clear();
        if (_directed) {
            for (auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) {
                    _acc.add(g, g, -val);
                    _acc.add(h, h, val);
                } else {
                    int t = _b[y];
                    _acc.add(g, t, -val);
                    _acc.add(h, t, val);
                }
            }
            for (auto& [y, val] : (*lower.m_in)[x]) {
                if (y == x)
                    continue;
                int t = _b[y];
                _acc.add(t, g, -val);
                _acc.add(t, h, val);
            }
        } else {
            for (auto& [y, val] : (*lower.m_out)[x]) {
                if (y == x) {
                    _acc.add(g, g, -val);
                    _acc.add(h, h, val);
                } else {
                    int t = _b[y];
                    _acc.add(g, t, t == g ? -2 * val : -val);
                    _acc.add(h, t, t == h ? 2 * val : val);
                }
            }
        }
        _wd.clear();
        add_weight(g, -1);
        add_weight(h, +1);
        double dS = process(0, out, apply);
        if (apply) {
            pop_node(x, g);
            _b[x] = h;
            push_node(x, h);
        }
        return dS;
    }

    // Assign an empty (weight-zero) node to a group; this changes no statistics.
    void set_empty_label(int x, int g) { _b[x] = g; }

    int random_neighbor(int x, const LowerView& lower, RRng& rng) const {
        double total = (*lower.ep)[x] + (_directed ? (*lower.em)[x] : 0);
        if (total <= 0)
            return -1;
        double u = rng.unif() * total;
        for (auto& [y, val] : (*lower.m_out)[x]) {
            u -= val;
            if (u < 0)
                return y;
        }
        if (_directed) {
            for (auto& [y, val] : (*lower.m_in)[x]) {
                u -= val;
                if (u < 0)
                    return y;
            }
        }
        return x;
    }

    int sample_group(int x, double c, const LowerView& lower, RRng& rng) const {
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
    std::vector<int> _b, _wr, _pos, _occ_pos, _occupied, _empty;
    std::vector<long> _ep, _em;
    std::vector<std::unordered_map<int, long>> _m_out, _m_in;
    std::vector<std::vector<int>> _members;
    DeltaAcc _acc, _pairs;
    std::vector<std::pair<int, int>> _wd;

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
        for (auto& p : _wd)
            if (p.first == g)
                w += p.second;
        return w;
    }

    double process(int dN, ChangeSet& out, bool apply) {
        out.clear();
        long B0 = B(), B1 = B0;
        for (auto& [g, dw] : _wd) {
            int w0 = _wr[g], w1 = w0 + dw;
            if (w0 == 0 && w1 > 0) {
                ++B1;
                out.weights.emplace_back(g, +1);
            } else if (w0 > 0 && w1 == 0) {
                --B1;
                out.weights.emplace_back(g, -1);
            }
        }
        for (const Delta& e : _acc.items())
            if (e.d != 0)
                out.edges.push_back(e);

        double dS = 0;
        if (!apply) {
            dS += partition_term(_N + dN, B1) - partition_term(_N, B0);
            dS += edges_term(B1) - edges_term(B0);
            for (auto& [g, dw] : _wd)
                dS += lgam(_wr[g] + 1) - lgam(_wr[g] + dw + 1);
            _pairs.clear();
            for (const Delta& e : out.edges)
                _pairs.add(e.x, e.y, e.d);
            for (auto& [g, dw] : _wd) {
                if (dw == 0)
                    continue;
                for (auto& [t, m] : _m_out[g])
                    _pairs.add(g, t, 0);
                if (_directed)
                    for (auto& [t, m] : _m_in[g])
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

    static void add_m(std::unordered_map<int, long>& row, int y, long d) {
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

    void push_node(int x, int r) {
        bool was_empty = _wr[r] == 0;
        add_member(x, r);
        if (was_empty) {
            auto it = std::find(_empty.rbegin(), _empty.rend(), r);
            if (it != _empty.rend())
                _empty.erase(std::next(it).base());
        }
    }

    void pop_node(int x, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[x]] = last;
        _pos[last] = _pos[x];
        vs.pop_back();
        if (--_wr[r] == 0) {
            int i = _occ_pos[r];
            int other = _occupied.back();
            _occupied[i] = other;
            _occ_pos[other] = i;
            _occupied.pop_back();
            _occ_pos[r] = -1;
            _empty.push_back(r);
        }
    }
};

class Nested {
public:
    typedef std::vector<std::vector<int>> Hierarchy;

    Nested(const Graph& graph, bool dc, bool dist, size_t q_cache_max, int levels)
        : _base(graph, dc, dist, q_cache_max, levels <= 1) {
        int cap = std::max(_base.N(), 1);
        for (int l = 1; l < levels; ++l)
            _up.emplace_back(cap, _base.directed(), _base.E(), l == levels - 1);
        _cs.resize(levels + 1);
    }

    int levels() const { return static_cast<int>(_up.size()) + 1; }
    SBM& base() { return _base; }
    BlockLevel& upper(int l) { return _up[l - 1]; }
    LowerView lower(int l) const {
        if (l == 1)
            return {&_base.m_out(), &_base.m_in(), &_base.ep(), &_base.em(), &_base.wr()};
        return _up[l - 2].view();
    }

    void build(const Hierarchy& bs) {
        _base.set_partition(bs[0]);
        for (size_t l = 1; l < bs.size(); ++l)
            _up[l - 1].build(bs[l], lower(l));
    }

    Hierarchy hierarchy() const {
        Hierarchy bs;
        bs.push_back(_base.blocks());
        for (auto& u : _up)
            bs.push_back(u.blocks());
        return bs;
    }

    double entropy() const {
        double S = _base.entropy();
        for (auto& u : _up)
            S += u.entropy();
        return S;
    }

    std::vector<double> level_entropies() const {
        std::vector<double> S{_base.entropy()};
        for (auto& u : _up)
            S.push_back(u.entropy());
        return S;
    }

    int group(int l, int v) const { return l == 0 ? _base.group(v) : _up[l - 1].group(v); }
    int size(int l, int r) const { return l == 0 ? _base.size(r) : _up[l - 1].size(r); }

    double virtual_move(int l, int v, int s) { return change(l, v, s, false); }
    void move(int l, int v, int s) { change(l, v, s, true); }

    // A new group t at level l takes the upper-level label of group r.
    void inherit(int l, int t, int r) {
        if (l + 1 < levels())
            _up[l].set_empty_label(t, _up[l].group(r));
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
            cs.clear();
            bool empties = _base.size(r) == 1, fills = _base.size(s) == 0;
            if (apply)
                _base.move(v, s);
            else
                dS += _base.virtual_move(v, s);
            if (levels() == 1)
                return dS;
            cs.edges = _base.deltas();
            if (empties)
                cs.weights.emplace_back(r, -1);
            if (fills)
                cs.weights.emplace_back(s, +1);
        } else {
            dS += _up[l - 1].node_move(v, s, lower(l), cs, apply);
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
    typedef Nested::Hierarchy Snapshot;
    LevelView(Nested& n, int l) : _n(n), _l(l) {
        int cap = n.base().N();
        if (l == 0) {
            _nodes.resize(cap);
            std::iota(_nodes.begin(), _nodes.end(), 0);
        } else {
            const auto& wr = *n.lower(l).wr;
            for (int x = 0; x < cap; ++x)
                if (wr[x] > 0)
                    _nodes.push_back(x);
        }
    }
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
    void move(int v, int s) { _n.move(_l, v, s); }
    int sample_group(int v, double c, RRng& rng) {
        if (_l == 0)
            return _n.base().sample_group(v, c, rng);
        return _n.upper(_l).sample_group(v, c, _n.lower(_l), rng);
    }
    double entropy() const { return _n.entropy(); }
    Snapshot save() const { return _n.hierarchy(); }
    void restore(const Snapshot& bs) { _n.build(bs); }
    void spread() {
        for (int v : _nodes) {
            int r = group(v);
            if (size(r) <= 1)
                continue;
            int t = _l == 0 ? _n.base().empty_group() : _n.upper(_l).empty_group();
            _n.inherit(_l, t, r);
            move(v, t);
        }
    }
    void collapse(int t) {
        for (int v : _nodes)
            move(v, t);
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
    size_t N() const { return _nodes.size(); }
    const std::vector<int>& nodes() const { return _nodes; }
    int B() const { return _s.B(); }
    int group(int v) const { return _s.group(v); }
    int size(int r) const { return _s.size(r); }
    const std::vector<int>& occupied() const { return _s.occupied(); }
    const std::vector<int>& members(int r) const { return _s.members(r); }
    double virtual_move(int v, int s) { return _s.virtual_move(v, s); }
    void move(int v, int s) { _s.move(v, s); }
    int sample_group(int v, double c, RRng& rng) { return _s.sample_group(v, c, rng); }
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

    double virtual_merge(int r, int s) {
        std::vector<int> vs = _s.members(r);
        double dS = 0;
        size_t done = 0;
        for (int v : vs) {
            dS += _s.virtual_move(v, s);
            _s.move(v, s);
            ++done;
        }
        for (size_t i = 0; i < done; ++i)
            _s.move(vs[i], r);
        return dS;
    }

    void merge(int r, int s) {
        std::vector<int> vs = _s.members(r);
        for (int v : vs)
            _s.move(v, s);
    }

    double merge_sweep(size_t B_target, bool strict, size_t B_max) {
        std::vector<int> rlist = _s.occupied();
        std::vector<std::pair<int, int>> pairs;
        std::vector<double> dSs;
        for (int r : rlist) {
            int best = r;
            double best_dS = inf;
            std::unordered_set<int> past;
            auto find_candidates = [&](bool allow_random) {
                for (size_t i = 0; i < _p.merge_sweeps; ++i) {
                    const auto& vs = _s.members(r);
                    int v = vs[_rng.index(vs.size())];
                    int s = _s.sample_group(v, allow_random ? _p.c : 0, _rng);
                    if (s != r && past.find(s) == past.end()) {
                        double dS = virtual_merge(r, s);
                        if (dS < best_dS) {
                            best = s;
                            best_dS = dS;
                        }
                        past.insert(s);
                    }
                }
            };
            find_candidates(false);
            if (best == r)
                find_candidates(true);
            if (best != r) {
                pairs.emplace_back(r, best);
                dSs.push_back(best_dS);
            }
        }

        std::vector<size_t> pis(pairs.size());
        std::iota(pis.begin(), pis.end(), 0);
        _rng.shuffle(pis);
        auto cmp = [&](size_t i, size_t j) { return dSs[i] > dSs[j]; };
        std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)>
            queue(cmp, std::move(pis));

        std::unordered_map<int, int> root;
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
            double ndS = virtual_merge(s, t);
            if (!queue.empty() && ndS > dSs[queue.top()]) {
                dSs[i] = ndS;
                queue.push(i);
                continue;
            }
            if (ndS > 0 && strict && static_cast<size_t>(_s.B()) < B_max)
                break;
            root[s] = t;
            merge(s, t);
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

    // stage_multilevel() over all vertices at beta = infinity
    double stage_multilevel() {
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
        for (auto& [B, c] : _cache) {
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
    p.verbose = control.containsElementNamed("verbose") && Rcpp::as<bool>(control["verbose"]);
    return p;
}

std::vector<int> zero_based(const Graph& graph, const Rcpp::IntegerVector& blocks) {
    size_t N = boost::num_vertices(graph.adjacency);
    if (static_cast<size_t>(blocks.size()) != N)
        Rcpp::stop("blocks must have one entry per vertex");
    // relabel to 0..B-1 preserving first-appearance order
    std::unordered_map<int, int> map;
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
    std::unordered_map<int, int> map;
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

// Applies random single-vertex moves and returns, for each move, the
// difference between the incremental and the recomputed entropy change.
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
        int v = rng.index(state.N());
        int s = (rng.unif() < 0.2 && state.empty_group() >= 0)
                    ? state.empty_group() : state.random_occupied(rng);
        double dS = state.virtual_move(v, s);
        state.move(v, s);
        double S_new = state.entropy();
        err[i] = dS - (S_new - S);
        S = S_new;
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_sbm_minimize(SEXP pointer, bool degree_corrected, bool distributed,
                            Rcpp::List control, bool refine, double q_max) {
    Graph& graph = get_graph(pointer);
    Params p = read_params(control);
    SBM state(graph, degree_corrected, distributed, q_cache_max(q_max));
    FlatView view(state);
    Multilevel<FlatView> ml(view, p);
    std::vector<int> b(state.N(), 0);
    state.set_partition(b);
    double S = ml.stage_multilevel();
    double S_stage = S;
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
    int cap = std::max(N, 1);
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        std::vector<int> b(l == 0 ? N : cap, 0);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            b[i] = x[i] - 1;
        out.push_back(b);
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
        std::unordered_map<int, int> map;
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
    if (!refine)
        return;
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
    }
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
        int v = view.nodes()[rng.index(view.N())];
        int s;
        int empty = l == 0 ? n.base().empty_group() : n.upper(l).empty_group();
        if (rng.unif() < 0.2 && empty >= 0) {
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
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_minimize(SEXP pointer, bool degree_corrected, bool distributed,
                               Rcpp::List control, bool refine, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = boost::num_vertices(graph.adjacency);
    int L = N > 1 ? static_cast<int>(std::ceil(std::log2(double(N)))) + 1 : 1;
    Params p = read_params(control);
    double c = control.containsElementNamed("c") ? Rcpp::as<double>(control["c"]) : 1.0;
    Nested n(graph, degree_corrected, distributed, q_cache_max(q_max), L);
    Nested::Hierarchy bs(L, std::vector<int>(std::max(N, 1), 0));
    bs[0].resize(N);
    n.build(bs);
    std::vector<double> trace;
    if (N > 1)
        nested_minimize(n, p, c, refine, trace);
    std::vector<double> S = n.level_entropies();
    return Rcpp::List::create(Rcpp::Named("bs") = external_hierarchy(n.hierarchy()),
                              Rcpp::Named("level_entropy") = S,
                              Rcpp::Named("trace") = trace);
}
