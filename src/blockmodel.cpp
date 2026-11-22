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
    SBM(const Graph& graph, bool dc, bool dist, size_t q_cache_max)
        : _directed(graph.directed), _dc(dc), _dist(dist) {
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
               lbinom(pairs + _E - 1, _E);
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
    bool _directed, _dc, _dist;
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

class Multilevel {
public:
    Multilevel(SBM& state, const Params& p) : _s(state), _p(p) {}

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
        std::vector<int> vs(N);
        std::iota(vs.begin(), vs.end(), 0);
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
        _cache[_s.B()] = {S, _s.blocks()};
        _trace_B.push_back(_s.B());
        _trace_S.push_back(S);
    }

    double get_cache(size_t B) {
        auto& c = _cache.at(B);
        _s.set_partition(c.second);
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
        std::vector<int> b(N);
        if (B_min == 1) {
            std::fill(b.begin(), b.end(), _s.group(_rng.index(N)));
            _s.set_partition(b);
            put_cache(_s.entropy());
        }
        std::iota(b.begin(), b.end(), 0);
        _s.set_partition(b);
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
    SBM& _s;
    Params _p;
    RRng _rng;
    std::map<size_t, std::pair<double, std::vector<int>>> _cache;
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
    Multilevel ml(state, p);
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
    Multilevel ml(state, p);
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
