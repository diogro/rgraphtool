// rgraphtool: shared dynamic-edge flat & nested SBM prior engine supporting
// exact O(1) edge multiplicity edits and O(deg(v)) vertex block moves.
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
// Copyright (C) 2026 Diogo Melo (R port)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#ifndef RGRAPHTOOL_DYN_SBM_H
#define RGRAPHTOOL_DYN_SBM_H

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>
#include "vendor/hash_map_wrap.hh"

namespace rgt_dyn {

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double LOG2 = 0.693147180559945309417232121458176568;

struct Rng {
    size_t index(size_t n) {
        return n <= 1 ? 0 : static_cast<size_t>(R_unif_index(static_cast<double>(n)));
    }
    double unif() { return unif_rand(); }
    template <class T>
    void shuffle(std::vector<T>& x) {
        for (size_t i = x.size(); i > 1; --i)
            std::swap(x[i - 1], x[index(i)]);
    }
};

inline uint64_t pair_key(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(b));
}

inline std::pair<int, int> unpack_key(uint64_t k) {
    return {static_cast<int>(k >> 32), static_cast<int>(k & 0xffffffffULL)};
}

inline double lgam(double x) { return std::lgamma(x); }

inline double lbinom_d(double n, double k) {
    if (k == 0.0 || n == k) return 0.0;
    if (n < k || k < 0.0) return 0.0;
    return lgam(n + 1.0) - lgam(k + 1.0) - lgam(n - k + 1.0);
}

inline double lbeta_d(double a, double b) {
    return lgam(a) + lgam(b) - lgam(a + b);
}

inline double log1mexp_d(double x) {
    if (x >= 0.0) return -INF;
    if (x < -LOG2) return std::log1p(-std::exp(x));
    return std::log(-std::expm1(x));
}

inline double log_sum_exp(double a, double b) {
    if (std::isinf(a) && a < 0) return b;
    if (std::isinf(b) && b < 0) return a;
    if (a < b) std::swap(a, b);
    return a + std::log1p(std::exp(b - a));
}

struct MathTab {
    static constexpr size_t MAX_TAB = 524288;
    std::vector<double> lgam_t, log_t;

    void ensure(size_t n) {
        size_t tgt = std::min(n + 1, MAX_TAB);
        if (tgt <= lgam_t.size()) return;
        size_t old = lgam_t.size();
        size_t cap = std::max<size_t>(4096, old);
        while (cap < tgt) cap = std::min(cap * 2, MAX_TAB);
        lgam_t.resize(cap);
        log_t.resize(cap);
        for (size_t i = old; i < cap; ++i) {
            lgam_t[i] = (i == 0) ? INF : std::lgamma(static_cast<double>(i));
            log_t[i] = (i == 0) ? -INF : std::log(static_cast<double>(i));
        }
    }
    inline double lgam1(int64_t x) const {
        if (x >= 0 && static_cast<size_t>(x + 1) < lgam_t.size())
            return lgam_t[static_cast<size_t>(x + 1)];
        return std::lgamma(static_cast<double>(x) + 1.0);
    }
    inline double logn(int64_t x) const {
        if (x >= 1 && static_cast<size_t>(x) < log_t.size())
            return log_t[static_cast<size_t>(x)];
        return std::log(static_cast<double>(x));
    }
    inline double lbinom(int64_t n, int64_t k) const {
        if (k == 0 || n == k) return 0.0;
        if (n < k || k < 0) return 0.0;
        return lgam1(n) - lgam1(k) - lgam1(n - k);
    }
};

inline MathTab& mtab() {
    static MathTab m;
    return m;
}

class LogQTab {
public:
    void ensure(size_t n_max) {
        if (n_max <= _n_max && !_table.empty()) return;
        _n_max = n_max;
        _table.assign(n_max * (n_max + 1) / 2 + 1, -INF);
        for (size_t n = 1; n <= n_max; ++n) {
            for (size_t k = 1; k <= n; ++k) {
                double a = (k == 1) ? -INF : get(n, k - 1);
                double b = get_raw(n - k, k);
                double hi = std::max(a, b), lo = std::min(a, b);
                at(n, k) = (hi == -INF) ? -INF : hi + std::log1p(std::exp(lo - hi));
            }
        }
    }
    double operator()(long n, long k) const {
        if (k > n) k = n;
        if (n == 0 && k == 0) return 0.0;
        if (n < 0 || k <= 0) return -INF;
        if (static_cast<size_t>(n) <= _n_max) return get(n, k);
        thread_local gt_hash_map<uint64_t, double> approx_cache;
        uint64_t key = (static_cast<uint64_t>(n) << 32) | static_cast<uint64_t>(k);
        auto it = approx_cache.find(key);
        if (it != approx_cache.end()) return it->second;
        double val = approx(static_cast<double>(n), static_cast<double>(k));
        if (approx_cache.size() >= 65536) approx_cache.clear();
        approx_cache.emplace(key, val);
        return val;
    }
    static double approx(double n, double k) {
        if (k < std::pow(n, 0.25))
            return lbinom_d(n - 1.0, k - 1.0) - lgam(k + 1.0);
        double u = k / std::sqrt(n);
        double v = u, delta = 1.0;
        while (delta > 1e-8) {
            double nv = u * std::sqrt(spence(std::exp(-v)));
            delta = std::abs(nv - v);
            v = nv;
        }
        double lf = std::log(v) - std::log1p(-std::exp(-v) * (1.0 + u * u / 2.0)) / 2.0 -
                    LOG2 * 1.5 - std::log(u) - std::log(M_PI);
        double g = 2.0 * v / u - u * std::log1p(-std::exp(-v));
        return lf - std::log(n) + std::sqrt(n) * g;
    }
    static double spence(double x) {
        if (x <= 0.0) return M_PI * M_PI / 6.0;
        if (x >= 1.0) return 0.0;
        if (x < 0.5) return M_PI * M_PI / 6.0 - std::log1p(-x) * std::log(x) - dilog(x);
        return dilog(1.0 - x);
    }
    static double dilog(double z) {
        if (z <= 0.0) return 0.0;
        if (z >= 1.0) return M_PI * M_PI / 6.0;
        if (z > 0.5) return M_PI * M_PI / 6.0 - std::log(z) * std::log1p(-z) - dilog(1.0 - z);
        double sum = 0.0, term = z;
        for (int k = 1; k < 200; ++k) {
            double add = term / (double(k) * k);
            sum += add;
            if (add < 1e-18 * sum) break;
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
        if (k > n) k = n;
        if (n == 0) return 0.0;
        return get(n, k);
    }
};

inline LogQTab& log_q() {
    static LogQTab q;
    return q;
}

struct Delta {
    int x, y;
    long d;
};

struct ChangeSet {
    std::vector<Delta> edges;
    std::vector<std::pair<int, int>> weights;
    void clear() { edges.clear(); weights.clear(); }
    bool empty() const { return edges.empty() && weights.empty(); }
};

class DeltaAcc {
public:
    void reset(bool directed) { _directed = directed; clear(); }
    void clear() { _v.clear(); _idx.clear(); }
    void add(int x, int y, long d) {
        if (!_directed && x > y) std::swap(x, y);
        uint64_t k = pair_key(x, y);
        auto it = _idx.find(k);
        if (it == _idx.end()) {
            _idx.emplace(k, _v.size());
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

class DynUpperLevel {
public:
    DynUpperLevel(int cap, bool directed, long E, bool top)
        : _cap(std::max(cap, 1)), _directed(directed), _top(top), _E(E),
          _b(_cap, 0), _wr(_cap, 0), _ep(_cap, 0), _em(_cap, 0),
          _m_out(_cap), _m_in(directed ? _cap : 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1), _members(_cap), _pos(_cap, 0) {
        _acc.reset(directed);
        _pairs.reset(directed);
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
    }

    int B() const { return static_cast<int>(_occupied.size()); }
    int N() const { return _N; }
    long E() const { return _E; }
    int group(int x) const { return _b[x]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<gt_hash_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<gt_hash_map<int, long>>& m_in() const { return _m_in; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    void set_E(long E) { _E = E; }
    void set_empty_label(int x, int g) { ensure_cap(std::max(x, g) + 1); _b[x] = g; }

    void ensure_cap(int need) {
        if (need <= _cap) return;
        int ncap = std::max(need, _cap * 2);
        int ocap = _cap;
        _cap = ncap;
        _b.resize(ncap, 0);
        _wr.resize(ncap, 0);
        _ep.resize(ncap, 0);
        _em.resize(ncap, 0);
        _m_out.resize(ncap);
        if (_directed) _m_in.resize(ncap);
        _occ_pos.resize(ncap, -1);
        _empty_pos.resize(ncap, -1);
        _members.resize(ncap);
        _pos.resize(ncap, 0);
        for (int r = ncap - 1; r >= ocap; --r) push_empty(r);
    }

    void build(const std::vector<int>& b, const std::vector<int>& locc,
               const std::vector<gt_hash_map<int, long>>& lm_out, long E_tot) {
        _E = E_tot;
        int max_id = 0;
        for (int x : locc) {
            if (x > max_id) max_id = x;
            int g = (x < static_cast<int>(b.size())) ? b[x] : 0;
            if (g > max_id) max_id = g;
        }
        ensure_cap(max_id + 1);
        for (int r : _occupied) {
            _wr[r] = 0; _ep[r] = 0; _em[r] = 0;
            _m_out[r].clear();
            if (_directed) _m_in[r].clear();
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        for (int x : locc) _b[x] = (x < static_cast<int>(b.size())) ? b[x] : 0;
        _N = 0;
        for (int x : locc) {
            int g = _b[x];
            if (_wr[g] == 0) take_empty(g);
            add_member(x, g);
            ++_N;
        }
        for (int x : locc) {
            int g = _b[x];
            for (const auto& [y, val] : lm_out[x]) {
                int h = _b[y];
                _m_out[g][h] += val;
                _ep[g] += val;
                if (_directed) {
                    _m_in[h][g] += val;
                    _em[h] += val;
                }
            }
        }
    }

    double partition_term(long Nv, long Bv) const {
        if (Nv == 0) return 0.0;
        const MathTab& M = mtab();
        return M.lbinom(Nv - 1, Bv - 1) + M.lgam1(Nv) + M.logn(Nv);
    }

    double edges_term(long Bv, long Ev) const {
        if (!_top || Ev <= 0 || Bv <= 0) return 0.0;
        int64_t pairs = _directed ? int64_t(Bv) * Bv : int64_t(Bv) * (Bv + 1) / 2;
        return mtab().lbinom(pairs + Ev - 1, Ev);
    }

    double dense_term(int x, int y, long m, long wx, long wy) const {
        if (m <= 0) return 0.0;
        int64_t nn, e = m;
        if (_directed || x != y) {
            nn = int64_t(wx) * wy;
        } else {
            nn = int64_t(wx) * (wx + 1) / 2;
            e = m >> 1;
        }
        return mtab().lbinom(nn + e - 1, e);
    }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0L : it->second;
    }

    double entropy() const {
        if (_N == 0) return 0.0;
        const MathTab& M = mtab();
        double S = partition_term(_N, B()) + edges_term(B(), _E);
        for (int r : _occupied) {
            S -= M.lgam1(_wr[r]);
            for (const auto& [t, m] : _m_out[r])
                if (_directed || t >= r)
                    S += dense_term(r, t, m, _wr[r], _wr[t]);
        }
        return S;
    }

    double edge_change(int gx, int gy, long dm_in, long dE, bool apply,
                       int& out_a, int& out_b, long& out_dm) {
        int a = _b[gx], b = _b[gy];
        long dm = dm_in;
        if (!_directed) {
            if (gx == gy) {
                a = b = _b[gx];
            } else if (a == b) {
                dm = 2 * dm_in;
            } else if (a > b) {
                std::swap(a, b);
            }
        }
        out_a = a; out_b = b; out_dm = dm;
        double dS = 0.0;
        if (dm != 0) {
            long m0 = get_m(a, b);
            dS += dense_term(a, b, m0 + dm, _wr[a], _wr[b]) -
                  dense_term(a, b, m0, _wr[a], _wr[b]);
        }
        if (_top && dE != 0) {
            dS += edges_term(B(), _E + dE) - edges_term(B(), _E);
        }
        if (apply) {
            if (dm != 0) {
                add_m(_m_out[a], b, dm);
                if (_directed) {
                    add_m(_m_in[b], a, dm);
                    _ep[a] += dm; _em[b] += dm;
                } else {
                    if (a != b) add_m(_m_out[b], a, dm);
                    _ep[a] += dm;
                    if (a != b) _ep[b] += dm;
                }
            }
            _E += dE;
        }
        return dS;
    }

    double lower_changes(const ChangeSet& in, ChangeSet& out, bool apply) {
        _acc.clear();
        for (const Delta& e : in.edges) {
            int gx = _b[e.x], gy = _b[e.y];
            if (_directed) _acc.add(gx, gy, e.d);
            else if (e.x == e.y) _acc.add(gx, gx, e.d);
            else if (gx == gy) _acc.add(gx, gx, 2 * e.d);
            else _acc.add(gx, gy, e.d);
        }
        _wd.clear();
        int dN = 0;
        for (const auto& [x, dw] : in.weights) {
            add_weight(_b[x], dw);
            dN += dw;
        }
        out.clear();
        for (const Delta& e : _acc.items())
            if (e.d != 0) out.edges.push_back(e);
        double dS = process(dN, out, apply);
        if (apply) {
            for (const auto& [x, dw] : in.weights) {
                if (dw > 0) push_node(x, _b[x]);
                else pop_node(x, _b[x]);
            }
        }
        return dS;
    }

    double node_move(int x, int h,
                     const std::vector<gt_hash_map<int, long>>& lm_out,
                     const std::vector<gt_hash_map<int, long>>& lm_in,
                     ChangeSet& out, bool apply) {
        int g = _b[x];
        out.clear();
        if (g == h) return 0.0;
        if (_directed) {
            gt_hash_map<int, long> co, ci;
            long lp = 0;
            for (const auto& [y, val] : lm_out[x]) {
                if (y == x) lp = val;
                else co[_b[y]] += val;
            }
            for (const auto& [y, val] : lm_in[x]) {
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
            gt_hash_map<int, long> co;
            long lp = 0;
            for (const auto& [y, val] : lm_out[x]) {
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

private:
    int _cap;
    bool _directed, _top;
    long _E;
    int _N = 0;
    std::vector<int> _b, _wr;
    std::vector<long> _ep, _em;
    std::vector<gt_hash_map<int, long>> _m_out, _m_in;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;
    DeltaAcc _acc, _pairs;
    std::vector<std::pair<int, int>> _wd;

    void add_weight(int g, int dw) {
        for (auto& p : _wd) if (p.first == g) { p.second += dw; return; }
        _wd.emplace_back(g, dw);
    }
    int new_weight(int g) const {
        int w = _wr[g];
        for (const auto& p : _wd) if (p.first == g) w += p.second;
        return w;
    }
    static void add_m(gt_hash_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) row.emplace(y, d);
        else { it->second += d; if (it->second == 0) row.erase(it); }
    }
    double process(int dN, ChangeSet& out, bool apply) {
        bool any_wd = false;
        long B0 = B(), B1 = B0;
        for (const auto& [g, dw] : _wd) {
            if (dw == 0) continue;
            any_wd = true;
            int w0 = _wr[g], w1 = w0 + dw;
            if (w0 == 0 && w1 > 0) { ++B1; out.weights.emplace_back(g, +1); }
            else if (w0 > 0 && w1 == 0) { --B1; out.weights.emplace_back(g, -1); }
        }
        if (dN == 0 && !any_wd && out.edges.empty()) return 0.0;
        const MathTab& M = mtab();
        double dS = 0.0;
        if (dN != 0 || B1 != B0) {
            dS += partition_term(_N + dN, B1) - partition_term(_N, B0);
            dS += edges_term(B1, _E) - edges_term(B0, _E);
        }
        for (const auto& [g, dw] : _wd)
            if (dw != 0) dS += M.lgam1(_wr[g]) - M.lgam1(_wr[g] + dw);
        _pairs.clear();
        for (const Delta& e : out.edges) _pairs.add(e.x, e.y, e.d);
        for (const auto& [g, dw] : _wd) {
            if (dw == 0) continue;
            for (const auto& [t, m] : _m_out[g]) _pairs.add(g, t, 0);
            if (_directed) for (const auto& [t, m] : _m_in[g]) _pairs.add(t, g, 0);
        }
        for (const Delta& p : _pairs.items()) {
            long m0 = get_m(p.x, p.y);
            dS += dense_term(p.x, p.y, m0 + p.d, new_weight(p.x), new_weight(p.y)) -
                  dense_term(p.x, p.y, m0, _wr[p.x], _wr[p.y]);
        }
        if (apply) {
            for (const Delta& e : out.edges) {
                add_m(_m_out[e.x], e.y, e.d);
                if (_directed) {
                    add_m(_m_in[e.y], e.x, e.d);
                    _ep[e.x] += e.d; _em[e.y] += e.d;
                } else {
                    if (e.x != e.y) add_m(_m_out[e.y], e.x, e.d);
                    _ep[e.x] += e.d;
                    if (e.x != e.y) _ep[e.y] += e.d;
                }
            }
            _N += dN;
        }
        return dS;
    }
    void push_empty(int r) { _empty_pos[r] = static_cast<int>(_empty.size()); _empty.push_back(r); }
    void take_empty(int s) {
        int i = _empty_pos[s];
        if (i >= 0) {
            int last = _empty.back();
            _empty[i] = last; _empty_pos[last] = i;
            _empty.pop_back(); _empty_pos[s] = -1;
        }
    }
    void remove_occupied(int r) {
        int i = _occ_pos[r], other = _occupied.back();
        _occupied[i] = other; _occ_pos[other] = i;
        _occupied.pop_back(); _occ_pos[r] = -1;
    }
    void add_member(int x, int r) {
        if (_wr[r] == 0) { _occ_pos[r] = static_cast<int>(_occupied.size()); _occupied.push_back(r); }
        _pos[x] = static_cast<int>(_members[r].size());
        _members[r].push_back(x);
        _wr[r]++;
    }
    void push_node(int x, int r) { if (_wr[r] == 0) take_empty(r); add_member(x, r); }
    void pop_node(int x, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[x]] = last; _pos[last] = _pos[x];
        vs.pop_back();
        if (--_wr[r] == 0) { remove_occupied(r); push_empty(r); }
    }
};

class DynSBM {
public:
    DynSBM(int N, bool directed, bool dc, bool dist, size_t q_max, bool edges_dl)
        : _N(N), _directed(directed), _dc(dc), _dist(dist), _edges_dl(edges_dl),
          _adj_out(N), _adj_in(directed ? N : 0),
          _kout(N, 0), _kin(N, 0), _loops(N, 0),
          _b(N, 0), _wr(std::max(N, 1), 0), _pos(N, 0),
          _occ_pos(std::max(N, 1), -1), _empty_pos(std::max(N, 1), -1),
          _ep(std::max(N, 1), 0), _em(std::max(N, 1), 0),
          _m_out(std::max(N, 1)), _m_in(directed ? std::max(N, 1) : 0),
          _h_out(std::max(N, 1)), _h_in(directed ? std::max(N, 1) : 0),
          _members(std::max(N, 1)) {
        mtab().ensure(static_cast<size_t>(std::max(N * 8, 4096)));
        if (dc && dist) log_q().ensure(std::min<size_t>(q_max, std::max<size_t>(N * 4, 256)));
        _empty.reserve(std::max(N, 1));
        for (int r = std::max(N, 1) - 1; r >= 0; --r) push_empty(r);
    }

    int N() const { return _N; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    bool directed() const { return _directed; }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<gt_hash_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<gt_hash_map<int, long>>& m_in() const { return _m_in; }
    const gt_hash_map<uint64_t, long>& edge_weights() const { return _w; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    uint64_t canon_key(int u, int v) const {
        if (!_directed && u > v) std::swap(u, v);
        return pair_key(u, v);
    }

    long get_w(int u, int v) const {
        auto it = _w.find(canon_key(u, v));
        return it == _w.end() ? 0L : it->second;
    }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0L : it->second;
    }

    void init_from_edges(const std::vector<std::tuple<int, int, long>>& edges,
                         const std::vector<int>& b_init) {
        _w.clear();
        for (int v = 0; v < _N; ++v) {
            _adj_out[v].clear();
            if (_directed) _adj_in[v].clear();
            _kout[v] = 0; _kin[v] = 0; _loops[v] = 0;
        }
        _E = 0;
        _const = 0.0;
        for (const auto& [u0, v0, m] : edges) {
            if (m <= 0) continue;
            int u = u0, v = v0;
            if (!_directed && u > v) std::swap(u, v);
            _w[pair_key(u, v)] += m;
        }
        for (const auto& [k, m] : _w) {
            auto [u, v] = unpack_key(k);
            _E += m;
            _adj_out[u][v] = m;
            if (_directed) {
                _adj_in[v][u] = m;
                _kout[u] += m; _kin[v] += m;
                if (u == v) _loops[u] += m;
            } else {
                if (u != v) _adj_out[v][u] = m;
                if (u == v) { _kout[u] += 2 * m; _loops[u] += m; }
                else { _kout[u] += m; _kout[v] += m; }
            }
        }
        mtab().ensure(static_cast<size_t>(std::max<int64_t>(2 * _E + _N + 64, 4096)));
        recompute_const();
        set_partition(b_init);
    }

    void set_partition(const std::vector<int>& b_in) {
        for (int r : _occupied) {
            _wr[r] = 0; _ep[r] = 0; _em[r] = 0;
            _m_out[r].clear();
            if (_directed) _m_in[r].clear();
            _h_out[r].clear();
            if (_directed) _h_in[r].clear();
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        for (int v = 0; v < _N; ++v) {
            int r = b_in[v];
            _b[v] = r;
            if (_wr[r] == 0) take_empty(r);
            add_member(v, r);
            _ep[r] += _kout[v];
            if (_directed) _em[r] += _kin[v];
            if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                if (_directed) _h_in[r][_kin[v]]++;
            }
        }
        for (const auto& [k, m] : _w) {
            auto [u, v] = unpack_key(k);
            int r = _b[u], s = _b[v];
            if (_directed) {
                _m_out[r][s] += m;
                _m_in[s][r] += m;
            } else if (r == s) {
                _m_out[r][r] += 2 * m;
            } else {
                _m_out[r][s] += m;
                _m_out[s][r] += m;
            }
        }
    }

    double entropy() const {
        if (_N == 0) return 0.0;
        const MathTab& M = mtab();
        double S = _const + global_term(B(), _E);
        for (int r : _occupied) {
            S += group_term(_wr[r], _ep[r], _em[r]);
            if (_dc && _dist) {
                for (const auto& [k, c] : _h_out[r]) S -= M.lgam1(c);
                if (_directed)
                    for (const auto& [k, c] : _h_in[r]) S -= M.lgam1(c);
            }
            for (const auto& [s, m] : _m_out[r]) {
                if (_directed || s >= r) S -= edge_term(r, s, m);
            }
        }
        return S;
    }

    double edge_delta(int u0, int v0, long dm, bool apply,
                      int& gx_out, int& gy_out, long& dm_m_out) {
        int u = u0, v = v0;
        if (!_directed && u > v) std::swap(u, v);
        uint64_t k = pair_key(u, v);
        auto wit = _w.find(k);
        long w0 = (wit == _w.end()) ? 0L : wit->second;
        long w1 = w0 + dm;
        if (w1 < 0) Rcpp::stop("edge multiplicity cannot become negative");
        if (dm == 0) { gx_out = _b[u]; gy_out = _b[v]; dm_m_out = 0; return 0.0; }

        mtab().ensure(static_cast<size_t>(std::max<int64_t>(2 * (_E + std::abs(dm)) + _N + 64, 4096)));
        const MathTab& M = mtab();

        double dS = 0.0;
        double dconst = M.lgam1(w1) - M.lgam1(w0);
        if (!_directed && u == v) dconst += static_cast<double>(dm) * LOG2;

        long kou0 = _kout[u], kov0 = _kout[v], kiv0 = _kin[v];
        long kou1 = kou0, kov1 = kov0, kiv1 = kiv0;
        if (_directed) {
            kou1 += dm; kiv1 += dm;
            if (u == v) { kov1 = kou1; }
            if (_dc) {
                dconst -= (M.lgam1(kou1) - M.lgam1(kou0)) + (M.lgam1(kiv1) - M.lgam1(kiv0));
            }
        } else {
            if (u == v) {
                kou1 += 2 * dm; kov1 = kou1;
                if (_dc) dconst -= M.lgam1(kou1) - M.lgam1(kou0);
            } else {
                kou1 += dm; kov1 += dm;
                if (_dc) {
                    dconst -= (M.lgam1(kou1) - M.lgam1(kou0)) + (M.lgam1(kov1) - M.lgam1(kov0));
                }
            }
        }
        dS += dconst;

        if (_edges_dl) {
            dS += global_term(B(), _E + dm) - global_term(B(), _E);
        }

        int r = _b[u], s = _b[v];
        int gx = r, gy = s;
        if (!_directed && gx > gy) std::swap(gx, gy);
        long dm_m = (!_directed && r == s) ? (2 * dm) : dm;
        gx_out = gx; gy_out = gy; dm_m_out = dm_m;
        long m0 = get_m(gx, gy);
        dS -= edge_term(gx, gy, m0 + dm_m) - edge_term(gx, gy, m0);

        if (r == s) {
            long dep = _directed ? dm : (2 * dm);
            long dem = _directed ? dm : 0;
            dS += group_term(_wr[r], _ep[r] + dep, _em[r] + dem) -
                  group_term(_wr[r], _ep[r], _em[r]);
            if (_dc && _dist) {
                if (_directed) {
                    dS += hist_change_1(_h_out[r], kou0, kou1);
                    dS += hist_change_1(_h_in[r], kiv0, kiv1);
                } else if (u == v) {
                    dS += hist_change_1(_h_out[r], kou0, kou1);
                } else {
                    dS += hist_change_2(_h_out[r], kou0, kou1, kov0, kov1);
                }
            }
        } else {
            long depr = dm, demr = 0;
            long deps = _directed ? 0 : dm, dems = _directed ? dm : 0;
            dS += group_term(_wr[r], _ep[r] + depr, _em[r] + demr) -
                  group_term(_wr[r], _ep[r], _em[r]);
            dS += group_term(_wr[s], _ep[s] + deps, _em[s] + dems) -
                  group_term(_wr[s], _ep[s], _em[s]);
            if (_dc && _dist) {
                dS += hist_change_1(_h_out[r], kou0, kou1);
                if (_directed) dS += hist_change_1(_h_in[s], kiv0, kiv1);
                else dS += hist_change_1(_h_out[s], kov0, kov1);
            }
        }

        if (apply) {
            _const += dconst;
            _E += dm;
            if (w1 == 0) {
                _w.erase(wit);
                _adj_out[u].erase(v);
                if (_directed) _adj_in[v].erase(u);
                else if (u != v) _adj_out[v].erase(u);
            } else {
                _w[k] = w1;
                _adj_out[u][v] = w1;
                if (_directed) _adj_in[v][u] = w1;
                else if (u != v) _adj_out[v][u] = w1;
            }
            if (u == v) _loops[u] += dm;
            if (_dc && _dist) {
                if (_directed) {
                    hist_apply(_h_out[r], kou0, kou1);
                    hist_apply(_h_in[s], kiv0, kiv1);
                } else if (u == v) {
                    hist_apply(_h_out[r], kou0, kou1);
                } else {
                    hist_apply(_h_out[r], kou0, kou1);
                    hist_apply(_h_out[s], kov0, kov1);
                }
            }
            _kout[u] = kou1;
            if (_directed) _kin[v] = kiv1;
            else if (u != v) _kout[v] = kov1;

            add_m(_m_out[gx], gy, dm_m);
            if (_directed) add_m(_m_in[gy], gx, dm_m);
            else if (gx != gy) add_m(_m_out[gy], gx, dm_m);

            if (r == s) {
                _ep[r] += (_directed ? dm : 2 * dm);
                if (_directed) _em[r] += dm;
            } else {
                _ep[r] += dm;
                if (_directed) _em[s] += dm;
                else _ep[s] += dm;
            }
        }
        return dS;
    }

    double move_vertex(int v, int s, bool apply, std::vector<Delta>& out_deltas) {
        int r = _b[v];
        out_deltas.clear();
        if (r == s) return 0.0;
        double dS = 0.0;
        int dB = (_wr[s] == 0 ? 1 : 0) - (_wr[r] == 1 ? 1 : 0);
        if (dB != 0) dS += global_term(B() + dB, _E) - global_term(B(), _E);
        long ko = _kout[v], ki = _kin[v];
        dS += group_term(_wr[r] - 1, _ep[r] - ko, _em[r] - ki) -
              group_term(_wr[r], _ep[r], _em[r]);
        dS += group_term(_wr[s] + 1, _ep[s] + ko, _em[s] + ki) -
              group_term(_wr[s], _ep[s], _em[s]);
        if (_dc && _dist) {
            dS += hist_move_delta(_h_out, r, s, ko);
            if (_directed) dS += hist_move_delta(_h_in, r, s, ki);
        }
        collect_vertex_deltas(v, r, s, out_deltas);
        for (const Delta& e : out_deltas) {
            long m = get_m(e.x, e.y);
            dS -= edge_term(e.x, e.y, m + e.d) - edge_term(e.x, e.y, m);
        }
        if (apply) {
            for (const Delta& e : out_deltas) {
                add_m(_m_out[e.x], e.y, e.d);
                if (_directed) add_m(_m_in[e.y], e.x, e.d);
                else if (e.x != e.y) add_m(_m_out[e.y], e.x, e.d);
            }
            _ep[r] -= ko; _em[r] -= ki;
            _ep[s] += ko; _em[s] += ki;
            if (_dc && _dist) {
                hist_move_apply(_h_out, r, s, ko);
                if (_directed) hist_move_apply(_h_in, r, s, ki);
            }
            remove_member(v, r);
            if (_wr[s] == 0) take_empty(s);
            add_member(v, s);
            _b[v] = s;
            if (_wr[r] == 0) push_empty(r);
        }
        return dS;
    }

    int random_occupied(Rng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    int sample_group(int v, double c, Rng& rng) const {
        if (_occupied.empty()) return 0;
        long kt = _kout[v] + (_directed ? _kin[v] : 0);
        if (std::isinf(c) || kt == 0) return random_occupied(rng);
        double u = rng.unif() * static_cast<double>(kt);
        int nbr = v;
        for (const auto& [w, m] : _adj_out[v]) {
            long wt = (!_directed && w == v) ? (2 * m) : m;
            u -= static_cast<double>(wt);
            if (u < 0.0) { nbr = w; break; }
        }
        if (u >= 0.0 && _directed) {
            for (const auto& [w, m] : _adj_in[v]) {
                u -= static_cast<double>(m);
                if (u < 0.0) { nbr = w; break; }
            }
        }
        int t = _b[nbr];
        long et = _ep[t] + (_directed ? _em[t] : 0);
        double Bo = static_cast<double>(_occupied.size());
        double p_rand = (c > 0.0) ? (c * Bo / (static_cast<double>(et) + c * Bo)) : 0.0;
        if (et <= 0 || (c > 0.0 && rng.unif() < p_rand)) return random_occupied(rng);
        double u2 = rng.unif() * static_cast<double>(et);
        for (const auto& [s, m] : _m_out[t]) {
            u2 -= static_cast<double>(m);
            if (u2 < 0.0) return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in[t]) {
                u2 -= static_cast<double>(m);
                if (u2 < 0.0) return s;
            }
        }
        return t;
    }

    double move_prob(int v, int s, double c) const {
        double Bo = static_cast<double>(_occupied.size());
        long kt = _kout[v] + (_directed ? _kin[v] : 0);
        if (kt == 0 || std::isinf(c)) return 1.0 / Bo;
        double prob = 0.0;
        auto add_nbr = [&](int nbr, long wt) {
            int t = _b[nbr];
            long et = _ep[t] + (_directed ? _em[t] : 0);
            long mts = get_m(t, s) + (_directed ? get_m(s, t) : 0);
            prob += (static_cast<double>(wt) / static_cast<double>(kt)) *
                    ((static_cast<double>(mts) + c) / (static_cast<double>(et) + c * Bo));
        };
        for (const auto& [w, m] : _adj_out[v]) {
            long wt = (!_directed && w == v) ? (2 * m) : m;
            add_nbr(w, wt);
        }
        if (_directed) {
            for (const auto& [w, m] : _adj_in[v]) add_nbr(w, m);
        }
        return prob;
    }

private:
    int _N;
    long _E = 0;
    bool _directed, _dc, _dist, _edges_dl;
    double _const = 0.0;
    gt_hash_map<uint64_t, long> _w;
    std::vector<gt_hash_map<int, long>> _adj_out, _adj_in;
    std::vector<long> _kout, _kin, _loops;

    std::vector<int> _b, _wr, _pos, _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<long> _ep, _em;
    std::vector<gt_hash_map<int, long>> _m_out, _m_in;
    std::vector<gt_hash_map<long, long>> _h_out, _h_in;
    std::vector<std::vector<int>> _members;

    void recompute_const() {
        const MathTab& M = mtab();
        double c = 0.0;
        for (const auto& [k, m] : _w) {
            auto [u, v] = unpack_key(k);
            c += M.lgam1(m);
            if (!_directed && u == v) c += static_cast<double>(m) * LOG2;
        }
        if (_dc) {
            for (int v = 0; v < _N; ++v) {
                c -= M.lgam1(_kout[v]);
                if (_directed) c -= M.lgam1(_kin[v]);
            }
        }
        _const = c;
    }

    double global_term(long Bv, long Ev) const {
        if (_N == 0) return 0.0;
        const MathTab& M = mtab();
        double S = M.lbinom(_N - 1, Bv - 1) + M.lgam1(_N) + M.logn(_N);
        if (_edges_dl && Ev > 0 && Bv > 0) {
            int64_t pairs = _directed ? int64_t(Bv) * Bv : int64_t(Bv) * (Bv + 1) / 2;
            S += M.lbinom(pairs + Ev - 1, Ev);
        }
        return S;
    }

    double group_term(long wr, long ep, long em) const {
        if (wr == 0) return 0.0;
        const MathTab& M = mtab();
        double S = -M.lgam1(wr);
        if (_dc) {
            S += M.lgam1(ep);
            if (_directed) S += M.lgam1(em);
            if (_dist) {
                const LogQTab& q = log_q();
                S += q(ep, wr) + M.lgam1(wr);
                if (_directed) S += q(em, wr) + M.lgam1(wr);
            } else {
                S += M.lbinom(wr + ep - 1, ep);
                if (_directed) S += M.lbinom(wr + em - 1, em);
            }
        } else {
            long deg = ep + (_directed ? em : 0);
            S += static_cast<double>(deg) * M.logn(wr);
        }
        return S;
    }

    double edge_term(int x, int y, long m) const {
        if (m == 0) return 0.0;
        const MathTab& M = mtab();
        if (_directed || x != y) return M.lgam1(m);
        long e = m >> 1;
        return M.lgam1(e) + static_cast<double>(e) * LOG2;
    }

    static long hget(const gt_hash_map<long, long>& h, long k) {
        auto it = h.find(k);
        return it == h.end() ? 0L : it->second;
    }

    static double hist_change_1(const gt_hash_map<long, long>& h, long k0, long k1) {
        if (k0 == k1) return 0.0;
        const MathTab& M = mtab();
        long c0 = hget(h, k0), c1 = hget(h, k1);
        return (M.lgam1(c0) - M.lgam1(c0 - 1)) - (M.lgam1(c1 + 1) - M.lgam1(c1));
    }

    static double hist_change_2(const gt_hash_map<long, long>& h,
                                long a0, long a1, long b0, long b1) {
        std::pair<long, int> ch[4] = {{a0, -1}, {a1, +1}, {b0, -1}, {b1, +1}};
        long keys[4]; int dlt[4]; int nk = 0;
        for (int i = 0; i < 4; ++i) {
            int f = -1;
            for (int j = 0; j < nk; ++j) if (keys[j] == ch[i].first) { f = j; break; }
            if (f >= 0) dlt[f] += ch[i].second;
            else { keys[nk] = ch[i].first; dlt[nk] = ch[i].second; ++nk; }
        }
        const MathTab& M = mtab();
        double dS = 0.0;
        for (int j = 0; j < nk; ++j) {
            if (dlt[j] == 0) continue;
            long c0 = hget(h, keys[j]);
            dS -= (M.lgam1(c0 + dlt[j]) - M.lgam1(c0));
        }
        return dS;
    }

    static void hist_apply(gt_hash_map<long, long>& h, long k0, long k1) {
        if (k0 == k1) return;
        auto it = h.find(k0);
        if (it != h.end() && --it->second == 0) h.erase(it);
        h[k1]++;
    }

    static double hist_move_delta(const std::vector<gt_hash_map<long, long>>& h,
                                  int r, int s, long k) {
        const MathTab& M = mtab();
        long cr = hget(h[r], k), cs = hget(h[s], k);
        return M.logn(cr) - M.logn(cs + 1);
    }

    static void hist_move_apply(std::vector<gt_hash_map<long, long>>& h,
                                int r, int s, long k) {
        auto it = h[r].find(k);
        if (it != h[r].end() && --it->second == 0) h[r].erase(it);
        h[s][k]++;
    }

    static void add_m(gt_hash_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) row.emplace(y, d);
        else { it->second += d; if (it->second == 0) row.erase(it); }
    }

    void collect_vertex_deltas(int v, int r, int s, std::vector<Delta>& deltas) const {
        deltas.clear();
        if (_directed) {
            gt_hash_map<int, long> co, ci;
            for (const auto& [u, m] : _adj_out[v]) if (u != v) co[_b[u]] += m;
            for (const auto& [u, m] : _adj_in[v])  if (u != v) ci[_b[u]] += m;
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k); return it == mp.end() ? 0L : it->second;
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
            gt_hash_map<int, long> co;
            for (const auto& [u, m] : _adj_out[v]) if (u != v) co[_b[u]] += m;
            auto get_c = [](const gt_hash_map<int, long>& mp, int k) -> long {
                auto it = mp.find(k); return it == mp.end() ? 0L : it->second;
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

    void push_empty(int r) { _empty_pos[r] = static_cast<int>(_empty.size()); _empty.push_back(r); }
    void take_empty(int s) {
        int i = _empty_pos[s];
        if (i >= 0) {
            int last = _empty.back();
            _empty[i] = last; _empty_pos[last] = i;
            _empty.pop_back(); _empty_pos[s] = -1;
        }
    }
    void remove_occupied(int r) {
        int i = _occ_pos[r], other = _occupied.back();
        _occupied[i] = other; _occ_pos[other] = i;
        _occupied.pop_back(); _occ_pos[r] = -1;
    }
    void add_member(int v, int r) {
        if (_wr[r] == 0) { _occ_pos[r] = static_cast<int>(_occupied.size()); _occupied.push_back(r); }
        _pos[v] = static_cast<int>(_members[r].size());
        _members[r].push_back(v);
        _wr[r]++;
    }
    void remove_member(int v, int r) {
        auto& vs = _members[r];
        int last = vs.back();
        vs[_pos[v]] = last; _pos[last] = _pos[v];
        vs.pop_back();
        if (--_wr[r] == 0) remove_occupied(r);
    }
};

class DynBlockPrior {
public:
    DynBlockPrior(int N, bool directed, bool dc, bool dist, size_t q_max, int levels)
        : _base(N, directed, dc, dist, q_max, levels <= 1) {
        _up.reserve(levels > 1 ? levels - 1 : 0);
        for (int l = 1; l < levels; ++l) {
            int init_cap = (l == 1) ? std::min(std::max(N, 1), 1024) : std::min(std::max(N, 1), 64);
            _up.emplace_back(init_cap, directed, 0, l == levels - 1);
        }
        _cs.resize(static_cast<size_t>(levels) + 1);
    }

    int levels() const { return static_cast<int>(_up.size()) + 1; }
    DynSBM& base() { return _base; }
    const DynSBM& base() const { return _base; }
    DynUpperLevel& upper(int l) { return _up[l - 1]; }
    const DynUpperLevel& upper(int l) const { return _up[l - 1]; }

    void init(const std::vector<std::tuple<int, int, long>>& edges,
              const std::vector<std::vector<int>>& bs) {
        _base.init_from_edges(edges, bs[0]);
        for (size_t l = 1; l < bs.size(); ++l) {
            const auto& locc = (l == 1) ? _base.occupied() : _up[l - 2].occupied();
            const auto& lm = (l == 1) ? _base.m_out() : _up[l - 2].m_out();
            _up[l - 1].build(bs[l], locc, lm, _base.E());
        }
    }

    double entropy() const {
        double S = _base.entropy();
        for (const auto& u : _up) S += u.entropy();
        return S;
    }

    double edge_delta(int u, int v, long dm, bool apply) {
        int gx = 0, gy = 0;
        long dm_m = 0;
        double dS = _base.edge_delta(u, v, dm, apply, gx, gy, dm_m);
        for (size_t l = 0; l < _up.size(); ++l) {
            int na = 0, nb = 0;
            long ndm = 0;
            dS += _up[l].edge_change(gx, gy, dm_m, dm, apply, na, nb, ndm);
            gx = na; gy = nb; dm_m = ndm;
        }
        return dS;
    }

    double vertex_move_dS(int v, int s) { return change_level0(v, s, false); }
    double vertex_move(int v, int s) { return change_level0(v, s, true); }

    void inherit(int l, int t, int r) {
        if (l + 1 < levels()) _up[l].set_empty_label(t, _up[l].group(r));
    }

    bool same_parent(int l, int r, int s) const {
        if (l + 1 >= levels()) return true;
        return _up[l].group(r) == _up[l].group(s);
    }

    double upper_move(int l, int x, int t, bool apply) {
        int r = _up[l - 1].group(x);
        if (r == t) return 0.0;
        bool empties = (_up[l - 1].size(r) == 1), fills = (_up[l - 1].size(t) == 0);
        const auto& lm_out = (l == 1) ? _base.m_out() : _up[l - 2].m_out();
        const auto& lm_in  = (l == 1) ? _base.m_in()  : _up[l - 2].m_in();
        double dS = _up[l - 1].node_move(x, t, lm_out, lm_in, _cs[l], apply);
        if (l + 1 < levels() && !empties && !fills && _up[l].group(r) == _up[l].group(t))
            return dS;
        for (int k = l + 1; k < levels(); ++k) {
            if (_cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k], apply);
        }
        return dS;
    }

    // Standard single-vertex + upper-level MCMC sweep over the current multigraph
    void mcmc_sweep(Rng& rng, double beta, double c, double d,
                    double& total_dS, size_t& nattempts, size_t& nmoves) {
        int N = _base.N();
        if (N <= 1) return;
        std::vector<int> vs(N);
        std::iota(vs.begin(), vs.end(), 0);
        rng.shuffle(vs);
        for (int v : vs) {
            ++nattempts;
            int r = _base.group(v);
            int B_cur = _base.B();
            bool can_create = (B_cur < N && d > 0.0 && _base.empty_group() >= 0);
            double p_new = can_create ? d : 0.0;
            bool propose_new = (can_create && rng.unif() < p_new);
            int s;
            double pfwd = 0.0;
            if (propose_new) {
                if (_base.size(r) == 1) continue;
                s = _base.empty_group();
                inherit(0, s, r);
                pfwd = p_new;
            } else {
                s = _base.sample_group(v, c, rng);
                if (s == r) continue;
                if (_base.size(r) == 1 && !same_parent(0, r, s)) continue;
                pfwd = (1.0 - p_new) * _base.move_prob(v, s, c);
            }
            bool empties = (_base.size(r) == 1);
            double dS = vertex_move_dS(v, s);
            if (std::isinf(beta)) {
                if (dS < 0.0) { vertex_move(v, s); total_dS += dS; ++nmoves; }
                continue;
            }
            if (pfwd <= 0.0) continue;
            vertex_move(v, s);
            double prev = 0.0;
            if (empties) {
                prev = (d > 0.0) ? d : 0.0;
            } else {
                bool rev_can = (_base.B() < N && d > 0.0);
                double p_new_rev = rev_can ? d : 0.0;
                prev = (1.0 - p_new_rev) * _base.move_prob(v, r, c);
            }
            bool accept = (prev > 0.0 &&
                           (rng.unif() < std::exp(-beta * dS + std::log(prev) - std::log(pfwd))));
            if (accept) { total_dS += dS; ++nmoves; }
            else vertex_move(v, r);
        }

        for (int l = 1; l < levels(); ++l) {
            const auto& nodes = (l == 1) ? _base.occupied() : _up[l - 2].occupied();
            if (nodes.size() <= 1) continue;
            std::vector<int> ord = nodes;
            rng.shuffle(ord);
            for (int x : ord) {
                int r = _up[l - 1].group(x);
                std::vector<int> cands = _up[l - 1].occupied();
                int emp = _up[l - 1].empty_group();
                if (emp >= 0 && _up[l - 1].size(r) > 1) cands.push_back(emp);
                if (cands.size() <= 1) continue;
                int t = cands[rng.index(cands.size())];
                if (t == r) continue;
                if (_up[l - 1].size(r) == 1 && !same_parent(l, r, t)) continue;
                if (_up[l - 1].size(t) == 0) inherit(l, t, r);
                ++nattempts;
                double pfwd = 1.0 / static_cast<double>(cands.size());
                double dS = upper_move(l, x, t, true);
                bool accept = false;
                if (std::isinf(beta)) accept = (dS < 0.0);
                else {
                    size_t rev_sz = _up[l - 1].occupied().size() +
                                    ((_up[l - 1].empty_group() >= 0 && _up[l - 1].size(t) > 1) ? 1 : 0);
                    double prev = rev_sz > 0 ? (1.0 / static_cast<double>(rev_sz)) : 0.0;
                    if (prev > 0.0) {
                        double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                        accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                    }
                }
                if (accept) { total_dS += dS; ++nmoves; }
                else upper_move(l, x, r, true);
            }
        }
    }

    Rcpp::List export_bs() const {
        Rcpp::List out(levels());
        std::vector<int> lower;
        for (int l = 0; l < levels(); ++l) {
            std::vector<int> items;
            if (l == 0) {
                items = _base.blocks();
            } else {
                const auto& b = _up[l - 1].blocks();
                items.reserve(lower.size());
                for (int g : lower) items.push_back(b[g]);
            }
            gt_hash_map<int, int> mp;
            std::vector<int> ord;
            Rcpp::IntegerVector x(items.size());
            for (size_t i = 0; i < items.size(); ++i) {
                auto it = mp.find(items[i]);
                if (it == mp.end()) {
                    it = mp.emplace(items[i], static_cast<int>(mp.size())).first;
                    ord.push_back(items[i]);
                }
                x[i] = it->second + 1;
            }
            out[l] = x;
            lower = std::move(ord);
        }
        return out;
    }

private:
    DynSBM _base;
    std::vector<DynUpperLevel> _up;
    std::vector<ChangeSet> _cs;

    double change_level0(int v, int s, bool apply) {
        int r = _base.group(v);
        if (r == s) return 0.0;
        bool empties = (_base.size(r) == 1), fills = (_base.size(s) == 0);
        _cs[0].clear();
        double dS = _base.move_vertex(v, s, apply, _cs[0].edges);
        if (_up.empty()) return dS;
        if (!empties && !fills && _up[0].group(r) == _up[0].group(s))
            return dS;
        if (empties) _cs[0].weights.emplace_back(r, -1);
        if (fills) _cs[0].weights.emplace_back(s, +1);
        for (int k = 1; k < levels(); ++k) {
            if (_cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k], apply);
        }
        return dS;
    }
};

}  // namespace rgt_dyn

#endif  // RGRAPHTOOL_DYN_SBM_H
