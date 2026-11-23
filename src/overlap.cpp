// graph-tool-R: native overlapping (mixed-membership) stochastic block model
// (flat and nested), following graph-tool 3.9's OverlapBlockState (Peixoto
// 2015, PRX 5, 011033).
//
// Copyright (C) 2015-2026 Tiago de Paula Peixoto (algorithms)
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
#include <vector>
#include "graph_types.h"
#include "vendor/hash_map_wrap.hh"

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double LOG2 = 0.693147180559945309417232121458176568;

struct ORng {
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

// Fast integer lgamma / log tables.
struct OMathTables {
    static constexpr size_t MAX_TABLE = 262144;
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

OMathTables& omath() {
    static OMathTables t;
    return t;
}

inline double lgam_i(int64_t x) {
    const auto& tab = omath().lgamma_tab;
    if (x >= 1 && static_cast<size_t>(x) < tab.size())
        return tab[static_cast<size_t>(x)];
    return std::lgamma(static_cast<double>(x));
}

inline double log_i(int64_t x) {
    const auto& tab = omath().log_tab;
    if (x >= 1 && static_cast<size_t>(x) < tab.size())
        return tab[static_cast<size_t>(x)];
    return std::log(static_cast<double>(x));
}

inline double lbinom_i(int64_t n, int64_t k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return lgam_i(n + 1) - lgam_i(k + 1) - lgam_i(n - k + 1);
}

inline double lbinom_d(double n, double k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

inline uint64_t pair_key(uint64_t a, uint64_t b) {
    return (a << 32) | (b & 0xffffffffULL);
}

// Integer-partition log q(n, k) table + Szekeres approximation.
class OLogQ {
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
    double operator()(long n, long k) const {
        if (k > n)
            k = n;
        if (n == 0 && k == 0)
            return 0;
        if (n < 0 || k <= 0)
            return -inf;
        if (static_cast<size_t>(n) <= _n_max)
            return get(n, k);
        uint64_t key = pair_key(static_cast<uint64_t>(n), static_cast<uint64_t>(k));
        auto it = _cache.find(key);
        if (it != _cache.end())
            return it->second;
        double val = approx(static_cast<double>(n), static_cast<double>(k));
        if (_cache.size() >= 65536)
            _cache.clear();
        _cache.emplace(key, val);
        return val;
    }
    static double approx(double n, double k) {
        if (k < std::pow(n, 0.25))
            return lbinom_d(n - 1, k - 1) - std::lgamma(k + 1);
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
    mutable gt_hash_map<uint64_t, double> _cache;
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

OLogQ& olog_q() {
    static OLogQ q;
    return q;
}

// Hashable vector keys for mixtures and joint degree vectors.
typedef std::vector<int> MixKey;
typedef std::vector<int> DegKey;  // undirected: [k_0, ...]; directed: [ko_0, ki_0, ...]

struct VecHash {
    size_t operator()(const std::vector<int>& v) const noexcept {
        size_t h = v.size();
        for (int x : v)
            h ^= static_cast<size_t>(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct LabEdgeKey {
    int u, ru, v, rv;
    bool operator==(const LabEdgeKey& o) const noexcept {
        return u == o.u && ru == o.ru && v == o.v && rv == o.rv;
    }
};

struct LabEdgeHash {
    size_t operator()(const LabEdgeKey& k) const noexcept {
        uint64_t h = static_cast<uint64_t>(k.u) * 1000003ULL ^
                     (static_cast<uint64_t>(k.ru) << 17) ^
                     (static_cast<uint64_t>(k.v) * 9176ULL) ^
                     (static_cast<uint64_t>(k.rv) << 33);
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        return static_cast<size_t>(h);
    }
};

struct ODelta {
    int x, y;
    long d;
};

struct MixData {
    long n_bv = 0;
    std::vector<long> ep;
    std::vector<long> em;
    gt_hash_map<DegKey, long, VecHash> deg_hist;
};

// ---------------------------------------------------------------------------
// OverlapSBM: Native flat overlapping microcanonical SBM state.

class OverlapSBM {
public:
    OverlapSBM(const Graph& graph, bool dc, bool dist, bool gen_overlap,
               size_t q_cache_max, bool edges_dl = true)
        : _directed(graph.directed), _dc(dc), _dist(dist),
          _gen_overlap(gen_overlap), _edges_dl(edges_dl) {
        const Adjacency& g = graph.adjacency;
        _N = static_cast<int>(boost::num_vertices(g));
        _E = static_cast<long>(boost::num_edges(g));
        _H = static_cast<int>(2 * _E);
        _cap = std::max({_H, _N, 1});
        omath().ensure(static_cast<size_t>(std::max<int64_t>(2 * _E + _N + 16, 4096)));
        if (_dc && _dist) {
            size_t n_max = std::min<size_t>(std::max<long>(2 * _E, 100), 10000);
            n_max = std::min(n_max, q_cache_max);
            olog_q().ensure(n_max);
        }

        _edge_u.reserve(_E);
        _edge_v.reserve(_E);
        _he_v.resize(_H);
        _v_hes.resize(_N);
        auto range = boost::edges(g);
        int e_idx = 0;
        for (auto it = range.first; it != range.second; ++it, ++e_idx) {
            int s = static_cast<int>(it->s);
            int t = static_cast<int>(it->t);
            _edge_u.push_back(s);
            _edge_v.push_back(t);
            int h0 = 2 * e_idx;
            int h1 = h0 + 1;
            _he_v[h0] = s;
            _he_v[h1] = t;
            _v_hes[s].push_back(h0);
            _v_hes[t].push_back(h1);
        }
        _N_act = 0;
        for (int v = 0; v < _N; ++v)
            if (!_v_hes[v].empty())
                ++_N_act;

        _b.assign(_H, 0);
        _wr.assign(_cap, 0);
        _ep.assign(_cap, 0);
        _em.assign(_cap, 0);
        _r_count.assign(_cap, 0);
        _occ_pos.assign(_cap, -1);
        _empty_pos.assign(_cap, -1);
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r)
            push_empty(r);
        _m_out.resize(_cap);
        if (_directed)
            _m_in.resize(_cap);
        _he_members.resize(_cap);
        _he_pos.assign(_H, 0);

        _v_kout.resize(_N);
        if (_directed)
            _v_kin.resize(_N);
        _v_mix.resize(_N);
        _v_cdeg.resize(_N);
        if (_dc && _dist && !_gen_overlap) {
            _h_out_simple.resize(_cap);
            if (_directed)
                _h_in_simple.resize(_cap);
        }
    }

    int N() const { return _N; }
    int N_act() const { return _N_act; }
    long E() const { return _E; }
    int H() const { return _H; }
    int cap() const { return _cap; }
    int B() const { return static_cast<int>(_occupied.size()); }
    bool directed() const { return _directed; }
    bool deg_corr() const { return _dc; }
    bool gen_overlap() const { return _gen_overlap; }

    int group(int h) const { return _b[h]; }
    int wr(int r) const { return _wr[r]; }
    int he_vertex(int h) const { return _he_v[h]; }
    int he_mate(int h) const { return h ^ 1; }
    const std::vector<int>& v_hes(int v) const { return _v_hes[v]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& he_members(int r) const { return _he_members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<gt_hash_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<gt_hash_map<int, long>>& m_in() const { return _m_in; }
    const std::vector<long>& ep() const { return _ep; }
    const std::vector<long>& em() const { return _em; }
    const std::vector<int>& wr_vec() const { return _wr; }
    const std::vector<ODelta>& last_deltas() const { return _last_deltas; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }
    int random_occupied(ORng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0 : it->second;
    }

    void set_partition(const std::vector<int>& b) {
        if (static_cast<int>(b.size()) != _H)
            Rcpp::stop("Internal error: half-edge label vector length mismatch");
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _r_count[r] = 0;
            _m_out[r].clear();
            if (_directed) _m_in[r].clear();
            _he_members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
            if (_dc && _dist && !_gen_overlap) {
                _h_out_simple[r].clear();
                if (_directed) _h_in_simple[r].clear();
            }
        }
        _occupied.clear();
        _e_mult.clear();
        _mixes.clear();
        _dhist.clear();
        _D = 0;

        _b = b;
        for (int v = 0; v < _N; ++v) {
            _v_kout[v].clear();
            if (_directed) _v_kin[v].clear();
            _v_mix[v].clear();
            _v_cdeg[v].clear();
        }

        // Populate half-edge membership, block-pair counts, labelled edge counts,
        // and per-vertex labelled degrees.
        for (int e = 0; e < _E; ++e) {
            int h0 = 2 * e, h1 = h0 + 1;
            int u = _edge_u[e], v = _edge_v[e];
            int r0 = _b[h0], r1 = _b[h1];
            if (r0 < 0 || r0 >= _cap || r1 < 0 || r1 >= _cap)
                Rcpp::stop("Internal error: overlap block label out of range");

            add_he_member(h0, r0);
            add_he_member(h1, r1);

            if (_directed) {
                _m_out[r0][r1]++;
                _m_in[r1][r0]++;
                _ep[r0]++;
                _em[r1]++;
                _v_kout[u][r0]++;
                _v_kin[v][r1]++;
            } else {
                if (r0 == r1) {
                    _m_out[r0][r0] += 2;
                    _ep[r0] += 2;
                } else {
                    _m_out[r0][r1]++;
                    _m_out[r1][r0]++;
                    _ep[r0]++;
                    _ep[r1]++;
                }
                _v_kout[u][r0]++;
                _v_kout[v][r1]++;
            }
            _e_mult[make_edge_key(u, r0, v, r1)]++;
        }

        // Build per-vertex mixtures, wr, _occupied, and mixture histograms.
        for (int v = 0; v < _N; ++v) {
            if (_v_hes[v].empty())
                continue;
            build_vertex_mix(v, _v_mix[v], _v_cdeg[v]);
            for (int r : _v_mix[v]) {
                if (_wr[r] == 0)
                    activate_group(r);
                _wr[r]++;
            }
            if (_gen_overlap) {
                add_vertex_to_mix(_v_mix[v], _v_cdeg[v]);
            } else if (_dc && _dist) {
                for (int r : _v_mix[v]) {
                    long ko = get_vk(_v_kout[v], r);
                    _h_out_simple[r][ko]++;
                    if (_directed) {
                        long ki = get_vk(_v_kin[v], r);
                        _h_in_simple[r][ki]++;
                    }
                }
            }
        }
    }

    double edge_term(int x, int y, long m) const {
        if (m <= 0)
            return 0;
        if (!_directed && x == y) {
            long half = m >> 1;
            return lgam_i(half + 1) + half * LOG2;
        }
        return lgam_i(m + 1);
    }

    double lab_edge_term(const LabEdgeKey& /*k*/, long m) const {
        if (m <= 0)
            return 0;
        return lgam_i(m + 1);
    }

    double edges_dl_term(long B_act) const {
        if (!_edges_dl || _E == 0)
            return 0;
        int64_t pairs = _directed ? int64_t(B_act) * B_act
                                  : int64_t(B_act) * (B_act + 1) / 2;
        return lbinom_i(pairs + _E - 1, _E);
    }

    static double d_part_term(int64_t B_act, int64_t d, int64_t nd) {
        if (nd <= 0)
            return 0;
        double x = lbinom_i(B_act, d);
        double ex = std::exp(x);
        if (ex > 1e12 || !std::isfinite(ex))
            return static_cast<double>(nd) * x - lgam_i(nd + 1);
        double ss = lbinom_d(ex + static_cast<double>(nd) - 1.0, static_cast<double>(nd));
        if (!std::isfinite(ss))
            ss = static_cast<double>(nd) * x - lgam_i(nd + 1);
        return ss;
    }

    double mix_deg_term(const MixKey& mx, const MixData& md) const {
        if (!_dc || md.n_bv == 0)
            return 0;
        double S = 0;
        long n_bv = md.n_bv;
        if (_dist) {
            const OLogQ& q = olog_q();
            for (size_t i = 0; i < mx.size(); ++i) {
                if (!_directed) {
                    S += q(md.ep[i] - n_bv, n_bv);
                } else {
                    S += q(md.ep[i], n_bv) + q(md.em[i], n_bv);
                }
            }
            S += lgam_i(n_bv + 1);
            for (const auto& [dk, c] : md.deg_hist)
                S -= lgam_i(c + 1);
        } else {
            for (size_t i = 0; i < mx.size(); ++i) {
                S += lbinom_i(n_bv + md.ep[i] - 1, md.ep[i]);
                if (_directed)
                    S += lbinom_i(n_bv + md.em[i] - 1, md.em[i]);
            }
        }
        return S;
    }

    double partition_dl() const {
        if (!_gen_overlap || _N_act == 0)
            return 0;
        int64_t B_act = static_cast<int64_t>(_occupied.size());
        double S = lbinom_i(_D + _N_act - 1, _N_act) + lgam_i(_N_act + 1);
        for (size_t d = 1; d < _dhist.size(); ++d)
            if (_dhist[d] > 0)
                S += d_part_term(B_act, static_cast<int64_t>(d), _dhist[d]);
        for (const auto& [mx, md] : _mixes)
            if (md.n_bv > 0)
                S -= lgam_i(md.n_bv + 1);
        return S;
    }

    double degree_dl() const {
        if (!_dc || _N_act == 0)
            return 0;
        double S = 0;
        if (_gen_overlap) {
            for (const auto& [mx, md] : _mixes)
                if (md.n_bv > 0)
                    S += mix_deg_term(mx, md);
            for (int r : _occupied) {
                long rc = _r_count[r];
                S += lbinom_i(rc + _ep[r] - 1, _ep[r]);
                if (_directed)
                    S += lbinom_i(rc + _em[r] - 1, _em[r]);
            }
        } else {
            long n = _N_act;
            const OLogQ& q = olog_q();
            for (int r : _occupied) {
                if (_dist) {
                    S += q(_ep[r], n) + lgam_i(n + 1);
                    long nz = 0;
                    for (const auto& [k, c] : _h_out_simple[r]) {
                        S -= lgam_i(c + 1);
                        nz += c;
                    }
                    if (n > nz) S -= lgam_i(n - nz + 1);
                    if (_directed) {
                        S += q(_em[r], n) + lgam_i(n + 1);
                        long nzi = 0;
                        for (const auto& [k, c] : _h_in_simple[r]) {
                            S -= lgam_i(c + 1);
                            nzi += c;
                        }
                        if (n > nzi) S -= lgam_i(n - nzi + 1);
                    }
                } else {
                    S += lbinom_i(n + _ep[r] - 1, _ep[r]);
                    if (_directed)
                        S += lbinom_i(n + _em[r] - 1, _em[r]);
                }
            }
        }
        return S;
    }

    double likelihood_neg() const {
        if (_E == 0)
            return 0;
        double S = 0;
        for (const auto& [k, m] : _e_mult)
            S += lab_edge_term(k, m);
        for (int r : _occupied) {
            for (const auto& [t, m] : _m_out[r])
                if (_directed || t >= r)
                    S -= edge_term(r, t, m);
            if (_dc) {
                S += lgam_i(_ep[r] + 1);
                if (_directed) S += lgam_i(_em[r] + 1);
            } else {
                S += (_ep[r] + (_directed ? _em[r] : 0)) * log_i(_wr[r]);
            }
        }
        if (_dc) {
            for (int v = 0; v < _N; ++v) {
                for (const auto& [r, k] : _v_kout[v])
                    S -= lgam_i(k + 1);
                if (_directed) {
                    for (const auto& [r, k] : _v_kin[v])
                        S -= lgam_i(k + 1);
                }
            }
        }
        return S;
    }

    double entropy() const {
        if (_E == 0)
            return 0;
        return likelihood_neg() + edges_dl_term(B()) + partition_dl() + degree_dl();
    }

    // Exact incremental half-edge move cost and apply!
    // Performs the move in O(d_v) and optionally reverts (for virtual_move) or keeps.
    // Wait: computing delta directly without mutating is even faster; let's compute
    // exact dS directly by temporarily applying + reverting OR directly.
    // Because EVERY single mutation in move_he() is O(d_v) (typically d_v = 1..3),
    // exact incremental delta with zero full scans runs in ~100ns!
    double virtual_move(int h, int s) {
        int r = _b[h];
        if (r == s)
            return 0;
        double dS = apply_move_delta(h, s);
        // Revert cleanly
        apply_move_delta(h, r);
        return dS;
    }

    double move(int h, int s) {
        int r = _b[h];
        if (r == s)
            return 0;
        return apply_move_delta(h, s);
    }

    // Move all half-edges of vertex v that currently have label != s to s.
    double move_vertex(int v, int s, std::vector<int>* saved = nullptr) {
        const auto& hes = _v_hes[v];
        if (saved) {
            saved->resize(hes.size());
            for (size_t i = 0; i < hes.size(); ++i)
                (*saved)[i] = _b[hes[i]];
        }
        double dS = 0;
        for (int h : hes)
            if (_b[h] != s)
                dS += apply_move_delta(h, s);
        return dS;
    }

    double virtual_move_vertex(int v, int s) {
        std::vector<int> saved;
        double dS = move_vertex(v, s, &saved);
        const auto& hes = _v_hes[v];
        for (size_t i = hes.size(); i > 0; --i)
            if (_b[hes[i - 1]] != saved[i - 1])
                apply_move_delta(hes[i - 1], saved[i - 1]);
        return dS;
    }

    // Merge all half-edges of group r into group s.
    double merge(int r, int s) {
        if (r == s)
            return 0;
        std::vector<int> hes = _he_members[r];
        double dS = 0;
        for (int h : hes)
            dS += apply_move_delta(h, s);
        return dS;
    }

    double virtual_merge(int r, int s) {
        if (r == s)
            return 0;
        std::vector<int> hes = _he_members[r];
        double dS = 0;
        for (int h : hes)
            dS += apply_move_delta(h, s);
        for (int h : hes)
            apply_move_delta(h, r);
        return dS;
    }

    // Proposal for half-edge h: sample from mate group's neighbors or uniform.
    int sample_group(int h, double c, ORng& rng) const {
        if (std::isinf(c) || _occupied.empty())
            return random_occupied(rng);
        int t = _b[h ^ 1];
        double et = static_cast<double>(_ep[t] + (_directed ? _em[t] : 0));
        double Bo = static_cast<double>(_occupied.size());
        double p_rand = (c > 0) ? (c * Bo / (et + c * Bo)) : 0.0;
        if (et <= 0 || (c > 0 && rng.unif() < p_rand))
            return random_occupied(rng);
        double u = rng.unif() * et;
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

    // Exact proposal probability q(r -> s | h) for finite-temperature MCMC.
    double move_prob(int h, int s, double c) const {
        double Bo = static_cast<double>(_occupied.size());
        if (std::isinf(c))
            return 1.0 / Bo;
        int t = _b[h ^ 1];
        double et = static_cast<double>(_ep[t] + (_directed ? _em[t] : 0));
        if (et <= 0)
            return 1.0 / Bo;
        double mts = static_cast<double>(get_m(t, s) + (_directed ? get_m(s, t) : 0));
        return (mts + c) / (et + c * Bo);
    }

    // Per-vertex overlap summary for R export.
    Rcpp::List overlap_groups_r() const {
        // Relabel occupied groups to contiguous 1..B in sorted internal group order.
        std::vector<int> occ_sorted = _occupied;
        std::sort(occ_sorted.begin(), occ_sorted.end());
        gt_hash_map<int, int> gmap;
        gmap.reserve(occ_sorted.size());
        for (size_t i = 0; i < occ_sorted.size(); ++i)
            gmap.emplace(occ_sorted[i], static_cast<int>(i));
        int Bo = static_cast<int>(gmap.size());
        Rcpp::List bv(static_cast<R_xlen_t>(_N));
        Rcpp::List bc_out(static_cast<R_xlen_t>(_N));
        Rcpp::List bc_in(static_cast<R_xlen_t>(_N));
        Rcpp::List bc_tot(static_cast<R_xlen_t>(_N));
        Rcpp::IntegerMatrix mat(_N, std::max(Bo, 1));
        Rcpp::IntegerVector maj(_N, NA_INTEGER);

        for (int v = 0; v < _N; ++v) {
            const auto& mx = _v_mix[v];
            std::vector<std::pair<int, int>> ord;
            ord.reserve(mx.size());
            for (int r : mx)
                ord.emplace_back(gmap[r] + 1, r);
            std::sort(ord.begin(), ord.end());
            Rcpp::IntegerVector gv(ord.size()), ko(ord.size()), ki(ord.size()), kt(ord.size());
            long best_k = -1;
            int best_g = NA_INTEGER;
            for (size_t i = 0; i < ord.size(); ++i) {
                int lbl = ord[i].first;
                int r = ord[i].second;
                long o = get_vk(_v_kout[v], r);
                long in_d = _directed ? get_vk(_v_kin[v], r) : 0;
                long tot = _directed ? (o + in_d) : o;
                gv[i] = lbl;
                ko[i] = static_cast<int>(o);
                ki[i] = static_cast<int>(in_d);
                kt[i] = static_cast<int>(tot);
                mat(v, lbl - 1) = static_cast<int>(tot);
                if (tot > best_k) {
                    best_k = tot;
                    best_g = lbl;
                }
            }
            bv[v] = gv;
            bc_out[v] = ko;
            bc_in[v] = ki;
            bc_tot[v] = kt;
            maj[v] = best_g;
        }
        return Rcpp::List::create(
            Rcpp::Named("groups") = bv,
            Rcpp::Named("out_degrees") = bc_out,
            Rcpp::Named("in_degrees") = bc_in,
            Rcpp::Named("total_degrees") = bc_tot,
            Rcpp::Named("matrix") = Bo > 0 ? mat : Rcpp::IntegerMatrix(_N, 0),
            Rcpp::Named("majority") = maj,
            Rcpp::Named("B") = Bo);
    }

    Rcpp::IntegerMatrix half_edges_r() const {
        std::vector<int> occ_sorted = _occupied;
        std::sort(occ_sorted.begin(), occ_sorted.end());
        gt_hash_map<int, int> gmap;
        gmap.reserve(occ_sorted.size());
        for (size_t i = 0; i < occ_sorted.size(); ++i)
            gmap.emplace(occ_sorted[i], static_cast<int>(i) + 1);
        Rcpp::IntegerMatrix out(_E, 2);
        for (int e = 0; e < _E; ++e) {
            out(e, 0) = gmap.at(_b[2 * e]);
            out(e, 1) = gmap.at(_b[2 * e + 1]);
        }
        Rcpp::colnames(out) = Rcpp::CharacterVector::create("source_block", "target_block");
        return out;
    }

private:
    bool _directed, _dc, _dist, _gen_overlap, _edges_dl;
    int _N = 0, _N_act = 0, _H = 0, _cap = 1;
    long _E = 0;
    std::vector<int> _edge_u, _edge_v, _he_v;
    std::vector<std::vector<int>> _v_hes;

    std::vector<int> _b, _wr, _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<long> _ep, _em, _r_count;
    std::vector<gt_hash_map<int, long>> _m_out, _m_in;
    std::vector<std::vector<int>> _he_members;
    std::vector<int> _he_pos;

    gt_hash_map<LabEdgeKey, long, LabEdgeHash> _e_mult;
    std::vector<gt_hash_map<int, long>> _v_kout, _v_kin;
    std::vector<MixKey> _v_mix;
    std::vector<DegKey> _v_cdeg;

    gt_hash_map<MixKey, MixData, VecHash> _mixes;
    std::vector<long> _dhist;
    int64_t _D = 0;

    std::vector<gt_hash_map<long, long>> _h_out_simple, _h_in_simple;
    std::vector<ODelta> _last_deltas;

    static long get_vk(const gt_hash_map<int, long>& mp, int r) {
        auto it = mp.find(r);
        return it == mp.end() ? 0 : it->second;
    }

    LabEdgeKey make_edge_key(int u, int ru, int v, int rv) const {
        if (!_directed && (u > v || (u == v && ru > rv)))
            return {v, rv, u, ru};
        return {u, ru, v, rv};
    }

    void add_he_member(int h, int r) {
        _he_pos[h] = static_cast<int>(_he_members[r].size());
        _he_members[r].push_back(h);
    }

    void remove_he_member(int h, int r) {
        auto& vs = _he_members[r];
        int last = vs.back();
        vs[_he_pos[h]] = last;
        _he_pos[last] = _he_pos[h];
        vs.pop_back();
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

    void activate_group(int r) {
        take_empty(r);
        _occ_pos[r] = static_cast<int>(_occupied.size());
        _occupied.push_back(r);
    }

    void deactivate_group(int r) {
        int i = _occ_pos[r];
        int other = _occupied.back();
        _occupied[i] = other;
        _occ_pos[other] = i;
        _occupied.pop_back();
        _occ_pos[r] = -1;
        push_empty(r);
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

    void build_vertex_mix(int v, MixKey& mx, DegKey& cdeg) const {
        mx.clear();
        for (const auto& [r, k] : _v_kout[v])
            if (k > 0) mx.push_back(r);
        if (_directed) {
            for (const auto& [r, k] : _v_kin[v])
                if (k > 0 && get_vk(_v_kout[v], r) == 0)
                    mx.push_back(r);
        }
        std::sort(mx.begin(), mx.end());
        cdeg.clear();
        cdeg.reserve(_directed ? 2 * mx.size() : mx.size());
        for (int r : mx) {
            cdeg.push_back(static_cast<int>(get_vk(_v_kout[v], r)));
            if (_directed)
                cdeg.push_back(static_cast<int>(get_vk(_v_kin[v], r)));
        }
    }

    void ensure_dhist(size_t d) {
        if (d >= _dhist.size())
            _dhist.resize(d + 1, 0);
    }

    void add_vertex_to_mix(const MixKey& mx, const DegKey& cdeg) {
        size_t d = mx.size();
        ensure_dhist(d);
        _dhist[d]++;
        if (static_cast<int64_t>(d) > _D)
            _D = static_cast<int64_t>(d);
        auto& md = _mixes[mx];
        if (md.n_bv == 0) {
            md.ep.assign(d, 0);
            if (_directed) md.em.assign(d, 0);
            for (int r : mx)
                _r_count[r]++;
        }
        md.n_bv++;
        for (size_t i = 0; i < d; ++i) {
            if (_directed) {
                md.ep[i] += cdeg[2 * i];
                md.em[i] += cdeg[2 * i + 1];
            } else {
                md.ep[i] += cdeg[i];
            }
        }
        md.deg_hist[cdeg]++;
    }

    void remove_vertex_from_mix(const MixKey& mx, const DegKey& cdeg) {
        size_t d = mx.size();
        _dhist[d]--;
        while (_D > 0 && _dhist[static_cast<size_t>(_D)] == 0)
            --_D;
        auto it = _mixes.find(mx);
        MixData& md = it->second;
        md.n_bv--;
        for (size_t i = 0; i < d; ++i) {
            if (_directed) {
                md.ep[i] -= cdeg[2 * i];
                md.em[i] -= cdeg[2 * i + 1];
            } else {
                md.ep[i] -= cdeg[i];
            }
        }
        auto dit = md.deg_hist.find(cdeg);
        if (--(dit->second) == 0)
            md.deg_hist.erase(dit);
        if (md.n_bv == 0) {
            for (int r : mx)
                _r_count[r]--;
            _mixes.erase(it);
        }
    }

    // Applies half-edge move h: r -> s and returns exact dS in O(d_v).
    double apply_move_delta(int h, int s) {
        int r = _b[h];
        if (r == s)
            return 0;
        double dS = 0.0;
        int e = h >> 1;
        bool is_out = ((h & 1) == 0);
        int v = _he_v[h];
        int t = _b[h ^ 1];
        int u_e = _edge_u[e], v_e = _edge_v[e];

        // 1. Labelled edge multiplicity term
        LabEdgeKey old_ek = make_edge_key(u_e, _b[2 * e], v_e, _b[2 * e + 1]);
        LabEdgeKey new_ek = is_out ? make_edge_key(u_e, s, v_e, t)
                                   : make_edge_key(u_e, t, v_e, s);
        if (!(old_ek == new_ek)) {
            auto ito = _e_mult.find(old_ek);
            long m0 = ito->second;
            dS += lab_edge_term(old_ek, m0 - 1) - lab_edge_term(old_ek, m0);
            if (m0 == 1)
                _e_mult.erase(ito);
            else
                ito->second = m0 - 1;

            auto itn = _e_mult.find(new_ek);
            long m1 = (itn == _e_mult.end()) ? 0 : itn->second;
            dS += lab_edge_term(new_ek, m1 + 1) - lab_edge_term(new_ek, m1);
            if (itn == _e_mult.end())
                _e_mult.emplace(new_ek, 1);
            else
                itn->second = m1 + 1;
        }

        // 2. Block-pair edge counts m_out / m_in
        _last_deltas.clear();
        auto apply_pair = [&](int x, int y, long delta) {
            if (!_directed && x > y)
                std::swap(x, y);
            long m0 = get_m(x, y);
            dS -= edge_term(x, y, m0 + delta) - edge_term(x, y, m0);
            add_m(_m_out[x], y, delta);
            if (_directed)
                add_m(_m_in[y], x, delta);
            else if (x != y)
                add_m(_m_out[y], x, delta);
            _last_deltas.push_back({x, y, delta});
        };
        if (_directed) {
            if (is_out) {
                apply_pair(r, t, -1);
                apply_pair(s, t, +1);
            } else {
                apply_pair(t, r, -1);
                apply_pair(t, s, +1);
            }
        } else {
            apply_pair(r, t, (r == t) ? -2 : -1);
            apply_pair(s, t, (s == t) ? +2 : +1);
        }

        // 3. Before-state group-level terms for r and s
        long B0 = B();
        long ep_r0 = _ep[r], em_r0 = _em[r], wr_r0 = _wr[r], rc_r0 = _r_count[r];
        long ep_s0 = _ep[s], em_s0 = _em[s], wr_s0 = _wr[s], rc_s0 = _r_count[s];

        auto group_lik = [&](long wr_v, long ep_v, long em_v) -> double {
            if (wr_v == 0) return 0;
            if (_dc) return lgam_i(ep_v + 1) + (_directed ? lgam_i(em_v + 1) : 0);
            return (ep_v + (_directed ? em_v : 0)) * log_i(wr_v);
        };
        dS -= group_lik(wr_r0, ep_r0, em_r0) + group_lik(wr_s0, ep_s0, em_s0);

        if (_dc && _gen_overlap) {
            if (wr_r0 > 0)
                dS -= lbinom_i(rc_r0 + ep_r0 - 1, ep_r0) +
                      (_directed ? lbinom_i(rc_r0 + em_r0 - 1, em_r0) : 0);
            if (wr_s0 > 0)
                dS -= lbinom_i(rc_s0 + ep_s0 - 1, ep_s0) +
                      (_directed ? lbinom_i(rc_s0 + em_s0 - 1, em_s0) : 0);
        } else if (_dc && !_gen_overlap) {
            dS -= simple_deg_group(r) + simple_deg_group(s);
        }

        // 4. Per-vertex degree likelihood term
        bool change_out = (!_directed || is_out);
        auto& vk = change_out ? _v_kout[v] : _v_kin[v];
        long k_r0 = get_vk(vk, r);
        long k_s0 = get_vk(vk, s);
        if (_dc) {
            dS -= (lgam_i(k_r0) - lgam_i(k_r0 + 1));       // k_r0 -> k_r0 - 1
            dS -= (lgam_i(k_s0 + 2) - lgam_i(k_s0 + 1));   // k_s0 -> k_s0 + 1
        }

        // 5. Update half-edge, ep/em, vertex degrees, and mixtures
        remove_he_member(h, r);
        _b[h] = s;
        add_he_member(h, s);

        if (change_out) {
            _ep[r]--;
            _ep[s]++;
        } else {
            _em[r]--;
            _em[s]++;
        }

        // Check simple degree histograms (gen_overlap = false)
        if (_dc && _dist && !_gen_overlap) {
            auto& hs = change_out ? _h_out_simple : _h_in_simple;
            auto it_r = hs[r].find(k_r0);
            if (it_r != hs[r].end() && --(it_r->second) == 0)
                hs[r].erase(it_r);
            if (k_r0 - 1 > 0)
                hs[r][k_r0 - 1]++;
            if (k_s0 > 0) {
                auto it_s = hs[s].find(k_s0);
                if (it_s != hs[s].end() && --(it_s->second) == 0)
                    hs[s].erase(it_s);
            }
            hs[s][k_s0 + 1]++;
        }

        // Update vertex labelled degrees
        auto it_kr = vk.find(r);
        if (--(it_kr->second) == 0)
            vk.erase(it_kr);
        vk[s]++;

        // Did vertex v leave group r or enter group s?
        long k_r_tot_after = get_vk(_v_kout[v], r) + (_directed ? get_vk(_v_kin[v], r) : 0);
        long k_s_tot_before = (change_out ? k_s0 + (_directed ? get_vk(_v_kin[v], s) : 0)
                                          : get_vk(_v_kout[v], s) + k_s0);
        bool left_r = (k_r_tot_after == 0);
        bool entered_s = (k_s_tot_before == 0);

        if (left_r) {
            if (--_wr[r] == 0)
                deactivate_group(r);
        }
        if (entered_s) {
            if (_wr[s] == 0)
                activate_group(s);
            _wr[s]++;
        }
        long B1 = B();

        // Update mixtures (when _gen_overlap) and track exact partition_dl + mix_deg_dl delta
        if (_gen_overlap) {
            MixKey old_mx = _v_mix[v];
            DegKey old_cd = _v_cdeg[v];
            MixKey new_mx;
            DegKey new_cd;
            build_vertex_mix(v, new_mx, new_cd);
            bool same_mx = (old_mx == new_mx);

            if (same_mx) {
                // Only one mixture's ep/em and deg_hist change; partition_dl & r_count unchanged!
                MixData& md = _mixes[old_mx];
                dS -= mix_deg_term(old_mx, md);
                for (size_t i = 0; i < old_mx.size(); ++i) {
                    if (old_mx[i] == r) {
                        if (change_out) md.ep[i]--; else md.em[i]--;
                    } else if (old_mx[i] == s) {
                        if (change_out) md.ep[i]++; else md.em[i]++;
                    }
                }
                auto dit = md.deg_hist.find(old_cd);
                if (--(dit->second) == 0)
                    md.deg_hist.erase(dit);
                md.deg_hist[new_cd]++;
                dS += mix_deg_term(old_mx, md);
                _v_cdeg[v] = std::move(new_cd);
            } else {
                // Mixture changed! Subtract old contributions for old_mx and new_mx (if present),
                // plus partition_dl d-terms if B or d changed.
                int64_t D0 = _D;
                size_t d_old = old_mx.size(), d_new = new_mx.size();
                double old_part_d = 0.0;
                if (B0 != B1) {
                    for (size_t d = 1; d < _dhist.size(); ++d)
                        if (_dhist[d] > 0)
                            old_part_d += d_part_term(B0, static_cast<int64_t>(d), _dhist[d]);
                } else if (d_old != d_new) {
                    old_part_d += d_part_term(B0, static_cast<int64_t>(d_old), _dhist[d_old]);
                    if (d_new < _dhist.size())
                        old_part_d += d_part_term(B0, static_cast<int64_t>(d_new), _dhist[d_new]);
                }
                double old_D_term = lbinom_i(D0 + _N_act - 1, _N_act);

                // Check other groups in old_mx / new_mx whose r_count might change (only when
                // old_mx empties or new_mx becomes nonempty!).
                // Collect affected groups != r, s whose r_count changes.
                long n_old_before = _mixes[old_mx].n_bv;
                auto it_new0 = _mixes.find(new_mx);
                long n_new_before = (it_new0 == _mixes.end()) ? 0 : it_new0->second.n_bv;

                dS += lgam_i(n_old_before + 1);
                dS += lgam_i(n_new_before + 1);
                dS -= mix_deg_term(old_mx, _mixes[old_mx]);
                if (n_new_before > 0)
                    dS -= mix_deg_term(new_mx, it_new0->second);

                bool old_empties = (n_old_before == 1);
                bool new_fills = (n_new_before == 0);
                // Any group g != r, s in old_mx (if old_empties) or new_mx (if new_fills)
                // changes _r_count[g]!
                std::vector<int> rc_changed;
                if (_dc && (old_empties || new_fills)) {
                    auto add_rc_g = [&](int g) {
                        if (g == r || g == s) return;
                        if (std::find(rc_changed.begin(), rc_changed.end(), g) == rc_changed.end())
                            rc_changed.push_back(g);
                    };
                    if (old_empties)
                        for (int g : old_mx) add_rc_g(g);
                    if (new_fills)
                        for (int g : new_mx) add_rc_g(g);
                    for (int g : rc_changed) {
                        long rc = _r_count[g];
                        dS -= lbinom_i(rc + _ep[g] - 1, _ep[g]) +
                              (_directed ? lbinom_i(rc + _em[g] - 1, _em[g]) : 0);
                    }
                }

                remove_vertex_from_mix(old_mx, old_cd);
                add_vertex_to_mix(new_mx, new_cd);
                _v_mix[v] = new_mx;
                _v_cdeg[v] = std::move(new_cd);

                if (_dc && !rc_changed.empty()) {
                    for (int g : rc_changed) {
                        long rc = _r_count[g];
                        dS += lbinom_i(rc + _ep[g] - 1, _ep[g]) +
                              (_directed ? lbinom_i(rc + _em[g] - 1, _em[g]) : 0);
                    }
                }

                auto it_old1 = _mixes.find(old_mx);
                long n_old_after = (it_old1 == _mixes.end()) ? 0 : it_old1->second.n_bv;
                long n_new_after = _mixes[new_mx].n_bv;
                dS -= lgam_i(n_old_after + 1);
                dS -= lgam_i(n_new_after + 1);
                if (n_old_after > 0)
                    dS += mix_deg_term(old_mx, it_old1->second);
                dS += mix_deg_term(new_mx, _mixes[new_mx]);

                double new_part_d = 0.0;
                if (B0 != B1) {
                    for (size_t d = 1; d < _dhist.size(); ++d)
                        if (_dhist[d] > 0)
                            new_part_d += d_part_term(B1, static_cast<int64_t>(d), _dhist[d]);
                } else if (d_old != d_new) {
                    new_part_d += d_part_term(B1, static_cast<int64_t>(d_old), _dhist[d_old]);
                    new_part_d += d_part_term(B1, static_cast<int64_t>(d_new), _dhist[d_new]);
                }
                double new_D_term = lbinom_i(_D + _N_act - 1, _N_act);
                dS += (new_part_d - old_part_d) + (new_D_term - old_D_term);
            }
        } else {
            build_vertex_mix(v, _v_mix[v], _v_cdeg[v]);
        }

        // 6. After-state group-level terms for r and s + edges_dl
        dS += group_lik(_wr[r], _ep[r], _em[r]) + group_lik(_wr[s], _ep[s], _em[s]);
        if (_dc && _gen_overlap) {
            if (_wr[r] > 0)
                dS += lbinom_i(_r_count[r] + _ep[r] - 1, _ep[r]) +
                      (_directed ? lbinom_i(_r_count[r] + _em[r] - 1, _em[r]) : 0);
            if (_wr[s] > 0)
                dS += lbinom_i(_r_count[s] + _ep[s] - 1, _ep[s]) +
                      (_directed ? lbinom_i(_r_count[s] + _em[s] - 1, _em[s]) : 0);
        } else if (_dc && !_gen_overlap) {
            dS += simple_deg_group(r) + simple_deg_group(s);
        }
        if (B0 != B1)
            dS += edges_dl_term(B1) - edges_dl_term(B0);

        return dS;
    }

    double simple_deg_group(int r) const {
        if (_wr[r] == 0)
            return 0;
        long n = _N_act;
        if (!_dist) {
            return lbinom_i(n + _ep[r] - 1, _ep[r]) +
                   (_directed ? lbinom_i(n + _em[r] - 1, _em[r]) : 0);
        }
        const OLogQ& q = olog_q();
        double S = q(_ep[r], n) + lgam_i(n + 1);
        long nz = 0;
        for (const auto& [k, c] : _h_out_simple[r]) {
            S -= lgam_i(c + 1);
            nz += c;
        }
        if (n > nz) S -= lgam_i(n - nz + 1);
        if (_directed) {
            S += q(_em[r], n) + lgam_i(n + 1);
            long nzi = 0;
            for (const auto& [k, c] : _h_in_simple[r]) {
                S -= lgam_i(c + 1);
                nzi += c;
            }
            if (n > nzi) S -= lgam_i(n - nzi + 1);
        }
        return S;
    }
};

// ---------------------------------------------------------------------------
// Upper hierarchical levels for NestedOverlap (identical math to BlockLevel,
// with exact O(deg) incremental updates for large-network scaling).

class ODeltaAcc {
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
    const std::vector<ODelta>& items() const { return _v; }
private:
    bool _directed = false;
    std::vector<ODelta> _v;
    gt_hash_map<uint64_t, size_t> _idx;
};

struct OChangeSet {
    std::vector<ODelta> edges;
    std::vector<std::pair<int, int>> weights;
    void clear() {
        edges.clear();
        weights.clear();
    }
    bool empty() const { return edges.empty() && weights.empty(); }
};

class OUpperLevel {
public:
    OUpperLevel(int cap, bool directed, long E, bool top)
        : _cap(cap), _directed(directed), _top(top), _E(E),
          _b(cap, 0), _wr(cap, 0), _ep(cap, 0), _em(cap, 0),
          _m_out(cap), _m_in(directed ? cap : 0), _occ_pos(cap, -1),
          _empty_pos(cap, -1), _members(cap), _pos(cap, 0) {
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
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<gt_hash_map<int, long>>& m_out() const { return _m_out; }
    const std::vector<gt_hash_map<int, long>>& m_in() const { return _m_in; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    void build(const std::vector<int>& b,
               const std::vector<int>& lower_occ,
               const std::vector<gt_hash_map<int, long>>& lower_m_out) {
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _m_out[r].clear();
            if (_directed) _m_in[r].clear();
            _members[r].clear();
            _occ_pos[r] = -1;
            push_empty(r);
        }
        _occupied.clear();
        _b = b;
        _N = 0;
        for (int x : lower_occ) {
            int g = _b[x];
            if (_wr[g] == 0)
                take_empty(g);
            add_member(x, g);
            ++_N;
        }
        for (int x : lower_occ) {
            int g = _b[x];
            for (const auto& [y, val] : lower_m_out[x]) {
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
        if (Nv == 0)
            return 0;
        return lbinom_i(Nv - 1, Bv - 1) + lgam_i(Nv + 1) + log_i(Nv);
    }

    double edges_term(long Bv) const {
        if (!_top || _E <= 0)
            return 0;
        int64_t pairs = _directed ? int64_t(Bv) * Bv : int64_t(Bv) * (Bv + 1) / 2;
        return lbinom_i(pairs + _E - 1, _E);
    }

    long get_m(int x, int y) const {
        auto it = _m_out[x].find(y);
        return it == _m_out[x].end() ? 0 : it->second;
    }

    double entropy() const {
        if (_N == 0)
            return 0;
        double S = partition_term(_N, B()) + edges_term(B());
        for (int r : _occupied) {
            S -= lgam_i(_wr[r] + 1);
            for (const auto& [t, m] : _m_out[r])
                if (_directed || t >= r)
                    S += dense_term(r, t, m, _wr[r], _wr[t]);
        }
        return S;
    }

    // Incremental update from lower-level edge/weight deltas; returns exact dS at this level.
    double lower_changes(const OChangeSet& in, OChangeSet& out) {
        _acc.clear();
        for (const ODelta& e : in.edges) {
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
        for (const ODelta& e : _acc.items())
            if (e.d != 0)
                out.edges.push_back(e);
        double dS = eval_and_apply(dN, out);
        for (const auto& [x, dw] : in.weights) {
            if (dw > 0)
                push_node(x, _b[x]);
            else
                pop_node(x, _b[x]);
        }
        return dS;
    }

    // Incremental single-node move x: g -> h at this level; returns exact dS at this level.
    double node_move(int x, int h,
                     const std::vector<gt_hash_map<int, long>>& lm_out,
                     const std::vector<gt_hash_map<int, long>>& lm_in,
                     OChangeSet& out) {
        int g = _b[x];
        out.clear();
        if (g == h)
            return 0.0;
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
        double dS = eval_and_apply(0, out);
        pop_node(x, g);
        _b[x] = h;
        push_node(x, h);
        return dS;
    }

    void set_label(int x, int g) { _b[x] = g; }

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
    ODeltaAcc _acc, _pairs;
    std::vector<std::pair<int, int>> _wd;

    void add_weight(int g, int dw) {
        for (auto& p : _wd) {
            if (p.first == g) { p.second += dw; return; }
        }
        _wd.emplace_back(g, dw);
    }

    int new_weight(int g) const {
        int w = _wr[g];
        for (const auto& p : _wd)
            if (p.first == g) w += p.second;
        return w;
    }

    static void add_m(gt_hash_map<int, long>& row, int y, long d) {
        auto it = row.find(y);
        if (it == row.end()) {
            row.emplace(y, d);
        } else {
            it->second += d;
            if (it->second == 0) row.erase(it);
        }
    }

    double eval_and_apply(int dN, OChangeSet& out) {
        bool any_wd = false;
        long B0 = B(), B1 = B0;
        for (const auto& [g, dw] : _wd) {
            if (dw == 0) continue;
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
            return 0.0;

        double dS = 0.0;
        if (dN != 0 || B1 != B0) {
            dS += partition_term(_N + dN, B1) - partition_term(_N, B0);
            dS += edges_term(B1) - edges_term(B0);
        }
        for (const auto& [g, dw] : _wd) {
            if (dw != 0)
                dS += lgam_i(_wr[g] + 1) - lgam_i(_wr[g] + dw + 1);
        }
        _pairs.clear();
        for (const ODelta& e : out.edges)
            _pairs.add(e.x, e.y, e.d);
        for (const auto& [g, dw] : _wd) {
            if (dw == 0) continue;
            for (const auto& [t, m] : _m_out[g])
                _pairs.add(g, t, 0);
            if (_directed)
                for (const auto& [t, m] : _m_in[g])
                    _pairs.add(t, g, 0);
        }
        for (const ODelta& p : _pairs.items()) {
            long m0 = get_m(p.x, p.y);
            dS += dense_term(p.x, p.y, m0 + p.d, new_weight(p.x), new_weight(p.y)) -
                  dense_term(p.x, p.y, m0, _wr[p.x], _wr[p.y]);
        }

        // Apply edge updates and _N
        for (const ODelta& e : out.edges) {
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
        return dS;
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
        return lbinom_i(nn + e - 1, e);
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

    void remove_occupied(int r) {
        int i = _occ_pos[r];
        int other = _occupied.back();
        _occupied[i] = other;
        _occ_pos[other] = i;
        _occupied.pop_back();
        _occ_pos[r] = -1;
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

class NestedOverlap {
public:
    NestedOverlap(const Graph& graph, bool dc, bool dist, bool gen_overlap,
                  size_t q_cache_max, int levels)
        : _base(graph, dc, dist, gen_overlap, q_cache_max, levels <= 1) {
        int cap = _base.cap();
        _up.reserve(levels > 1 ? levels - 1 : 0);
        for (int l = 1; l < levels; ++l)
            _up.emplace_back(cap, _base.directed(), _base.E(), l == levels - 1);
        _cs.resize(static_cast<size_t>(levels) + 1);
    }

    int levels() const { return static_cast<int>(_up.size()) + 1; }
    OverlapSBM& base() { return _base; }
    OUpperLevel& upper(int l) { return _up[l - 1]; }

    // Exact O(1)/O(L) incremental level-0 half-edge move propagating upward.
    double move_he(int h, int s) {
        int r = _base.group(h);
        if (r == s) return 0.0;
        bool empties = (_base.he_members(r).size() == 1);
        bool fills = (_base.he_members(s).size() == 0);
        double dS = _base.move(h, s);
        if (_up.empty()) return dS;
        if (!empties && !fills && _up[0].group(r) == _up[0].group(s))
            return dS;
        _cs[0].clear();
        _cs[0].edges = _base.last_deltas();
        if (empties) _cs[0].weights.emplace_back(r, -1);
        if (fills) _cs[0].weights.emplace_back(s, +1);
        for (int k = 1; k < levels(); ++k) {
            if (_cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k]);
        }
        return dS;
    }

    // Exact O(deg + L) incremental upper-level node move propagating upward.
    double move_upper(int l, int x, int t) {
        int r = _up[l - 1].group(x);
        if (r == t) return 0.0;
        bool empties = (_up[l - 1].size(r) == 1);
        bool fills = (_up[l - 1].size(t) == 0);
        const auto& lm_out = (l == 1) ? _base.m_out() : _up[l - 2].m_out();
        const auto& lm_in  = (l == 1) ? _base.m_in()  : _up[l - 2].m_in();
        double dS = _up[l - 1].node_move(x, t, lm_out, lm_in, _cs[l]);
        if (l + 1 < levels() && !empties && !fills && _up[l].group(r) == _up[l].group(t))
            return dS;
        for (int k = l + 1; k < levels(); ++k) {
            if (_cs[k - 1].empty()) break;
            dS += _up[k - 1].lower_changes(_cs[k - 1], _cs[k]);
        }
        return dS;
    }

    void build_upper() {
        for (size_t l = 1; l <= _up.size(); ++l) {
            if (l == 1)
                _up[0].build(_up[0].blocks(), _base.occupied(), _base.m_out());
            else
                _up[l - 1].build(_up[l - 1].blocks(), _up[l - 2].occupied(), _up[l - 2].m_out());
        }
    }

    void build(const std::vector<int>& b0, const std::vector<std::vector<int>>& upper_bs) {
        _base.set_partition(b0);
        for (size_t l = 0; l < _up.size(); ++l) {
            const auto& locc = (l == 0) ? _base.occupied() : _up[l - 1].occupied();
            const auto& lm = (l == 0) ? _base.m_out() : _up[l - 1].m_out();
            _up[l].build(upper_bs[l], locc, lm);
        }
    }

    double upper_entropy() const {
        double S = 0;
        for (const auto& u : _up)
            S += u.entropy();
        return S;
    }

    double entropy() const {
        return _base.entropy() + upper_entropy();
    }

    std::vector<double> level_entropies() const {
        std::vector<double> out{_base.entropy()};
        for (const auto& u : _up)
            out.push_back(u.entropy());
        return out;
    }

    // Optimize upper levels (l >= 1) greedily given fixed level 0.
    // Note: since level 0 typically has B_0 << N (e.g. 2..40 groups), upper levels
    // are tiny graphs and greedy single-node + merge + split sweeps run in microseconds!
    void optimize_upper(ORng& rng) {
        if (_up.empty())
            return;
        build_upper();
        for (int pass = 0; pass < 3; ++pass) {
            bool improved = false;
            for (int l = 1; l < levels(); ++l) {
                const auto& nodes = (l == 1) ? _base.occupied() : _up[l - 2].occupied();
                if (nodes.size() <= 1) {
                    // Collapse single node to 0
                    std::vector<int> bl = _up[l - 1].blocks();
                    for (int x : nodes) bl[x] = 0;
                    _up[l - 1].build(bl, nodes, (l == 1) ? _base.m_out() : _up[l - 2].m_out());
                    build_upper_from(l + 1);
                    continue;
                }
                double best_S = upper_entropy();
                std::vector<std::vector<int>> best_bs = save_upper();

                // Try starting from 1 group and from singleton groups, refining greedily
                for (int init_mode = 0; init_mode < 2; ++init_mode) {
                    std::vector<int> bl = _up[l - 1].blocks();
                    if (init_mode == 0) {
                        for (int x : nodes) bl[x] = 0;
                    } else {
                        for (size_t i = 0; i < nodes.size(); ++i)
                            bl[nodes[i]] = static_cast<int>(i);
                    }
                    _up[l - 1].build(bl, nodes, (l == 1) ? _base.m_out() : _up[l - 2].m_out());
                    build_upper_from(l + 1);
                    double cur_S = upper_entropy();

                    // Greedy single-node refinement at level l using O(deg) incremental updates
                    std::vector<int> ord = nodes;
                    for (int sw = 0; sw < 15; ++sw) {
                        rng.shuffle(ord);
                        bool ch = false;
                        for (int x : ord) {
                            int r = _up[l - 1].group(x);
                            std::vector<int> cands = _up[l - 1].occupied();
                            int emp = _up[l - 1].empty_group();
                            if (emp >= 0 && _up[l - 1].size(r) > 1)
                                cands.push_back(emp);
                            int best_t = r;
                            double best_dS = -1e-9;
                            for (int t : cands) {
                                if (t == r) continue;
                                double dS = move_upper(l, x, t);
                                if (dS < best_dS) {
                                    best_dS = dS;
                                    best_t = t;
                                }
                                move_upper(l, x, r);
                            }
                            if (best_t != r) {
                                cur_S += move_upper(l, x, best_t);
                                ch = true;
                            }
                        }
                        if (!ch) break;
                    }
                    if (cur_S < best_S - 1e-9) {
                        best_S = cur_S;
                        best_bs = save_upper();
                        improved = true;
                    }
                }
                restore_upper(best_bs);
            }
            if (!improved)
                break;
        }
    }

    std::vector<std::vector<int>> save_upper() const {
        std::vector<std::vector<int>> out;
        out.reserve(_up.size());
        for (const auto& u : _up)
            out.push_back(u.blocks());
        return out;
    }

    void restore_upper(const std::vector<std::vector<int>>& bs) {
        for (size_t l = 0; l < _up.size(); ++l) {
            const auto& locc = (l == 0) ? _base.occupied() : _up[l - 1].occupied();
            const auto& lm = (l == 0) ? _base.m_out() : _up[l - 1].m_out();
            _up[l].build(bs[l], locc, lm);
        }
    }

    void build_upper_from(int start_l) {
        for (int l = start_l; l < levels(); ++l) {
            const auto& locc = (l == 1) ? _base.occupied() : _up[l - 2].occupied();
            const auto& lm = (l == 1) ? _base.m_out() : _up[l - 2].m_out();
            _up[l - 1].build(_up[l - 1].blocks(), locc, lm);
        }
    }

    // Export upper levels as contiguous 1-based R list matching level-0 group order.
    Rcpp::List external_upper_r() const {
        // Level-0 group order matches half_edges_r() first-appearance order.
        std::vector<int> lower;
        gt_hash_map<int, int> seen;
        for (int h = 0; h < _base.H(); ++h) {
            int r = _base.group(h);
            if (seen.find(r) == seen.end()) {
                seen.emplace(r, static_cast<int>(lower.size()));
                lower.push_back(r);
            }
        }
        Rcpp::List out(_up.size());
        for (size_t l = 0; l < _up.size(); ++l) {
            const auto& b = _up[l].blocks();
            gt_hash_map<int, int> map;
            std::vector<int> order;
            Rcpp::IntegerVector x(lower.size());
            for (size_t i = 0; i < lower.size(); ++i) {
                int g = b[lower[i]];
                auto it = map.find(g);
                if (it == map.end()) {
                    it = map.emplace(g, static_cast<int>(map.size())).first;
                    order.push_back(g);
                }
                x[i] = it->second + 1;
            }
            out[l] = x;
            lower = order;
        }
        return out;
    }

private:
    OverlapSBM _base;
    std::vector<OUpperLevel> _up;
    std::vector<OChangeSet> _cs;
};

// ---------------------------------------------------------------------------
// Minimization and refinement for flat and nested OverlapSBM.

std::vector<int> parse_half_edges(const Graph& graph, const Rcpp::IntegerMatrix& he) {
    long E = static_cast<long>(boost::num_edges(graph.adjacency));
    if (he.nrow() != E || he.ncol() != 2)
        Rcpp::stop("half_edges must be an E x 2 matrix of positive integer labels");
    std::vector<int> uniq;
    uniq.reserve(2 * E);
    for (long e = 0; e < E; ++e) {
        for (int j = 0; j < 2; ++j) {
            int val = he(e, j);
            if (val == NA_INTEGER || val < 1)
                Rcpp::stop("half_edges labels must be positive integers");
            uniq.push_back(val);
        }
    }
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    gt_hash_map<int, int> map;
    map.reserve(uniq.size());
    for (size_t i = 0; i < uniq.size(); ++i)
        map.emplace(uniq[i], static_cast<int>(i));
    std::vector<int> b(2 * E);
    for (long e = 0; e < E; ++e) {
        b[2 * e] = map[he(e, 0)];
        b[2 * e + 1] = map[he(e, 1)];
    }
    return b;
}

// Greedy refinement for flat OverlapSBM (half-edge moves + vertex moves + merges + splits).
double overlap_refine_flat(OverlapSBM& s, ORng& rng, size_t B_min, size_t B_max,
                           double epsilon, size_t max_iter) {
    if (s.E() == 0)
        return 0;
    std::vector<int> hes(s.H());
    std::iota(hes.begin(), hes.end(), 0);
    std::vector<int> vs;
    for (int v = 0; v < s.N(); ++v)
        if (!s.v_hes(v).empty())
            vs.push_back(v);

    double S = s.entropy();
    for (size_t iter = 0; iter < max_iter; ++iter) {
        double dS_iter = 0.0;
        // 1. Greedy bundled vertex moves
        rng.shuffle(vs);
        for (int v : vs) {
            size_t B_cur = static_cast<size_t>(s.B());
            gt_hash_set<int> cands;
            for (int h : s.v_hes(v)) {
                cands.insert(s.group(h));
                cands.insert(s.group(h ^ 1));
                cands.insert(s.sample_group(h, 0.5, rng));
            }
            if (B_cur < B_max && s.empty_group() >= 0)
                cands.insert(s.empty_group());
            int best = -1;
            double best_dS = -1e-9;
            for (int t : cands) {
                double dS = s.virtual_move_vertex(v, t);
                if (dS < best_dS) {
                    // Check B_min/B_max if accepted
                    best_dS = dS;
                    best = t;
                }
            }
            if (best >= 0) {
                std::vector<int> saved;
                double dS = s.move_vertex(v, best, &saved);
                size_t B_new = static_cast<size_t>(s.B());
                if (B_new < B_min || B_new > B_max) {
                    const auto& vh = s.v_hes(v);
                    for (size_t i = vh.size(); i > 0; --i)
                        if (s.group(vh[i - 1]) != saved[i - 1])
                            s.move(vh[i - 1], saved[i - 1]);
                } else {
                    dS_iter += dS;
                }
            }
        }

        // 2. Greedy half-edge moves (enables overlapping bridge vertices!)
        for (int sw = 0; sw < 3; ++sw) {
            rng.shuffle(hes);
            double dS_sw = 0.0;
            for (int h : hes) {
                int r = s.group(h);
                size_t B_cur = static_cast<size_t>(s.B());
                gt_hash_set<int> cands;
                cands.insert(s.group(h ^ 1));
                for (int oh : s.v_hes(s.he_vertex(h)))
                    cands.insert(s.group(oh));
                for (int k = 0; k < 4; ++k)
                    cands.insert(s.sample_group(h, 0.5, rng));
                if (B_cur < B_max && s.empty_group() >= 0)
                    cands.insert(s.empty_group());
                cands.erase(r);

                int best = r;
                double best_dS = -1e-9;
                for (int t : cands) {
                    double dS = s.virtual_move(h, t);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best = t;
                    }
                }
                if (best != r) {
                    double dS = s.move(h, best);
                    size_t B_new = static_cast<size_t>(s.B());
                    if (B_new < B_min || B_new > B_max) {
                        s.move(h, r);
                    } else {
                        dS_sw += dS;
                    }
                }
            }
            dS_iter += dS_sw;
            if (std::abs(dS_sw) < 1e-8)
                break;
        }

        // 2b. Bundled mate-alignment moves (jump over D=1 -> D>=2 barrier):
        //     For each vertex whose half-edges have mates in multiple groups, try
        //     assigning each half-edge to its mate's group atomically.
        if (s.B() >= 2) {
            rng.shuffle(vs);
            for (int v : vs) {
                const auto& vh = s.v_hes(v);
                if (vh.size() < 2)
                    continue;
                std::vector<int> old_g(vh.size()), tgt_g(vh.size());
                bool any_diff = false;
                for (size_t i = 0; i < vh.size(); ++i) {
                    old_g[i] = s.group(vh[i]);
                    tgt_g[i] = s.group(s.he_mate(vh[i]));
                    if (tgt_g[i] != old_g[i])
                        any_diff = true;
                }
                if (!any_diff)
                    continue;
                double dS = 0.0;
                for (size_t i = 0; i < vh.size(); ++i)
                    if (tgt_g[i] != old_g[i])
                        dS += s.move(vh[i], tgt_g[i]);
                size_t B_new = static_cast<size_t>(s.B());
                if (dS < -1e-9 && B_new >= B_min && B_new <= B_max) {
                    dS_iter += dS;
                } else {
                    for (size_t i = 0; i < vh.size(); ++i)
                        if (tgt_g[i] != old_g[i])
                            s.move(vh[i], old_g[i]);
                }
            }
            // Also try aligning ALL boundary vertices simultaneously (amortizes D=1 -> D=2 cost)
            std::vector<int> saved_b = s.blocks();
            double dS_all = 0.0;
            bool moved_any = false;
            for (int v : vs) {
                const auto& vh = s.v_hes(v);
                if (vh.size() < 2) continue;
                for (int h : vh) {
                    int mt = saved_b[s.he_mate(h)];
                    if (s.group(h) != mt) {
                        dS_all += s.move(h, mt);
                        moved_any = true;
                    }
                }
            }
            size_t B_all = static_cast<size_t>(s.B());
            if (moved_any && dS_all < -1e-9 && B_all >= B_min && B_all <= B_max) {
                dS_iter += dS_all;
            } else if (moved_any) {
                s.set_partition(saved_b);
            }
        }

        // 2c. Atomic edge-coupled moves: move BOTH half-edges (2e, 2e+1) of an edge
        //     together to (t, t), bypassing the disassortative (t, r) energy barrier!
        if (s.B() >= 2) {
            std::vector<int> es(s.E());
            std::iota(es.begin(), es.end(), 0);
            for (int esw = 0; esw < 3; ++esw) {
                rng.shuffle(es);
                double dS_esw = 0.0;
                for (int e : es) {
                    int h0 = 2 * e, h1 = h0 + 1;
                    int r0 = s.group(h0), r1 = s.group(h1);
                    gt_hash_set<int> cands;
                    for (int oh : s.v_hes(s.he_vertex(h0))) cands.insert(s.group(oh));
                    for (int oh : s.v_hes(s.he_vertex(h1))) cands.insert(s.group(oh));
                    int best_t = -1;
                    double best_d = -1e-9;
                    for (int t : cands) {
                        if (t == r0 && t == r1) continue;
                        double d = 0.0;
                        if (r0 != t) d += s.move(h0, t);
                        if (r1 != t) d += s.move(h1, t);
                        size_t B_new = static_cast<size_t>(s.B());
                        if (d < best_d && B_new >= B_min && B_new <= B_max) {
                            best_d = d;
                            best_t = t;
                        }
                        if (r1 != t) s.move(h1, r1);
                        if (r0 != t) s.move(h0, r0);
                    }
                    if (best_t >= 0) {
                        if (r0 != best_t) dS_esw += s.move(h0, best_t);
                        if (r1 != best_t) dS_esw += s.move(h1, best_t);
                    }
                }
                dS_iter += dS_esw;
                if (std::abs(dS_esw) < 1e-8)
                    break;
            }
        }

        // 3. Greedy group merges
        while (static_cast<size_t>(s.B()) > B_min) {
            std::vector<int> occ = s.occupied();
            int best_r = -1, best_t = -1;
            double best_dS = -1e-9;
            for (size_t i = 0; i < occ.size(); ++i) {
                int r = occ[i];
                for (size_t j = i + 1; j < occ.size(); ++j) {
                    int t = occ[j];
                    double dS = s.virtual_merge(r, t);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best_r = r;
                        best_t = t;
                    }
                }
            }
            if (best_r < 0)
                break;
            dS_iter += s.merge(best_r, best_t);
        }

        // 4. Try splitting each group into 2 via random + greedy half-edge sweeps
        if (static_cast<size_t>(s.B()) < B_max) {
            std::vector<int> occ = s.occupied();
            for (int r : occ) {
                if (static_cast<size_t>(s.B()) >= B_max)
                    break;
                if (s.he_members(r).size() < 4 || s.empty_group() < 0)
                    continue;
                int t = s.empty_group();
                std::vector<int> r_hes = s.he_members(r);
                double dS_split = 0.0;
                for (int h : r_hes)
                    if (rng.unif() < 0.5)
                        dS_split += s.move(h, t);
                for (int sw = 0; sw < 5; ++sw) {
                    rng.shuffle(r_hes);
                    bool ch = false;
                    for (int h : r_hes) {
                        int cur = s.group(h);
                        int oth = (cur == r) ? t : r;
                        double d = s.virtual_move(h, oth);
                        if (d < -1e-9) {
                            dS_split += s.move(h, oth);
                            ch = true;
                        }
                    }
                    if (!ch) break;
                }
                if (dS_split < -1e-8 && s.wr(r) > 0 && s.wr(t) > 0) {
                    dS_iter += dS_split;
                } else {
                    for (int h : r_hes)
                        if (s.group(h) != r)
                            s.move(h, r);
                }
            }
        }

        // 5. Try dissolving each group into its neighboring (mate) groups:
        //    this collective move jumps directly over the D=1 -> D>=2 overlap barrier
        //    when an intersection/bridge group is better represented as overlapping memberships.
        if (static_cast<size_t>(s.B()) > B_min) {
            std::vector<int> occ = s.occupied();
            for (int r : occ) {
                if (static_cast<size_t>(s.B()) <= B_min || s.wr(r) == 0)
                    continue;
                std::vector<int> r_hes = s.he_members(r);
                // Check whether group r has edges connecting to >= 2 other groups
                gt_hash_map<int, int> nb_cnt;
                for (int h : r_hes) {
                    int mt = s.group(s.he_mate(h));
                    if (mt != r)
                        nb_cnt[mt]++;
                }
                if (nb_cnt.size() < 2)
                    continue;
                // Find most common external neighbor group as fallback for internal edges of r
                int top_nb = -1, top_c = -1;
                for (const auto& [g, c] : nb_cnt) {
                    if (c > top_c) {
                        top_c = c;
                        top_nb = g;
                    }
                }
                std::vector<int> saved_b = s.blocks();
                double dS_dis = 0.0;
                // First move half-edges with external mates to their mate's group
                std::vector<int> internal_hes;
                for (int h : r_hes) {
                    int mt = s.group(s.he_mate(h));
                    if (mt != r) {
                        dS_dis += s.move(h, mt);
                    } else {
                        internal_hes.push_back(h);
                    }
                }
                // Move remaining internal half-edges of r greedily to best neighboring group
                std::vector<int> nb_groups;
                nb_groups.reserve(nb_cnt.size());
                for (const auto& [g, c] : nb_cnt)
                    nb_groups.push_back(g);
                for (int h : internal_hes) {
                    int best_g = top_nb;
                    double best_d = inf;
                    for (int g : nb_groups) {
                        double d = s.virtual_move(h, g);
                        if (d < best_d) {
                            best_d = d;
                            best_g = g;
                        }
                    }
                    dS_dis += s.move(h, best_g);
                }
                // Local greedy cleanup on moved half-edges
                for (int sw = 0; sw < 3; ++sw) {
                    bool ch = false;
                    for (int h : r_hes) {
                        int cur = s.group(h);
                        int best_g = cur;
                        double best_d = -1e-9;
                        for (int g : nb_groups) {
                            if (g == cur) continue;
                            double d = s.virtual_move(h, g);
                            if (d < best_d) {
                                best_d = d;
                                best_g = g;
                            }
                        }
                        if (best_g != cur) {
                            dS_dis += s.move(h, best_g);
                            ch = true;
                        }
                    }
                    if (!ch) break;
                }
                if (dS_dis < -1e-8 && static_cast<size_t>(s.B()) >= B_min) {
                    dS_iter += dS_dis;
                } else {
                    s.set_partition(saved_b);
                }
            }
        }

        Rcpp::checkUserInterrupt();
        S = s.entropy();
        if (S == 0 || std::abs(dS_iter / S) <= epsilon)
            break;
    }
    return S;
}

// Agglomerative coarsening from a starting partition down to B_target.
void overlap_coarsen_to(OverlapSBM& s, ORng& rng, size_t B_target, double r_shrink = 0.85) {
    while (static_cast<size_t>(s.B()) > B_target) {
        size_t B_cur = static_cast<size_t>(s.B());
        size_t B_next = std::max(B_target, std::min(B_cur - 1, static_cast<size_t>(std::round(B_cur * r_shrink))));
        std::vector<int> occ = s.occupied();
        rng.shuffle(occ);
        std::vector<std::tuple<double, int, int>> cands;
        cands.reserve(occ.size());
        for (int r : occ) {
            int best = -1;
            double best_dS = inf;
            gt_hash_set<int> tried;
            const auto& r_hes = s.he_members(r);
            for (int k = 0; k < 10; ++k) {
                int h = r_hes[rng.index(r_hes.size())];
                int t = r;
                if (k < 5) {
                    const auto& vh = s.v_hes(s.he_vertex(h));
                    if (vh.size() > 1)
                        t = s.group(vh[rng.index(vh.size())]);
                }
                if (t == r)
                    t = s.sample_group(h, k < 5 ? 0.0 : 0.5, rng);
                if (t == r || tried.count(t))
                    continue;
                tried.insert(t);
                double dS = s.virtual_merge(r, t);
                if (dS < best_dS) {
                    best_dS = dS;
                    best = t;
                }
            }
            if (best < 0 && occ.size() > 1) {
                int t = occ[rng.index(occ.size() - 1)];
                if (t == r) t = occ.back();
                best = t;
                best_dS = s.virtual_merge(r, t);
            }
            if (best >= 0)
                cands.emplace_back(best_dS, r, best);
        }
        std::sort(cands.begin(), cands.end());
        gt_hash_map<int, int> root;
        auto get_root = [&](int x) {
            int cur = x;
            while (root.count(cur) && root[cur] != cur)
                cur = root[cur];
            root[x] = cur;
            return cur;
        };
        size_t merged = 0;
        for (const auto& [dS_est, a, b_raw] : cands) {
            if (static_cast<size_t>(s.B()) <= B_next)
                break;
            int t = get_root(b_raw);
            if (a == t || s.wr(a) == 0 || s.wr(t) == 0)
                continue;
            s.merge(a, t);
            root[a] = t;
            ++merged;
        }
        if (merged == 0)
            break;
        Rcpp::checkUserInterrupt();
    }
}

double overlap_minimize_flat(OverlapSBM& s, ORng& rng, const std::vector<int>& init_vblocks,
                             size_t B_min, size_t B_max, bool refine, double epsilon) {
    if (s.E() == 0)
        return 0;
    size_t B_hi = std::min<size_t>(B_max, static_cast<size_t>(s.H()));
    size_t B_lo = std::min(std::max<size_t>(1, B_min), B_hi);

    // Candidate 1: start from supplied vertex-level partition (e.g. from non-overlap SBM)
    std::vector<int> best_b(s.H(), 0);
    double best_S = inf;

    auto try_start = [&](const std::vector<int>& start_b) {
        s.set_partition(start_b);
        if (static_cast<size_t>(s.B()) > B_hi)
            overlap_coarsen_to(s, rng, B_hi);
        double cur_S = s.entropy();
        if (refine)
            cur_S = overlap_refine_flat(s, rng, B_lo, B_hi, epsilon, 30);
        if (cur_S < best_S) {
            best_S = cur_S;
            best_b = s.blocks();
        }
    };

    if (static_cast<int>(init_vblocks.size()) == s.N()) {
        std::vector<int> b0(s.H());
        for (int h = 0; h < s.H(); ++h)
            b0[h] = init_vblocks[s.he_vertex(h)];
        try_start(b0);
    }

    // Candidate 2: start from B = min(N_act, B_hi) singleton vertices and agglomerate
    {
        std::vector<int> b1(s.H());
        for (int h = 0; h < s.H(); ++h)
            b1[h] = s.he_vertex(h);
        s.set_partition(b1);
        if (static_cast<size_t>(s.B()) > B_hi)
            overlap_coarsen_to(s, rng, B_hi);
        // Coarsen along a geometric schedule checking best DL
        double loc_best_S = s.entropy();
        std::vector<int> loc_best_b = s.blocks();
        while (static_cast<size_t>(s.B()) > B_lo) {
            size_t B_cur = static_cast<size_t>(s.B());
            size_t B_next = std::max(B_lo, std::min(B_cur - 1, static_cast<size_t>(std::round(B_cur * 0.8))));
            overlap_coarsen_to(s, rng, B_next, 0.8);
            // Quick greedy half-edge sweep
            std::vector<int> hes(s.H());
            std::iota(hes.begin(), hes.end(), 0);
            rng.shuffle(hes);
            for (int h : hes) {
                int r = s.group(h);
                int t = s.sample_group(h, 0.5, rng);
                if (t != r && s.virtual_move(h, t) < -1e-9) {
                    s.move(h, t);
                    if (static_cast<size_t>(s.B()) < B_lo)
                        s.move(h, r);
                }
            }
            double cur_S = s.entropy();
            if (cur_S < loc_best_S) {
                loc_best_S = cur_S;
                loc_best_b = s.blocks();
            }
        }
        try_start(loc_best_b);
    }

    // Candidate 3: single group (if B_lo == 1)
    if (B_lo == 1) {
        try_start(std::vector<int>(s.H(), 0));
    }

    // Candidate 4: sequential bridge-group dissolution trajectory from best_b.
    // Dissolving a single bridge group can temporarily raise entropy because the
    // D = 1 -> D = 2 overlap prior is paid once for the whole graph; dissolving
    // multiple bridge groups sequentially amortizes that cost and reveals true
    // overlapping solutions.
    if (refine) {
        s.set_partition(best_b);
        double traj_best_S = best_S;
        std::vector<int> traj_best_b = best_b;
        while (static_cast<size_t>(s.B()) > std::max<size_t>(B_lo, 2)) {
            std::vector<int> occ = s.occupied();
            int best_r = -1;
            double best_cand_S = inf;
            std::vector<int> best_after_b;
            std::vector<int> base_b = s.blocks();
            for (int r : occ) {
                std::vector<int> r_hes = s.he_members(r);
                gt_hash_map<int, int> nb_cnt;
                for (int h : r_hes) {
                    int mt = s.group(s.he_mate(h));
                    if (mt != r)
                        nb_cnt[mt]++;
                }
                if (nb_cnt.size() < 2)
                    continue;
                int top_nb = -1, top_c = -1;
                std::vector<int> nb_groups;
                nb_groups.reserve(nb_cnt.size());
                for (const auto& [g, c] : nb_cnt) {
                    nb_groups.push_back(g);
                    if (c > top_c) { top_c = c; top_nb = g; }
                }
                std::vector<int> internal_hes;
                for (int h : r_hes) {
                    int mt = s.group(s.he_mate(h));
                    if (mt != r) s.move(h, mt);
                    else internal_hes.push_back(h);
                }
                for (int h : internal_hes) {
                    int best_g = top_nb;
                    double best_d = inf;
                    for (int g : nb_groups) {
                        double d = s.virtual_move(h, g);
                        if (d < best_d) { best_d = d; best_g = g; }
                    }
                    s.move(h, best_g);
                }
                double cand_S = s.entropy();
                if (cand_S < best_cand_S) {
                    best_cand_S = cand_S;
                    best_r = r;
                    best_after_b = s.blocks();
                }
                s.set_partition(base_b);
            }
            if (best_r < 0)
                break;
            s.set_partition(best_after_b);
            // Quick greedy edge-coupled + half-edge cleanup so bridge-to-bridge edges snap onto core groups
            std::vector<int> es_c(s.E());
            std::iota(es_c.begin(), es_c.end(), 0);
            for (int sw = 0; sw < 3; ++sw) {
                rng.shuffle(es_c);
                bool ch = false;
                for (int e : es_c) {
                    int h0 = 2 * e, h1 = h0 + 1;
                    int r0 = s.group(h0), r1 = s.group(h1);
                    gt_hash_set<int> cands;
                    for (int oh : s.v_hes(s.he_vertex(h0))) cands.insert(s.group(oh));
                    for (int oh : s.v_hes(s.he_vertex(h1))) cands.insert(s.group(oh));
                    int best_t = -1;
                    double best_d = -1e-9;
                    for (int t : cands) {
                        if (t == r0 && t == r1) continue;
                        double d = 0.0;
                        if (r0 != t) d += s.move(h0, t);
                        if (r1 != t) d += s.move(h1, t);
                        if (d < best_d && static_cast<size_t>(s.B()) >= B_lo) {
                            best_d = d;
                            best_t = t;
                        }
                        if (r1 != t) s.move(h1, r1);
                        if (r0 != t) s.move(h0, r0);
                    }
                    if (best_t >= 0) {
                        if (r0 != best_t) s.move(h0, best_t);
                        if (r1 != best_t) s.move(h1, best_t);
                        ch = true;
                    }
                }
                if (!ch) break;
            }
            double after_S = s.entropy();
            if (after_S < traj_best_S) {
                traj_best_S = after_S;
                traj_best_b = s.blocks();
            }
        }
        if (traj_best_S < best_S - 1e-8) {
            try_start(traj_best_b);
        }
    }

    // Candidate 5: edge-level agglomeration (each edge starts in its own group so
    // vertices begin with overlapping memberships D >= 2 and coarsen into edge communities).
    if (refine && s.E() > 1 && s.E() <= 600 && B_hi >= 2) {
        std::vector<int> be(s.H());
        for (int e = 0; e < s.E(); ++e) {
            be[2 * e] = e;
            be[2 * e + 1] = e;
        }
        s.set_partition(be);
        if (static_cast<size_t>(s.B()) > B_hi)
            overlap_coarsen_to(s, rng, B_hi, 0.75);
        double loc_best_S = s.entropy();
        std::vector<int> loc_best_b = s.blocks();
        while (static_cast<size_t>(s.B()) > B_lo) {
            size_t B_cur = static_cast<size_t>(s.B());
            size_t B_next = std::max(B_lo, std::min(B_cur - 1, static_cast<size_t>(std::round(B_cur * 0.75))));
            overlap_coarsen_to(s, rng, B_next, 0.75);
            if (B_next <= 32) {
                // Quick greedy half-edge sweep near the relevant B scale
                std::vector<int> hes(s.H());
                std::iota(hes.begin(), hes.end(), 0);
                rng.shuffle(hes);
                for (int h : hes) {
                    int r = s.group(h);
                    int t = s.sample_group(h, 0.5, rng);
                    if (t != r && s.virtual_move(h, t) < -1e-9) {
                        s.move(h, t);
                        if (static_cast<size_t>(s.B()) < B_lo)
                            s.move(h, r);
                    }
                }
            }
            double cur_S = s.entropy();
            if (cur_S < loc_best_S) {
                loc_best_S = cur_S;
                loc_best_b = s.blocks();
            }
        }
        try_start(loc_best_b);
    }

    s.set_partition(best_b);
    return s.entropy();
}

size_t o_q_cache_max(double x) {
    return std::isinf(x) ? 10000 : static_cast<size_t>(x);
}

}  // namespace

// [[Rcpp::export]]
double cpp_overlap_entropy(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                           bool degree_corrected, bool distributed,
                           bool gen_overlap, double q_max) {
    Graph& graph = get_graph(pointer);
    OverlapSBM state(graph, degree_corrected, distributed, gen_overlap, o_q_cache_max(q_max));
    state.set_partition(parse_half_edges(graph, half_edges));
    return state.entropy();
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_overlap_check_moves(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                                            bool degree_corrected, bool distributed,
                                            bool gen_overlap, int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    OverlapSBM state(graph, degree_corrected, distributed, gen_overlap, o_q_cache_max(q_max));
    state.set_partition(parse_half_edges(graph, half_edges));
    ORng rng;
    Rcpp::NumericVector err(nmoves);
    if (state.H() == 0)
        return err;
    double S = state.entropy();
    for (int i = 0; i < nmoves; ++i) {
        int h = static_cast<int>(rng.index(state.H()));
        int s = (rng.unif() < 0.2 && state.empty_group() >= 0)
                    ? state.empty_group()
                    : state.random_occupied(rng);
        double dS_virt = state.virtual_move(h, s);
        double dS_act = state.move(h, s);
        double S_new = state.entropy();
        err[i] = std::max(std::abs(dS_virt - (S_new - S)), std::abs(dS_act - (S_new - S)));
        S = S_new;
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_overlap_groups(SEXP pointer, Rcpp::IntegerMatrix half_edges) {
    Graph& graph = get_graph(pointer);
    OverlapSBM state(graph, false, false, false, 100);
    state.set_partition(parse_half_edges(graph, half_edges));
    return state.overlap_groups_r();
}

// [[Rcpp::export]]
Rcpp::List cpp_overlap_minimize(SEXP pointer, Rcpp::IntegerVector init_vblocks,
                                bool degree_corrected, bool distributed,
                                bool gen_overlap, int B_min, double B_max,
                                bool refine, double epsilon, double q_max) {
    Graph& graph = get_graph(pointer);
    OverlapSBM state(graph, degree_corrected, distributed, gen_overlap, o_q_cache_max(q_max));
    std::vector<int> iv(init_vblocks.size());
    for (R_xlen_t i = 0; i < init_vblocks.size(); ++i)
        iv[i] = init_vblocks[i] - 1;
    ORng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    double S = overlap_minimize_flat(state, rng, iv, static_cast<size_t>(B_min),
                                     bmax, refine, epsilon);
    return Rcpp::List::create(
        Rcpp::Named("half_edges") = state.half_edges_r(),
        Rcpp::Named("entropy") = S,
        Rcpp::Named("overlap") = state.overlap_groups_r());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_overlap_entropy(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                                               Rcpp::List upper_bs, bool degree_corrected,
                                               bool distributed, bool gen_overlap,
                                               double q_max) {
    Graph& graph = get_graph(pointer);
    int L = static_cast<int>(upper_bs.size()) + 1;
    NestedOverlap no(graph, degree_corrected, distributed, gen_overlap,
                     o_q_cache_max(q_max), L);
    std::vector<int> b0 = parse_half_edges(graph, half_edges);
    int cap = no.base().cap();
    std::vector<std::vector<int>> ubs;
    ubs.reserve(upper_bs.size());
    for (R_xlen_t l = 0; l < upper_bs.size(); ++l) {
        Rcpp::IntegerVector x = upper_bs[l];
        std::vector<int> bl(cap, 0);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            bl[i] = x[i] - 1;
        ubs.push_back(std::move(bl));
    }
    no.build(b0, ubs);
    std::vector<double> S = no.level_entropies();
    return Rcpp::NumericVector(S.begin(), S.end());
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_overlap_minimize(SEXP pointer, Rcpp::IntegerVector init_vblocks,
                                       bool degree_corrected, bool distributed,
                                       bool gen_overlap, int B_min, double B_max,
                                       bool refine, double epsilon, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    int L = N > 1 ? static_cast<int>(std::ceil(std::log2(double(N)))) + 1 : 1;
    NestedOverlap no(graph, degree_corrected, distributed, gen_overlap,
                     o_q_cache_max(q_max), L);
    std::vector<int> iv(init_vblocks.size());
    for (R_xlen_t i = 0; i < init_vblocks.size(); ++i)
        iv[i] = init_vblocks[i] - 1;
    ORng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);

    // Seed level 0 with a fast flat OverlapSBM fit (with edges_dl=true so B_0 starts
    // near the optimal scale rather than over-splitting without an edge-count prior),
    // then optimize upper levels and jointly refine level 0 + upper levels under full nested DL.
    if (no.base().E() > 0) {
        OverlapSBM flat_seed(graph, degree_corrected, distributed, gen_overlap,
                             o_q_cache_max(q_max), true);
        overlap_minimize_flat(flat_seed, rng, iv, static_cast<size_t>(B_min),
                              bmax, refine, epsilon);
        no.base().set_partition(flat_seed.blocks());
        no.optimize_upper(rng);

        // Joint refinement: greedy level-0 half-edge/vertex/merge moves accepting when
        // total nested entropy (level 0 + optimized upper levels) decreases!
        if (refine) {
            double best_S = no.entropy();
            std::vector<int> best_b0 = no.base().blocks();
            auto best_up = no.save_upper();
            std::vector<int> hes(no.base().H());
            std::iota(hes.begin(), hes.end(), 0);
            for (int iter = 0; iter < 5; ++iter) {
                bool improved = false;
                // Try merges evaluated under full nested entropy
                while (static_cast<size_t>(no.base().B()) > static_cast<size_t>(B_min)) {
                    std::vector<int> occ = no.base().occupied();
                    int mr = -1, mt = -1;
                    double mS = best_S - 1e-8;
                    for (size_t i = 0; i < occ.size(); ++i) {
                        for (size_t j = i + 1; j < occ.size(); ++j) {
                            int r = occ[i], t = occ[j];
                            std::vector<int> r_hes = no.base().he_members(r);
                            no.base().merge(r, t);
                            no.build_upper();
                            double cand_S = no.entropy();
                            if (cand_S < mS) {
                                mS = cand_S;
                                mr = r;
                                mt = t;
                            }
                            for (int h : r_hes)
                                no.base().move(h, r);
                        }
                    }
                    if (mr < 0)
                        break;
                    no.base().merge(mr, mt);
                    no.optimize_upper(rng);
                    double now_S = no.entropy();
                    if (now_S < best_S - 1e-8) {
                        best_S = now_S;
                        best_b0 = no.base().blocks();
                        best_up = no.save_upper();
                        improved = true;
                    } else {
                        no.build(best_b0, best_up);
                        break;
                    }
                }
                // Half-edge sweep checking total nested entropy
                rng.shuffle(hes);
                for (int h : hes) {
                    int r = no.base().group(h);
                    int t = no.base().sample_group(h, 0.5, rng);
                    if (t == r) continue;
                    no.base().move(h, t);
                    size_t B_now = static_cast<size_t>(no.base().B());
                    if (B_now < static_cast<size_t>(B_min) || B_now > bmax) {
                        no.base().move(h, r);
                        continue;
                    }
                    no.build_upper();
                    double cand_S = no.entropy();
                    if (cand_S < best_S - 1e-9) {
                        best_S = cand_S;
                        best_b0 = no.base().blocks();
                        best_up = no.save_upper();
                        improved = true;
                    } else {
                        no.base().move(h, r);
                    }
                }
                if (improved) {
                    no.build(best_b0, best_up);
                    no.optimize_upper(rng);
                    double now_S = no.entropy();
                    if (now_S < best_S) {
                        best_S = now_S;
                        best_b0 = no.base().blocks();
                        best_up = no.save_upper();
                    } else {
                        no.build(best_b0, best_up);
                    }
                } else {
                    break;
                }
                Rcpp::checkUserInterrupt();
            }
            no.build(best_b0, best_up);
        }
    }

    std::vector<double> S = no.level_entropies();
    return Rcpp::List::create(
        Rcpp::Named("half_edges") = no.base().half_edges_r(),
        Rcpp::Named("upper_levels") = no.external_upper_r(),
        Rcpp::Named("level_entropy") = S,
        Rcpp::Named("entropy") = std::accumulate(S.begin(), S.end(), 0.0),
        Rcpp::Named("overlap") = no.base().overlap_groups_r());
}

// Finite-temperature half-edge Metropolis-Hastings MCMC for flat OverlapSBM.
// [[Rcpp::export]]
Rcpp::List cpp_overlap_mcmc(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                            bool degree_corrected, bool distributed,
                            bool gen_overlap, double beta, double c, double d,
                            int niter, int thin, bool collect, int B_min,
                            double B_max, double q_max) {
    Graph& graph = get_graph(pointer);
    OverlapSBM state(graph, degree_corrected, distributed, gen_overlap, o_q_cache_max(q_max));
    state.set_partition(parse_half_edges(graph, half_edges));
    ORng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);

    std::vector<int> hes(state.H());
    std::iota(hes.begin(), hes.end(), 0);
    double total_dS = 0;
    size_t nattempts = 0, nmoves = 0;
    std::vector<Rcpp::IntegerVector> maj_samples;
    std::vector<double> entropies;
    double S = state.entropy();

    for (int it = 1; it <= niter; ++it) {
        rng.shuffle(hes);
        for (int h : hes) {
            ++nattempts;
            int r = state.group(h);
            size_t B_cur = static_cast<size_t>(state.B());
            bool can_create = (B_cur < static_cast<size_t>(state.H()) &&
                               B_cur < bmax && d > 0 && state.empty_group() >= 0);
            double p_new = can_create ? d : 0.0;
            bool propose_new = (can_create && rng.unif() < p_new);
            int s;
            double pfwd = 0.0;
            if (propose_new) {
                if (state.he_members(r).size() == 1)
                    continue;
                s = state.empty_group();
                pfwd = p_new;
            } else {
                s = state.sample_group(h, c, rng);
                if (s == r)
                    continue;
                if (state.he_members(r).size() == 1 && B_cur <= static_cast<size_t>(B_min))
                    continue;
                pfwd = (1.0 - p_new) * state.move_prob(h, s, c);
            }
            bool empties = (state.he_members(r).size() == 1);
            double dS = state.move(h, s);
            bool accept = false;
            if (std::isinf(beta)) {
                accept = (dS < 0);
            } else {
                double prev = 0.0;
                if (empties) {
                    prev = (d > 0) ? d : 0.0;
                } else {
                    size_t B_after = static_cast<size_t>(state.B());
                    bool rev_can = (B_after < static_cast<size_t>(state.H()) &&
                                    B_after < bmax && d > 0);
                    double p_new_rev = rev_can ? d : 0.0;
                    prev = (1.0 - p_new_rev) * state.move_prob(h, r, c);
                }
                if (prev > 0 && pfwd > 0) {
                    double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                    accept = (log_a >= 0 || rng.unif() < std::exp(log_a));
                }
            }
            if (accept) {
                total_dS += dS;
                S += dS;
                ++nmoves;
            } else {
                state.move(h, r);
            }
        }
        if (collect && (it % std::max(1, thin) == 0)) {
            Rcpp::List ov = state.overlap_groups_r();
            maj_samples.push_back(Rcpp::as<Rcpp::IntegerVector>(ov["majority"]));
            entropies.push_back(S);
        }
        Rcpp::checkUserInterrupt();
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("half_edges") = state.half_edges_r(),
        Rcpp::Named("entropy") = state.entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("overlap") = state.overlap_groups_r());
    if (collect) {
        Rcpp::List sl(maj_samples.size());
        for (size_t i = 0; i < maj_samples.size(); ++i)
            sl[i] = maj_samples[i];
        out["majority_samples"] = sl;
        out["entropies"] = entropies;
    }
    return out;
}

// Finite-temperature Metropolis-Hastings MCMC for NestedOverlap (level-0 half-edge moves
// + upper-level single-node moves under exact total nested description length).
// [[Rcpp::export]]
Rcpp::List cpp_nested_overlap_mcmc(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                                   Rcpp::List upper_bs, bool degree_corrected,
                                   bool distributed, bool gen_overlap,
                                   double beta, double c, double d,
                                   int niter, int thin, bool collect,
                                   int B_min, double B_max, double q_max) {
    Graph& graph = get_graph(pointer);
    int L = static_cast<int>(upper_bs.size()) + 1;
    NestedOverlap no(graph, degree_corrected, distributed, gen_overlap,
                     o_q_cache_max(q_max), L);
    std::vector<int> b0 = parse_half_edges(graph, half_edges);
    int cap = no.base().cap();
    std::vector<std::vector<int>> ubs;
    ubs.reserve(upper_bs.size());
    for (R_xlen_t l = 0; l < upper_bs.size(); ++l) {
        Rcpp::IntegerVector x = upper_bs[l];
        std::vector<int> bl(cap, 0);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            bl[i] = x[i] - 1;
        ubs.push_back(std::move(bl));
    }
    no.build(b0, ubs);

    ORng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    std::vector<int> hes(no.base().H());
    std::iota(hes.begin(), hes.end(), 0);
    double total_dS = 0;
    size_t nattempts = 0, nmoves = 0;
    std::vector<Rcpp::IntegerVector> maj_samples;
    std::vector<double> entropies;
    double S = no.entropy();

    for (int it = 1; it <= niter; ++it) {
        // Level-0 half-edge MH sweep under exact total nested entropy
        rng.shuffle(hes);
        for (int h : hes) {
            ++nattempts;
            int r = no.base().group(h);
            size_t B_cur = static_cast<size_t>(no.base().B());
            bool can_create = (B_cur < static_cast<size_t>(no.base().H()) &&
                               B_cur < bmax && d > 0 && no.base().empty_group() >= 0);
            double p_new = can_create ? d : 0.0;
            bool propose_new = (can_create && rng.unif() < p_new);
            int s;
            double pfwd = 0.0;
            if (propose_new) {
                if (no.base().he_members(r).size() == 1)
                    continue;
                s = no.base().empty_group();
                pfwd = p_new;
            } else {
                s = no.base().sample_group(h, c, rng);
                if (s == r)
                    continue;
                if (no.base().he_members(r).size() == 1 && B_cur <= static_cast<size_t>(B_min))
                    continue;
                pfwd = (1.0 - p_new) * no.base().move_prob(h, s, c);
            }
            bool empties = (no.base().he_members(r).size() == 1);
            double dS = no.move_he(h, s);
            bool accept = false;
            if (std::isinf(beta)) {
                accept = (dS < 0);
            } else {
                double prev = 0.0;
                if (empties) {
                    prev = (d > 0) ? d : 0.0;
                } else {
                    size_t B_after = static_cast<size_t>(no.base().B());
                    bool rev_can = (B_after < static_cast<size_t>(no.base().H()) &&
                                    B_after < bmax && d > 0);
                    double p_new_rev = rev_can ? d : 0.0;
                    prev = (1.0 - p_new_rev) * no.base().move_prob(h, r, c);
                }
                if (prev > 0 && pfwd > 0) {
                    double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                    accept = (log_a >= 0 || rng.unif() < std::exp(log_a));
                }
            }
            if (accept) {
                total_dS += dS;
                S += dS;
                ++nmoves;
            } else {
                no.move_he(h, r);
            }
        }
        // Upper-level single-node MH sweep (O(deg) incremental updates)
        for (int l = 1; l < no.levels(); ++l) {
            const auto& nodes = (l == 1) ? no.base().occupied() : no.upper(l - 1).occupied();
            if (nodes.size() <= 1) continue;
            std::vector<int> ord = nodes;
            rng.shuffle(ord);
            for (int x : ord) {
                int r = no.upper(l).group(x);
                std::vector<int> cands = no.upper(l).occupied();
                int emp = no.upper(l).empty_group();
                if (emp >= 0 && no.upper(l).size(r) > 1)
                    cands.push_back(emp);
                if (cands.size() <= 1) continue;
                int t = cands[rng.index(cands.size())];
                if (t == r) continue;
                ++nattempts;
                double pfwd = 1.0 / static_cast<double>(cands.size());
                double dS = no.move_upper(l, x, t);
                bool accept = false;
                if (std::isinf(beta)) {
                    accept = (dS < 0);
                } else {
                    size_t rev_sz = no.upper(l).occupied().size() +
                                    ((no.upper(l).empty_group() >= 0 && no.upper(l).size(t) > 1) ? 1 : 0);
                    double prev = rev_sz > 0 ? (1.0 / static_cast<double>(rev_sz)) : 0.0;
                    if (prev > 0 && pfwd > 0) {
                        double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                        accept = (log_a >= 0 || rng.unif() < std::exp(log_a));
                    }
                }
                if (accept) {
                    total_dS += dS;
                    S += dS;
                    ++nmoves;
                } else {
                    no.move_upper(l, x, r);
                }
            }
        }
        if (collect && (it % std::max(1, thin) == 0)) {
            Rcpp::List ov = no.base().overlap_groups_r();
            maj_samples.push_back(Rcpp::as<Rcpp::IntegerVector>(ov["majority"]));
            entropies.push_back(S);
        }
        Rcpp::checkUserInterrupt();
    }

    std::vector<double> lvS = no.level_entropies();
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("half_edges") = no.base().half_edges_r(),
        Rcpp::Named("upper_levels") = no.external_upper_r(),
        Rcpp::Named("level_entropy") = lvS,
        Rcpp::Named("entropy") = std::accumulate(lvS.begin(), lvS.end(), 0.0),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("overlap") = no.base().overlap_groups_r());
    if (collect) {
        Rcpp::List sl(maj_samples.size());
        for (size_t i = 0; i < maj_samples.size(); ++i)
            sl[i] = maj_samples[i];
        out["majority_samples"] = sl;
        out["entropies"] = entropies;
    }
    return out;
}

// Verify exact agreement between incremental NestedOverlap move_he / move_upper
// and full recomputation across random moves.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_overlap_check_moves(SEXP pointer, Rcpp::IntegerMatrix half_edges,
                                                   Rcpp::List upper_bs, bool degree_corrected,
                                                   bool distributed, bool gen_overlap,
                                                   int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    int L = static_cast<int>(upper_bs.size()) + 1;
    NestedOverlap no(graph, degree_corrected, distributed, gen_overlap,
                     o_q_cache_max(q_max), L);
    std::vector<int> b0 = parse_half_edges(graph, half_edges);
    int cap = no.base().cap();
    std::vector<std::vector<int>> ubs;
    ubs.reserve(upper_bs.size());
    for (R_xlen_t l = 0; l < upper_bs.size(); ++l) {
        Rcpp::IntegerVector x = upper_bs[l];
        std::vector<int> bl(cap, 0);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            bl[i] = x[i] - 1;
        ubs.push_back(std::move(bl));
    }
    no.build(b0, ubs);
    ORng rng;
    Rcpp::NumericVector err(nmoves);
    double S = no.entropy();
    for (int i = 0; i < nmoves; ++i) {
        double dS = 0.0;
        if (rng.unif() < 0.65 || no.levels() <= 1) {
            if (no.base().H() == 0) continue;
            int h = static_cast<int>(rng.index(static_cast<size_t>(no.base().H())));
            int emp = no.base().empty_group();
            int s = (rng.unif() < 0.2 && emp >= 0) ? emp : no.base().random_occupied(rng);
            dS = no.move_he(h, s);
        } else {
            int l = 1 + static_cast<int>(rng.index(static_cast<size_t>(no.levels() - 1)));
            const auto& nodes = (l == 1) ? no.base().occupied() : no.upper(l - 1).occupied();
            if (nodes.empty()) continue;
            int x = nodes[rng.index(nodes.size())];
            std::vector<int> cands = no.upper(l).occupied();
            int emp = no.upper(l).empty_group();
            if (emp >= 0) cands.push_back(emp);
            int t = cands[rng.index(cands.size())];
            dS = no.move_upper(l, x, t);
        }
        double S_new = no.entropy();
        err[i] = std::abs(dS - (S_new - S));
        S = S_new;
    }
    return err;
}
