// graph-tool-R: native C++ flat and hierarchical weighted stochastic block models
// (Peixoto 2018, Phys. Rev. E 97, 012306).
//
// Supported edge covariate models (single or multiple simultaneous covariates):
//   - "real-exponential"    (microcanonical or Gamma(alpha, beta) conjugate prior)
//   - "real-normal"         (microcanonical or Normal-Inverse-Chi-Squared prior)
//   - "discrete-geometric"  (microcanonical or Beta(alpha, beta) conjugate prior)
//   - "discrete-poisson"    (microcanonical or Gamma(alpha, beta) conjugate prior)
//   - "discrete-binomial"   (microcanonical or Beta(alpha, beta) with N trials)
//
// Copyright (C) 2026 Diogo Melo.
// Algorithms and entropy formulas follow graph-tool by Tiago de Paula Peixoto.
// LGPL-3.0-or-later.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include "graph_types.h"
#include "parallel.h"
#include "weight_cov.h"
#include "logq.h"

namespace {

using namespace rgt_cov;

template <class K, class V, class H = boost::hash<K>>
using wt_hash_map = boost::unordered_flat_map<K, V, H>;

template <class K, class H = boost::hash<K>>
using wt_hash_set = boost::unordered_flat_set<K, H>;

constexpr double W_LOG2 = 0.693147180559945309417;

// Precomputed integer lgamma and log lookup tables.
struct WMath {
    std::vector<double> lg, lg1, ln;
    explicit WMath(size_t max_n = 262144) { ensure(max_n); }
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
    inline double lgam(int64_t n) const {
        if (n <= 1) return 0.0;
        if (static_cast<size_t>(n) < lg.size()) return lg[n];
        return std::lgamma(static_cast<double>(n));
    }
    inline double logi(int64_t n) const {
        if (n <= 1) return 0.0;
        if (static_cast<size_t>(n) < ln.size()) return ln[n];
        return std::log(static_cast<double>(n));
    }
    inline double lbinom(int64_t n, int64_t k) const {
        if (k <= 0 || k >= n) return 0.0;
        return lgam1(n) - lgam1(k) - lgam1(n - k);
    }
};

WMath& wmath() {
    static WMath m;
    return m;
}

// Shared log q(n, k) table (src/logq.h).
using WLogQ = rgt_logq::LogQ;
inline WLogQ& wlog_q() { return rgt_logq::shared(); }

struct WRng {
    double unif() { return R::unif_rand(); }
    size_t index(size_t n) {
        size_t i = static_cast<size_t>(R::unif_rand() * static_cast<double>(n));
        return i < n ? i : n - 1;
    }
    template <class T>
    void shuffle(std::vector<T>& a) {
        for (size_t i = a.size(); i > 1; --i)
            std::swap(a[i - 1], a[index(i)]);
    }
};

struct WPairStat {
    long m_topo = 0;  // actual undirected/directed edge multiplicity (e_rr = 2*m_topo for undir r==r)
    long m_cov = 0;   // covariate observation count (== m_topo at level 0; # lower non-empty pairs at l >= 1)
    WVec w;

    WPairStat() = default;
    explicit WPairStat(int K) : w(K) {}
};

// Edge representation inside a WeightedLevel
struct WGraphEdge {
    int u, v;
    long m_topo;      // 1 at level 0; m_topo of lower pair at l >= 1
    long m_cov;       // 1 at all levels
    WVec w;           // w.x(k) = edge weight; w.x2(k) = w^2
};

struct WAdjRef {
    int nbr;
    int eidx;
    bool out;  // true if edge is v -> nbr (or undirected)
};

// Pack (a, b) into uint64_t key
inline uint64_t pack_pair(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(b));
}

inline uint64_t canon_pair(int a, int b, bool directed) {
    if (!directed && a > b) std::swap(a, b);
    return pack_pair(a, b);
}

inline std::pair<int, int> unpack_pair(uint64_t key) {
    return {static_cast<int>(static_cast<uint32_t>(key >> 32)),
            static_cast<int>(static_cast<uint32_t>(key & 0xffffffffULL))};
}

// One level of a WeightedSBM (base level l=0 or upper dense level l>=1).
class WeightedLevel {
public:
    WeightedLevel(int N, bool directed, bool dc, bool dist, bool is_base,
                  bool is_top, std::vector<CovSpec> covs)
        : _N(N), _directed(directed), _dc(dc), _dist(dist),
          _is_base(is_base), _is_top(is_top), _covs(std::move(covs)),
          _K(static_cast<int>(_covs.size())),
          _cap(std::max(2, N)),
          _adj(N), _kout(N, 0), _kin(directed ? N : 0, 0),
          _b(N, 0), _wr(_cap, 0), _ep(_cap, 0), _em(directed ? _cap : 0, 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(N, 0),
          _m_out_topo(_cap), _m_in_topo(directed ? _cap : 0),
          _recsum(_K, 0.0), _recx2(_K, 0.0), _recdx(_K, 0.0),
          _h_out(dc && dist ? _cap : 0),
          _h_in(dc && dist && directed ? _cap : 0) {
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        for (int k = 0; k < _K; ++k)
            if (_covs[k].type == CovType::REAL_NORMAL)
                _normal_idx.push_back(k);
    }

    void set_edges(std::vector<WGraphEdge> edges) {
        _edges = std::move(edges);
        _E = 0;
        _adj.assign(_N, {});
        _kout.assign(_N, 0);
        if (_directed) _kin.assign(_N, 0);
        _recsum.assign(_K, 0.0);
        for (int ei = 0; ei < static_cast<int>(_edges.size()); ++ei) {
            const auto& e = _edges[ei];
            _E += e.m_topo;
            for (int k = 0; k < _K; ++k) _recsum[k] += e.w.x(k);
            if (_directed) {
                _adj[e.u].push_back({e.v, ei, true});
                _adj[e.v].push_back({e.u, ei, false});
                _kout[e.u] += e.m_topo;
                _kin[e.v] += e.m_topo;
            } else {
                _adj[e.u].push_back({e.v, ei, true});
                _kout[e.u] += e.m_topo;
                if (e.u != e.v) {
                    _adj[e.v].push_back({e.u, ei, true});
                    _kout[e.v] += e.m_topo;
                } else {
                    _kout[e.u] += e.m_topo;
                }
            }
        }
        // Compute base-level constant likelihood factor (parallel edges, degrees, cov consts)
        _base_const = 0.0;
        if (_is_base && _E > 0) {
            const WMath& M = wmath();
            wt_hash_map<uint64_t, long> emult;
            emult.reserve(_edges.size());
            for (const auto& e : _edges) {
                emult[canon_pair(e.u, e.v, _directed)] += e.m_topo;
            }
            for (const auto& [key, mult] : emult) {
                auto [u, v] = unpack_pair(key);
                _base_const += M.lgam1(mult);
                if (!_directed && u == v) _base_const += static_cast<double>(mult) * W_LOG2;
            }
            if (_dc) {
                for (int v = 0; v < _N; ++v) {
                    _base_const -= M.lgam1(_kout[v]);
                    if (_directed) _base_const -= M.lgam1(_kin[v]);
                }
            }
            for (int k = 0; k < _K; ++k) _base_const += _covs[k].const_term;
        }
        wmath().ensure(static_cast<size_t>(_N + _E + 64));
        if (_dc && _dist) wlog_q().init(rgt_logq::table_size(static_cast<double>(_E), 5000));
    }

    int N() const { return _N; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    int K() const { return _K; }
    bool directed() const { return _directed; }
    bool is_base() const { return _is_base; }
    bool is_top() const { return _is_top; }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    long B_E() const { return _B_E; }
    long B_E_D() const { return _B_E_D; }
    double recdx(int k) const { return _recdx[k]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<CovSpec>& covs() const { return _covs; }
    const std::vector<int>& normal_idx() const { return _normal_idx; }
    const wt_hash_map<uint64_t, WPairStat>& pairs() const { return _pairs; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }

    template <class Rng>
    int random_occupied(Rng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    void set_partition(const std::vector<int>& b) {
        if (static_cast<int>(b.size()) != _N)
            Rcpp::stop("Internal error: partition length mismatch");
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            if (_directed) _em[r] = 0;
            _members[r].clear();
            _occ_pos[r] = -1;
            _m_out_topo[r].clear();
            if (_directed) _m_in_topo[r].clear();
            if (_dc && _dist) {
                _h_out[r].clear();
                if (_directed) _h_in[r].clear();
            }
            push_empty(r);
        }
        _occupied.clear();
        _pairs.clear();
        _B_E = 0;
        _B_E_D = 0;
        _recx2.assign(_K, 0.0);
        _recdx.assign(_K, 0.0);

        _b = b;
        for (int v = 0; v < _N; ++v) {
            int r = _b[v];
            if (r < 0 || r >= _cap) Rcpp::stop("Block label out of range");
            add_member(v, r);
            if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                if (_directed) _h_in[r][_kin[v]]++;
            }
        }

        for (const auto& e : _edges) {
            int r = _b[e.u], s = _b[e.v];
            uint64_t key = canon_pair(r, s, _directed);
            auto it = _pairs.find(key);
            if (it == _pairs.end()) {
                it = _pairs.emplace(key, WPairStat(_K)).first;
            }
            WPairStat& ps = it->second;
            ps.m_topo += e.m_topo;
            ps.m_cov += e.m_cov;
            for (int k = 0; k < _K; ++k) {
                ps.w.x(k) += e.w.x(k);
                ps.w.x2(k) += e.w.x2(k);
            }
            if (_directed) {
                _m_out_topo[r][s] += e.m_topo;
                _m_in_topo[s][r] += e.m_topo;
                _ep[r] += e.m_topo;
                _em[s] += e.m_topo;
            } else {
                if (r == s) {
                    _m_out_topo[r][r] += 2 * e.m_topo;
                    _ep[r] += 2 * e.m_topo;
                } else {
                    _m_out_topo[r][s] += e.m_topo;
                    _m_out_topo[s][r] += e.m_topo;
                    _ep[r] += e.m_topo;
                    _ep[s] += e.m_topo;
                }
            }
        }
        refresh_cov_globals();
    }

    void refresh_cov_globals() {
        _B_E = 0;
        _B_E_D = 0;
        _recx2.assign(_K, 0.0);
        _recdx.assign(_K, 0.0);
        for (const auto& [key, ps] : _pairs) {
            if (ps.m_cov > 0) {
                ++_B_E;
                if (ps.m_cov > 1) ++_B_E_D;
                for (int k : _normal_idx) {
                    double x = ps.w.x(k);
                    _recx2[k] += x * x;
                    if (ps.m_cov > 1) {
                        double dx = ps.w.x2(k) - x * (x / static_cast<double>(ps.m_cov));
                        if (dx > 0.0) _recdx[k] += dx;
                    }
                }
            }
        }
    }

    // Topology entropy (exact match to BlockLevel in blockmodel.cpp)
    double topo_entropy() const {
        if (_N == 0) return 0.0;
        const WMath& M = wmath();
        long B_act = B();
        double S = M.lbinom(_N - 1, B_act - 1) + M.lgam1(_N) + M.logi(_N);
        if (_is_top && _E > 0) {
            int64_t pairs = _directed ? int64_t(B_act) * B_act
                                      : int64_t(B_act) * (B_act + 1) / 2;
            S += M.lbinom(pairs + _E - 1, _E);
        }
        if (_is_base) {
            S += _base_const;
            for (int r : _occupied) {
                S -= M.lgam1(_wr[r]);
                if (_dc) {
                    S += M.lgam1(_ep[r]);
                    if (_directed) S += M.lgam1(_em[r]);
                    S += deg_dl_group(r);
                } else {
                    S += static_cast<double>(_ep[r] + (_directed ? _em[r] : 0)) * M.logi(_wr[r]);
                }
            }
            for (const auto& [key, ps] : _pairs) {
                if (ps.m_topo <= 0) continue;
                auto [r, s] = unpack_pair(key);
                S -= base_eterm(r, s, ps.m_topo);
            }
        } else {
            for (int r : _occupied) S -= M.lgam1(_wr[r]);
            for (const auto& [key, ps] : _pairs) {
                if (ps.m_topo <= 0) continue;
                auto [r, s] = unpack_pair(key);
                S += dense_eterm(r, s, ps.m_topo, _wr[r], _wr[s]);
            }
        }
        return S;
    }

    // Covariate entropy (local + rec_dl if is_top; excludes nested Lrecdx global term)
    double cov_entropy(bool include_rec = true, bool include_dl = true) const {
        if (_K == 0 || _E == 0) return 0.0;
        double S = 0.0;
        if (include_rec) {
            for (const auto& [key, ps] : _pairs) {
                if (ps.m_cov <= 0) continue;
                double m = static_cast<double>(ps.m_cov);
                for (int k = 0; k < _K; ++k)
                    S += cov_pair_neg_log_P(_covs[k], m, ps.w.x(k), ps.w.x2(k));
            }
        }
        for (int k = 0; k < _K; ++k) {
            const CovSpec& c = _covs[k];
            if (c.type == CovType::REAL_NORMAL) {
                if (c.micro) {
                    if (include_dl && _is_top)
                        S += -signed_w_log_P(static_cast<double>(_B_E), _recsum[k],
                                             _recx2[k], c.p0, c.p1, c.p2, c.p3, c.eps);
                    S += -positive_w_log_P(static_cast<double>(_B_E_D), _recdx[k],
                                           c.p2, c.p3, c.eps);
                }
            } else if (include_dl && _is_top && c.micro) {
                if (c.type == CovType::REAL_EXPONENTIAL) {
                    S += -positive_w_log_P(static_cast<double>(_B_E), _recsum[k],
                                           c.p0, c.p1, c.eps);
                } else {
                    // geometric, poisson, and binomial all use microcanonical geometric prior on recsum
                    S += -geometric_w_log_P(static_cast<double>(_B_E), _recsum[k],
                                            NAN, NAN);
                }
            }
        }
        return S;
    }

    double entropy(bool include_rec = true, bool include_dl = true) const {
        return topo_entropy() + cov_entropy(include_rec, include_dl);
    }

    struct PairDelta {
        uint64_t key;
        long dm_topo = 0;
        long dm_cov = 0;
        WVec dw;
        explicit PairDelta(int K = 0) : key(0), dw(K) {}
    };

    struct MoveBuf {
        std::vector<PairDelta> pdeltas;
        wt_hash_map<uint64_t, size_t> scratch_idx;
        std::vector<uint64_t> old_keys;
        long d_ep_r = 0, d_ep_s = 0, d_em_r = 0, d_em_s = 0;
        long d_B_E = 0, d_B_E_D = 0;
        std::vector<double> d_recx2, d_recdx;
    };

    MoveBuf& default_buf() const { return _scratch_buf; }

    // Exact O(deg(v)) incremental virtual_move(v, s) for flat WeightedLevel.
    double virtual_move(int v, int s) const {
        int r = _b[v];
        if (r == s) return 0.0;
        return compute_move_delta(v, r, s, &_scratch_buf);
    }

    // Thread-safe read-only virtual_move(v, s) using caller-supplied MoveBuf.
    double virtual_move_const(int v, int s, MoveBuf& buf) const {
        int r = _b[v];
        if (r == s) return 0.0;
        return compute_move_delta(v, r, s, &buf);
    }

    // Exact O(deg(v)) incremental move(v, s) returning exact dS.
    double move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double dS = compute_move_delta(v, r, s, &_scratch_buf);
        apply_move(v, r, s, _scratch_buf);
        return dS;
    }

    // Apply the move most recently evaluated by virtual_move(v, s).
    void apply_last_virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s) return;
        apply_move(v, r, s, _scratch_buf);
    }

    // Direct group-level virtual_merge(r, t): move all vertices in r into t.
    double virtual_merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        for (int v : verts) move(v, r);
        return dS;
    }

    // Thread-safe read-only virtual_merge(r, t) for parallel coarsening/merge sweeps.
    double virtual_merge_const(int r, int t, MoveBuf& buf) const {
        if (r == t || _wr[r] == 0 || _wr[t] == 0) return 0.0;
        const WMath& M = wmath();
        auto& pds = buf.pdeltas;
        auto& sidx = buf.scratch_idx;
        auto& old_keys = buf.old_keys;
        pds.clear();
        sidx.clear();
        old_keys.clear();

        auto find_or_add = [&](uint64_t key) -> PairDelta& {
            auto it = sidx.find(key);
            if (it == sidx.end()) {
                it = sidx.emplace(key, pds.size()).first;
                pds.emplace_back(_K);
                pds.back().key = key;
            }
            return pds[it->second];
        };

        // Enumerate unique non-empty pairs incident on group r via _m_out_topo / _m_in_topo.
        for (const auto& [h, _m] : _m_out_topo[r]) {
            (void)_m;
            uint64_t ok = canon_pair(r, h, _directed);
            if (sidx.emplace(ok, 0).second) old_keys.push_back(ok);
        }
        if (_directed) {
            for (const auto& [h, _m] : _m_in_topo[r]) {
                (void)_m;
                uint64_t ok = pack_pair(h, r);
                if (sidx.emplace(ok, 0).second) old_keys.push_back(ok);
            }
        }
        sidx.clear();

        auto add_pd = [&](uint64_t key, long dmt, long dmc, const WVec& w, double sign) {
            PairDelta& pd = find_or_add(key);
            pd.dm_topo += dmt;
            pd.dm_cov += dmc;
            for (int k = 0; k < _K; ++k) {
                pd.dw.x(k) += sign * w.x(k);
                pd.dw.x2(k) += sign * w.x2(k);
            }
        };
        for (uint64_t ok : old_keys) {
            auto pit = _pairs.find(ok);
            if (pit == _pairs.end()) continue;
            const WPairStat& ps = pit->second;
            auto [a, b_g] = unpack_pair(ok);
            int na = (a == r) ? t : a;
            int nb = (b_g == r) ? t : b_g;
            uint64_t nk = canon_pair(na, nb, _directed);
            add_pd(ok, -ps.m_topo, -ps.m_cov, ps.w, -1.0);
            add_pd(nk, +ps.m_topo, +ps.m_cov, ps.w, +1.0);
        }

        long wr_r0 = _wr[r], wr_r1 = 0;
        long wr_t0 = _wr[t], wr_t1 = wr_t0 + wr_r0;
        long B0 = B(), B1 = B0 - 1;

        double dS = 0.0;
        dS += M.lbinom(_N - 1, B1 - 1) - M.lbinom(_N - 1, B0 - 1);
        if (_is_top && _E > 0) {
            int64_t p0 = _directed ? int64_t(B0) * B0 : int64_t(B0) * (B0 + 1) / 2;
            int64_t p1 = _directed ? int64_t(B1) * B1 : int64_t(B1) * (B1 + 1) / 2;
            dS += M.lbinom(p1 + _E - 1, _E) - M.lbinom(p0 + _E - 1, _E);
        }
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_t1) - M.lgam1(wr_t0));

        if (_is_base) {
            long ep_r0 = _ep[r], ep_t0 = _ep[t], ep_t1 = ep_t0 + ep_r0;
            long em_r0 = _directed ? _em[r] : 0;
            long em_t0 = _directed ? _em[t] : 0;
            long em_t1 = em_t0 + em_r0;
            if (_dc) {
                dS += M.lgam1(ep_t1) - M.lgam1(ep_t0) - M.lgam1(ep_r0);
                if (_directed) dS += M.lgam1(em_t1) - M.lgam1(em_t0) - M.lgam1(em_r0);
                if (!_dist) {
                    dS += M.lbinom(wr_t1 + ep_t1 - 1, ep_t1) -
                          M.lbinom(wr_t0 + ep_t0 - 1, ep_t0) -
                          M.lbinom(wr_r0 + ep_r0 - 1, ep_r0);
                    if (_directed) {
                        dS += M.lbinom(wr_t1 + em_t1 - 1, em_t1) -
                              M.lbinom(wr_t0 + em_t0 - 1, em_t0) -
                              M.lbinom(wr_r0 + em_r0 - 1, em_r0);
                    }
                } else {
                    const WLogQ& q = wlog_q();
                    dS += (q(ep_t1, wr_t1) + M.lgam1(wr_t1)) -
                          (q(ep_t0, wr_t0) + M.lgam1(wr_t0)) -
                          (q(ep_r0, wr_r0) + M.lgam1(wr_r0));
                    for (const auto& [dk, cr] : _h_out[r]) {
                        auto ito = _h_out[t].find(dk);
                        long ct0 = (ito == _h_out[t].end()) ? 0 : ito->second;
                        dS += M.lgam1(cr) - (M.lgam1(ct0 + cr) - M.lgam1(ct0));
                    }
                    if (_directed) {
                        dS += (q(em_t1, wr_t1) + M.lgam1(wr_t1)) -
                              (q(em_t0, wr_t0) + M.lgam1(wr_t0)) -
                              (q(em_r0, wr_r0) + M.lgam1(wr_r0));
                        for (const auto& [dk, cr] : _h_in[r]) {
                            auto iti = _h_in[t].find(dk);
                            long ct0 = (iti == _h_in[t].end()) ? 0 : iti->second;
                            dS += M.lgam1(cr) - (M.lgam1(ct0 + cr) - M.lgam1(ct0));
                        }
                    }
                }
            } else {
                long dr0 = ep_r0 + em_r0, dt0 = ep_t0 + em_t0;
                dS += static_cast<double>(dr0 + dt0) * M.logi(wr_t1) -
                      static_cast<double>(dr0) * M.logi(wr_r0) -
                      static_cast<double>(dt0) * M.logi(wr_t0);
            }
        } else {
            // Dense upper level: also touch any existing pair incident on t whose wr[t] grows.
            for (const auto& [h, m2] : _m_out_topo[t]) {
                (void)m2;
                find_or_add(canon_pair(t, h, _directed));
            }
            if (_directed) {
                for (const auto& [h, m2] : _m_in_topo[t]) {
                    (void)m2;
                    find_or_add(pack_pair(h, t));
                }
            }
        }

        auto wr_after = [&](int g) -> long {
            if (g == r) return wr_r1;
            if (g == t) return wr_t1;
            return _wr[g];
        };
        dS += eval_pair_cov_deltas(buf, wr_after);
        return dS;
    }

    double merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        return dS;
    }

    // Fast neighbor-based proposal (matching blockmodel.cpp)
    template <class Rng>
    int sample_group(int v, double c, Rng& rng) const {
        if (_occupied.size() <= 1) return _b[v];
        if (std::isinf(c) || _adj[v].empty()) return random_occupied(rng);
        const auto& adj = _adj[v];
        int u = adj[rng.index(adj.size())].nbr;
        int t = _b[u];
        long kt = _ep[t] + (_directed ? _em[t] : 0);
        double B_act = static_cast<double>(_occupied.size());
        if (kt == 0 || rng.unif() < (c * B_act) / (static_cast<double>(kt) + c * B_act))
            return random_occupied(rng);
        long total = kt;
        long pick = static_cast<long>(rng.index(static_cast<size_t>(total)));
        if (pick < _ep[t]) {
            for (const auto& [g, m] : _m_out_topo[t]) {
                if (pick < m) return g;
                pick -= m;
            }
        } else if (_directed) {
            pick -= _ep[t];
            for (const auto& [g, m] : _m_in_topo[t]) {
                if (pick < m) return g;
                pick -= m;
            }
        }
        return random_occupied(rng);
    }

    double move_prob(int v, int s, double c) const {
        double B_act = static_cast<double>(_occupied.size());
        if (B_act <= 0.0) return 0.0;
        if (std::isinf(c) || _adj[v].empty()) return 1.0 / B_act;
        const auto& adj = _adj[v];
        double prob = 0.0;
        for (const auto& a : adj) {
            int t = _b[a.nbr];
            long kt = _ep[t] + (_directed ? _em[t] : 0);
            double m_ts = 0.0;
            auto ito = _m_out_topo[t].find(s);
            if (ito != _m_out_topo[t].end()) m_ts += static_cast<double>(ito->second);
            if (_directed) {
                auto iti = _m_in_topo[t].find(s);
                if (iti != _m_in_topo[t].end()) m_ts += static_cast<double>(iti->second);
            }
            prob += (m_ts + c) / (static_cast<double>(kt) + c * B_act);
        }
        return prob / static_cast<double>(adj.size());
    }

    // Export 1-based contiguous block labels
    Rcpp::IntegerVector external_blocks_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        wt_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i) gmap.emplace(occ[i], static_cast<int>(i) + 1);
        Rcpp::IntegerVector out(_N);
        for (int v = 0; v < _N; ++v) out[v] = gmap[_b[v]];
        return out;
    }

private:
    int _N;
    bool _directed, _dc, _dist, _is_base, _is_top;
    std::vector<CovSpec> _covs;
    int _K, _cap;
    long _E = 0;
    double _base_const = 0.0;

    std::vector<WGraphEdge> _edges;
    std::vector<std::vector<WAdjRef>> _adj;
    std::vector<long> _kout, _kin;

    std::vector<int> _b, _wr;
    std::vector<long> _ep, _em;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;

    wt_hash_map<uint64_t, WPairStat> _pairs;
    std::vector<wt_hash_map<int, long>> _m_out_topo, _m_in_topo;

    long _B_E = 0, _B_E_D = 0;
    std::vector<double> _recsum, _recx2, _recdx;
    std::vector<int> _normal_idx;
    std::vector<wt_hash_map<long, long>> _h_out, _h_in;
    mutable MoveBuf _scratch_buf;

    inline double base_eterm(int r, int s, long m_topo) const {
        if (m_topo <= 0) return 0.0;
        const WMath& M = wmath();
        if (!_directed && r == s)
            return M.lgam1(m_topo) + static_cast<double>(m_topo) * W_LOG2;
        return M.lgam1(m_topo);
    }

    inline double dense_eterm(int r, int s, long m_topo, long wr_r, long wr_s) const {
        if (m_topo <= 0) return 0.0;
        int64_t nn = (_directed || r != s) ? (int64_t(wr_r) * wr_s)
                                           : (int64_t(wr_r) * (wr_r + 1) / 2);
        return wmath().lbinom(nn + m_topo - 1, m_topo);
    }

    double deg_dl_group(int r) const {
        if (!_dc || _wr[r] == 0) return 0.0;
        const WMath& M = wmath();
        long n = _wr[r];
        if (!_dist) {
            double s = M.lbinom(n + _ep[r] - 1, _ep[r]);
            if (_directed) s += M.lbinom(n + _em[r] - 1, _em[r]);
            return s;
        }
        const WLogQ& q = wlog_q();
        double s = q(_ep[r], n) + M.lgam1(n);
        for (const auto& [dk, cnt] : _h_out[r])
            s -= M.lgam1(cnt);
        if (_directed) {
            s += q(_em[r], n) + M.lgam1(n);
            for (const auto& [dk, cnt] : _h_in[r])
                s -= M.lgam1(cnt);
        }
        return s;
    }

    double deg_dl_delta_group(int r, long n0, long n1, long ep0, long ep1,
                              long em0, long em1, long ko_v, long ki_v, long dcnt) const {
        if (!_dc) return 0.0;
        const WMath& M = wmath();
        if (!_dist) {
            double s0 = (n0 > 0) ? (M.lbinom(n0 + ep0 - 1, ep0) +
                                    (_directed ? M.lbinom(n0 + em0 - 1, em0) : 0.0)) : 0.0;
            double s1 = (n1 > 0) ? (M.lbinom(n1 + ep1 - 1, ep1) +
                                    (_directed ? M.lbinom(n1 + em1 - 1, em1) : 0.0)) : 0.0;
            return s1 - s0;
        }
        const WLogQ& q = wlog_q();
        double s0 = (n0 > 0) ? (q(ep0, n0) + M.lgam1(n0)) : 0.0;
        double s1 = (n1 > 0) ? (q(ep1, n1) + M.lgam1(n1)) : 0.0;
        auto ito = _h_out[r].find(ko_v);
        long co0 = (ito == _h_out[r].end()) ? 0 : ito->second;
        double d = (s1 - s0) - (M.lgam1(co0 + dcnt) - M.lgam1(co0));
        if (_directed) {
            double si0 = (n0 > 0) ? (q(em0, n0) + M.lgam1(n0)) : 0.0;
            double si1 = (n1 > 0) ? (q(em1, n1) + M.lgam1(n1)) : 0.0;
            auto iti = _h_in[r].find(ki_v);
            long ci0 = (iti == _h_in[r].end()) ? 0 : iti->second;
            d += (si1 - si0) - (M.lgam1(ci0 + dcnt) - M.lgam1(ci0));
        }
        return d;
    }

    template <class WrFn>
    double eval_pair_cov_deltas(MoveBuf& buf, WrFn&& wr_after) const {
        double dS = 0.0;
        long d_B_E = 0, d_B_E_D = 0;
        auto& d_recx2 = buf.d_recx2;
        auto& d_recdx = buf.d_recdx;
        d_recx2.assign(_K, 0.0);
        d_recdx.assign(_K, 0.0);

        for (const auto& pd : buf.pdeltas) {
            auto [a, b_g] = unpack_pair(pd.key);
            auto it = _pairs.find(pd.key);
            long mt0 = (it == _pairs.end()) ? 0 : it->second.m_topo;
            long mc0 = (it == _pairs.end()) ? 0 : it->second.m_cov;
            long mt1 = mt0 + pd.dm_topo;
            long mc1 = mc0 + pd.dm_cov;

            if (_is_base) {
                if (mt0 != mt1)
                    dS -= base_eterm(a, b_g, mt1) - base_eterm(a, b_g, mt0);
            } else {
                dS += dense_eterm(a, b_g, mt1, wr_after(a), wr_after(b_g)) -
                      dense_eterm(a, b_g, mt0, _wr[a], _wr[b_g]);
            }

            if (_K > 0 && (mc0 > 0 || mc1 > 0)) {
                d_B_E += (mc1 > 0 ? 1 : 0) - (mc0 > 0 ? 1 : 0);
                d_B_E_D += (mc1 > 1 ? 1 : 0) - (mc0 > 1 ? 1 : 0);
                for (int k = 0; k < _K; ++k) {
                    double x0 = (it == _pairs.end()) ? 0.0 : it->second.w.x(k);
                    double x20 = (it == _pairs.end()) ? 0.0 : it->second.w.x2(k);
                    double x1 = (mc1 > 0) ? (x0 + pd.dw.x(k)) : 0.0;
                    double x21 = (mc1 > 0) ? (x20 + pd.dw.x2(k)) : 0.0;
                    dS += cov_pair_neg_log_P(_covs[k], static_cast<double>(mc1), x1, x21) -
                          cov_pair_neg_log_P(_covs[k], static_cast<double>(mc0), x0, x20);
                    if (_covs[k].type == CovType::REAL_NORMAL) {
                        d_recx2[k] += (mc1 > 0 ? x1 * x1 : 0.0) - (mc0 > 0 ? x0 * x0 : 0.0);
                        double dx0 = (mc0 > 1) ? std::max(0.0, x20 - x0 * (x0 / static_cast<double>(mc0))) : 0.0;
                        double dx1 = (mc1 > 1) ? std::max(0.0, x21 - x1 * (x1 / static_cast<double>(mc1))) : 0.0;
                        d_recdx[k] += dx1 - dx0;
                    }
                }
            }
        }

        if (_K > 0) {
            long be0 = _B_E, be1 = be0 + d_B_E;
            long bed0 = _B_E_D, bed1 = bed0 + d_B_E_D;
            for (int k = 0; k < _K; ++k) {
                const CovSpec& c = _covs[k];
                if (c.type == CovType::REAL_NORMAL) {
                    if (c.micro) {
                        if (_is_top) {
                            dS += -signed_w_log_P(static_cast<double>(be1), _recsum[k],
                                                  std::max(0.0, _recx2[k] + d_recx2[k]),
                                                  c.p0, c.p1, c.p2, c.p3, c.eps) -
                                  (-signed_w_log_P(static_cast<double>(be0), _recsum[k],
                                                   _recx2[k], c.p0, c.p1, c.p2, c.p3, c.eps));
                        }
                        dS += -positive_w_log_P(static_cast<double>(bed1),
                                                std::max(0.0, _recdx[k] + d_recdx[k]),
                                                c.p2, c.p3, c.eps) -
                              (-positive_w_log_P(static_cast<double>(bed0), _recdx[k],
                                                 c.p2, c.p3, c.eps));
                    }
                } else if (_is_top && c.micro && be0 != be1) {
                    if (c.type == CovType::REAL_EXPONENTIAL) {
                        dS += -positive_w_log_P(static_cast<double>(be1), _recsum[k],
                                                c.p0, c.p1, c.eps) -
                              (-positive_w_log_P(static_cast<double>(be0), _recsum[k],
                                                 c.p0, c.p1, c.eps));
                    } else {
                        dS += -geometric_w_log_P(static_cast<double>(be1), _recsum[k], NAN, NAN) -
                              (-geometric_w_log_P(static_cast<double>(be0), _recsum[k], NAN, NAN));
                    }
                }
            }
        }

        buf.d_B_E = d_B_E;
        buf.d_B_E_D = d_B_E_D;
        return dS;
    }

    double compute_move_delta(int v, int r, int s, MoveBuf* out) const {
        const WMath& M = wmath();
        MoveBuf& buf = out ? *out : _scratch_buf;
        auto& pds = buf.pdeltas;
        auto& sidx = buf.scratch_idx;
        pds.clear();
        bool use_map = false;
        auto find_or_add = [&](uint64_t key) -> PairDelta& {
            if (!use_map) {
                for (size_t i = 0; i < pds.size(); ++i) {
                    if (pds[i].key == key) return pds[i];
                }
                if (pds.size() >= 24) {
                    use_map = true;
                    sidx.clear();
                    for (size_t i = 0; i < pds.size(); ++i) sidx.emplace(pds[i].key, i);
                } else {
                    pds.emplace_back(_K);
                    pds.back().key = key;
                    return pds.back();
                }
            }
            auto it = sidx.find(key);
            if (it == sidx.end()) {
                it = sidx.emplace(key, pds.size()).first;
                pds.emplace_back(_K);
                pds.back().key = key;
            }
            return pds[it->second];
        };
        auto add_pd = [&](uint64_t key, long dmt, long dmc, const WVec& w, double sign) {
            PairDelta& pd = find_or_add(key);
            pd.dm_topo += dmt;
            pd.dm_cov += dmc;
            for (int k = 0; k < _K; ++k) {
                pd.dw.x(k) += sign * w.x(k);
                pd.dw.x2(k) += sign * w.x2(k);
            }
        };

        long d_ep_r = 0, d_ep_s = 0, d_em_r = 0, d_em_s = 0;
        for (const auto& a : _adj[v]) {
            const auto& e = _edges[a.eidx];
            if (e.u == e.v) {
                // Self-loop on v (appears once in _adj[v] for undir, twice for dir!)
                if (_directed && !a.out) continue;  // process directed self-loop once
                add_pd(canon_pair(r, r, _directed), -e.m_topo, -e.m_cov, e.w, -1.0);
                add_pd(canon_pair(s, s, _directed), +e.m_topo, +e.m_cov, e.w, +1.0);
                if (_directed) {
                    d_ep_r -= e.m_topo; d_em_r -= e.m_topo;
                    d_ep_s += e.m_topo; d_em_s += e.m_topo;
                } else {
                    d_ep_r -= 2 * e.m_topo;
                    d_ep_s += 2 * e.m_topo;
                }
            } else {
                int t = _b[a.nbr];
                if (_directed) {
                    if (a.out) {
                        add_pd(pack_pair(r, t), -e.m_topo, -e.m_cov, e.w, -1.0);
                        add_pd(pack_pair(s, t), +e.m_topo, +e.m_cov, e.w, +1.0);
                        d_ep_r -= e.m_topo;
                        d_ep_s += e.m_topo;
                    } else {
                        add_pd(pack_pair(t, r), -e.m_topo, -e.m_cov, e.w, -1.0);
                        add_pd(pack_pair(t, s), +e.m_topo, +e.m_cov, e.w, +1.0);
                        d_em_r -= e.m_topo;
                        d_em_s += e.m_topo;
                    }
                } else {
                    add_pd(canon_pair(r, t, false), -e.m_topo, -e.m_cov, e.w, -1.0);
                    add_pd(canon_pair(s, t, false), +e.m_topo, +e.m_cov, e.w, +1.0);
                    d_ep_r -= e.m_topo;
                    d_ep_s += e.m_topo;
                }
            }
        }

        long wr_r0 = _wr[r], wr_r1 = wr_r0 - 1;
        long wr_s0 = _wr[s], wr_s1 = wr_s0 + 1;
        long B0 = B();
        long B1 = B0 - (wr_r1 == 0 ? 1 : 0) + (wr_s0 == 0 ? 1 : 0);

        double dS = 0.0;
        // Partition DL + edges_dl
        if (B0 != B1) {
            dS += M.lbinom(_N - 1, B1 - 1) - M.lbinom(_N - 1, B0 - 1);
            if (_is_top && _E > 0) {
                int64_t p0 = _directed ? int64_t(B0) * B0 : int64_t(B0) * (B0 + 1) / 2;
                int64_t p1 = _directed ? int64_t(B1) * B1 : int64_t(B1) * (B1 + 1) / 2;
                dS += M.lbinom(p1 + _E - 1, _E) - M.lbinom(p0 + _E - 1, _E);
            }
        }
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_s1) - M.lgam1(wr_s0));

        if (_is_base) {
            long ep_r0 = _ep[r], ep_r1 = ep_r0 + d_ep_r;
            long ep_s0 = _ep[s], ep_s1 = ep_s0 + d_ep_s;
            long em_r0 = _directed ? _em[r] : 0, em_r1 = em_r0 + d_em_r;
            long em_s0 = _directed ? _em[s] : 0, em_s1 = em_s0 + d_em_s;
            if (_dc) {
                dS += (M.lgam1(ep_r1) - M.lgam1(ep_r0)) + (M.lgam1(ep_s1) - M.lgam1(ep_s0));
                if (_directed)
                    dS += (M.lgam1(em_r1) - M.lgam1(em_r0)) + (M.lgam1(em_s1) - M.lgam1(em_s0));
                long ko_v = _kout[v];
                long ki_v = _directed ? _kin[v] : 0;
                dS += deg_dl_delta_group(r, wr_r0, wr_r1, ep_r0, ep_r1, em_r0, em_r1, ko_v, ki_v, -1);
                dS += deg_dl_delta_group(s, wr_s0, wr_s1, ep_s0, ep_s1, em_s0, em_s1, ko_v, ki_v, +1);
            } else {
                long dr0 = ep_r0 + em_r0, dr1 = ep_r1 + em_r1;
                long ds0 = ep_s0 + em_s0, ds1 = ep_s1 + em_s1;
                dS += (wr_r1 > 0 ? dr1 * M.logi(wr_r1) : 0.0) - dr0 * M.logi(wr_r0);
                dS += ds1 * M.logi(wr_s1) - (wr_s0 > 0 ? ds0 * M.logi(wr_s0) : 0.0);
            }
        } else {
            // Dense upper level: pairs involving r or s whose wr changes even if m_topo didn't!
            // Collect all existing pairs incident on r or s that are NOT yet in pds
            auto touch_dense = [&](int g) {
                for (const auto& [h, m2] : _m_out_topo[g]) {
                    find_or_add(canon_pair(g, h, _directed));
                }
                if (_directed) {
                    for (const auto& [h, m2] : _m_in_topo[g]) {
                        find_or_add(pack_pair(h, g));
                    }
                }
            };
            touch_dense(r);
            touch_dense(s);
        }

        auto wr_after = [&](int g) -> long {
            if (g == r) return wr_r1;
            if (g == s) return wr_s1;
            return _wr[g];
        };
        dS += eval_pair_cov_deltas(buf, wr_after);

        buf.d_ep_r = d_ep_r;
        buf.d_ep_s = d_ep_s;
        buf.d_em_r = d_em_r;
        buf.d_em_s = d_em_s;
        return dS;
    }

    void apply_move(int v, int r, int s, const MoveBuf& buf) {
        if (_dc && _dist) {
            auto ito = _h_out[r].find(_kout[v]);
            if (ito != _h_out[r].end() && --ito->second == 0) _h_out[r].erase(ito);
            _h_out[s][_kout[v]]++;
            if (_directed) {
                auto iti = _h_in[r].find(_kin[v]);
                if (iti != _h_in[r].end() && --iti->second == 0) _h_in[r].erase(iti);
                _h_in[s][_kin[v]]++;
            }
        }
        remove_member(v, r);
        _b[v] = s;
        add_member(v, s);

        _ep[r] += buf.d_ep_r;
        _ep[s] += buf.d_ep_s;
        if (_directed) {
            _em[r] += buf.d_em_r;
            _em[s] += buf.d_em_s;
        }
        _B_E += buf.d_B_E;
        _B_E_D += buf.d_B_E_D;
        for (int k : _normal_idx) {
            _recx2[k] = std::max(0.0, _recx2[k] + buf.d_recx2[k]);
            _recdx[k] = std::max(0.0, _recdx[k] + buf.d_recdx[k]);
        }

        for (const auto& pd : buf.pdeltas) {
            auto [a, b_g] = unpack_pair(pd.key);
            auto it = _pairs.find(pd.key);
            if (it == _pairs.end()) {
                if (pd.dm_topo <= 0 && pd.dm_cov <= 0) continue;
                it = _pairs.emplace(pd.key, WPairStat(_K)).first;
            }
            WPairStat& ps = it->second;
            ps.m_topo += pd.dm_topo;
            ps.m_cov += pd.dm_cov;
            for (int k = 0; k < _K; ++k) {
                ps.w.x(k) += pd.dw.x(k);
                ps.w.x2(k) += pd.dw.x2(k);
            }
            if (ps.m_topo <= 0 && ps.m_cov <= 0)
                _pairs.erase(it);

            // Update fast adjacency count maps for proposals
            auto add_m = [](wt_hash_map<int, long>& mp, int key_g, long d) {
                if (d == 0) return;
                auto mit = mp.find(key_g);
                if (mit == mp.end()) {
                    if (d > 0) mp.emplace(key_g, d);
                } else {
                    mit->second += d;
                    if (mit->second <= 0) mp.erase(mit);
                }
            };
            if (_directed) {
                add_m(_m_out_topo[a], b_g, pd.dm_topo);
                add_m(_m_in_topo[b_g], a, pd.dm_topo);
            } else {
                if (a == b_g) {
                    add_m(_m_out_topo[a], a, 2 * pd.dm_topo);
                } else {
                    add_m(_m_out_topo[a], b_g, pd.dm_topo);
                    add_m(_m_out_topo[b_g], a, pd.dm_topo);
                }
            }
        }
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

    void add_member(int v, int r) {
        if (_wr[r] == 0) {
            take_empty(r);
            _occ_pos[r] = static_cast<int>(_occupied.size());
            _occupied.push_back(r);
        }
        _pos[v] = static_cast<int>(_members[r].size());
        _members[r].push_back(v);
        _wr[r]++;
    }

    void remove_member(int v, int r) {
        int p = _pos[v];
        int last = _members[r].back();
        _members[r][p] = last;
        _pos[last] = p;
        _members[r].pop_back();
        if (--_wr[r] == 0) {
            int op = _occ_pos[r];
            int lo = _occupied.back();
            _occupied[op] = lo;
            _occ_pos[lo] = op;
            _occupied.pop_back();
            _occ_pos[r] = -1;
            push_empty(r);
        }
    }
};

std::vector<WGraphEdge> build_base_edges(const Graph& graph,
                                         const Rcpp::NumericMatrix& weights) {
    int K = weights.ncol();
    size_t E = boost::num_edges(graph.adjacency);
    if (static_cast<size_t>(weights.nrow()) != E)
        Rcpp::stop("weights must have nrow(weights) == graph_num_edges(graph)");
    std::vector<WGraphEdge> edges;
    edges.reserve(E);
    auto range = boost::edges(graph.adjacency);
    size_t idx = 0;
    for (auto it = range.first; it != range.second; ++it, ++idx) {
        WGraphEdge ge;
        ge.u = static_cast<int>(it->s);
        ge.v = static_cast<int>(it->t);
        ge.m_topo = 1;
        ge.m_cov = 1;
        ge.w.zero(K);
        for (int k = 0; k < K; ++k) {
            double val = weights(static_cast<R_xlen_t>(idx), k);
            ge.w.x(k) = val;
            ge.w.x2(k) = val * val;
        }
        edges.push_back(std::move(ge));
    }
    return edges;
}

// Build upper-level edge list from lower level's non-empty block pairs (reindexed 0..B_lower-1).
std::vector<WGraphEdge> build_upper_edges(const WeightedLevel& lower,
                                          wt_hash_map<int, int>& lower_to_node) {
    std::vector<int> occ = lower.occupied();
    std::sort(occ.begin(), occ.end());
    lower_to_node.clear();
    lower_to_node.reserve(occ.size());
    for (size_t i = 0; i < occ.size(); ++i)
        lower_to_node.emplace(occ[i], static_cast<int>(i));

    int K = lower.K();
    std::vector<WGraphEdge> edges;
    edges.reserve(lower.pairs().size());
    for (const auto& [key, ps] : lower.pairs()) {
        if (ps.m_topo <= 0) continue;
        auto [a, b] = unpack_pair(key);
        WGraphEdge ge;
        ge.u = lower_to_node.at(a);
        ge.v = lower_to_node.at(b);
        ge.m_topo = ps.m_topo;
        ge.m_cov = 1;  // each non-empty lower block pair is 1 covariate observation at l >= 1
        ge.w.zero(K);
        for (int k = 0; k < K; ++k) {
            double x = ps.w.x(k);
            ge.w.x(k) = x;
            ge.w.x2(k) = x * x;
        }
        edges.push_back(std::move(ge));
    }
    return edges;
}

// Parse 1-based R partition into 0-based sorted-label vector
std::vector<int> parse_w_blocks(Rcpp::IntegerVector blocks, int N) {
    if (blocks.size() != N) Rcpp::stop("blocks must have length N");
    std::vector<int> uniq;
    uniq.reserve(N);
    for (int i = 0; i < N; ++i) {
        if (blocks[i] == NA_INTEGER || blocks[i] < 1)
            Rcpp::stop("blocks must contain positive integer labels");
        uniq.push_back(blocks[i]);
    }
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    wt_hash_map<int, int> mp;
    for (size_t i = 0; i < uniq.size(); ++i) mp.emplace(uniq[i], static_cast<int>(i));
    std::vector<int> out(N);
    for (int i = 0; i < N; ++i) out[i] = mp[blocks[i]];
    return out;
}

// ---------------------------------------------------------------------------
// Flat WeightedSBM minimization and refinement (agglomerative + greedy sweeps)

double weighted_refine_flat(WeightedLevel& s, WRng& rng, size_t B_min,
                            size_t B_max, double epsilon, size_t max_iter = 30,
                            int nthreads = 0) {
    if (s.N() <= 1 || s.E() == 0) return s.entropy();
    int nw = rgt_parallel::resolve_threads(nthreads);
    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    double S = s.entropy();

    for (size_t iter = 0; iter < max_iter; ++iter) {
        // 1. Greedy single-node sweeps (batch-parallel candidate screening for N >= 512)
        for (int sw = 0; sw < 10; ++sw) {
            rng.shuffle(vs);
            double dS_sw = 0.0;
            if (nw > 1 && s.N() >= 512 && s.B() > 1) {
                size_t N = static_cast<size_t>(s.N());
                size_t batch = std::min(N, std::max<size_t>(512, static_cast<size_t>(nw) * 256));
                std::vector<WeightedLevel::MoveBuf> bufs(static_cast<size_t>(nw));
                std::vector<wt_hash_set<int>> cand_sets(static_cast<size_t>(nw));
                std::vector<int> cand_t;
                for (size_t base = 0; base < N; base += batch) {
                    size_t len = std::min(batch, N - base);
                    cand_t.assign(len, -1);
                    uint64_t batch_seed = static_cast<uint64_t>(rng.index(1ULL << 30)) ^
                                          (static_cast<uint64_t>(rng.index(1ULL << 30)) << 30) ^
                                          0x9e3779b97f4a7c15ULL;
                    const size_t B_snap = static_cast<size_t>(s.B());
                    rgt_parallel::parallel_for(len, nw, [&](int tid, size_t lo, size_t hi) {
                        rgt_parallel::ThreadRng trng(
                            batch_seed ^ (static_cast<uint64_t>(lo + 1) * 0xbf58476d1ce4e5b9ULL));
                        auto& buf = bufs[static_cast<size_t>(tid)];
                        auto& cands = cand_sets[static_cast<size_t>(tid)];
                        for (size_t idx = lo; idx < hi; ++idx) {
                            int v = vs[base + idx];
                            int r = s.group(v);
                            cands.clear();
                            for (int k = 0; k < 6; ++k) cands.insert(s.sample_group(v, 0.5, trng));
                            if (B_snap <= 20) {
                                for (int g : s.occupied()) cands.insert(g);
                            }
                            if (B_snap < B_max && s.size(r) > 1 && s.empty_group() >= 0)
                                cands.insert(s.empty_group());
                            cands.erase(r);
                            int best = -1;
                            double best_dS = -1e-9;
                            for (int t : cands) {
                                if (s.size(r) == 1 && B_snap <= B_min && s.size(t) > 0) continue;
                                double d = s.virtual_move_const(v, t, buf);
                                if (d < best_dS) { best_dS = d; best = t; }
                            }
                            cand_t[idx] = best;
                        }
                    }, 32);

                    for (size_t idx = 0; idx < len; ++idx) {
                        int t = cand_t[idx];
                        if (t < 0) continue;
                        int v = vs[base + idx];
                        int r = s.group(v);
                        if (t == r) continue;
                        size_t B_cur = static_cast<size_t>(s.B());
                        if (s.size(r) == 1 && s.size(t) > 0 && B_cur <= B_min) continue;
                        if (s.size(t) == 0 && (B_cur >= B_max && s.size(r) > 1)) continue;
                        double d = s.virtual_move(v, t);
                        if (d < -1e-9) dS_sw += s.move(v, t);
                    }
                }
            } else {
                for (int v : vs) {
                    int r = s.group(v);
                    size_t B_cur = static_cast<size_t>(s.B());
                    wt_hash_set<int> cands;
                    for (int k = 0; k < 6; ++k) cands.insert(s.sample_group(v, 0.5, rng));
                    if (B_cur <= 20) {
                        for (int g : s.occupied()) cands.insert(g);
                    }
                    if (B_cur < B_max && s.size(r) > 1 && s.empty_group() >= 0)
                        cands.insert(s.empty_group());
                    cands.erase(r);
                    int best = r;
                    double best_dS = -1e-9;
                    for (int t : cands) {
                        if (s.size(r) == 1 && B_cur <= B_min && s.size(t) > 0) continue;
                        double d = s.virtual_move(v, t);
                        if (d < best_dS) { best_dS = d; best = t; }
                    }
                    if (best != r) dS_sw += s.move(v, best);
                }
            }
            if (std::abs(dS_sw) < 1e-8) break;
        }

        // 2. Strict greedy merges (parallel candidate evaluation for occ.size() >= 32)
        while (static_cast<size_t>(s.B()) > B_min) {
            std::vector<int> occ = s.occupied();
            int mr = -1, mt = -1;
            double best_dS = -1e-9;
            if (occ.size() <= 40) {
                for (size_t i = 0; i < occ.size(); ++i) {
                    for (size_t j = i + 1; j < occ.size(); ++j) {
                        double d = s.virtual_merge_const(occ[i], occ[j], s.default_buf());
                        if (d < best_dS) { best_dS = d; mr = occ[i]; mt = occ[j]; }
                    }
                }
            } else if (nw > 1) {
                size_t B_occ = occ.size();
                std::vector<int> cand_t(B_occ, -1);
                std::vector<double> cand_val(B_occ, 0.0);
                std::vector<WeightedLevel::MoveBuf> bufs(static_cast<size_t>(nw));
                std::vector<wt_hash_set<int>> tried_sets(static_cast<size_t>(nw));
                uint64_t m_seed = static_cast<uint64_t>(rng.index(1ULL << 30)) ^
                                  (static_cast<uint64_t>(rng.index(1ULL << 30)) << 30) ^
                                  0xbf58476d1ce4e5b9ULL;
                rgt_parallel::parallel_for(B_occ, nw, [&](int tid, size_t lo, size_t hi) {
                    rgt_parallel::ThreadRng trng(
                        m_seed ^ (static_cast<uint64_t>(lo + 1) * 0x94d049bb133111ebULL));
                    auto& buf = bufs[static_cast<size_t>(tid)];
                    auto& tried = tried_sets[static_cast<size_t>(tid)];
                    for (size_t idx = lo; idx < hi; ++idx) {
                        int r = occ[idx];
                        const auto& mem = s.members(r);
                        tried.clear();
                        int best_loc = -1;
                        double best_loc_dS = -1e-9;
                        for (int k = 0; k < 8; ++k) {
                            int v = mem[trng.index(mem.size())];
                            int t = s.sample_group(v, k < 4 ? 0.0 : 0.5, trng);
                            if (t == r || tried.count(t)) continue;
                            tried.insert(t);
                            double d = s.virtual_merge_const(r, t, buf);
                            if (d < best_loc_dS) { best_loc_dS = d; best_loc = t; }
                        }
                        if (best_loc >= 0) {
                            cand_t[idx] = best_loc;
                            cand_val[idx] = best_loc_dS;
                        }
                    }
                }, 4);
                for (size_t idx = 0; idx < B_occ; ++idx) {
                    if (cand_t[idx] >= 0 && cand_val[idx] < best_dS) {
                        best_dS = cand_val[idx];
                        mr = occ[idx];
                        mt = cand_t[idx];
                    }
                }
            } else {
                for (int r : occ) {
                    const auto& mem = s.members(r);
                    wt_hash_set<int> tried;
                    for (int k = 0; k < 8; ++k) {
                        int v = mem[rng.index(mem.size())];
                        int t = s.sample_group(v, k < 4 ? 0.0 : 0.5, rng);
                        if (t == r || tried.count(t)) continue;
                        tried.insert(t);
                        double d = s.virtual_merge(r, t);
                        if (d < best_dS) { best_dS = d; mr = r; mt = t; }
                    }
                }
            }
            if (mr < 0) break;
            s.merge(mr, mt);
        }

        // 3. Try splitting each group into 2
        if (static_cast<size_t>(s.B()) < B_max) {
            std::vector<int> occ = s.occupied();
            for (int r : occ) {
                if (static_cast<size_t>(s.B()) >= B_max) break;
                if (s.size(r) < 2 || s.empty_group() < 0) continue;
                int t = s.empty_group();
                std::vector<int> mem = s.members(r);
                double dS_sp = 0.0;
                for (int v : mem)
                    if (rng.unif() < 0.5) dS_sp += s.move(v, t);
                if (s.size(r) == 0 || s.size(t) == 0) {
                    for (int v : mem) if (s.group(v) != r) s.move(v, r);
                    continue;
                }
                for (int sw = 0; sw < 5; ++sw) {
                    rng.shuffle(mem);
                    bool ch = false;
                    for (int v : mem) {
                        int cur = s.group(v);
                        int oth = (cur == r) ? t : r;
                        if (s.size(cur) <= 1) continue;
                        double d = s.virtual_move(v, oth);
                        if (d < -1e-9) { dS_sp += s.move(v, oth); ch = true; }
                    }
                    if (!ch) break;
                }
                if (!(dS_sp < -1e-8 && s.size(r) > 0 && s.size(t) > 0)) {
                    for (int v : mem) if (s.group(v) != r) s.move(v, r);
                }
            }
        }

        s.refresh_cov_globals();
        Rcpp::checkUserInterrupt();
        double S_new = s.entropy();
        if (S_new == 0.0 || std::abs((S_new - S) / std::max(1.0, std::abs(S))) <= epsilon) {
            S = S_new;
            break;
        }
        S = S_new;
    }
    return S;
}

void weighted_coarsen_to(WeightedLevel& s, WRng& rng, size_t B_target,
                         double shrink = 0.8, int nthreads = 0) {
    int nw = rgt_parallel::resolve_threads(nthreads);
    while (static_cast<size_t>(s.B()) > B_target) {
        size_t B_cur = static_cast<size_t>(s.B());
        size_t B_next = std::max(B_target, std::min(B_cur - 1, static_cast<size_t>(std::round(B_cur * shrink))));
        std::vector<int> occ = s.occupied();
        rng.shuffle(occ);
        std::vector<std::tuple<double, int, int>> cands;
        cands.reserve(occ.size());
        if (nw > 1 && occ.size() >= 32) {
            size_t B_occ = occ.size();
            std::vector<int> best_tgt(B_occ, -1);
            std::vector<double> best_val(B_occ, W_INF);
            std::vector<WeightedLevel::MoveBuf> bufs(static_cast<size_t>(nw));
            std::vector<wt_hash_set<int>> tried_sets(static_cast<size_t>(nw));
            uint64_t c_seed = static_cast<uint64_t>(rng.index(1ULL << 30)) ^
                              (static_cast<uint64_t>(rng.index(1ULL << 30)) << 30) ^
                              0xbf58476d1ce4e5b9ULL;
            rgt_parallel::parallel_for(B_occ, nw, [&](int tid, size_t lo, size_t hi) {
                rgt_parallel::ThreadRng trng(
                    c_seed ^ (static_cast<uint64_t>(lo + 1) * 0x94d049bb133111ebULL));
                auto& buf = bufs[static_cast<size_t>(tid)];
                auto& tried = tried_sets[static_cast<size_t>(tid)];
                for (size_t idx = lo; idx < hi; ++idx) {
                    int r = occ[idx];
                    const auto& mem = s.members(r);
                    int best = -1;
                    double best_dS = W_INF;
                    tried.clear();
                    for (int k = 0; k < 10; ++k) {
                        int v = mem[trng.index(mem.size())];
                        int t = s.sample_group(v, k < 5 ? 0.0 : 0.5, trng);
                        if (t == r || tried.count(t)) continue;
                        tried.insert(t);
                        double d = s.virtual_merge_const(r, t, buf);
                        if (d < best_dS) { best_dS = d; best = t; }
                    }
                    if (best < 0 && B_occ > 1) {
                        int t = occ[trng.index(B_occ - 1)];
                        if (t == r) t = occ.back();
                        best = t;
                        best_dS = s.virtual_merge_const(r, t, buf);
                    }
                    if (best >= 0) {
                        best_tgt[idx] = best;
                        best_val[idx] = best_dS;
                    }
                }
            }, 8);
            for (size_t idx = 0; idx < B_occ; ++idx) {
                if (best_tgt[idx] >= 0) cands.emplace_back(best_val[idx], occ[idx], best_tgt[idx]);
            }
        } else {
            for (int r : occ) {
                const auto& mem = s.members(r);
                int best = -1;
                double best_dS = W_INF;
                wt_hash_set<int> tried;
                for (int k = 0; k < 10; ++k) {
                    int v = mem[rng.index(mem.size())];
                    int t = s.sample_group(v, k < 5 ? 0.0 : 0.5, rng);
                    if (t == r || tried.count(t)) continue;
                    tried.insert(t);
                    double d = s.virtual_merge(r, t);
                    if (d < best_dS) { best_dS = d; best = t; }
                }
                if (best < 0 && occ.size() > 1) {
                    int t = occ[rng.index(occ.size() - 1)];
                    if (t == r) t = occ.back();
                    best = t;
                    best_dS = s.virtual_merge(r, t);
                }
                if (best >= 0) cands.emplace_back(best_dS, r, best);
            }
        }
        std::sort(cands.begin(), cands.end());
        wt_hash_map<int, int> root;
        auto get_root = [&](int x) {
            int cur = x;
            while (root.count(cur) && root[cur] != cur) cur = root[cur];
            root[x] = cur;
            return cur;
        };
        size_t merged = 0;
        for (const auto& [dS_est, a, b_raw] : cands) {
            if (static_cast<size_t>(s.B()) <= B_next) break;
            int t = get_root(b_raw);
            if (a == t || s.size(a) == 0 || s.size(t) == 0) continue;
            s.merge(a, t);
            root[a] = t;
            ++merged;
        }
        if (merged == 0) break;
        Rcpp::checkUserInterrupt();
    }
}

double weighted_minimize_flat(WeightedLevel& s, WRng& rng, size_t B_min,
                              size_t B_max, bool refine, double epsilon,
                              int nthreads = 0) {
    if (s.N() <= 1 || s.E() == 0) {
        s.set_partition(std::vector<int>(s.N(), 0));
        return s.entropy();
    }
    int nw = rgt_parallel::resolve_threads(nthreads);
    size_t B_hi = std::min<size_t>(B_max, static_cast<size_t>(s.N()));
    size_t B_lo = std::min(std::max<size_t>(1, B_min), B_hi);

    // Start from singleton partition and coarsen along geometric schedule
    std::vector<int> b_init(s.N());
    std::iota(b_init.begin(), b_init.end(), 0);
    s.set_partition(b_init);
    if (static_cast<size_t>(s.B()) > B_hi)
        weighted_coarsen_to(s, rng, B_hi, 0.8, nw);

    double best_S = s.entropy();
    std::vector<int> best_b = s.blocks();

    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    std::vector<WeightedLevel::MoveBuf> sw_bufs(static_cast<size_t>(std::max(1, nw)));
    while (static_cast<size_t>(s.B()) > B_lo) {
        size_t B_cur = static_cast<size_t>(s.B());
        size_t B_next = std::max(B_lo, std::min(B_cur - 1, static_cast<size_t>(std::round(B_cur * 0.8))));
        weighted_coarsen_to(s, rng, B_next, 0.8, nw);
        // Quick greedy single-node sweep at each coarsening stop
        rng.shuffle(vs);
        if (nw > 1 && s.N() >= 512 && s.B() > 1) {
            size_t N = static_cast<size_t>(s.N());
            size_t batch = std::min(N, std::max<size_t>(512, static_cast<size_t>(nw) * 256));
            std::vector<int> cand_t;
            for (size_t base = 0; base < N; base += batch) {
                size_t len = std::min(batch, N - base);
                cand_t.assign(len, -1);
                uint64_t b_seed = static_cast<uint64_t>(rng.index(1ULL << 30)) ^
                                  (static_cast<uint64_t>(rng.index(1ULL << 30)) << 30) ^
                                  0x9e3779b97f4a7c15ULL;
                rgt_parallel::parallel_for(len, nw, [&](int tid, size_t lo, size_t hi) {
                    rgt_parallel::ThreadRng trng(
                        b_seed ^ (static_cast<uint64_t>(lo + 1) * 0xbf58476d1ce4e5b9ULL));
                    auto& buf = sw_bufs[static_cast<size_t>(tid)];
                    for (size_t idx = lo; idx < hi; ++idx) {
                        int v = vs[base + idx];
                        int r = s.group(v);
                        int t = s.sample_group(v, 0.5, trng);
                        if (t != r && s.virtual_move_const(v, t, buf) < -1e-9) {
                            cand_t[idx] = t;
                        }
                    }
                }, 64);
                for (size_t idx = 0; idx < len; ++idx) {
                    int t = cand_t[idx];
                    if (t < 0) continue;
                    int v = vs[base + idx];
                    int r = s.group(v);
                    if (t == r || (s.size(r) <= 1 && static_cast<size_t>(s.B()) <= B_lo)) continue;
                    if (s.virtual_move(v, t) < -1e-9) s.move(v, t);
                }
            }
        } else {
            for (int v : vs) {
                int r = s.group(v);
                if (s.size(r) <= 1 && static_cast<size_t>(s.B()) <= B_lo) continue;
                int t = s.sample_group(v, 0.5, rng);
                if (t != r && s.virtual_move(v, t) < -1e-9) s.move(v, t);
            }
        }
        s.refresh_cov_globals();
        double cur_S = s.entropy();
        if (cur_S < best_S) {
            best_S = cur_S;
            best_b = s.blocks();
        }
    }

    // Refine best agglomerative state
    s.set_partition(best_b);
    if (refine) {
        best_S = weighted_refine_flat(s, rng, B_lo, B_hi, epsilon, 30, nw);
        best_b = s.blocks();
    }

    // Also try single-group start (with split refinement) when B_lo == 1
    if (B_lo == 1) {
        s.set_partition(std::vector<int>(s.N(), 0));
        double s1 = refine ? weighted_refine_flat(s, rng, B_lo, B_hi, epsilon, 15, nw) : s.entropy();
        if (s1 < best_S) {
            best_S = s1;
            best_b = s.blocks();
        }
    }

    s.set_partition(best_b);
    return s.entropy();
}

// ---------------------------------------------------------------------------
// Hierarchical NestedWeightedSBM (exact match to graph-tool 3.9 NestedBlockState
// with WeightedBlockState base_state, including Lrecdx for real-normal).

class NestedWeightedSBM {
public:
    NestedWeightedSBM(const Graph& graph, const Rcpp::NumericMatrix& weights,
                      const std::vector<CovSpec>& base_covs,
                      bool dc, bool dist)
        : _directed(graph.directed), _dc(dc), _dist(dist),
          _base_covs(base_covs), _up_covs(make_upper_covs(base_covs)),
          _base_N(static_cast<int>(boost::num_vertices(graph.adjacency))),
          _base_edges(build_base_edges(graph, weights)) {}

    // Build all L levels from 0-based level partitions bs[0..L-1]
    void build(const std::vector<std::vector<int>>& bs) {
        _levels.clear();
        int L = static_cast<int>(bs.size());
        if (L < 1) Rcpp::stop("Nested weighted state requires at least 1 level");
        _levels.reserve(L);
        for (int l = 0; l < L; ++l) {
            bool is_base = (l == 0);
            bool is_top = (l == L - 1);
            if (is_base) {
                _levels.emplace_back(_base_N, _directed, _dc, _dist, true, is_top, _base_covs);
                _levels.back().set_edges(_base_edges);
                _levels.back().set_partition(bs[0]);
            } else {
                wt_hash_map<int, int> lmap;
                std::vector<WGraphEdge> uedges = build_upper_edges(_levels[l - 1], lmap);
                int Nl = static_cast<int>(lmap.size());
                _levels.emplace_back(Nl, _directed, false, false, false, is_top, _up_covs);
                _levels.back().set_edges(std::move(uedges));
                if (static_cast<int>(bs[l].size()) != Nl)
                    Rcpp::stop("Upper level %d partition length (%d) does not match lower group count (%d)",
                               l, static_cast<int>(bs[l].size()), Nl);
                _levels.back().set_partition(bs[l]);
            }
        }
    }

    int levels() const { return static_cast<int>(_levels.size()); }
    WeightedLevel& level(int l) { return _levels[l]; }
    const WeightedLevel& level(int l) const { return _levels[l]; }

    // Exact per-level entropies (including graph-tool 3.9's Lrecdx nested normal term!)
    std::vector<double> level_entropies() const {
        int L = levels();
        std::vector<double> out(L, 0.0);
        // Precompute Lrecdx across all levels (active only when L > 1 and REAL_NORMAL exists)
        long L_active = 0;
        std::vector<double> L_recdx(_base_covs.size(), 0.0);
        bool has_normal = !_levels[0].normal_idx().empty();
        if (L > 1 && has_normal) {
            for (int l = 0; l < L; ++l) {
                long bed = _levels[l].B_E_D();
                if (bed > 0) ++L_active;
                for (int k : _levels[0].normal_idx())
                    L_recdx[k] += _levels[l].recdx(k) * static_cast<double>(bed);
            }
        }
        for (int l = 0; l < L; ++l) {
            double s = _levels[l].entropy(true, l == L - 1);
            if (L > 1 && has_normal) {
                for (int k : _levels[0].normal_idx()) {
                    s += -w_safelog(static_cast<double>(_levels[l].B_E_D()));
                    if (l == 0) {
                        const CovSpec& c = _base_covs[k];
                        s += -positive_w_log_P(static_cast<double>(L_active), L_recdx[k],
                                               c.p2, c.p3, c.eps);
                    }
                }
            }
            out[l] = s;
        }
        return out;
    }

    double entropy() const {
        auto lv = level_entropies();
        return std::accumulate(lv.begin(), lv.end(), 0.0);
    }

    // Reconstruct contiguous 1-based R list of partitions
    Rcpp::List external_bs_r() const {
        int L = levels();
        Rcpp::List out(L);
        for (int l = 0; l < L; ++l) out[l] = _levels[l].external_blocks_r();
        return out;
    }

    // Build upper levels greedily from a fixed base partition b0 (and trim redundant top 1-block levels)
    void optimize_from_base(const std::vector<int>& b0, WRng& rng, bool refine,
                            double epsilon, int nthreads = 0) {
        std::vector<std::vector<int>> bs;
        bs.push_back(b0);
        // Repeatedly minimize upper level on lower block graph until B == 1
        WeightedLevel cur_base(_base_N, _directed, _dc, _dist, true, false, _base_covs);
        cur_base.set_edges(_base_edges);
        cur_base.set_partition(b0);

        WeightedLevel* prev = &cur_base;
        std::vector<WeightedLevel> built_ups;
        built_ups.reserve(16);
        while (prev->B() > 1) {
            wt_hash_map<int, int> lmap;
            std::vector<WGraphEdge> uedges = build_upper_edges(*prev, lmap);
            int Nl = static_cast<int>(lmap.size());
            // Minimize upper level with edges_dl=true (as tentative top)
            built_ups.emplace_back(Nl, _directed, false, false, false, true, _up_covs);
            built_ups.back().set_edges(std::move(uedges));
            weighted_minimize_flat(built_ups.back(), rng, 1, static_cast<size_t>(Nl),
                                   refine, epsilon, nthreads);
            // Map labels to contiguous 0..B_l-1
            Rcpp::IntegerVector eb = built_ups.back().external_blocks_r();
            std::vector<int> bl(Nl);
            for (int i = 0; i < Nl; ++i) bl[i] = eb[i] - 1;
            int Bl = built_ups.back().B();
            if (Bl == Nl) {
                // Coarsen to 1 group at top so hierarchy terminates cleanly
                bs.emplace_back(Nl, 0);
                break;
            }
            bs.push_back(bl);
            built_ups.back().set_partition(bl);
            prev = &built_ups.back();
        }
        if (bs.size() == 1) {
            // Ensure at least 2 levels for nested state (top single-group level)
            bs.emplace_back(cur_base.B(), 0);
        }
        build(bs);
    }

    // Joint greedy refinement across levels (level 0 single-node & merge moves evaluated
    // under exact total nested entropy + upper-level re-optimization).
    void refine_joint(WRng& rng, size_t B_min, size_t B_max, double epsilon) {
        if (_base_N <= 1 || _base_edges.empty()) return;
        double best_S = entropy();
        std::vector<std::vector<int>> best_bs = save_bs();
        std::vector<int> vs(_base_N);
        std::iota(vs.begin(), vs.end(), 0);

        for (int pass = 0; pass < 4; ++pass) {
            bool improved = false;
            // 1. Level-0 greedy sweeps using fast base level proposal + exact nested check
            //    (Rebuild upper levels only for candidate moves that pass a loose base filter)
            WeightedLevel& lv0 = _levels[0];
            rng.shuffle(vs);
            for (int v : vs) {
                int r = lv0.group(v);
                size_t B_cur = static_cast<size_t>(lv0.B());
                wt_hash_set<int> cands;
                for (int k = 0; k < 4; ++k) cands.insert(lv0.sample_group(v, 0.5, rng));
                cands.erase(r);
                int best_t = r;
                double best_cand_S = best_S - 1e-8;
                for (int t : cands) {
                    if (lv0.size(r) == 1 && B_cur <= B_min) continue;
                    // Fast base-level filter: skip obviously terrible moves
                    double d0 = lv0.virtual_move(v, t);
                    if (d0 > 25.0) continue;
                    int old_B0 = lv0.B();
                    lv0.move(v, t);
                    bool changed_B0 = (lv0.B() != old_B0);
                    double cand_S = W_INF;
                    if (!changed_B0) {
                        // Rebuild upper levels keeping current upper assignments!
                        rebuild_upper_keep_bs();
                        cand_S = entropy();
                    } else {
                        // Number of level-0 blocks changed: rebuild upper via optimize
                        std::vector<int> b0_cur = lv0.blocks();
                        Rcpp::IntegerVector eb0 = lv0.external_blocks_r();
                        std::vector<int> b0_c(_base_N);
                        for (int i = 0; i < _base_N; ++i) b0_c[i] = eb0[i] - 1;
                        NestedWeightedSBM tmp = *this;
                        tmp.optimize_from_base(b0_c, rng, false, epsilon);
                        cand_S = tmp.entropy();
                    }
                    lv0.move(v, r);
                    if (!changed_B0) rebuild_upper_keep_bs();
                    if (cand_S < best_cand_S) {
                        best_cand_S = cand_S;
                        best_t = t;
                    }
                }
                if (best_t != r) {
                    int old_B0 = lv0.B();
                    lv0.move(v, best_t);
                    if (lv0.B() != old_B0) {
                        Rcpp::IntegerVector eb0 = lv0.external_blocks_r();
                        std::vector<int> b0_c(_base_N);
                        for (int i = 0; i < _base_N; ++i) b0_c[i] = eb0[i] - 1;
                        optimize_from_base(b0_c, rng, true, epsilon);
                    } else {
                        rebuild_upper_keep_bs();
                    }
                    double now_S = entropy();
                    if (now_S < best_S - 1e-8) {
                        best_S = now_S;
                        best_bs = save_bs();
                        improved = true;
                    } else {
                        build(best_bs);
                    }
                }
            }
            // 2. Try level-0 merges under full nested entropy
            while (static_cast<size_t>(_levels[0].B()) > B_min) {
                std::vector<int> occ = _levels[0].occupied();
                if (occ.size() > 30) break;
                int mr = -1;
                double mS = best_S - 1e-8;
                std::vector<std::vector<int>> m_bs;
                for (size_t i = 0; i < occ.size(); ++i) {
                    for (size_t j = i + 1; j < occ.size(); ++j) {
                        std::vector<int> b0_cand = _levels[0].blocks();
                        for (int& g : b0_cand) if (g == occ[i]) g = occ[j];
                        // Compress labels
                        std::vector<int> uq = b0_cand;
                        std::sort(uq.begin(), uq.end());
                        uq.erase(std::unique(uq.begin(), uq.end()), uq.end());
                        wt_hash_map<int, int> mp;
                        for (size_t q = 0; q < uq.size(); ++q) mp.emplace(uq[q], static_cast<int>(q));
                        for (int& g : b0_cand) g = mp[g];
                        NestedWeightedSBM tmp = *this;
                        tmp.optimize_from_base(b0_cand, rng, true, epsilon);
                        double cand_S = tmp.entropy();
                        if (cand_S < mS) {
                            mS = cand_S;
                            mr = occ[i];
                            m_bs = tmp.save_bs();
                        }
                    }
                }
                if (mr < 0) break;
                build(m_bs);
                best_S = mS;
                best_bs = m_bs;
                improved = true;
            }
            if (!improved) break;
            Rcpp::checkUserInterrupt();
        }
        build(best_bs);
    }

    std::vector<std::vector<int>> save_bs() const {
        std::vector<std::vector<int>> out;
        out.reserve(_levels.size());
        for (const auto& lv : _levels) {
            Rcpp::IntegerVector eb = lv.external_blocks_r();
            std::vector<int> v(eb.size());
            for (R_xlen_t i = 0; i < eb.size(); ++i) v[i] = eb[i] - 1;
            out.push_back(std::move(v));
        }
        return out;
    }

    void rebuild_upper_keep_bs() {
        for (int l = 1; l < levels(); ++l) {
            std::vector<int> bl = _levels[l].blocks();
            wt_hash_map<int, int> lmap;
            std::vector<WGraphEdge> uedges = build_upper_edges(_levels[l - 1], lmap);
            _levels[l].set_edges(std::move(uedges));
            _levels[l].set_partition(bl);
        }
    }

private:
    bool _directed, _dc, _dist;
    std::vector<CovSpec> _base_covs, _up_covs;
    int _base_N;
    std::vector<WGraphEdge> _base_edges;
    std::vector<WeightedLevel> _levels;
};

// ---------------------------------------------------------------------------
// Finite-temperature MCMC (single-node + merge-split) for flat WeightedLevel

void weighted_flat_mcmc_sweep(WeightedLevel& s, WRng& rng, double beta,
                              double c, double d, bool multiflip,
                              size_t B_min, size_t B_max,
                              double& total_dS, size_t& nattempts,
                              size_t& nmoves) {
    if (s.N() <= 1 || s.E() == 0) return;
    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    rng.shuffle(vs);

    // 1. Single-node Metropolis-Hastings sweep (exact detailed balance!)
    for (int v : vs) {
        ++nattempts;
        int r = s.group(v);
        size_t B_cur = static_cast<size_t>(s.B());
        bool can_create = (B_cur < static_cast<size_t>(s.N()) &&
                           B_cur < B_max && d > 0.0 && s.empty_group() >= 0);
        double p_new = can_create ? d : 0.0;
        bool propose_new = (can_create && rng.unif() < p_new);
        int t;
        double pfwd = 0.0;
        if (propose_new) {
            if (s.size(r) == 1) continue;
            t = s.empty_group();
            pfwd = p_new;
        } else {
            t = s.sample_group(v, c, rng);
            if (t == r) continue;
            if (s.size(r) == 1 && B_cur <= B_min) continue;
            pfwd = (1.0 - p_new) * s.move_prob(v, t, c);
        }
        bool empties = (s.size(r) == 1);
        double dS = s.virtual_move(v, t);
        if (std::isinf(beta)) {
            if (dS < 0.0) {
                s.apply_last_virtual_move(v, t);
                total_dS += dS;
                ++nmoves;
            }
            continue;
        }
        if (pfwd <= 0.0) continue;
        double u_rand = rng.unif();
        double log_u = (u_rand > 0.0) ? std::log(u_rand) : -W_INF;
        // Since prev <= 1.0, log_a <= -beta * dS - log(pfwd): exact early reject without mutating state!
        if (log_u >= -beta * dS - std::log(pfwd)) continue;
        s.apply_last_virtual_move(v, t);
        double prev = 0.0;
        if (empties) {
            prev = (d > 0.0) ? d : 0.0;
        } else {
            size_t B_after = static_cast<size_t>(s.B());
            bool rev_can = (B_after < static_cast<size_t>(s.N()) &&
                            B_after < B_max && d > 0.0);
            double p_new_rev = rev_can ? d : 0.0;
            prev = (1.0 - p_new_rev) * s.move_prob(v, r, c);
        }
        bool accept = (prev > 0.0 && log_u < -beta * dS + std::log(prev) - std::log(pfwd));
        if (accept) {
            total_dS += dS;
            ++nmoves;
        } else {
            s.move(v, r);
        }
    }

    // 2. Merge-split moves (auxiliary Gibbs 2-way split / merge, exact MH)
    if (multiflip && s.B() >= 1) {
        // A state-independent number of merge-split attempts per sweep (as in
        // blockmodel.cpp: psplit + pmerge = 2). Using s.B() made the number of
        // pi-invariant steps depend on the current state, which breaks
        // invariance (states with more groups got more merge attempts).
        const size_t n_ms = 2;
        for (size_t step = 0; step < n_ms; ++step) {
            bool try_split = (rng.unif() < 0.5);
            if (try_split) {
                if (static_cast<size_t>(s.B()) >= B_max || s.empty_group() < 0) continue;
                int r = s.random_occupied(rng);
                if (s.size(r) < 2) continue;
                ++nattempts;
                int t = s.empty_group();
                std::vector<int> mem = s.members(r);
                size_t nr = mem.size();
                size_t B_before = static_cast<size_t>(s.B());
                // Uniform random non-trivial 2-way split (excluding empty subsets)
                double dS = 0.0;
                size_t cnt_t = 0;
                for (int v : mem) {
                    if (rng.unif() < 0.5) {
                        dS += s.move(v, t);
                        ++cnt_t;
                    }
                }
                if (cnt_t == 0 || cnt_t == nr) {
                    for (int v : mem) if (s.group(v) != r) s.move(v, r);
                    continue;
                }
                // Exact proposal probabilities of this scheme (the 1/2 for
                // choosing split vs merge cancels):
                //  - split: pick r with 1/B, send each member to t with 1/2;
                //    trivial outcomes are rejected (self-transitions), and an
                //    unordered split {S, S^c} arises from 2 assignments, so
                //    q_split = (1/B) 2^(1 - n_r);
                //  - reverse merge: r and t drawn independently from B + 1
                //    groups (r == t rejected), so an unordered pair has
                //    q_merge = 2 / (B + 1)^2.
                const double LN2 = std::log(2.0);
                double log_pfwd = -std::log(static_cast<double>(B_before)) +
                                  (1.0 - static_cast<double>(nr)) * LN2;
                double log_prev = LN2 - 2.0 * std::log(static_cast<double>(B_before + 1));
                bool accept = false;
                if (std::isinf(beta)) {
                    accept = (dS < 0.0);
                } else {
                    double log_a = -beta * dS + log_prev - log_pfwd;
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) {
                    total_dS += dS;
                    ++nmoves;
                } else {
                    for (int v : mem) if (s.group(v) != r) s.move(v, r);
                }
            } else {
                if (static_cast<size_t>(s.B()) <= B_min || s.B() < 2) continue;
                ++nattempts;
                size_t B_before = static_cast<size_t>(s.B());
                int r = s.random_occupied(rng);
                int t = s.random_occupied(rng);
                if (r == t) continue;
                std::vector<int> r_mem = s.members(r);
                size_t nr_tot = r_mem.size() + static_cast<size_t>(s.size(t));
                double dS = s.merge(r, t);
                // Merge of an unordered pair drawn as two independent picks:
                // 2 / B^2; reverse split of the merged group: (1/(B-1)) 2^(1-n).
                const double LN2 = std::log(2.0);
                double log_pfwd = LN2 - 2.0 * std::log(static_cast<double>(B_before));
                double log_prev = -std::log(static_cast<double>(B_before - 1)) +
                                  (1.0 - static_cast<double>(nr_tot)) * LN2;
                bool accept = false;
                if (std::isinf(beta)) {
                    accept = (dS < 0.0);
                } else {
                    double log_a = -beta * dS + log_prev - log_pfwd;
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) {
                    total_dS += dS;
                    ++nmoves;
                } else {
                    for (int v : r_mem) s.move(v, r);
                }
            }
        }
    }
}

size_t w_q_cache_max(double x) {
    return (!std::isfinite(x) || x <= 0) ? 10000 : std::min<size_t>(static_cast<size_t>(x), 10000);
}

}  // namespace

// ---------------------------------------------------------------------------
// Rcpp Exports

// [[Rcpp::export]]
double cpp_weighted_entropy(SEXP pointer, Rcpp::IntegerVector blocks,
                            Rcpp::NumericMatrix weights,
                            Rcpp::CharacterVector types, Rcpp::List params,
                            bool degree_corrected, bool distributed,
                            bool include_rec, bool include_dl, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    WeightedLevel lv(N, graph.directed, degree_corrected, distributed, true, true, covs);
    lv.set_edges(build_base_edges(graph, weights));
    lv.set_partition(parse_w_blocks(blocks, N));
    return lv.entropy(include_rec, include_dl);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_weighted_check_moves(SEXP pointer,
                                             Rcpp::IntegerVector blocks,
                                             Rcpp::NumericMatrix weights,
                                             Rcpp::CharacterVector types,
                                             Rcpp::List params,
                                             bool degree_corrected,
                                             bool distributed, int nmoves,
                                             double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    WeightedLevel lv(N, graph.directed, degree_corrected, distributed, true, true, covs);
    lv.set_edges(build_base_edges(graph, weights));
    lv.set_partition(parse_w_blocks(blocks, N));
    WRng rng;
    Rcpp::NumericVector err(nmoves);
    if (N <= 1) return err;
    double S = lv.entropy();
    for (int i = 0; i < nmoves; ++i) {
        int v = static_cast<int>(rng.index(N));
        int s = (rng.unif() < 0.2 && lv.empty_group() >= 0)
                    ? lv.empty_group()
                    : lv.random_occupied(rng);
        double dS_v = lv.virtual_move(v, s);
        double dS_a = lv.move(v, s);
        double S_new = lv.entropy();
        double max_e = std::max(std::abs(dS_v - (S_new - S)), std::abs(dS_a - (S_new - S)));
        if (lv.B() >= 2 && (i % 5 == 0)) {
            int mr = lv.random_occupied(rng);
            int mt = lv.random_occupied(rng);
            if (mr != mt) {
                WeightedLevel::MoveBuf mb;
                double mc = lv.virtual_merge_const(mr, mt, mb);
                double ma = lv.virtual_merge(mr, mt);
                max_e = std::max(max_e, std::abs(mc - ma));
            }
        }
        err[i] = max_e;
        S = S_new;
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_weighted_minimize(SEXP pointer, Rcpp::NumericMatrix weights,
                                 Rcpp::CharacterVector types, Rcpp::List params,
                                 bool degree_corrected, bool distributed,
                                 int B_min, double B_max, bool refine,
                                 double epsilon, double q_max, int nthreads = 0) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    WeightedLevel lv(N, graph.directed, degree_corrected, distributed, true, true, covs);
    lv.set_edges(build_base_edges(graph, weights));
    WRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    double S = weighted_minimize_flat(lv, rng, static_cast<size_t>(B_min), bmax,
                                      refine, epsilon, nthreads);
    return Rcpp::List::create(
        Rcpp::Named("blocks") = lv.external_blocks_r(),
        Rcpp::Named("entropy") = S,
        Rcpp::Named("B") = lv.B());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_weighted_entropy(SEXP pointer, Rcpp::List bs,
                                                Rcpp::NumericMatrix weights,
                                                Rcpp::CharacterVector types,
                                                Rcpp::List params,
                                                bool degree_corrected,
                                                bool distributed, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    NestedWeightedSBM nw(graph, weights, covs, degree_corrected, distributed);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_w_blocks(x, cur_n);
        int max_lbl = 0;
        for (int v : bl) if (v + 1 > max_lbl) max_lbl = v + 1;
        cur_n = max_lbl;
        b_vec.push_back(std::move(bl));
    }
    nw.build(b_vec);
    auto lvS = nw.level_entropies();
    return Rcpp::NumericVector(lvS.begin(), lvS.end());
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_weighted_minimize(SEXP pointer, Rcpp::NumericMatrix weights,
                                        Rcpp::CharacterVector types,
                                        Rcpp::List params, bool degree_corrected,
                                        bool distributed, int B_min, double B_max,
                                        bool refine, double epsilon, double q_max,
                                        int nthreads = 0) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    WRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);

    // Seed level 0 via fast flat weighted minimization, build optimal hierarchy,
    // and jointly refine under exact total nested weighted DL!
    WeightedLevel seed(N, graph.directed, degree_corrected, distributed, true, true, covs);
    seed.set_edges(build_base_edges(graph, weights));
    weighted_minimize_flat(seed, rng, static_cast<size_t>(B_min), bmax, refine, epsilon, nthreads);
    Rcpp::IntegerVector eb0 = seed.external_blocks_r();
    std::vector<int> b0(N);
    for (int i = 0; i < N; ++i) b0[i] = eb0[i] - 1;

    NestedWeightedSBM nw(graph, weights, covs, degree_corrected, distributed);
    nw.optimize_from_base(b0, rng, refine, epsilon, nthreads);
    if (refine)
        nw.refine_joint(rng, static_cast<size_t>(B_min), bmax, epsilon);

    auto lvS = nw.level_entropies();
    return Rcpp::List::create(
        Rcpp::Named("bs") = nw.external_bs_r(),
        Rcpp::Named("level_entropy") = lvS,
        Rcpp::Named("entropy") = std::accumulate(lvS.begin(), lvS.end(), 0.0));
}

// [[Rcpp::export]]
Rcpp::List cpp_weighted_mcmc(SEXP pointer, Rcpp::IntegerVector blocks,
                             Rcpp::NumericMatrix weights,
                             Rcpp::CharacterVector types, Rcpp::List params,
                             bool degree_corrected, bool distributed,
                             Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                             double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    WeightedLevel lv(N, graph.directed, degree_corrected, distributed, true, true, covs);
    lv.set_edges(build_base_edges(graph, weights));
    lv.set_partition(parse_w_blocks(blocks, N));
    WRng rng;

    double beta = Rcpp::as<double>(mcmc_opts["beta"]);
    double c = Rcpp::as<double>(mcmc_opts["c"]);
    double d = Rcpp::as<double>(mcmc_opts["d"]);
    int niter = static_cast<int>(Rcpp::as<double>(mcmc_opts["niter"]));
    bool multiflip = Rcpp::as<bool>(mcmc_opts["multiflip"]);
    size_t B_min = static_cast<size_t>(Rcpp::as<double>(mcmc_opts["B_min"]));
    double bmax_d = Rcpp::as<double>(mcmc_opts["B_max"]);
    size_t B_max = std::isinf(bmax_d) ? std::numeric_limits<size_t>::max()
                                      : static_cast<size_t>(bmax_d);

    std::string mode = Rcpp::as<std::string>(eq_opts["mode"]);
    double total_dS = 0.0;
    size_t nattempts = 0, nmoves = 0;
    int steps_done = 0;

    if (mode == "sweep") {
        for (int it = 0; it < niter; ++it) {
            weighted_flat_mcmc_sweep(lv, rng, beta, c, d, multiflip, B_min, B_max,
                                     total_dS, nattempts, nmoves);
            if ((it + 1) % 25 == 0) lv.refresh_cov_globals();
            ++steps_done;
            Rcpp::checkUserInterrupt();
        }
        return Rcpp::List::create(
            Rcpp::Named("blocks") = lv.external_blocks_r(),
            Rcpp::Named("entropy") = lv.entropy(),
            Rcpp::Named("dS") = total_dS,
            Rcpp::Named("nattempts") = static_cast<double>(nattempts),
            Rcpp::Named("nmoves") = static_cast<double>(nmoves),
            Rcpp::Named("niter") = steps_done);
    }

    // Equilibrate / collect / anneal
    int wait = static_cast<int>(Rcpp::as<double>(eq_opts["wait"]));
    int nbreaks = static_cast<int>(Rcpp::as<double>(eq_opts["nbreaks"]));
    double max_niter_d = Rcpp::as<double>(eq_opts["max_niter"]);
    int max_niter = std::isinf(max_niter_d) ? std::numeric_limits<int>::max()
                                            : static_cast<int>(max_niter_d);
    int force_niter = static_cast<int>(Rcpp::as<double>(eq_opts["force_niter"]));
    double epsilon = Rcpp::as<double>(eq_opts["epsilon"]);
    bool keep_hist = Rcpp::as<bool>(eq_opts["history"]);
    bool collect = Rcpp::as<bool>(eq_opts["collect"]);
    int thin = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["thin"])));

    std::vector<double> hist_S;
    std::vector<Rcpp::IntegerVector> samples;
    std::vector<double> sample_entropies;

    auto run_eq_loop = [&](double cur_beta, int max_it, int f_it) {
        double best_S = lv.entropy();
        double worst_S = best_S;
        int no_break = 0;
        int breaks_left = nbreaks;
        for (int step = 1; step <= max_it; ++step) {
            for (int sub = 0; sub < std::max(1, niter); ++sub) {
                weighted_flat_mcmc_sweep(lv, rng, cur_beta, c, d, multiflip,
                                         B_min, B_max, total_dS, nattempts, nmoves);
            }
            if (step % 25 == 0) lv.refresh_cov_globals();
            double S_now = lv.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_now);
            if (collect && (step % thin == 0)) {
                samples.push_back(lv.external_blocks_r());
                sample_entropies.push_back(S_now);
            }
            if (f_it > 0) {
                if (step >= f_it) break;
            } else {
                bool broke = false;
                if (S_now < best_S - epsilon) { best_S = S_now; broke = true; }
                if (S_now > worst_S + epsilon) { worst_S = S_now; broke = true; }
                if (broke) {
                    no_break = 0;
                } else if (++no_break >= wait) {
                    if (--breaks_left <= 0) break;
                    best_S = worst_S = S_now;
                    no_break = 0;
                }
            }
            Rcpp::checkUserInterrupt();
        }
    };

    if (mode == "anneal") {
        double b0 = Rcpp::as<double>(eq_opts["beta_start"]);
        double b1 = Rcpp::as<double>(eq_opts["beta_end"]);
        int an = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["anneal_niter"])));
        for (int k = 0; k < an; ++k) {
            double frac = (an == 1) ? 1.0 : static_cast<double>(k) / static_cast<double>(an - 1);
            double cur_beta = b0 * std::pow(b1 / b0, frac);
            run_eq_loop(cur_beta, max_niter, force_niter);
        }
    } else {
        int lim = (force_niter > 0) ? force_niter : max_niter;
        run_eq_loop(beta, lim, force_niter);
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("blocks") = lv.external_blocks_r(),
        Rcpp::Named("entropy") = lv.entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("niter") = steps_done);
    if (keep_hist) out["history"] = hist_S;
    if (collect) {
        Rcpp::List sl(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) sl[i] = samples[i];
        out["samples"] = sl;
        out["sample_entropies"] = sample_entropies;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_weighted_mcmc(SEXP pointer, Rcpp::List bs,
                                    Rcpp::NumericMatrix weights,
                                    Rcpp::CharacterVector types,
                                    Rcpp::List params, bool degree_corrected,
                                    bool distributed, Rcpp::List mcmc_opts,
                                    Rcpp::List eq_opts, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) wlog_q().init(w_q_cache_max(q_max));
    auto covs = parse_covs(weights, types, params);
    NestedWeightedSBM nw(graph, weights, covs, degree_corrected, distributed);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_w_blocks(x, cur_n);
        int max_lbl = 0;
        for (int v : bl) if (v + 1 > max_lbl) max_lbl = v + 1;
        cur_n = max_lbl;
        b_vec.push_back(std::move(bl));
    }
    nw.build(b_vec);
    WRng rng;

    double beta = Rcpp::as<double>(mcmc_opts["beta"]);
    double c = Rcpp::as<double>(mcmc_opts["c"]);
    int niter = static_cast<int>(Rcpp::as<double>(mcmc_opts["niter"]));
    size_t B_min = static_cast<size_t>(Rcpp::as<double>(mcmc_opts["B_min"]));
    std::string mode = Rcpp::as<std::string>(eq_opts["mode"]);

    double total_dS = 0.0;
    size_t nattempts = 0, nmoves = 0;
    int steps_done = 0;
    double S_cur = nw.entropy();

    // Probability of a B-changing proposal per node visit.
    double d_bc = mcmc_opts.containsElementNamed("d") ? Rcpp::as<double>(mcmc_opts["d"]) : 0.01;
    if (!(d_bc > 0.0)) d_bc = 0.0;

    // B-changing moves at every level, edited on the contiguous level
    // vectors (level l + 1 node i is level-l group i) and evaluated by
    // rebuilding the hierarchy:
    //  - split-off: v (in a group r of size >= 2, parent P) moves to a new
    //    group under P;
    //  - absorb: a singleton v (group t, parent P) joins one of the K groups
    //    s != t with parent P, chosen uniformly; t disappears.
    // These are exact reverses of each other (P always keeps a member, so
    // nothing cascades upward): q_rev / q_fwd = 1/K for a split-off, K for an
    // absorb.
    // At the top level the groups share one implicit parent.
    auto bchange_move = [&](int l, int v, double cur_beta) {
        std::vector<std::vector<int>> bs_old = nw.save_bs();
        std::vector<std::vector<int>> bs_new = bs_old;
        std::vector<int>& bl = bs_new[l];
        const bool has_up = (l + 1 < nw.levels());
        std::vector<int> no_up;
        std::vector<int>& bu = has_up ? bs_new[l + 1] : no_up;
        int r = bl[v];
        int B_l = has_up ? static_cast<int>(bu.size())
                         : (*std::max_element(bl.begin(), bl.end()) + 1);
        auto parent = [&](int g) { return has_up ? bu[g] : 0; };
        int nr = 0;
        for (int x : bl) if (x == r) ++nr;
        int P = parent(r);
        double log_q_ratio;
        if (nr >= 2) {
            bl[v] = B_l;
            if (has_up) bu.push_back(P);
            int K = 0;  // groups with parent P other than the new one (includes r)
            for (int g = 0; g < B_l; ++g) if (parent(g) == P) ++K;
            log_q_ratio = -std::log(static_cast<double>(K));
        } else {
            if (l == 0 && static_cast<size_t>(B_l) <= B_min) return;
            std::vector<int> cand;
            for (int s = 0; s < B_l; ++s) if (s != r && parent(s) == P) cand.push_back(s);
            if (cand.empty()) return;
            int K = static_cast<int>(cand.size());
            int tgt = cand[rng.index(cand.size())];
            bl[v] = tgt;
            for (int& x : bl) if (x > r) --x;
            if (has_up) bu.erase(bu.begin() + r);
            log_q_ratio = std::log(static_cast<double>(K));
        }
        ++nattempts;
        double S_old = S_cur;
        nw.build(bs_new);
        double S_new = nw.entropy();
        double dS = S_new - S_old;
        bool accept;
        if (std::isinf(cur_beta)) {
            accept = (dS < 0.0);
        } else {
            double log_a = -cur_beta * dS + log_q_ratio;
            accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
        }
        if (accept) {
            total_dS += dS;
            S_cur = S_new;
            ++nmoves;
        } else {
            nw.build(bs_old);
        }
    };

    auto one_nested_sweep = [&](double cur_beta) {
        // Level-by-level MH sweeps evaluated under exact total nested weighted entropy.
        for (int l = 0; l < nw.levels(); ++l) {
            if (nw.level(l).N() <= 1) continue;
            const bool can_bc = d_bc > 0.0;
            std::vector<int> ord(nw.level(l).N());
            std::iota(ord.begin(), ord.end(), 0);
            rng.shuffle(ord);
            for (int v : ord) {
                if (can_bc && rng.unif() < d_bc) {
                    bchange_move(l, v, cur_beta);
                    continue;
                }
                WeightedLevel& lv = nw.level(l);  // re-fetch: build() reallocates levels
                if (lv.B() <= 1) continue;
                int r = lv.group(v);
                if (lv.size(r) <= 1) continue;  // B-preserving move; B changes above
                int t = lv.sample_group(v, c, rng);
                if (t == r) continue;
                ++nattempts;
                double pfwd = lv.move_prob(v, t, c);
                lv.move(v, t);
                if (l + 1 < nw.levels()) nw.rebuild_upper_keep_bs();
                double S_new = nw.entropy();
                double dS = S_new - S_cur;
                bool accept = false;
                if (std::isinf(cur_beta)) {
                    accept = (dS < 0.0);
                } else {
                    double prev = nw.level(l).move_prob(v, r, c);
                    if (prev > 0.0 && pfwd > 0.0) {
                        double log_a = -cur_beta * dS + std::log(prev) - std::log(pfwd);
                        accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                    }
                }
                if (accept) {
                    total_dS += dS;
                    S_cur = S_new;
                    ++nmoves;
                } else {
                    nw.level(l).move(v, r);
                    if (l + 1 < nw.levels()) nw.rebuild_upper_keep_bs();
                }
            }
        }
    };

    if (mode == "sweep") {
        for (int it = 0; it < niter; ++it) {
            one_nested_sweep(beta);
            ++steps_done;
            Rcpp::checkUserInterrupt();
        }
        return Rcpp::List::create(
            Rcpp::Named("bs") = nw.external_bs_r(),
            Rcpp::Named("entropy") = nw.entropy(),
            Rcpp::Named("dS") = total_dS,
            Rcpp::Named("nattempts") = static_cast<double>(nattempts),
            Rcpp::Named("nmoves") = static_cast<double>(nmoves),
            Rcpp::Named("niter") = steps_done);
    }

    int wait = static_cast<int>(Rcpp::as<double>(eq_opts["wait"]));
    int nbreaks = static_cast<int>(Rcpp::as<double>(eq_opts["nbreaks"]));
    double max_niter_d = Rcpp::as<double>(eq_opts["max_niter"]);
    int max_niter = std::isinf(max_niter_d) ? std::numeric_limits<int>::max()
                                            : static_cast<int>(max_niter_d);
    int force_niter = static_cast<int>(Rcpp::as<double>(eq_opts["force_niter"]));
    double epsilon = Rcpp::as<double>(eq_opts["epsilon"]);
    bool keep_hist = Rcpp::as<bool>(eq_opts["history"]);
    bool collect = Rcpp::as<bool>(eq_opts["collect"]);
    int thin = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["thin"])));

    std::vector<double> hist_S;
    std::vector<Rcpp::List> samples;
    std::vector<double> sample_entropies;

    auto run_eq = [&](double cur_beta, int max_it, int f_it) {
        double best_S = S_cur, worst_S = S_cur;
        int no_break = 0, breaks_left = nbreaks;
        for (int step = 1; step <= max_it; ++step) {
            for (int sub = 0; sub < std::max(1, niter); ++sub)
                one_nested_sweep(cur_beta);
            S_cur = nw.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_cur);
            if (collect && (step % thin == 0)) {
                samples.push_back(nw.external_bs_r());
                sample_entropies.push_back(S_cur);
            }
            if (f_it > 0) {
                if (step >= f_it) break;
            } else {
                bool broke = false;
                if (S_cur < best_S - epsilon) { best_S = S_cur; broke = true; }
                if (S_cur > worst_S + epsilon) { worst_S = S_cur; broke = true; }
                if (broke) no_break = 0;
                else if (++no_break >= wait) {
                    if (--breaks_left <= 0) break;
                    best_S = worst_S = S_cur;
                    no_break = 0;
                }
            }
            Rcpp::checkUserInterrupt();
        }
    };

    if (mode == "anneal") {
        double b0 = Rcpp::as<double>(eq_opts["beta_start"]);
        double b1 = Rcpp::as<double>(eq_opts["beta_end"]);
        int an = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["anneal_niter"])));
        for (int k = 0; k < an; ++k) {
            double frac = (an == 1) ? 1.0 : static_cast<double>(k) / static_cast<double>(an - 1);
            run_eq(b0 * std::pow(b1 / b0, frac), max_niter, force_niter);
        }
    } else {
        int lim = (force_niter > 0) ? force_niter : max_niter;
        run_eq(beta, lim, force_niter);
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("bs") = nw.external_bs_r(),
        Rcpp::Named("entropy") = nw.entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("niter") = steps_done);
    if (keep_hist) out["history"] = hist_S;
    if (collect) {
        Rcpp::List sl(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) sl[i] = samples[i];
        out["samples"] = sl;
        out["sample_entropies"] = sample_entropies;
    }
    return out;
}
