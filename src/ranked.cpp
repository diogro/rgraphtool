// Native C++ flat and nested Ranked (ordered) stochastic block model engine
// matching graph-tool 3.9's RankedBlockState and
// NestedBlockState(base_state=RankedBlockState) to machine precision.
//
// References:
//   Peixoto, T. P. (2022). Ordered community detection in directed networks.
//   Physical Review E 106, 024305.
//
// Adapted from graph-tool (https://graph-tool.skewed.de), LGPL-3.0-or-later.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph_types.h"
#include "logq.h"

namespace {

constexpr double RK_INF = std::numeric_limits<double>::infinity();
constexpr double RK_LOG2 = 0.693147180559945309417232121458176568;

template <class K, class V>
using rk_hash_map = std::unordered_map<K, V>;
template <class K>
using rk_hash_set = std::unordered_set<K>;

struct RKMath {
    std::vector<double> lgam;
    std::vector<double> logi;

    void ensure(size_t n) {
        if (n < lgam.size()) return;
        size_t old = lgam.size();
        size_t sz = std::max<size_t>(n + 64, old * 2 + 64);
        lgam.resize(sz);
        logi.resize(sz);
        for (size_t i = old; i < sz; ++i) {
            lgam[i] = std::lgamma(static_cast<double>(i) + 1.0);
            logi[i] = (i > 0) ? std::log(static_cast<double>(i)) : 0.0;
        }
    }

    inline double lgam1(long n) const {
        if (n <= 1) return 0.0;
        if (static_cast<size_t>(n) < lgam.size()) return lgam[n];
        return std::lgamma(static_cast<double>(n) + 1.0);
    }

    inline double logn(long n) const {
        if (n <= 1) return 0.0;
        if (static_cast<size_t>(n) < logi.size()) return logi[n];
        return std::log(static_cast<double>(n));
    }

    inline double lbinom(int64_t n, int64_t k) const {
        if (k < 0 || n < 0 || k > n) return 0.0;
        if (k == 0 || k == n) return 0.0;
        if (static_cast<size_t>(n) < lgam.size())
            return lgam[n] - lgam[k] - lgam[n - k];
        return std::lgamma(static_cast<double>(n) + 1.0) -
               std::lgamma(static_cast<double>(k) + 1.0) -
               std::lgamma(static_cast<double>(n - k) + 1.0);
    }
};

RKMath& rkmath() {
    static RKMath m;
    return m;
}
// Shared log q(n, k) table (src/logq.h).
using RKLogQ = rgt_logq::LogQ;
inline RKLogQ& rklog_q() { return rgt_logq::shared(); }

size_t rk_q_cache_max(double q_max) {
    if (!std::isfinite(q_max) || q_max <= 0.0) return 10000;
    return std::min<size_t>(static_cast<size_t>(q_max), 10000);
}

struct RKRng {
    inline double unif() { return R::runif(0.0, 1.0); }
    inline size_t index(size_t n) {
        size_t i = static_cast<size_t>(std::floor(unif() * static_cast<double>(n)));
        return (i >= n) ? (n - 1) : i;
    }
    template <class Vec>
    void shuffle(Vec& v) {
        for (size_t i = v.size(); i > 1; --i) std::swap(v[i - 1], v[index(i)]);
    }
};

inline uint64_t rk_pack(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(b));
}

inline uint64_t rk_undir_pack(int a, int b) {
    if (a > b) std::swap(a, b);
    return rk_pack(a, b);
}

inline std::pair<int, int> rk_unpack(uint64_t key) {
    return {static_cast<int>(static_cast<uint32_t>(key >> 32)),
            static_cast<int>(static_cast<uint32_t>(key & 0xffffffffULL))};
}

// Exact graph-tool 3.9 distance-model helper functions (spec_base.hh lines 309-379)
inline long rk_floor_half(long m) { return m / 2; }
inline long rk_ceil_half(long m) { return (m + 1) / 2; }

inline long rk_edp(long ers, long mrs) {
    if (mrs <= 0) return 0;
    long fl = rk_floor_half(mrs);
    long ce = rk_ceil_half(mrs);
    long val = 0;
    if (ers <= fl) val += ers;
    if (ers >= ce) val += (mrs - ers);
    if (2 * ers == mrs) val -= ers;
    return val;
}

inline long rk_edm(long ers, long mrs) {
    if (mrs <= 0) return 0;
    long fl = rk_floor_half(mrs);
    long ce = rk_ceil_half(mrs);
    long val = 0;
    if (ers <= fl) val += (fl - ers);
    if (ers >= ce) val += (ers - ce);
    return val;
}

inline double rk_Sd(long edp, long edm) {
    const RKMath& M = rkmath();
    return M.lgam1(edp + edm + 1) - M.lgam1(edp) - M.lgam1(edm);
}

inline double rk_pair_orient_term(long ers, long mrs, bool distance, bool same_rank) {
    if (mrs <= 0) return 0.0;
    const RKMath& M = rkmath();
    if (!distance || same_rank) {
        return -M.lbinom(mrs, ers);
    }
    long fl = rk_floor_half(mrs);
    long ce = rk_ceil_half(mrs);
    double s = 0.0;
    if (ers < fl) s -= M.lbinom(fl, ers);
    if (ers > ce) s -= M.lbinom(mrs - ce, ers - ce);
    return s;
}

// ---------------------------------------------------------------------------
// Base-level directed Ranked SBM state (matches graph-tool 3.9 RankedBlockState)

class RankedState {
public:
    RankedState(const Graph& graph, bool distance, bool unique, long M_pos,
                bool dc, bool dist, bool is_top)
        : _N(static_cast<int>(boost::num_vertices(graph.adjacency))),
          _distance(distance), _unique(unique),
          _M(M_pos > 0 ? M_pos : std::max<long>(1, _N)),
          _dc(dc), _dist(dist), _is_top(is_top),
          _cap(std::max(2, _N)),
          _adj_out(_N), _adj_in(_N), _kout(_N, 0), _kin(_N, 0),
          _b(_N, 0), _u(_cap, 0), _wr(_cap, 0), _ep(_cap, 0), _em(_cap, 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(_N, 0),
          _m_out(_cap), _m_in(_cap),
          _h_out(dc && dist ? _cap : 0),
          _h_in(dc && dist ? _cap : 0) {
        if (!graph.directed)
            Rcpp::stop("RankedBlockState can only be used with directed graphs");
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        for (int r = 0; r < _cap; ++r) _u[r] = r;

        _E = 0;
        rk_hash_map<uint64_t, long> emult;
        for (auto [it, end] = boost::edges(graph.adjacency); it != end; ++it) {
            int u = static_cast<int>(boost::source(*it, graph.adjacency));
            int v = static_cast<int>(boost::target(*it, graph.adjacency));
            _edges.push_back({u, v});
            _adj_out[u].push_back(v);
            _adj_in[v].push_back(u);
            _kout[u]++;
            _kin[v]++;
            ++_E;
            emult[rk_pack(u, v)]++;
        }

        rkmath().ensure(static_cast<size_t>(_N + _M + 2 * _E + 64));
        if (_dc && _dist)
            rklog_q().init(rgt_logq::table_size(static_cast<double>(_E), 5000));

        const RKMath& M = rkmath();
        _base_const = 0.0;
        for (const auto& [pk, mult] : emult) _base_const += M.lgam1(mult);
        if (_dc) {
            for (int v = 0; v < _N; ++v)
                _base_const -= M.lgam1(_kout[v]) + M.lgam1(_kin[v]);
        }
    }

    int N() const { return _N; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    bool distance() const { return _distance; }
    bool unique() const { return _unique; }
    long M() const { return _M; }
    bool is_top() const { return _is_top; }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    long rank_of(int r) const { return _u[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<long>& ranks() const { return _u; }
    const std::array<long, 4>& Es() const { return _Es; }
    const rk_hash_map<uint64_t, long>& dir_pairs() const { return _dir_pairs; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }
    int random_occupied(RKRng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    long get_mrs(int r, int s) const {
        auto it = _dir_pairs.find(rk_pack(r, s));
        return (it == _dir_pairs.end()) ? 0 : it->second;
    }

    void set_partition(const std::vector<int>& b, const std::vector<long>* u_in = nullptr) {
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _em[r] = 0;
            _members[r].clear();
            _occ_pos[r] = -1;
            _m_out[r].clear();
            _m_in[r].clear();
            if (_dc && _dist) {
                _h_out[r].clear();
                _h_in[r].clear();
            }
        }
        _occupied.clear();
        _empty.clear();
        _empty_pos.assign(_cap, -1);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        _dir_pairs.clear();
        _Es = {0, 0, 0, 0};
        _ed.clear();
        _Dp = 0;
        _Dm = 0;

        for (int r = 0; r < _cap; ++r) _u[r] = r;
        if (u_in != nullptr) {
            for (size_t i = 0; i < u_in->size() && i < static_cast<size_t>(_cap); ++i)
                _u[i] = (*u_in)[i];
        }
        // Exact graph-tool 3.9 rule (spec_base.hh line 206):
        // when distance=TRUE and unique=TRUE, group ranks are initialized to u[r] = r
        if (_distance && _unique) {
            for (int r = 0; r < _cap; ++r) _u[r] = r;
        }

        for (int v = 0; v < _N; ++v) {
            int r = b[v];
            _b[v] = r;
            add_member(v, r);
            _ep[r] += _kout[v];
            _em[r] += _kin[v];
            if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                _h_in[r][_kin[v]]++;
            }
        }

        for (const auto& e : _edges) {
            int r = _b[e.first], s = _b[e.second];
            _dir_pairs[rk_pack(r, s)]++;
            _m_out[r][s]++;
            _m_in[s][r]++;
            _Es[stream_dir(r, s)]++;
            if (r == s) _Es[3]++;
        }

        if (_distance) rebuild_distance_stats();
    }

    // Recompute distance-based stats (_ed, _Dp, _Dm) from _dir_pairs
    void rebuild_distance_stats() {
        _ed.clear();
        _Dp = 0;
        _Dm = 0;
        rk_hash_set<uint64_t> seen;
        for (const auto& [pk, ers_raw] : _dir_pairs) {
            if (ers_raw <= 0) continue;
            auto [r, s] = rk_unpack(pk);
            if (r == s) continue;
            uint64_t upk = rk_undir_pack(r, s);
            if (!seen.insert(upk).second) continue;
            long d = std::abs(_u[r] - _u[s]);
            if (d == 0) continue;
            // Canonical orientation: _u[a] < _u[b]
            int a = (_u[r] < _u[s]) ? r : s;
            int b_g = (_u[r] < _u[s]) ? s : r;
            long eab = get_mrs(a, b_g);
            long eba = get_mrs(b_g, a);
            long m = eab + eba;
            auto& entry = _ed[d];
            entry.first += rk_edp(eab, m);
            entry.second += rk_edm(eab, m);
            if (eab < rk_floor_half(m)) _Dp++;
            if (eab > rk_ceil_half(m)) _Dm++;
        }
    }

    // Exact graph-tool 3.9 RankedBlockState description length
    double entropy(bool include_label_dl = true) const {
        const RKMath& M = rkmath();
        long Bn = B();
        double S = _base_const;

        // 1. Base directed BlockState (with edges_dl = FALSE)
        if (_N > 0) {
            S += M.logn(_N) + M.lbinom(_N - 1, Bn - 1) + M.lgam1(_N);
            for (int r : _occupied) S -= M.lgam1(_wr[r]);
        }
        for (int r : _occupied) {
            if (_dc) {
                S += M.lgam1(_ep[r]) + M.lgam1(_em[r]) + deg_dl_group(r);
            } else {
                S += static_cast<double>(_ep[r] + _em[r]) * M.logn(_wr[r]);
            }
        }
        for (const auto& [pk, m] : _dir_pairs) {
            if (m > 0) S -= M.lgam1(m);
        }

        // 2. Ranked edges_dl replacement (always active!)
        if (_is_top) {
            int64_t E_tot = _Es[0] + _Es[1] + _Es[2];
            int64_t P_undir = int64_t(Bn) * (Bn + 1) / 2;
            if (E_tot > 0) S += M.lbinom(P_undir + E_tot - 1, E_tot);
        }
        if (_distance) {
            S += M.lbinom(_N, Bn);
        } else {
            S += rk_Sd(_Es[0], _Es[2]);
        }
        S += static_cast<double>(_Es[1] - _Es[3]) * RK_LOG2;

        // Pair orientation terms (iterate each unordered pair {r, s}, r != s once)
        rk_hash_set<uint64_t> seen;
        for (const auto& [pk, m_raw] : _dir_pairs) {
            if (m_raw <= 0) continue;
            auto [r, s] = rk_unpack(pk);
            if (r == s) continue;
            uint64_t upk = rk_undir_pack(r, s);
            if (!seen.insert(upk).second) continue;
            bool same = (_u[r] == _u[s]);
            int a = (!same && _u[r] > _u[s]) ? s : r;
            int b_g = (!same && _u[r] > _u[s]) ? r : s;
            long eab = get_mrs(a, b_g);
            long eba = get_mrs(b_g, a);
            S += rk_pair_orient_term(eab, eab + eba, _distance, same);
        }

        if (_distance) {
            for (const auto& [d, ep_em] : _ed)
                S += rk_Sd(ep_em.first, ep_em.second);
            S += rk_Sd(_Dp, _Dm);
        }

        // 3. Label DL (if !unique)
        if (!_unique && include_label_dl) {
            S += M.lbinom(_M + Bn - 1, Bn);
        }
        return S;
    }

    // Exact O(deg(v)) incremental virtual_move(v, s)
    double virtual_move(int v, int s) const {
        int r = _b[v];
        if (r == s) return 0.0;
        return compute_move_delta(v, r, s, &_scratch_buf);
    }

    double move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double dS = compute_move_delta(v, r, s, &_scratch_buf);
        apply_move(v, r, s, _scratch_buf);
        return dS;
    }

    void apply_last_virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s) return;
        apply_move(v, r, s, _scratch_buf);
    }

    double virtual_merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        for (int v : verts) move(v, r);
        return dS;
    }

    double merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        return dS;
    }

    // Exact O(deg_bg(r)) incremental rank relabeling: set _u[r] = new_u
    double relabel_group_dS(int r, long new_u) const {
        long old_u = _u[r];
        if (old_u == new_u || _wr[r] == 0) return 0.0;

        // Collect all distinct neighbor blocks t != r connected to r in either direction
        rk_hash_set<int> nbrs;
        for (const auto& [t, m] : _m_out[r]) if (t != r && m > 0) nbrs.insert(t);
        for (const auto& [t, m] : _m_in[r]) if (t != r && m > 0) nbrs.insert(t);
        if (nbrs.empty()) return 0.0;

        std::array<long, 4> dE = {0, 0, 0, 0};
        double dS = 0.0;
        rk_hash_map<long, std::pair<long, long>> ed_d;
        long dDp = 0, dDm = 0;

        for (int t : nbrs) {
            long ert = get_mrs(r, t);
            long etr = get_mrs(t, r);
            long m = ert + etr;
            long ut = _u[t];

            int sd_rt0 = udir(old_u, ut), sd_tr0 = udir(ut, old_u);
            int sd_rt1 = udir(new_u, ut), sd_tr1 = udir(ut, new_u);
            dE[sd_rt0] -= ert; dE[sd_tr0] -= etr;
            dE[sd_rt1] += ert; dE[sd_tr1] += etr;

            if (_distance) {
                long d0 = std::abs(old_u - ut);
                long d1 = std::abs(new_u - ut);
                // Old pair orientation contribution
                bool same0 = (d0 == 0);
                long eab0 = (!same0 && old_u > ut) ? etr : ert;
                dS -= rk_pair_orient_term(eab0, m, true, same0);
                if (!same0) {
                    ed_d[d0].first -= rk_edp(eab0, m);
                    ed_d[d0].second -= rk_edm(eab0, m);
                    if (eab0 < rk_floor_half(m)) dDp--;
                    if (eab0 > rk_ceil_half(m)) dDm--;
                }
                // New pair orientation contribution
                bool same1 = (d1 == 0);
                long eab1 = (!same1 && new_u > ut) ? etr : ert;
                dS += rk_pair_orient_term(eab1, m, true, same1);
                if (!same1) {
                    ed_d[d1].first += rk_edp(eab1, m);
                    ed_d[d1].second += rk_edm(eab1, m);
                    if (eab1 < rk_floor_half(m)) dDp++;
                    if (eab1 > rk_ceil_half(m)) dDm++;
                }
            }
        }

        if (!_distance) {
            dS += rk_Sd(_Es[0] + dE[0], _Es[2] + dE[2]) - rk_Sd(_Es[0], _Es[2]);
        } else {
            for (const auto& [d, dd] : ed_d) {
                auto it = _ed.find(d);
                long ep0 = (it == _ed.end()) ? 0 : it->second.first;
                long em0 = (it == _ed.end()) ? 0 : it->second.second;
                dS += rk_Sd(ep0 + dd.first, em0 + dd.second) - rk_Sd(ep0, em0);
            }
            dS += rk_Sd(_Dp + dDp, _Dm + dDm) - rk_Sd(_Dp, _Dm);
        }
        dS += static_cast<double>(dE[1]) * RK_LOG2;
        return dS;
    }

    double relabel_group(int r, long new_u) {
        long old_u = _u[r];
        if (old_u == new_u || _wr[r] == 0) return 0.0;
        double dS = relabel_group_dS(r, new_u);

        rk_hash_set<int> nbrs;
        for (const auto& [t, m] : _m_out[r]) if (t != r && m > 0) nbrs.insert(t);
        for (const auto& [t, m] : _m_in[r]) if (t != r && m > 0) nbrs.insert(t);
        for (int t : nbrs) {
            long ert = get_mrs(r, t);
            long etr = get_mrs(t, r);
            long m = ert + etr;
            long ut = _u[t];
            _Es[udir(old_u, ut)] -= ert;
            _Es[udir(ut, old_u)] -= etr;
            _Es[udir(new_u, ut)] += ert;
            _Es[udir(ut, new_u)] += etr;
            if (_distance) {
                long d0 = std::abs(old_u - ut);
                long d1 = std::abs(new_u - ut);
                if (d0 > 0) {
                    long eab0 = (old_u > ut) ? etr : ert;
                    auto& e0 = _ed[d0];
                    e0.first -= rk_edp(eab0, m);
                    e0.second -= rk_edm(eab0, m);
                    if (e0.first == 0 && e0.second == 0) _ed.erase(d0);
                    if (eab0 < rk_floor_half(m)) _Dp--;
                    if (eab0 > rk_ceil_half(m)) _Dm--;
                }
                if (d1 > 0) {
                    long eab1 = (new_u > ut) ? etr : ert;
                    auto& e1 = _ed[d1];
                    e1.first += rk_edp(eab1, m);
                    e1.second += rk_edm(eab1, m);
                    if (eab1 < rk_floor_half(m)) _Dp++;
                    if (eab1 > rk_ceil_half(m)) _Dm++;
                }
            }
        }
        _u[r] = new_u;
        return dS;
    }

    // Swap ranks of two occupied groups r and s (valid even when unique=TRUE!)
    double swap_group_ranks_dS(int r, int s) {
        if (r == s || _u[r] == _u[s]) return 0.0;
        long ur = _u[r], us = _u[s];
        double dS1 = relabel_group(r, us);
        double dS2 = relabel_group(s, ur);
        relabel_group(s, us);
        relabel_group(r, ur);
        return dS1 + dS2;
    }

    double swap_group_ranks(int r, int s) {
        if (r == s || _u[r] == _u[s]) return 0.0;
        long ur = _u[r], us = _u[s];
        double dS1 = relabel_group(r, us);
        double dS2 = relabel_group(s, ur);
        return dS1 + dS2;
    }

    // Peixoto neighbor-based ergodic proposal
    int sample_group(int v, double c, RKRng& rng) const {
        if (_occupied.empty()) return 0;
        size_t dout = _adj_out[v].size(), din = _adj_in[v].size();
        size_t deg = dout + din;
        if (deg == 0 || std::isinf(c))
            return _occupied[rng.index(_occupied.size())];
        size_t idx = rng.index(deg);
        int nbr = (idx < dout) ? _adj_out[v][idx] : _adj_in[v][idx - dout];
        int t = _b[nbr];
        long kt = _ep[t] + _em[t];
        double B_act = static_cast<double>(_occupied.size());
        double p_rand = (c * B_act) / (static_cast<double>(kt) + c * B_act);
        if (rng.unif() < p_rand)
            return _occupied[rng.index(_occupied.size())];
        double u = rng.unif() * static_cast<double>(kt);
        double acc = 0.0;
        for (const auto& [s, m] : _m_out[t]) {
            acc += static_cast<double>(m);
            if (u <= acc) return s;
        }
        for (const auto& [s, m] : _m_in[t]) {
            acc += static_cast<double>(m);
            if (u <= acc) return s;
        }
        return _occupied[rng.index(_occupied.size())];
    }

    double move_prob(int v, int s, double c) const {
        if (_occupied.empty()) return 0.0;
        size_t dout = _adj_out[v].size(), din = _adj_in[v].size();
        size_t deg = dout + din;
        double B_act = static_cast<double>(_occupied.size());
        if (deg == 0 || std::isinf(c)) return 1.0 / B_act;
        double prob = 0.0;
        auto add_nbr = [&](int nbr) {
            int t = _b[nbr];
            long kt = _ep[t] + _em[t];
            long m_ts = 0;
            auto it1 = _m_out[t].find(s);
            if (it1 != _m_out[t].end()) m_ts += it1->second;
            auto it2 = _m_in[t].find(s);
            if (it2 != _m_in[t].end()) m_ts += it2->second;
            prob += (static_cast<double>(m_ts) + c) / (static_cast<double>(kt) + c * B_act);
        };
        for (int nbr : _adj_out[v]) add_nbr(nbr);
        for (int nbr : _adj_in[v]) add_nbr(nbr);
        return prob / static_cast<double>(deg);
    }

    // Export contiguous 1-based blocks and matching per-block ranks u[1..B]
    std::pair<Rcpp::IntegerVector, Rcpp::NumericVector> external_blocks_and_u() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        rk_hash_map<int, int> gmap;
        Rcpp::NumericVector u_out(occ.size());
        for (size_t i = 0; i < occ.size(); ++i) {
            gmap.emplace(occ[i], static_cast<int>(i) + 1);
            u_out[i] = static_cast<double>(_u[occ[i]]);
        }
        Rcpp::IntegerVector b_out(_N);
        for (int v = 0; v < _N; ++v) b_out[v] = gmap[_b[v]];
        return {b_out, u_out};
    }

    // Export undirected block graph edges (1-based lower blocks 0..B-1) for upper nested levels
    int coarsen_undirected(std::vector<std::tuple<int, int, long>>& out_edges) const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        rk_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i)
            gmap.emplace(occ[i], static_cast<int>(i));

        rk_hash_map<uint64_t, long> undir_m;
        for (const auto& [pk, m] : _dir_pairs) {
            if (m <= 0) continue;
            auto [r, s] = rk_unpack(pk);
            int gr = gmap[r], gs = gmap[s];
            undir_m[rk_undir_pack(gr, gs)] += m;
        }
        out_edges.clear();
        out_edges.reserve(undir_m.size());
        for (const auto& [upk, m] : undir_m) {
            auto [a, b_g] = rk_unpack(upk);
            out_edges.emplace_back(a, b_g, m);
        }
        return static_cast<int>(occ.size());
    }

private:
    struct PairDelta {
        uint64_t key = 0;  // directed pack (a, b)
        long dm = 0;
    };

    struct MoveBuf {
        std::vector<PairDelta> pdeltas;
        std::array<long, 4> dE = {0, 0, 0, 0};
        rk_hash_map<long, std::pair<long, long>> ed_d;
        long dDp = 0, dDm = 0;
    };

    int _N;
    bool _distance, _unique;
    long _M;
    bool _dc, _dist, _is_top;
    int _cap;
    long _E = 0;
    double _base_const = 0.0;

    std::vector<std::pair<int, int>> _edges;
    std::vector<std::vector<int>> _adj_out, _adj_in;
    std::vector<long> _kout, _kin;

    std::vector<int> _b;
    std::vector<long> _u;
    std::vector<int> _wr;
    std::vector<long> _ep, _em;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;

    rk_hash_map<uint64_t, long> _dir_pairs;
    std::vector<rk_hash_map<int, long>> _m_out, _m_in;
    std::array<long, 4> _Es = {0, 0, 0, 0};
    rk_hash_map<long, std::pair<long, long>> _ed;
    long _Dp = 0, _Dm = 0;

    std::vector<rk_hash_map<long, long>> _h_out, _h_in;
    mutable MoveBuf _scratch_buf;
    mutable rk_hash_map<uint64_t, size_t> _scratch_idx;

    static inline int udir(long x, long y) {
        if (x < y) return 0;   // upstream
        if (x > y) return 2;   // downstream
        return 1;              // lateral
    }

    inline int stream_dir(int r, int s) const { return udir(_u[r], _u[s]); }

    double deg_dl_group(int r) const {
        if (!_dc || _wr[r] == 0) return 0.0;
        const RKMath& M = rkmath();
        long n = _wr[r];
        if (!_dist) {
            return M.lbinom(n + _ep[r] - 1, _ep[r]) + M.lbinom(n + _em[r] - 1, _em[r]);
        }
        const RKLogQ& q = rklog_q();
        double s = q(_ep[r], n) + M.lgam1(n) + q(_em[r], n) + M.lgam1(n);
        for (const auto& [dk, cnt] : _h_out[r]) s -= M.lgam1(cnt);
        for (const auto& [dk, cnt] : _h_in[r]) s -= M.lgam1(cnt);
        return s;
    }

    double deg_dl_delta_group(int r, long n0, long n1, long ep0, long ep1,
                              long em0, long em1, long ko_v, long ki_v, long dcnt) const {
        if (!_dc) return 0.0;
        const RKMath& M = rkmath();
        if (!_dist) {
            double s0 = (n0 > 0) ? (M.lbinom(n0 + ep0 - 1, ep0) + M.lbinom(n0 + em0 - 1, em0)) : 0.0;
            double s1 = (n1 > 0) ? (M.lbinom(n1 + ep1 - 1, ep1) + M.lbinom(n1 + em1 - 1, em1)) : 0.0;
            return s1 - s0;
        }
        const RKLogQ& q = rklog_q();
        double s0 = (n0 > 0) ? (q(ep0, n0) + M.lgam1(n0) + q(em0, n0) + M.lgam1(n0)) : 0.0;
        double s1 = (n1 > 0) ? (q(ep1, n1) + M.lgam1(n1) + q(em1, n1) + M.lgam1(n1)) : 0.0;
        auto ito = _h_out[r].find(ko_v);
        long co0 = (ito == _h_out[r].end()) ? 0 : ito->second;
        auto iti = _h_in[r].find(ki_v);
        long ci0 = (iti == _h_in[r].end()) ? 0 : iti->second;
        return (s1 - s0) - (M.lgam1(co0 + dcnt) - M.lgam1(co0))
                         - (M.lgam1(ci0 + dcnt) - M.lgam1(ci0));
    }

    double compute_move_delta(int v, int r, int s, MoveBuf* out) const {
        const RKMath& M = rkmath();
        MoveBuf& buf = out ? *out : _scratch_buf;
        auto& pds = buf.pdeltas;
        pds.clear();
        buf.dE = {0, 0, 0, 0};
        buf.ed_d.clear();
        buf.dDp = 0;
        buf.dDm = 0;

        bool use_map = false;
        auto add_pd = [&](uint64_t key, long dm) {
            if (!use_map) {
                for (size_t i = 0; i < pds.size(); ++i) {
                    if (pds[i].key == key) { pds[i].dm += dm; return; }
                }
                if (pds.size() >= 24) {
                    use_map = true;
                    _scratch_idx.clear();
                    for (size_t i = 0; i < pds.size(); ++i) _scratch_idx.emplace(pds[i].key, i);
                } else {
                    pds.push_back({key, dm});
                    return;
                }
            }
            auto it = _scratch_idx.find(key);
            if (it == _scratch_idx.end()) {
                _scratch_idx.emplace(key, pds.size());
                pds.push_back({key, dm});
            } else {
                pds[it->second].dm += dm;
            }
        };

        for (int u : _adj_out[v]) {
            if (u == v) {
                add_pd(rk_pack(r, r), -1);
                add_pd(rk_pack(s, s), +1);
            } else {
                int t = _b[u];
                add_pd(rk_pack(r, t), -1);
                add_pd(rk_pack(s, t), +1);
            }
        }
        for (int u : _adj_in[v]) {
            if (u == v) continue;  // self-loop already handled above
            int t = _b[u];
            add_pd(rk_pack(t, r), -1);
            add_pd(rk_pack(t, s), +1);
        }

        auto get_dm = [&](int a, int b_g) -> long {
            uint64_t k = rk_pack(a, b_g);
            if (!use_map) {
                for (const auto& pd : pds) if (pd.key == k) return pd.dm;
                return 0;
            }
            auto it = _scratch_idx.find(k);
            return (it == _scratch_idx.end()) ? 0 : pds[it->second].dm;
        };

        long wr_r0 = _wr[r], wr_r1 = wr_r0 - 1;
        long wr_s0 = _wr[s], wr_s1 = wr_s0 + 1;
        long B0 = B();
        long B1 = B0 - (wr_r1 == 0 ? 1 : 0) + (wr_s0 == 0 ? 1 : 0);

        double dS = 0.0;
        // Partition DL + top-level undirected edges DL + label DL + distance N-choose-B
        if (B0 != B1) {
            dS += M.lbinom(_N - 1, B1 - 1) - M.lbinom(_N - 1, B0 - 1);
            if (_is_top && _E > 0) {
                int64_t p0 = int64_t(B0) * (B0 + 1) / 2;
                int64_t p1 = int64_t(B1) * (B1 + 1) / 2;
                dS += M.lbinom(p1 + _E - 1, _E) - M.lbinom(p0 + _E - 1, _E);
            }
            if (_distance) dS += M.lbinom(_N, B1) - M.lbinom(_N, B0);
            if (!_unique) dS += M.lbinom(_M + B1 - 1, B1) - M.lbinom(_M + B0 - 1, B0);
        }
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_s1) - M.lgam1(wr_s0));

        // Vertex/group degree + degree DL
        long ko = _kout[v], ki = _kin[v];
        long ep_r0 = _ep[r], ep_r1 = ep_r0 - ko;
        long ep_s0 = _ep[s], ep_s1 = ep_s0 + ko;
        long em_r0 = _em[r], em_r1 = em_r0 - ki;
        long em_s0 = _em[s], em_s1 = em_s0 + ki;
        if (_dc) {
            dS += (M.lgam1(ep_r1) - M.lgam1(ep_r0)) + (M.lgam1(ep_s1) - M.lgam1(ep_s0));
            dS += (M.lgam1(em_r1) - M.lgam1(em_r0)) + (M.lgam1(em_s1) - M.lgam1(em_s0));
            dS += deg_dl_delta_group(r, wr_r0, wr_r1, ep_r0, ep_r1, em_r0, em_r1, ko, ki, -1);
            dS += deg_dl_delta_group(s, wr_s0, wr_s1, ep_s0, ep_s1, em_s0, em_s1, ko, ki, +1);
        } else {
            long dr0 = ep_r0 + em_r0, dr1 = ep_r1 + em_r1;
            long ds0 = ep_s0 + em_s0, ds1 = ep_s1 + em_s1;
            dS += (wr_r1 > 0 ? dr1 * M.logn(wr_r1) : 0.0) - dr0 * M.logn(wr_r0);
            dS += ds1 * M.logn(wr_s1) - (wr_s0 > 0 ? ds0 * M.logn(wr_s0) : 0.0);
        }

        // Directed base eterm deltas + _Es deltas + unordered pair orientation deltas
        rk_hash_set<uint64_t> seen_undir;
        for (const auto& pd : pds) {
            if (pd.dm == 0) continue;
            auto [a, b_g] = rk_unpack(pd.key);
            long m0 = get_mrs(a, b_g);
            long m1 = m0 + pd.dm;
            dS -= M.lgam1(m1) - M.lgam1(m0);

            buf.dE[stream_dir(a, b_g)] += pd.dm;
            if (a == b_g) {
                buf.dE[3] += pd.dm;
                continue;
            }

            uint64_t upk = rk_undir_pack(a, b_g);
            if (!seen_undir.insert(upk).second) continue;

            bool same = (_u[a] == _u[b_g]);
            int p = (!same && _u[a] > _u[b_g]) ? b_g : a;
            int q = (!same && _u[a] > _u[b_g]) ? a : b_g;
            long epq0 = get_mrs(p, q), eqp0 = get_mrs(q, p);
            long dm_pq = get_dm(p, q), dm_qp = get_dm(q, p);
            long epq1 = epq0 + dm_pq, eqp1 = eqp0 + dm_qp;
            long m_tot0 = epq0 + eqp0, m_tot1 = epq1 + eqp1;

            dS += rk_pair_orient_term(epq1, m_tot1, _distance, same) -
                  rk_pair_orient_term(epq0, m_tot0, _distance, same);

            if (_distance && !same) {
                long d = std::abs(_u[p] - _u[q]);
                auto& dd = buf.ed_d[d];
                dd.first += rk_edp(epq1, m_tot1) - rk_edp(epq0, m_tot0);
                dd.second += rk_edm(epq1, m_tot1) - rk_edm(epq0, m_tot0);
                buf.dDp += (epq1 < rk_floor_half(m_tot1) ? 1 : 0) -
                           (epq0 < rk_floor_half(m_tot0) ? 1 : 0);
                buf.dDm += (epq1 > rk_ceil_half(m_tot1) ? 1 : 0) -
                           (epq0 > rk_ceil_half(m_tot0) ? 1 : 0);
            }
        }

        if (!_distance) {
            dS += rk_Sd(_Es[0] + buf.dE[0], _Es[2] + buf.dE[2]) - rk_Sd(_Es[0], _Es[2]);
        } else {
            for (const auto& [d, dd] : buf.ed_d) {
                auto it = _ed.find(d);
                long ep0 = (it == _ed.end()) ? 0 : it->second.first;
                long em0 = (it == _ed.end()) ? 0 : it->second.second;
                dS += rk_Sd(ep0 + dd.first, em0 + dd.second) - rk_Sd(ep0, em0);
            }
            dS += rk_Sd(_Dp + buf.dDp, _Dm + buf.dDm) - rk_Sd(_Dp, _Dm);
        }
        dS += static_cast<double>(buf.dE[1] - buf.dE[3]) * RK_LOG2;
        return dS;
    }

    void apply_move(int v, int r, int s, const MoveBuf& buf) {
        if (_dc && _dist) {
            auto ito = _h_out[r].find(_kout[v]);
            if (ito != _h_out[r].end() && --ito->second == 0) _h_out[r].erase(ito);
            _h_out[s][_kout[v]]++;
            auto iti = _h_in[r].find(_kin[v]);
            if (iti != _h_in[r].end() && --iti->second == 0) _h_in[r].erase(iti);
            _h_in[s][_kin[v]]++;
        }
        remove_member(v, r);
        _b[v] = s;
        add_member(v, s);
        _ep[r] -= _kout[v]; _ep[s] += _kout[v];
        _em[r] -= _kin[v];  _em[s] += _kin[v];
        for (int i = 0; i < 4; ++i) _Es[i] += buf.dE[i];
        if (_distance) {
            for (const auto& [d, dd] : buf.ed_d) {
                auto& entry = _ed[d];
                entry.first += dd.first;
                entry.second += dd.second;
                if (entry.first == 0 && entry.second == 0) _ed.erase(d);
            }
            _Dp += buf.dDp;
            _Dm += buf.dDm;
        }
        auto add_m = [](rk_hash_map<int, long>& mp, int key_g, long d) {
            if (d == 0) return;
            auto mit = mp.find(key_g);
            if (mit == mp.end()) {
                if (d > 0) mp.emplace(key_g, d);
            } else {
                mit->second += d;
                if (mit->second <= 0) mp.erase(mit);
            }
        };
        for (const auto& pd : buf.pdeltas) {
            if (pd.dm == 0) continue;
            auto [a, b_g] = rk_unpack(pd.key);
            auto it = _dir_pairs.find(pd.key);
            if (it == _dir_pairs.end()) {
                if (pd.dm > 0) _dir_pairs.emplace(pd.key, pd.dm);
            } else {
                it->second += pd.dm;
                if (it->second <= 0) _dir_pairs.erase(it);
            }
            add_m(_m_out[a], b_g, pd.dm);
            add_m(_m_in[b_g], a, pd.dm);
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

// ---------------------------------------------------------------------------
// Dense undirected upper level for NestedRankedSBM (levels l >= 1)

class DenseUndirLevel {
public:
    DenseUndirLevel(int N, bool is_top)
        : _N(N), _is_top(is_top), _cap(std::max(2, N)),
          _adj(N), _b(N, 0), _wr(_cap, 0), _ep(_cap, 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(N, 0), _m_out(_cap) {
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
    }

    void set_edges(std::vector<std::tuple<int, int, long>> edges) {
        _edges = std::move(edges);
        _E = 0;
        _adj.assign(_N, {});
        for (int ei = 0; ei < static_cast<int>(_edges.size()); ++ei) {
            auto [u, v, m] = _edges[ei];
            _E += m;
            _adj[u].push_back({v, ei});
            if (u != v) _adj[v].push_back({u, ei});
        }
        rkmath().ensure(static_cast<size_t>(_N + _E + 64));
    }

    int N() const { return _N; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }
    int random_occupied(RKRng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    void set_partition(const std::vector<int>& b) {
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            _members[r].clear();
            _occ_pos[r] = -1;
            _m_out[r].clear();
        }
        _occupied.clear();
        _empty.clear();
        _empty_pos.assign(_cap, -1);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        _pairs.clear();
        for (int v = 0; v < _N; ++v) {
            _b[v] = b[v];
            add_member(v, b[v]);
        }
        for (const auto& [u, v, m] : _edges) {
            int r = _b[u], s = _b[v];
            _pairs[rk_undir_pack(r, s)] += m;
            _ep[r] += m;
            _ep[s] += m;
            if (r == s) _m_out[r][r] += 2 * m;
            else { _m_out[r][s] += m; _m_out[s][r] += m; }
        }
    }

    double entropy() const {
        const RKMath& M = rkmath();
        long Bn = B();
        double S = 0.0;
        if (_N > 0) {
            S += M.logn(_N) + M.lbinom(_N - 1, Bn - 1) + M.lgam1(_N);
            for (int r : _occupied) S -= M.lgam1(_wr[r]);
        }
        if (_is_top && _E > 0) {
            int64_t P = int64_t(Bn) * (Bn + 1) / 2;
            S += M.lbinom(P + _E - 1, _E);
        }
        for (const auto& [pk, m] : _pairs) {
            if (m <= 0) continue;
            auto [r, s] = rk_unpack(pk);
            int64_t nn = (r == s) ? (int64_t(_wr[r]) * (_wr[r] + 1) / 2)
                                  : (int64_t(_wr[r]) * _wr[s]);
            S += M.lbinom(nn + m - 1, m);
        }
        return S;
    }

    double move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double s0 = entropy();
        // Simple exact update (upper levels are small!)
        for (const auto& [nbr, ei] : _adj[v]) {
            auto [eu, ev, m] = _edges[ei];
            if (eu == ev) {
                add_pair(r, r, -m);
                add_pair(s, s, +m);
                _ep[r] -= 2 * m;
                _ep[s] += 2 * m;
            } else {
                int t = _b[nbr];
                add_pair(r, t, -m);
                add_pair(s, t, +m);
                _ep[r] -= m;
                _ep[s] += m;
            }
        }
        remove_member(v, r);
        _b[v] = s;
        add_member(v, s);
        return entropy() - s0;
    }

    double virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double dS = move(v, s);
        move(v, r);
        return dS;
    }

    double merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        return dS;
    }

    double virtual_merge(int r, int t) {
        if (r == t || _wr[r] == 0) return 0.0;
        std::vector<int> verts = _members[r];
        double dS = 0.0;
        for (int v : verts) dS += move(v, t);
        for (int v : verts) move(v, r);
        return dS;
    }

    int sample_group(int v, double c, RKRng& rng) const {
        if (_occupied.empty()) return 0;
        if (_adj[v].empty() || std::isinf(c))
            return _occupied[rng.index(_occupied.size())];
        int nbr = _adj[v][rng.index(_adj[v].size())].first;
        int t = _b[nbr];
        long kt = _ep[t];
        double B_act = static_cast<double>(_occupied.size());
        if (rng.unif() < (c * B_act) / (static_cast<double>(kt) + c * B_act))
            return _occupied[rng.index(_occupied.size())];
        double u = rng.unif() * static_cast<double>(kt);
        double acc = 0.0;
        for (const auto& [s, m] : _m_out[t]) {
            acc += static_cast<double>(m);
            if (u <= acc) return s;
        }
        return _occupied[rng.index(_occupied.size())];
    }

    double move_prob(int v, int s, double c) const {
        if (_occupied.empty()) return 0.0;
        double B_act = static_cast<double>(_occupied.size());
        if (_adj[v].empty() || std::isinf(c)) return 1.0 / B_act;
        double prob = 0.0;
        for (const auto& [nbr, ei] : _adj[v]) {
            (void)ei;
            int t = _b[nbr];
            long kt = _ep[t];
            auto it = _m_out[t].find(s);
            long m_ts = (it == _m_out[t].end()) ? 0 : it->second;
            prob += (static_cast<double>(m_ts) + c) / (static_cast<double>(kt) + c * B_act);
        }
        return prob / static_cast<double>(_adj[v].size());
    }

    int coarsen_undirected(std::vector<std::tuple<int, int, long>>& out_edges) const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        rk_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i) gmap.emplace(occ[i], static_cast<int>(i));
        out_edges.clear();
        for (const auto& [pk, m] : _pairs) {
            if (m <= 0) continue;
            auto [r, s] = rk_unpack(pk);
            out_edges.emplace_back(gmap[r], gmap[s], m);
        }
        return static_cast<int>(occ.size());
    }

    Rcpp::IntegerVector external_blocks_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        rk_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i) gmap.emplace(occ[i], static_cast<int>(i) + 1);
        Rcpp::IntegerVector out(_N);
        for (int v = 0; v < _N; ++v) out[v] = gmap[_b[v]];
        return out;
    }

private:
    int _N;
    bool _is_top;
    int _cap;
    long _E = 0;
    std::vector<std::tuple<int, int, long>> _edges;
    std::vector<std::vector<std::pair<int, int>>> _adj;
    std::vector<int> _b, _wr;
    std::vector<long> _ep;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;
    rk_hash_map<uint64_t, long> _pairs;
    std::vector<rk_hash_map<int, long>> _m_out;

    void add_pair(int a, int b_g, long dm) {
        uint64_t pk = rk_undir_pack(a, b_g);
        auto it = _pairs.find(pk);
        if (it == _pairs.end()) {
            if (dm > 0) _pairs.emplace(pk, dm);
        } else {
            it->second += dm;
            if (it->second <= 0) _pairs.erase(it);
        }
        auto add_m = [](rk_hash_map<int, long>& mp, int k, long d) {
            auto mit = mp.find(k);
            if (mit == mp.end()) { if (d > 0) mp.emplace(k, d); }
            else { mit->second += d; if (mit->second <= 0) mp.erase(mit); }
        };
        if (a == b_g) add_m(_m_out[a], a, 2 * dm);
        else { add_m(_m_out[a], b_g, dm); add_m(_m_out[b_g], a, dm); }
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

// ---------------------------------------------------------------------------
// Rank optimization & label MCMC sweeps on RankedState

void ranked_optimize_u(RankedState& st, RKRng& rng, int max_passes = 8) {
    std::vector<int> occ = st.occupied();
    if (occ.size() <= 1) return;

    // 1. Greedy pairwise swaps (works for BOTH unique=TRUE and unique=FALSE!)
    for (int pass = 0; pass < max_passes; ++pass) {
        bool improved = false;
        rng.shuffle(occ);
        for (size_t i = 0; i < occ.size(); ++i) {
            int r = occ[i];
            if (!st.unique()) {
                // Try setting u[r] to candidate ranks: existing ranks, +/-1 neighbors, and endpoints
                rk_hash_set<long> cand_u;
                long min_u = st.rank_of(occ[0]), max_u = min_u;
                for (int g : occ) {
                    long ug = st.rank_of(g);
                    cand_u.insert(ug);
                    if (ug < min_u) min_u = ug;
                    if (ug > max_u) max_u = ug;
                }
                if (min_u > 0) cand_u.insert(min_u - 1);
                if (max_u + 1 < st.M()) cand_u.insert(max_u + 1);
                for (int g : occ) {
                    long ug = st.rank_of(g);
                    if (ug > 0) cand_u.insert(ug - 1);
                    if (ug + 1 < st.M()) cand_u.insert(ug + 1);
                }
                long best_u = st.rank_of(r);
                double best_dS = -1e-11;
                for (long nu : cand_u) {
                    if (nu == st.rank_of(r) || nu < 0 || nu >= st.M()) continue;
                    double dS = st.relabel_group_dS(r, nu);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best_u = nu;
                    }
                }
                if (best_u != st.rank_of(r)) {
                    st.relabel_group(r, best_u);
                    improved = true;
                }
            }
            // Also try swapping rank with another occupied group
            if (occ.size() <= 30) {
                int best_s = -1;
                double best_dS = -1e-11;
                for (size_t j = 0; j < occ.size(); ++j) {
                    if (occ[j] == r) continue;
                    double dS = st.swap_group_ranks_dS(r, occ[j]);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best_s = occ[j];
                    }
                }
                if (best_s >= 0) {
                    st.swap_group_ranks(r, best_s);
                    improved = true;
                }
            } else {
                for (int k = 0; k < 8; ++k) {
                    int s = occ[rng.index(occ.size())];
                    if (s == r) continue;
                    double dS = st.swap_group_ranks_dS(r, s);
                    if (dS < -1e-11) {
                        st.swap_group_ranks(r, s);
                        improved = true;
                    }
                }
            }
        }
        if (!improved) break;
    }
    // Canonicalize orientation so upstream >= downstream (pure symmetry when distance=FALSE)
    if (!st.distance() && st.Es()[2] > st.Es()[0]) {
        long max_u = 0;
        for (int g : occ) max_u = std::max(max_u, st.rank_of(g));
        // Reflect all ranks: u' = max_u - u (exact symmetry: swaps E[0] and E[2]!)
        std::vector<long> new_u(occ.size());
        for (size_t i = 0; i < occ.size(); ++i) new_u[i] = max_u - st.rank_of(occ[i]);
        // Apply via temporary negative offset to avoid collisions if needed
        for (size_t i = 0; i < occ.size(); ++i) st.relabel_group(occ[i], -1000000 - static_cast<long>(i));
        for (size_t i = 0; i < occ.size(); ++i) st.relabel_group(occ[i], new_u[i]);
    }
}

void ranked_label_mcmc_sweep(RankedState& st, RKRng& rng, double beta,
                             double& total_dS, size_t& nattempts, size_t& nmoves) {
    std::vector<int> occ = st.occupied();
    if (occ.size() <= 1) return;
    rng.shuffle(occ);

    for (int r : occ) {
        ++nattempts;
        if (st.unique()) {
            // Symmetric swap proposal between r and another occupied group s
            int s = occ[rng.index(occ.size())];
            if (s == r) continue;
            double dS = st.swap_group_ranks_dS(r, s);
            bool accept = std::isinf(beta)
                              ? (dS < 0.0)
                              : (dS <= 0.0 || rng.unif() < std::exp(-beta * dS));
            if (accept) {
                st.swap_group_ranks(r, s);
                total_dS += dS;
                ++nmoves;
            }
        } else {
            // Symmetric proposal over {0, ..., M - 1} or swap with another group
            if (rng.unif() < 0.5) {
                int s = occ[rng.index(occ.size())];
                if (s == r) continue;
                double dS = st.swap_group_ranks_dS(r, s);
                bool accept = std::isinf(beta)
                                  ? (dS < 0.0)
                                  : (dS <= 0.0 || rng.unif() < std::exp(-beta * dS));
                if (accept) {
                    st.swap_group_ranks(r, s);
                    total_dS += dS;
                    ++nmoves;
                }
            } else {
                // Uniformly pick among {0, ..., M-1} (symmetric proposal -> Hastings ratio = 1!)
                long nu = static_cast<long>(rng.index(static_cast<size_t>(st.M())));
                if (nu == st.rank_of(r)) continue;
                double dS = st.relabel_group_dS(r, nu);
                bool accept = std::isinf(beta)
                                  ? (dS < 0.0)
                                  : (dS <= 0.0 || rng.unif() < std::exp(-beta * dS));
                if (accept) {
                    st.relabel_group(r, nu);
                    total_dS += dS;
                    ++nmoves;
                }
            }
        }
    }
}

void ranked_greedy_sweeps(RankedState& st, RKRng& rng, size_t B_min, size_t B_max,
                          int max_sweeps, double epsilon) {
    if (st.N() <= 1) return;
    std::vector<int> vs(st.N());
    std::iota(vs.begin(), vs.end(), 0);
    for (int sw = 0; sw < max_sweeps; ++sw) {
        double tot_dS = 0.0;
        rng.shuffle(vs);
        for (int v : vs) {
            int r = st.group(v);
            size_t B_cur = static_cast<size_t>(st.B());
            double best_dS = -1e-11;
            int best_s = r;
            rk_hash_set<int> cands;
            const auto& occ = st.occupied();
            if (occ.size() <= 14) {
                for (int g : occ) cands.insert(g);
            } else {
                for (int k = 0; k < 8; ++k) cands.insert(st.sample_group(v, 0.5, rng));
                for (int k = 0; k < 4; ++k) cands.insert(st.random_occupied(rng));
            }
            if (B_cur < B_max && B_cur < static_cast<size_t>(st.N()) &&
                st.size(r) > 1 && st.empty_group() >= 0) {
                cands.insert(st.empty_group());
            }
            cands.erase(r);
            for (int s : cands) {
                if (st.size(r) == 1 && st.size(s) > 0 && B_cur <= B_min) continue;
                double dS = st.virtual_move(v, s);
                if (dS < best_dS) { best_dS = dS; best_s = s; }
            }
            if (best_s != r) tot_dS += st.move(v, best_s);
        }
        if (sw % 2 == 0) ranked_optimize_u(st, rng, 3);
        if (std::abs(tot_dS) <= epsilon) break;
    }
}

void ranked_flat_mcmc_sweep(RankedState& s, RKRng& rng, double beta, double c,
                            double d, bool multiflip, size_t B_min, size_t B_max,
                            double& total_dS, size_t& nattempts, size_t& nmoves);

double ranked_minimize_flat(RankedState& st, RKRng& rng, size_t B_min, size_t B_max,
                            bool refine, double epsilon,
                            std::vector<std::pair<std::vector<int>, std::vector<long>>>* cps = nullptr) {
    int N = st.N();
    if (N <= 1) {
        std::vector<int> b(N, 0);
        st.set_partition(b);
        if (cps) cps->push_back({b, std::vector<long>(1, 0)});
        return st.entropy();
    }
    B_min = std::max<size_t>(1, std::min<size_t>(B_min, static_cast<size_t>(N)));
    B_max = std::max<size_t>(B_min, std::min<size_t>(B_max, static_cast<size_t>(N)));

    auto record_cp = [&]() {
        if (!cps) return;
        auto [eb, eu] = st.external_blocks_and_u();
        std::vector<int> b(N);
        for (int i = 0; i < N; ++i) b[i] = eb[i] - 1;
        std::vector<long> u(eu.size());
        for (int i = 0; i < eu.size(); ++i) u[i] = static_cast<long>(eu[i]);
        cps->push_back({std::move(b), std::move(u)});
    };

    std::vector<int> best_b(N, 0);
    std::vector<long> best_u(N, 0);
    double best_S = RK_INF;

    // Two-stage fast minimization: brief thermal mixing (beta=1 -> 2 -> Inf)
    // crosses shallow rank/block barriers cleanly, landing at deep minima!
    std::vector<size_t> init_Bs = {B_min};
    size_t mid_B = std::min<size_t>(B_max, std::max<size_t>(B_min, std::min<size_t>(N, 16)));
    if (mid_B != B_min) init_Bs.push_back(mid_B);

    for (size_t b0_k : init_Bs) {
        std::vector<int> b_init(N);
        for (int v = 0; v < N; ++v) b_init[v] = static_cast<int>(v % b0_k);
        st.set_partition(b_init);
        double tot_dS = 0.0;
        size_t na = 0, nm = 0;

        // Warm thermal sweeps (beta = 1 -> 2 -> 4) to escape local traps (~0.03s)
        int warm_sw = refine ? 120 : 40;
        for (int sw = 0; sw < warm_sw; ++sw) {
            double b_cur = (sw < warm_sw / 2) ? 1.0 : ((sw < 3 * warm_sw / 4) ? 2.0 : 4.0);
            ranked_flat_mcmc_sweep(st, rng, b_cur, 0.5, 0.05, true, B_min, B_max, tot_dS, na, nm);
            if (sw % 10 == 0) ranked_optimize_u(st, rng, 3);
        }

        // Greedy quench at beta = Inf
        double cur_S = st.entropy();
        int no_imp = 0;
        int max_sw = refine ? 40 : 15;
        for (int sw = 0; sw < max_sw; ++sw) {
            ranked_flat_mcmc_sweep(st, rng, RK_INF, 0.5, 0.05, true, B_min, B_max, tot_dS, na, nm);
            if (sw % 4 == 0) ranked_optimize_u(st, rng, 4);
            double new_S = st.entropy();
            if (new_S < cur_S - epsilon) {
                cur_S = new_S;
                no_imp = 0;
            } else if (++no_imp >= 5) {
                break;
            }
        }
        if (refine) {
            ranked_greedy_sweeps(st, rng, B_min, B_max, 6, epsilon);
            ranked_optimize_u(st, rng, 8);
            cur_S = st.entropy();
        }
        if (cps) record_cp();
        if (cur_S < best_S) {
            best_S = cur_S;
            best_b = st.blocks();
            best_u = st.ranks();
        }
    }

    st.set_partition(best_b, &best_u);
    if (!st.distance() && st.B() > 0) {
        // Compress distinct ranks to contiguous 0..K-1 (exact ordinal symmetry)
        // and orient so upstream >= downstream.
        std::vector<int> occ = st.occupied();
        std::vector<long> vals;
        vals.reserve(occ.size());
        for (int g : occ) vals.push_back(st.rank_of(g));
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        long K = static_cast<long>(vals.size());
        bool flip = (st.Es()[2] > st.Es()[0]);
        std::vector<long> new_u(occ.size());
        for (size_t i = 0; i < occ.size(); ++i) {
            long idx = static_cast<long>(
                std::lower_bound(vals.begin(), vals.end(), st.rank_of(occ[i])) - vals.begin());
            new_u[i] = flip ? (K - 1 - idx) : idx;
        }
        for (size_t i = 0; i < occ.size(); ++i)
            st.relabel_group(occ[i], -1000000 - static_cast<long>(i));
        for (size_t i = 0; i < occ.size(); ++i)
            st.relabel_group(occ[i], new_u[i]);
    }
    return st.entropy();
}

void ranked_flat_mcmc_sweep(RankedState& s, RKRng& rng, double beta, double c,
                            double d, bool multiflip, size_t B_min, size_t B_max,
                            double& total_dS, size_t& nattempts, size_t& nmoves) {
    if (s.N() <= 1 || s.E() == 0) return;
    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    rng.shuffle(vs);

    // 1. Single-node vertex moves
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
            if (dS < 0.0) { s.apply_last_virtual_move(v, t); total_dS += dS; ++nmoves; }
            continue;
        }
        if (pfwd <= 0.0) continue;
        double u_rand = rng.unif();
        double log_u = (u_rand > 0.0) ? std::log(u_rand) : -RK_INF;
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
        if (accept) { total_dS += dS; ++nmoves; }
        else s.move(v, r);
    }

    // 2. Group rank label sweep
    ranked_label_mcmc_sweep(s, rng, beta, total_dS, nattempts, nmoves);

    // 3. Merge-split moves
    if (multiflip && s.B() >= 1) {
        size_t n_ms = std::max<size_t>(1, static_cast<size_t>(s.B()));
        for (size_t step = 0; step < n_ms; ++step) {
            bool try_split = (rng.unif() < 0.5);
            if (try_split) {
                if (static_cast<size_t>(s.B()) >= B_max || s.empty_group() < 0) continue;
                int r = s.random_occupied(rng);
                if (s.size(r) < 2) continue;
                ++nattempts;
                int t = s.empty_group();
                long old_ut = s.rank_of(t);
                if (!s.unique()) {
                    long ur = s.rank_of(r);
                    long init_ut = ur;
                    if (std::isinf(beta) && s.M() > 1) {
                        bool go_up = (rng.unif() < 0.5);
                        if (go_up) init_ut = (ur + 1 < s.M()) ? (ur + 1) : (ur - 1);
                        else       init_ut = (ur > 0) ? (ur - 1) : (ur + 1);
                    }
                    s.relabel_group(t, init_ut);
                }
                std::vector<int> mem = s.members(r);
                size_t nr = mem.size();
                size_t B_before = static_cast<size_t>(s.B());
                double dS = 0.0, log_q = 0.0;
                int cnt_t = 0;
                for (int v : mem) {
                    if (rng.unif() < 0.5) { dS += s.move(v, t); ++cnt_t; }
                }
                if (cnt_t == 0) { dS += s.move(mem[0], t); cnt_t = 1; }
                else if (cnt_t == static_cast<int>(nr)) { dS += s.move(mem[0], r); --cnt_t; }
                int n_gsw = std::isinf(beta) ? 8 : 3;
                for (int gsw = 0; gsw < n_gsw; ++gsw) {
                    rng.shuffle(mem);
                    for (int v : mem) {
                        int cur = s.group(v);
                        int oth = (cur == r) ? t : r;
                        if (s.size(cur) <= 1) continue;
                        double d_mv = s.virtual_move(v, oth);
                        double p_oth;
                        if (std::isinf(beta)) p_oth = (d_mv < 0.0) ? 1.0 : 0.0;
                        else {
                            double z = -beta * d_mv;
                            p_oth = (z > 35.0) ? 1.0 : (z < -35.0 ? 0.0 : 1.0 / (1.0 + std::exp(-z)));
                        }
                        if (rng.unif() < p_oth) {
                            s.apply_last_virtual_move(v, oth);
                            dS += d_mv;
                            if (gsw == n_gsw - 1) log_q += std::log(std::max(1e-300, p_oth));
                        } else {
                            if (gsw == n_gsw - 1) log_q += std::log(std::max(1e-300, 1.0 - p_oth));
                        }
                    }
                }
                if (std::isinf(beta) && !s.unique() && s.M() > 1) {
                    long ur = s.rank_of(r);
                    long best_u = s.rank_of(t);
                    double best_du = 0.0;
                    long cands[3] = {std::max<long>(0, ur - 1), ur, std::min<long>(s.M() - 1, ur + 1)};
                    for (long cu : cands) {
                        if (cu == s.rank_of(t)) continue;
                        double du = s.relabel_group_dS(t, cu);
                        if (du < best_du - 1e-11) { best_du = du; best_u = cu; }
                    }
                    if (best_u != s.rank_of(t)) dS += s.relabel_group(t, best_u);
                }
                bool accept = false;
                if (std::isinf(beta)) accept = (dS < 0.0);
                else {
                    size_t B_after = static_cast<size_t>(s.B());
                    double log_pf = -std::log(static_cast<double>(B_before)) + log_q;
                    double log_pr = -std::log(static_cast<double>(B_after) * (B_after - 1) / 2.0);
                    double log_a = -beta * dS + log_pr - log_pf;
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) { total_dS += dS; ++nmoves; }
                else {
                    std::vector<int> mem_t = s.members(t);
                    for (int v : mem_t) s.move(v, r);
                    if (!s.unique()) s.relabel_group(t, old_ut);
                }
            } else {
                if (static_cast<size_t>(s.B()) <= B_min || s.B() < 2) continue;
                ++nattempts;
                int r = s.random_occupied(rng);
                int t = s.random_occupied(rng);
                if (r == t) continue;
                size_t B_before = static_cast<size_t>(s.B());
                std::vector<int> mem_r = s.members(r);
                double dS = s.merge(r, t);
                bool accept = false;
                if (std::isinf(beta)) accept = (dS < 0.0);
                else {
                    size_t B_after = static_cast<size_t>(s.B());
                    double log_pf = -std::log(static_cast<double>(B_before) * (B_before - 1) / 2.0);
                    double log_pr = -std::log(static_cast<double>(B_after)) -
                                    static_cast<double>(mem_r.size()) * RK_LOG2;
                    double log_a = -beta * dS + log_pr - log_pf;
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) { total_dS += dS; ++nmoves; }
                else { for (int v : mem_r) s.move(v, r); }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Nested Ranked SBM (level 0 = RankedState(is_top=false), levels 1..L-1 = DenseUndirLevel)

class NestedRankedSBM {
public:
    NestedRankedSBM(const Graph& graph, bool distance, bool unique, long M_pos,
                    bool dc, bool dist)
        : _base(graph, distance, unique, M_pos, dc, dist, false) {}

    void build(const std::vector<std::vector<int>>& bs, const std::vector<long>* u_in = nullptr) {
        if (bs.size() < 2)
            Rcpp::stop("Nested ranked state requires at least 2 hierarchy levels");
        _base.set_partition(bs[0], u_in);
        rebuild_upper_with_bs(bs);
    }

    void rebuild_upper_with_bs(const std::vector<std::vector<int>>& bs) {
        _upper.clear();
        int L = static_cast<int>(bs.size());
        std::vector<std::tuple<int, int, long>> cur_edges;
        int cur_N = _base.coarsen_undirected(cur_edges);
        for (int l = 1; l < L; ++l) {
            bool is_top = (l == L - 1);
            DenseUndirLevel ulv(cur_N, is_top);
            ulv.set_edges(std::move(cur_edges));
            ulv.set_partition(bs[l]);
            if (l + 1 < L) cur_N = ulv.coarsen_undirected(cur_edges);
            _upper.push_back(std::move(ulv));
        }
    }

    void rebuild_upper_keep_bs() {
        std::vector<std::vector<int>> bs = save_bs();
        rebuild_upper_with_bs(bs);
    }

    int num_levels() const { return 1 + static_cast<int>(_upper.size()); }
    RankedState& base() { return _base; }
    const RankedState& base() const { return _base; }
    DenseUndirLevel& upper(int idx) { return _upper[idx]; }

    std::vector<double> level_entropies() const {
        std::vector<double> out(1 + _upper.size());
        out[0] = _base.entropy(true);
        for (size_t i = 0; i < _upper.size(); ++i) out[i + 1] = _upper[i].entropy();
        return out;
    }

    double entropy() const {
        double S = _base.entropy(true);
        for (const auto& u : _upper) S += u.entropy();
        return S;
    }

    std::vector<std::vector<int>> save_bs() const {
        std::vector<std::vector<int>> bs;
        bs.reserve(1 + _upper.size());
        auto [eb0, eu0] = _base.external_blocks_and_u();
        (void)eu0;
        std::vector<int> b0(eb0.size());
        for (int i = 0; i < eb0.size(); ++i) b0[i] = eb0[i] - 1;
        bs.push_back(std::move(b0));
        for (const auto& u : _upper) {
            Rcpp::IntegerVector eb = u.external_blocks_r();
            std::vector<int> bl(eb.size());
            for (int i = 0; i < eb.size(); ++i) bl[i] = eb[i] - 1;
            bs.push_back(std::move(bl));
        }
        return bs;
    }

    Rcpp::List external_bs_r() const {
        Rcpp::List out(1 + _upper.size());
        out[0] = _base.external_blocks_and_u().first;
        for (size_t i = 0; i < _upper.size(); ++i) out[i + 1] = _upper[i].external_blocks_r();
        return out;
    }

    Rcpp::NumericVector external_u_r() const {
        return _base.external_blocks_and_u().second;
    }

    // Build greedy optimal undirected hierarchy above base partition b0 + ranks u0
    void optimize_from_base(const std::vector<int>& b0, const std::vector<long>& u0,
                            RKRng& rng, double epsilon) {
        _base.set_partition(b0, &u0);
        _upper.clear();
        std::vector<std::tuple<int, int, long>> cur_edges;
        int cur_N = _base.coarsen_undirected(cur_edges);
        if (cur_N <= 1) {
            DenseUndirLevel ulv(cur_N, true);
            ulv.set_edges(std::move(cur_edges));
            ulv.set_partition(std::vector<int>(cur_N, 0));
            _upper.push_back(std::move(ulv));
            return;
        }
        while (cur_N > 1) {
            DenseUndirLevel ulv(cur_N, true);
            ulv.set_edges(cur_edges);
            // Fast greedy minimization of dense undirected upper level (starts from B=1 + singletons)
            std::vector<int> best_b(cur_N, 0);
            ulv.set_partition(best_b);
            double best_S = ulv.entropy();
            if (cur_N > 2) {
                // Start from min(cur_N, 8) groups and run fast greedy vertex moves + merges
                int k0 = std::min(cur_N, 8);
                std::vector<int> b_init(cur_N);
                for (int i = 0; i < cur_N; ++i) b_init[i] = i % k0;
                ulv.set_partition(b_init);
                std::vector<int> uvs(cur_N);
                std::iota(uvs.begin(), uvs.end(), 0);
                for (int sw = 0; sw < 12; ++sw) {
                    bool moved = false;
                    rng.shuffle(uvs);
                    for (int v : uvs) {
                        int r = ulv.group(v);
                        int best_t = r;
                        double best_d = -1e-10;
                        std::vector<int> occ = ulv.occupied();
                        for (int t : occ) {
                            if (t == r) continue;
                            double d = ulv.virtual_move(v, t);
                            if (d < best_d) { best_d = d; best_t = t; }
                        }
                        if (ulv.size(r) > 1 && ulv.empty_group() >= 0) {
                            int e_g = ulv.empty_group();
                            double d = ulv.virtual_move(v, e_g);
                            if (d < best_d) { best_d = d; best_t = e_g; }
                        }
                        if (best_t != r) { ulv.move(v, best_t); moved = true; }
                    }
                    if (!moved) break;
                }
                // Try pairwise merges while improving
                while (ulv.B() > 1) {
                    std::vector<int> occ = ulv.occupied();
                    int mr = -1, mt = -1;
                    double mdS = -1e-10;
                    for (size_t i = 0; i < occ.size(); ++i) {
                        for (size_t j = i + 1; j < occ.size(); ++j) {
                            double dS = ulv.virtual_merge(occ[i], occ[j]);
                            if (dS < mdS) { mdS = dS; mr = occ[i]; mt = occ[j]; }
                        }
                    }
                    if (mr < 0) break;
                    ulv.merge(mr, mt);
                }
                double s_try = ulv.entropy();
                if (s_try < best_S - 1e-10) { best_S = s_try; best_b = ulv.blocks(); }
            }
            ulv.set_partition(best_b);
            if (ulv.B() >= cur_N) {
                ulv.set_partition(std::vector<int>(cur_N, 0));
            }
            if (ulv.B() == 1) {
                _upper.push_back(std::move(ulv));
                break;
            }
            std::vector<int> bl_ext;
            {
                Rcpp::IntegerVector eb = ulv.external_blocks_r();
                bl_ext.resize(eb.size());
                for (int i = 0; i < eb.size(); ++i) bl_ext[i] = eb[i] - 1;
            }
            DenseUndirLevel ulv_mid(cur_N, false);
            ulv_mid.set_edges(std::move(cur_edges));
            ulv_mid.set_partition(bl_ext);
            cur_N = ulv_mid.coarsen_undirected(cur_edges);
            _upper.push_back(std::move(ulv_mid));
        }
        (void)rng;
        (void)epsilon;
    }

    void refine_joint(RKRng& rng, size_t B_min, size_t B_max, double epsilon) {
        double best_S = entropy();
        std::vector<int> vs(_base.N());
        std::iota(vs.begin(), vs.end(), 0);
        for (int pass = 0; pass < 3; ++pass) {
            bool improved = false;
            rng.shuffle(vs);
            for (int v : vs) {
                int r = _base.group(v);
                size_t B_cur = static_cast<size_t>(_base.B());
                rk_hash_set<int> cands;
                for (int k = 0; k < 4; ++k) cands.insert(_base.sample_group(v, 0.5, rng));
                cands.erase(r);
                int best_t = r;
                double best_cS = best_S - 1e-8;
                for (int t : cands) {
                    if (_base.size(r) == 1 && B_cur <= B_min) continue;
                    double d0 = _base.virtual_move(v, t);
                    if (d0 > 20.0) continue;
                    if (_base.size(r) == 1 || _base.size(t) == 0) continue;  // fixed-B fast step
                    _base.apply_last_virtual_move(v, t);
                    rebuild_upper_keep_bs();
                    double cS = entropy();
                    _base.move(v, r);
                    rebuild_upper_keep_bs();
                    if (cS < best_cS) { best_cS = cS; best_t = t; }
                }
                if (best_t != r) {
                    _base.move(v, best_t);
                    rebuild_upper_keep_bs();
                    best_S = entropy();
                    improved = true;
                }
            }
            ranked_optimize_u(_base, rng, 6);
            double s_u = entropy();
            if (s_u < best_S - 1e-8) { best_S = s_u; improved = true; }
            if (!improved) break;
        }
        (void)B_max;
        (void)epsilon;
    }

    void mcmc_sweep(RKRng& rng, double beta, double c, double d,
                    size_t B_min, size_t B_max, double& total_dS,
                    size_t& nattempts, size_t& nmoves) {
        double cur_S = entropy();
        // 1. Level-0 single-node moves (fixed-B and empty/create)
        std::vector<int> vs(_base.N());
        std::iota(vs.begin(), vs.end(), 0);
        rng.shuffle(vs);
        for (int v : vs) {
            ++nattempts;
            int r = _base.group(v);
            size_t B_cur = static_cast<size_t>(_base.B());
            bool can_create = (B_cur < static_cast<size_t>(_base.N()) &&
                               B_cur < B_max && d > 0.0 && _base.empty_group() >= 0);
            double p_new = can_create ? d : 0.0;
            bool propose_new = (can_create && rng.unif() < p_new);
            int t;
            double pfwd = 0.0;
            if (propose_new) {
                if (_base.size(r) == 1) continue;
                t = _base.empty_group();
                pfwd = p_new;
            } else {
                t = _base.sample_group(v, c, rng);
                if (t == r) continue;
                if (_base.size(r) == 1 && B_cur <= B_min) continue;
                pfwd = (1.0 - p_new) * _base.move_prob(v, t, c);
            }
            if (pfwd <= 0.0) continue;
            double d0 = _base.virtual_move(v, t);
            if (!std::isinf(beta) && -beta * d0 < -40.0) continue;

            bool empties = (_base.size(r) == 1);
            bool creates = (_base.size(t) == 0);
            if (!empties && !creates) {
                _base.apply_last_virtual_move(v, t);
                double prev = (1.0 - p_new) * _base.move_prob(v, r, c);
                rebuild_upper_keep_bs();
                double new_S = entropy();
                double dS = new_S - cur_S;
                bool accept = std::isinf(beta)
                                  ? (dS < 0.0)
                                  : (prev > 0.0 && (dS <= 0.0 || rng.unif() < std::exp(-beta * dS + std::log(prev) - std::log(pfwd))));
                if (accept) {
                    cur_S = new_S;
                    total_dS += dS;
                    ++nmoves;
                } else {
                    _base.move(v, r);
                    rebuild_upper_keep_bs();
                }
            }
        }

        // 2. Group rank label sweep (doesn't change undirected block graph -> O(deg_bg)!)
        ranked_label_mcmc_sweep(_base, rng, beta, total_dS, nattempts, nmoves);
        cur_S = entropy();

        // 3. Upper-level undirected moves
        for (size_t ui = 0; ui < _upper.size(); ++ui) {
            DenseUndirLevel& ulv = _upper[ui];
            if (ulv.N() <= 1 || (ui + 1 == _upper.size() && ulv.B() == 1)) continue;
            std::vector<int> uvs(ulv.N());
            std::iota(uvs.begin(), uvs.end(), 0);
            rng.shuffle(uvs);
            double cl = c * std::pow(2.0, static_cast<double>(ui + 1));
            for (int v : uvs) {
                ++nattempts;
                int r = ulv.group(v);
                int t = ulv.sample_group(v, cl, rng);
                if (t == r || ulv.size(r) == 1) continue;
                double pfwd = ulv.move_prob(v, t, cl);
                if (pfwd <= 0.0) continue;
                ulv.move(v, t);
                double prev = ulv.move_prob(v, r, cl);
                if (ui + 1 < _upper.size()) rebuild_upper_keep_bs();
                double new_S = entropy();
                double dS = new_S - cur_S;
                bool accept = std::isinf(beta)
                                  ? (dS < 0.0)
                                  : (prev > 0.0 && (dS <= 0.0 || rng.unif() < std::exp(-beta * dS + std::log(prev) - std::log(pfwd))));
                if (accept) {
                    cur_S = new_S;
                    total_dS += dS;
                    ++nmoves;
                } else {
                    _upper[ui].move(v, r);
                    if (ui + 1 < _upper.size()) rebuild_upper_keep_bs();
                }
            }
        }
    }

private:
    RankedState _base;
    std::vector<DenseUndirLevel> _upper;
};

std::vector<int> parse_rk_blocks(const Rcpp::IntegerVector& blocks, int N) {
    if (blocks.size() != N)
        Rcpp::stop("blocks must have length equal to number of vertices");
    std::vector<int> uniq;
    uniq.reserve(N);
    for (int i = 0; i < N; ++i) {
        int val = blocks[i];
        if (val == NA_INTEGER || val < 1)
            Rcpp::stop("blocks must contain positive integer labels");
        uniq.push_back(val);
    }
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    rk_hash_map<int, int> remap;
    for (size_t i = 0; i < uniq.size(); ++i)
        remap.emplace(uniq[i], static_cast<int>(i));
    std::vector<int> out(N);
    for (int i = 0; i < N; ++i) out[i] = remap[blocks[i]];
    return out;
}

std::vector<long> parse_rk_u(const Rcpp::NumericVector& u_r, int B_exp) {
    std::vector<long> out(std::max<int>(B_exp, static_cast<int>(u_r.size())), 0);
    for (int i = 0; i < static_cast<int>(u_r.size()); ++i) {
        if (!std::isfinite(u_r[i])) Rcpp::stop("u must contain finite integer ranks");
        out[i] = static_cast<long>(std::llround(u_r[i]));
    }
    return out;
}

}  // namespace

// [[Rcpp::export]]
double cpp_ranked_entropy(SEXP pointer, Rcpp::IntegerVector blocks,
                          Rcpp::NumericVector u, bool distance, bool unique,
                          double M_pos, bool degree_corrected, bool distributed,
                          bool include_label_dl, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    RankedState st(graph, distance, unique, static_cast<long>(M_pos),
                   degree_corrected, distributed, true);
    auto b = parse_rk_blocks(blocks, N);
    auto uv = parse_rk_u(u, N);
    st.set_partition(b, u.size() > 0 ? &uv : nullptr);
    return st.entropy(include_label_dl);
}

// [[Rcpp::export]]
Rcpp::List cpp_ranked_stats(SEXP pointer, Rcpp::IntegerVector blocks,
                            Rcpp::NumericVector u, bool distance, bool unique,
                            double M_pos, bool degree_corrected, bool distributed,
                            double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    RankedState st(graph, distance, unique, static_cast<long>(M_pos),
                   degree_corrected, distributed, true);
    auto b = parse_rk_blocks(blocks, N);
    auto uv = parse_rk_u(u, N);
    st.set_partition(b, u.size() > 0 ? &uv : nullptr);
    auto es = st.Es();
    auto [eb, eu] = st.external_blocks_and_u();
    return Rcpp::List::create(
        Rcpp::Named("upstream") = static_cast<double>(es[0]),
        Rcpp::Named("lateral") = static_cast<double>(es[1]),
        Rcpp::Named("downstream") = static_cast<double>(es[2]),
        Rcpp::Named("internal") = static_cast<double>(es[3]),
        Rcpp::Named("u") = eu);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_ranked_check_moves(SEXP pointer, Rcpp::IntegerVector blocks,
                                           Rcpp::NumericVector u, bool distance,
                                           bool unique, double M_pos,
                                           bool degree_corrected, bool distributed,
                                           int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    RankedState st(graph, distance, unique, static_cast<long>(M_pos),
                   degree_corrected, distributed, true);
    auto b = parse_rk_blocks(blocks, N);
    auto uv = parse_rk_u(u, N);
    st.set_partition(b, u.size() > 0 ? &uv : nullptr);
    RKRng rng;
    Rcpp::NumericVector err(nmoves);
    if (N <= 1) return err;
    double S = st.entropy();
    for (int i = 0; i < nmoves; ++i) {
        if (i % 3 == 0 && !unique && st.B() > 1) {
            // Check rank relabeling delta exactness too!
            int r = st.random_occupied(rng);
            long nu = static_cast<long>(rng.index(static_cast<size_t>(st.M())));
            double dS_v = st.relabel_group_dS(r, nu);
            double dS_a = st.relabel_group(r, nu);
            double S_new = st.entropy();
            err[i] = std::max(std::abs(dS_v - (S_new - S)), std::abs(dS_a - (S_new - S)));
            S = S_new;
        } else {
            int v = static_cast<int>(rng.index(N));
            int s = (rng.unif() < 0.2 && st.empty_group() >= 0)
                        ? st.empty_group()
                        : st.random_occupied(rng);
            double dS_v = st.virtual_move(v, s);
            double dS_a = st.move(v, s);
            double S_new = st.entropy();
            err[i] = std::max(std::abs(dS_v - (S_new - S)), std::abs(dS_a - (S_new - S)));
            S = S_new;
        }
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_ranked_minimize(SEXP pointer, bool distance, bool unique,
                               double M_pos, bool degree_corrected,
                               bool distributed, int B_min, double B_max,
                               bool refine, double epsilon, double q_max) {
    Graph& graph = get_graph(pointer);
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    RankedState st(graph, distance, unique, static_cast<long>(M_pos),
                   degree_corrected, distributed, true);
    RKRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    double S = ranked_minimize_flat(st, rng, static_cast<size_t>(B_min), bmax,
                                    refine, epsilon);
    auto [eb, eu] = st.external_blocks_and_u();
    return Rcpp::List::create(
        Rcpp::Named("blocks") = eb,
        Rcpp::Named("u") = eu,
        Rcpp::Named("entropy") = S,
        Rcpp::Named("B") = st.B());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_ranked_entropy(SEXP pointer, Rcpp::List bs,
                                              Rcpp::NumericVector u, bool distance,
                                              bool unique, double M_pos,
                                              bool degree_corrected, bool distributed,
                                              double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    NestedRankedSBM nr(graph, distance, unique, static_cast<long>(M_pos),
                       degree_corrected, distributed);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_rk_blocks(x, cur_n);
        cur_n = *std::max_element(bl.begin(), bl.end()) + 1;
        b_vec.push_back(std::move(bl));
    }
    auto uv = parse_rk_u(u, N);
    nr.build(b_vec, u.size() > 0 ? &uv : nullptr);
    auto lvS = nr.level_entropies();
    return Rcpp::NumericVector(lvS.begin(), lvS.end());
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_ranked_minimize(SEXP pointer, bool distance, bool unique,
                                      double M_pos, bool degree_corrected,
                                      bool distributed, int B_min, double B_max,
                                      bool refine, double epsilon, double q_max) {
    Graph& graph = get_graph(pointer);
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    // Use is_top = true on lv0 so base search respects proper B regularization
    RankedState lv0(graph, distance, unique, static_cast<long>(M_pos),
                    degree_corrected, distributed, true);
    RKRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    std::vector<std::pair<std::vector<int>, std::vector<long>>> cps;
    ranked_minimize_flat(lv0, rng, static_cast<size_t>(B_min), bmax, refine, epsilon, &cps);

    NestedRankedSBM best_nr(graph, distance, unique, static_cast<long>(M_pos),
                            degree_corrected, distributed);
    double best_S = RK_INF;
    for (const auto& [b0, u0] : cps) {
        NestedRankedSBM cand(graph, distance, unique, static_cast<long>(M_pos),
                             degree_corrected, distributed);
        cand.optimize_from_base(b0, u0, rng, epsilon);
        double s_cand = cand.entropy();
        if (s_cand < best_S) {
            best_S = s_cand;
            best_nr = std::move(cand);
        }
    }
    if (refine)
        best_nr.refine_joint(rng, static_cast<size_t>(B_min), bmax, epsilon);

    auto lvS = best_nr.level_entropies();
    return Rcpp::List::create(
        Rcpp::Named("bs") = best_nr.external_bs_r(),
        Rcpp::Named("u") = best_nr.external_u_r(),
        Rcpp::Named("level_entropy") = lvS,
        Rcpp::Named("entropy") = std::accumulate(lvS.begin(), lvS.end(), 0.0));
}

// [[Rcpp::export]]
Rcpp::List cpp_ranked_mcmc(SEXP pointer, Rcpp::IntegerVector blocks,
                           Rcpp::NumericVector u, bool distance, bool unique,
                           double M_pos, bool degree_corrected, bool distributed,
                           Rcpp::List mcmc_opts, Rcpp::List eq_opts, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    RankedState st(graph, distance, unique, static_cast<long>(M_pos),
                   degree_corrected, distributed, true);
    auto b = parse_rk_blocks(blocks, N);
    auto uv = parse_rk_u(u, N);
    st.set_partition(b, u.size() > 0 ? &uv : nullptr);
    RKRng rng;

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
            ranked_flat_mcmc_sweep(st, rng, beta, c, d, multiflip, B_min, B_max,
                                   total_dS, nattempts, nmoves);
            ++steps_done;
            Rcpp::checkUserInterrupt();
        }
        auto [eb, eu] = st.external_blocks_and_u();
        return Rcpp::List::create(
            Rcpp::Named("blocks") = eb,
            Rcpp::Named("u") = eu,
            Rcpp::Named("entropy") = st.entropy(),
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
    std::vector<Rcpp::IntegerVector> samples;
    std::vector<double> sample_entropies;

    auto run_eq_loop = [&](double cur_beta, int max_it, int f_it) {
        double best_S = st.entropy();
        double worst_S = best_S;
        int no_break = 0;
        int breaks_left = nbreaks;
        for (int step = 1; step <= max_it; ++step) {
            for (int sub = 0; sub < std::max(1, niter); ++sub) {
                ranked_flat_mcmc_sweep(st, rng, cur_beta, c, d, multiflip,
                                       B_min, B_max, total_dS, nattempts, nmoves);
            }
            double S_now = st.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_now);
            if (collect && (step % thin == 0)) {
                samples.push_back(st.external_blocks_and_u().first);
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
        double beta_lo = Rcpp::as<double>(eq_opts["beta_start"]);
        double beta_hi = Rcpp::as<double>(eq_opts["beta_end"]);
        int nsteps = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["anneal_niter"])));
        for (int i = 0; i < nsteps; ++i) {
            double frac = (nsteps == 1) ? 1.0 : static_cast<double>(i) / (nsteps - 1);
            double cur_beta = beta_lo * std::pow(beta_hi / beta_lo, frac);
            run_eq_loop(cur_beta, max_niter, force_niter);
        }
    } else {
        run_eq_loop(beta, max_niter, force_niter);
    }

    auto [eb, eu] = st.external_blocks_and_u();
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("blocks") = eb,
        Rcpp::Named("u") = eu,
        Rcpp::Named("entropy") = st.entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("niter") = steps_done);
    if (keep_hist) out["history"] = Rcpp::NumericVector(hist_S.begin(), hist_S.end());
    if (collect) {
        Rcpp::List s_list(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) s_list[i] = samples[i];
        out["samples"] = s_list;
        out["sample_entropies"] = Rcpp::NumericVector(sample_entropies.begin(), sample_entropies.end());
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_ranked_mcmc(SEXP pointer, Rcpp::List bs, Rcpp::NumericVector u,
                                  bool distance, bool unique, double M_pos,
                                  bool degree_corrected, bool distributed,
                                  Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                                  double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) rklog_q().init(rk_q_cache_max(q_max));
    NestedRankedSBM nr(graph, distance, unique, static_cast<long>(M_pos),
                       degree_corrected, distributed);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_rk_blocks(x, cur_n);
        cur_n = *std::max_element(bl.begin(), bl.end()) + 1;
        b_vec.push_back(std::move(bl));
    }
    auto uv = parse_rk_u(u, N);
    nr.build(b_vec, u.size() > 0 ? &uv : nullptr);
    RKRng rng;

    double beta = Rcpp::as<double>(mcmc_opts["beta"]);
    double c = Rcpp::as<double>(mcmc_opts["c"]);
    double d = Rcpp::as<double>(mcmc_opts["d"]);
    int niter = static_cast<int>(Rcpp::as<double>(mcmc_opts["niter"]));
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
            nr.mcmc_sweep(rng, beta, c, d, B_min, B_max, total_dS, nattempts, nmoves);
            ++steps_done;
            Rcpp::checkUserInterrupt();
        }
        return Rcpp::List::create(
            Rcpp::Named("bs") = nr.external_bs_r(),
            Rcpp::Named("u") = nr.external_u_r(),
            Rcpp::Named("entropy") = nr.entropy(),
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

    auto run_eq_loop = [&](double cur_beta, int max_it, int f_it) {
        double best_S = nr.entropy();
        double worst_S = best_S;
        int no_break = 0;
        int breaks_left = nbreaks;
        for (int step = 1; step <= max_it; ++step) {
            for (int sub = 0; sub < std::max(1, niter); ++sub) {
                nr.mcmc_sweep(rng, cur_beta, c, d, B_min, B_max, total_dS, nattempts, nmoves);
            }
            double S_now = nr.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_now);
            if (collect && (step % thin == 0)) {
                samples.push_back(nr.external_bs_r());
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
        double beta_lo = Rcpp::as<double>(eq_opts["beta_start"]);
        double beta_hi = Rcpp::as<double>(eq_opts["beta_end"]);
        int nsteps = std::max(1, static_cast<int>(Rcpp::as<double>(eq_opts["anneal_niter"])));
        for (int i = 0; i < nsteps; ++i) {
            double frac = (nsteps == 1) ? 1.0 : static_cast<double>(i) / (nsteps - 1);
            double cur_beta = beta_lo * std::pow(beta_hi / beta_lo, frac);
            run_eq_loop(cur_beta, max_niter, force_niter);
        }
    } else {
        run_eq_loop(beta, max_niter, force_niter);
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("bs") = nr.external_bs_r(),
        Rcpp::Named("u") = nr.external_u_r(),
        Rcpp::Named("entropy") = nr.entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves),
        Rcpp::Named("niter") = steps_done);
    if (keep_hist) out["history"] = Rcpp::NumericVector(hist_S.begin(), hist_S.end());
    if (collect) {
        Rcpp::List s_list(samples.size());
        for (size_t i = 0; i < samples.size(); ++i) s_list[i] = samples[i];
        out["samples"] = s_list;
        out["sample_entropies"] = Rcpp::NumericVector(sample_entropies.begin(), sample_entropies.end());
    }
    return out;
}
