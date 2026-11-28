// Native C++ flat and nested layered stochastic block model engine
// matching graph-tool 3.9's LayeredBlockState (both independent=FALSE
// edge-covariate model and independent=TRUE conditionally independent layers)
// and NestedBlockState(base_state=LayeredBlockState) to machine precision.
//
// References:
//   Peixoto, T. P. (2015). Inferring the mesoscale structure of layered,
//   edge-valued, and time-varying networks. Physical Review E 92, 042807.
//
// Adapted from graph-tool (https://graph-tool.skewed.de), LGPL-3.0-or-later.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph_types.h"
#include "parallel.h"
#include "weight_cov.h"
#include "logq.h"

namespace {

using namespace rgt_cov;

constexpr double L_INF = std::numeric_limits<double>::infinity();
constexpr double L_LOG2 = 0.693147180559945309417232121458176568;

template <class K, class V>
using lt_hash_map = std::unordered_map<K, V>;
template <class K>
using lt_hash_set = std::unordered_set<K>;

// ---------------------------------------------------------------------------
// Fast math + Szekeres log_q cache (exact match to blockmodel.cpp / graph-tool 3.9)

struct LMath {
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

LMath& lmath() {
    static LMath m;
    return m;
}
// Shared log q(n, k) table (src/logq.h).
using LLogQ = rgt_logq::LogQ;
inline LLogQ& llog_q() { return rgt_logq::shared(); }

size_t l_q_cache_max(double q_max) {
    if (!std::isfinite(q_max) || q_max <= 0.0) return 10000;
    return std::min<size_t>(static_cast<size_t>(q_max), 10000);
}

// RNG backed by R's RNG state
struct LRng {
    inline double unif() { return R::runif(0.0, 1.0); }
    inline size_t index(size_t n) {
        size_t i = static_cast<size_t>(std::floor(unif() * static_cast<double>(n)));
        return (i >= n) ? (n - 1) : i;
    }
    template <class Vec>
    void shuffle(Vec& v) {
        for (size_t i = v.size(); i > 1; --i) {
            std::swap(v[i - 1], v[index(i)]);
        }
    }
};

// ---------------------------------------------------------------------------
// Sparse sorted layer count map (small vector of (layer, count))

struct LayerCountVec {
    std::vector<std::pair<int, long>> items;  // sorted by layer

    bool empty() const { return items.empty(); }
    size_t size() const { return items.size(); }
    void clear() { items.clear(); }

    long get(int l) const {
        for (const auto& kv : items) {
            if (kv.first == l) return kv.second;
            if (kv.first > l) break;
        }
        return 0;
    }

    void add(int l, long d) {
        if (d == 0) return;
        for (auto it = items.begin(); it != items.end(); ++it) {
            if (it->first == l) {
                it->second += d;
                if (it->second == 0) items.erase(it);
                return;
            }
            if (it->first > l) {
                items.insert(it, {l, d});
                return;
            }
        }
        items.push_back({l, d});
    }

    void add_all(const LayerCountVec& other, long sign = 1) {
        for (const auto& [l, c] : other.items) add(l, sign * c);
    }
};

// Sparse sorted per-layer covariate accumulators
struct LayerWItem {
    int l;
    long m_cov = 0;
    WVec w;
};

struct LayerWVec {
    std::vector<LayerWItem> items;  // sorted by l

    bool empty() const { return items.empty(); }
    void clear() { items.clear(); }

    const LayerWItem* find_const(int l) const {
        for (const auto& it : items) {
            if (it.l == l) return &it;
            if (it.l > l) break;
        }
        return nullptr;
    }

    void add_cov(int l, long dm_cov, const WVec& dw, double sign, int K,
                 bool erase_zero_m = false) {
        for (auto it = items.begin(); it != items.end(); ++it) {
            if (it->l == l) {
                it->m_cov += dm_cov;
                for (int k = 0; k < K; ++k) {
                    it->w.x(k) += sign * dw.x(k);
                    it->w.x2(k) += sign * dw.x2(k);
                }
                if (erase_zero_m && it->m_cov <= 0) items.erase(it);
                return;
            }
            if (it->l > l) {
                if (erase_zero_m && dm_cov <= 0) return;
                LayerWItem nw;
                nw.l = l;
                nw.m_cov = dm_cov;
                nw.w.zero(K);
                for (int k = 0; k < K; ++k) {
                    nw.w.x(k) = sign * dw.x(k);
                    nw.w.x2(k) = sign * dw.x2(k);
                }
                items.insert(it, std::move(nw));
                return;
            }
        }
        if (erase_zero_m && dm_cov <= 0) return;
        LayerWItem nw;
        nw.l = l;
        nw.m_cov = dm_cov;
        nw.w.zero(K);
        for (int k = 0; k < K; ++k) {
            nw.w.x(k) = sign * dw.x(k);
            nw.w.x2(k) = sign * dw.x2(k);
        }
        items.push_back(std::move(nw));
    }
};

struct LPairStat {
    long m_tot = 0;
    long m_cov = 0;      // used when _K > 0 && !_independent
    LayerCountVec ml;
    WVec w;              // used when _K > 0 && !_independent
    LayerWVec wl;        // used when _K > 0 && _independent
};

struct LGraphEdge {
    int u, v;
    long m_tot;          // sum of ml counts (1 at level 0)
    int single_layer;    // >= 0 at level 0 (fast path!), -1 at upper levels
    LayerCountVec ml;    // used when single_layer < 0
    long m_cov = 1;      // used when K > 0 && (single_layer >= 0 || !independent)
    WVec w;              // used when K > 0 && (single_layer >= 0 || !independent)
    LayerWVec wl;        // used when K > 0 && independent && single_layer < 0
};

struct LAdjRef {
    int nbr;
    int eidx;
    bool out;
};

inline uint64_t lpack_pair(int a, int b) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(b));
}

inline uint64_t lcanon_pair(int a, int b, bool directed) {
    if (!directed && a > b) std::swap(a, b);
    return lpack_pair(a, b);
}

inline std::pair<int, int> lunpack_pair(uint64_t key) {
    return {static_cast<int>(static_cast<uint32_t>(key >> 32)),
            static_cast<int>(static_cast<uint32_t>(key & 0xffffffffULL))};
}

// ---------------------------------------------------------------------------
// One level of a LayeredSBM (base level l=0 or upper dense level l>=1)

class LayeredLevel {
public:
    LayeredLevel(int N, int C, bool directed, bool independent, bool dc,
                 bool dist, bool is_base, bool is_top,
                 std::vector<CovSpec> covs = {})
        : _N(N), _C(C), _directed(directed), _independent(independent),
          _dc(dc), _dist(dist), _is_base(is_base), _is_top(is_top),
          _covs(std::move(covs)),
          _K(static_cast<int>(_covs.size())),
          _num_cov_strata((_K == 0) ? 0 : (independent ? C : 1)),
          _cap(std::max(2, N)),
          _adj(N), _kout(N, 0), _kin(directed ? N : 0, 0),
          _v_layers(independent ? N : 0),
          _v_kout_l(independent && is_base ? N : 0),
          _v_kin_l(independent && is_base && directed ? N : 0),
          _E_l(C, 0),
          _b(N, 0), _wr(_cap, 0), _ep(_cap, 0), _em(directed ? _cap : 0, 0),
          _wr_l(independent ? _cap : 0),
          _ep_l(independent && is_base ? _cap : 0),
          _em_l(independent && is_base && directed ? _cap : 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(N, 0),
          _m_out_tot(_cap), _m_in_tot(directed ? _cap : 0),
          _h_out(!independent && is_base && dc && dist ? _cap : 0),
          _h_in(!independent && is_base && dc && dist && directed ? _cap : 0),
          _h_out_l(independent && is_base && dc && dist ? C : 0),
          _h_in_l(independent && is_base && dc && dist && directed ? C : 0),
          _B_E_c(_num_cov_strata, 0),
          _B_E_D_c(_num_cov_strata, 0),
          _recsum_c(_num_cov_strata, std::vector<double>(_K, 0.0)),
          _recx2_c(_num_cov_strata, std::vector<double>(_K, 0.0)),
          _recdx_c(_num_cov_strata, std::vector<double>(_K, 0.0)) {
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        if (independent && is_base && dc && dist) {
            for (int l = 0; l < C; ++l) {
                _h_out_l[l].resize(_cap);
                if (directed) _h_in_l[l].resize(_cap);
            }
        }
        for (int k = 0; k < _K; ++k) {
            if (_covs[k].type == CovType::REAL_NORMAL) _normal_idx.push_back(k);
        }
    }

    // Set edges for level 0 (or upper levels; for upper levels with independent=TRUE,
    // also pass per-vertex active layers from the lower level's occupied groups!)
    void set_edges(std::vector<LGraphEdge> edges,
                   const std::vector<std::vector<int>>* upper_v_layers = nullptr) {
        _edges = std::move(edges);
        _E = 0;
        _E_l.assign(_C, 0);
        _adj.assign(_N, {});
        _kout.assign(_N, 0);
        if (_directed) _kin.assign(_N, 0);
        if (_independent) {
            _v_layers.assign(_N, {});
            if (_is_base) {
                _v_kout_l.assign(_N, {});
                if (_directed) _v_kin_l.assign(_N, {});
            }
        }
        for (int st = 0; st < _num_cov_strata; ++st) _recsum_c[st].assign(_K, 0.0);

        for (int ei = 0; ei < static_cast<int>(_edges.size()); ++ei) {
            const auto& e = _edges[ei];
            _E += e.m_tot;
            if (e.single_layer >= 0) {
                _E_l[e.single_layer] += e.m_tot;
            } else {
                for (const auto& [l, c] : e.ml.items) _E_l[l] += c;
            }
            if (_K > 0) {
                if (e.single_layer >= 0) {
                    int st = _independent ? e.single_layer : 0;
                    for (int k = 0; k < _K; ++k) _recsum_c[st][k] += e.w.x(k);
                } else if (!_independent) {
                    for (int k = 0; k < _K; ++k) _recsum_c[0][k] += e.w.x(k);
                } else {
                    for (const auto& wi : e.wl.items) {
                        for (int k = 0; k < _K; ++k) _recsum_c[wi.l][k] += wi.w.x(k);
                    }
                }
            }
            if (_directed) {
                _adj[e.u].push_back({e.v, ei, true});
                _adj[e.v].push_back({e.u, ei, false});
                _kout[e.u] += e.m_tot;
                _kin[e.v] += e.m_tot;
                if (_independent && _is_base) {
                    if (e.single_layer >= 0) {
                        _v_kout_l[e.u].add(e.single_layer, e.m_tot);
                        _v_kin_l[e.v].add(e.single_layer, e.m_tot);
                    } else {
                        _v_kout_l[e.u].add_all(e.ml, +1);
                        _v_kin_l[e.v].add_all(e.ml, +1);
                    }
                }
            } else {
                _adj[e.u].push_back({e.v, ei, true});
                _kout[e.u] += e.m_tot;
                if (e.u != e.v) {
                    _adj[e.v].push_back({e.u, ei, true});
                    _kout[e.v] += e.m_tot;
                } else {
                    _kout[e.u] += e.m_tot;
                }
                if (_independent && _is_base) {
                    long mult = (e.u == e.v) ? 2 : 1;
                    if (e.single_layer >= 0) {
                        _v_kout_l[e.u].add(e.single_layer, mult * e.m_tot);
                        if (e.u != e.v) _v_kout_l[e.v].add(e.single_layer, e.m_tot);
                    } else {
                        _v_kout_l[e.u].add_all(e.ml, mult);
                        if (e.u != e.v) _v_kout_l[e.v].add_all(e.ml, +1);
                    }
                }
            }
        }

        // Count active layers
        _C_active = 0;
        for (int l = 0; l < _C; ++l)
            if (_E_l[l] > 0) _C_active++;

        // Populate _v_layers for independent = TRUE
        if (_independent) {
            if (_is_base) {
                for (int v = 0; v < _N; ++v) {
                    std::vector<int> ls;
                    for (const auto& [l, c] : _v_kout_l[v].items)
                        if (c > 0) ls.push_back(l);
                    if (_directed) {
                        for (const auto& [l, c] : _v_kin_l[v].items)
                            if (c > 0) ls.push_back(l);
                        std::sort(ls.begin(), ls.end());
                        ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
                    }
                    _v_layers[v] = std::move(ls);
                }
            } else if (upper_v_layers != nullptr) {
                _v_layers = *upper_v_layers;
            } else {
                // Fallback: derive from incident edges
                for (int v = 0; v < _N; ++v) {
                    std::vector<int> ls;
                    for (const auto& a : _adj[v]) {
                        const auto& e = _edges[a.eidx];
                        if (e.single_layer >= 0) ls.push_back(e.single_layer);
                        else for (const auto& [l, c] : e.ml.items) if (c > 0) ls.push_back(l);
                    }
                    std::sort(ls.begin(), ls.end());
                    ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
                    _v_layers[v] = std::move(ls);
                }
            }
        }

        // Base-level constant likelihood factor (per-layer parallel edges & degrees)
        _base_const = 0.0;
        const LMath& M = lmath();
        lmath().ensure(static_cast<size_t>(_N + _E + 2 * _C + 64));
        if (_is_base && _E > 0) {
            // Parallel edges within each layer (u, v, l)
            // Pack ((u, v), l) via map from canon_pair(u, v) -> LayerCountVec
            lt_hash_map<uint64_t, LayerCountVec> emult;
            emult.reserve(_edges.size());
            for (const auto& e : _edges) {
                uint64_t pk = lcanon_pair(e.u, e.v, _directed);
                if (e.single_layer >= 0) emult[pk].add(e.single_layer, e.m_tot);
                else emult[pk].add_all(e.ml, +1);
            }
            for (const auto& [pk, lcv] : emult) {
                auto [u, v] = lunpack_pair(pk);
                for (const auto& [l, mult] : lcv.items) {
                    if (mult <= 0) continue;
                    _base_const += M.lgam1(mult);
                    if (!_directed && u == v)
                        _base_const += static_cast<double>(mult) * L_LOG2;
                }
            }
            if (_dc) {
                if (!_independent) {
                    for (int v = 0; v < _N; ++v) {
                        _base_const -= M.lgam1(_kout[v]);
                        if (_directed) _base_const -= M.lgam1(_kin[v]);
                    }
                } else {
                    for (int v = 0; v < _N; ++v) {
                        for (const auto& [l, k] : _v_kout_l[v].items)
                            _base_const -= M.lgam1(k);
                        if (_directed) {
                            for (const auto& [l, k] : _v_kin_l[v].items)
                                _base_const -= M.lgam1(k);
                        }
                    }
                }
            }
            for (int k = 0; k < _K; ++k) _base_const += _covs[k].const_term;
        }
        if (_is_base && _dc && _dist) {
            llog_q().init(rgt_logq::table_size(static_cast<double>(_E), 5000));
        }
    }

    int N() const { return _N; }
    int C() const { return _C; }
    int C_active() const { return _C_active; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    int K() const { return _K; }
    int num_cov_strata() const { return _num_cov_strata; }
    bool directed() const { return _directed; }
    bool independent() const { return _independent; }
    bool is_base() const { return _is_base; }
    bool is_top() const { return _is_top; }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    long B_E_D_c(int st) const { return _B_E_D_c[st]; }
    double recdx_c(int st, int k) const { return _recdx_c[st][k]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    const std::vector<CovSpec>& covs() const { return _covs; }
    const std::vector<int>& normal_idx() const { return _normal_idx; }
    const lt_hash_map<uint64_t, LPairStat>& pairs() const { return _pairs; }
    const LayerCountVec& wr_l(int r) const { return _wr_l[r]; }
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
            if (_independent) {
                _wr_l[r].clear();
                if (_is_base) {
                    _ep_l[r].clear();
                    if (_directed) _em_l[r].clear();
                }
            }
            _members[r].clear();
            _occ_pos[r] = -1;
            _m_out_tot[r].clear();
            if (_directed) _m_in_tot[r].clear();
            if (!_independent && _is_base && _dc && _dist) {
                _h_out[r].clear();
                if (_directed) _h_in[r].clear();
            }
            if (_independent && _is_base && _dc && _dist) {
                for (int l = 0; l < _C; ++l) {
                    _h_out_l[l][r].clear();
                    if (_directed) _h_in_l[l][r].clear();
                }
            }
        }
        _occupied.clear();
        _empty.clear();
        _empty_pos.assign(_cap, -1);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        _pairs.clear();

        for (int v = 0; v < _N; ++v) {
            int r = b[v];
            if (r < 0 || r >= _cap)
                Rcpp::stop("Internal error: invalid group index");
            _b[v] = r;
            add_member(v, r);
            if (_independent) {
                for (int l : _v_layers[v]) _wr_l[r].add(l, 1);
            }
            if (_is_base && _dc && _dist) {
                if (!_independent) {
                    _h_out[r][_kout[v]]++;
                    if (_directed) _h_in[r][_kin[v]]++;
                } else {
                    for (int l : _v_layers[v]) {
                        _h_out_l[l][r][_v_kout_l[v].get(l)]++;
                        if (_directed) _h_in_l[l][r][_v_kin_l[v].get(l)]++;
                    }
                }
            }
        }

        for (const auto& e : _edges) {
            int r = _b[e.u], s = _b[e.v];
            uint64_t pk = lcanon_pair(r, s, _directed);
            LPairStat& ps = _pairs[pk];
            ps.m_tot += e.m_tot;
            if (e.single_layer >= 0) ps.ml.add(e.single_layer, e.m_tot);
            else ps.ml.add_all(e.ml, +1);

            if (_K > 0) {
                if (!_independent) {
                    ps.m_cov += e.m_cov;
                    if (ps.w.k == 0) ps.w.zero(_K);
                    for (int k = 0; k < _K; ++k) {
                        ps.w.x(k) += e.w.x(k);
                        ps.w.x2(k) += e.w.x2(k);
                    }
                } else if (e.single_layer >= 0) {
                    ps.wl.add_cov(e.single_layer, e.m_cov, e.w, +1.0, _K);
                } else {
                    for (const auto& wi : e.wl.items) {
                        ps.wl.add_cov(wi.l, wi.m_cov, wi.w, +1.0, _K);
                    }
                }
            }

            if (_directed) {
                _ep[r] += e.m_tot;
                _em[s] += e.m_tot;
                _m_out_tot[r][s] += e.m_tot;
                _m_in_tot[s][r] += e.m_tot;
                if (_independent && _is_base) {
                    if (e.single_layer >= 0) {
                        _ep_l[r].add(e.single_layer, e.m_tot);
                        _em_l[s].add(e.single_layer, e.m_tot);
                    } else {
                        _ep_l[r].add_all(e.ml, +1);
                        _em_l[s].add_all(e.ml, +1);
                    }
                }
            } else {
                _ep[r] += e.m_tot;
                _ep[s] += e.m_tot;
                if (r == s) {
                    _m_out_tot[r][r] += 2 * e.m_tot;
                } else {
                    _m_out_tot[r][s] += e.m_tot;
                    _m_out_tot[s][r] += e.m_tot;
                }
                if (_independent && _is_base) {
                    if (e.single_layer >= 0) {
                        _ep_l[r].add(e.single_layer, e.m_tot);
                        _ep_l[s].add(e.single_layer, e.m_tot);
                    } else {
                        _ep_l[r].add_all(e.ml, +1);
                        _ep_l[s].add_all(e.ml, +1);
                    }
                }
            }
        }
        refresh_cov_globals();
    }

    void refresh_cov_globals() {
        if (_K == 0) return;
        for (int st = 0; st < _num_cov_strata; ++st) {
            _B_E_c[st] = 0;
            _B_E_D_c[st] = 0;
            _recx2_c[st].assign(_K, 0.0);
            _recdx_c[st].assign(_K, 0.0);
        }
        for (const auto& [pk, ps] : _pairs) {
            (void)pk;
            if (!_independent) {
                if (ps.m_cov > 0) {
                    ++_B_E_c[0];
                    if (ps.m_cov > 1) ++_B_E_D_c[0];
                    for (int k : _normal_idx) {
                        double x = ps.w.x(k);
                        _recx2_c[0][k] += x * x;
                        if (ps.m_cov > 1) {
                            double dx = ps.w.x2(k) - x * (x / static_cast<double>(ps.m_cov));
                            if (dx > 0.0) _recdx_c[0][k] += dx;
                        }
                    }
                }
            } else {
                for (const auto& wi : ps.wl.items) {
                    if (wi.m_cov <= 0) continue;
                    int st = wi.l;
                    ++_B_E_c[st];
                    if (wi.m_cov > 1) ++_B_E_D_c[st];
                    for (int k : _normal_idx) {
                        double x = wi.w.x(k);
                        _recx2_c[st][k] += x * x;
                        if (wi.m_cov > 1) {
                            double dx = wi.w.x2(k) - x * (x / static_cast<double>(wi.m_cov));
                            if (dx > 0.0) _recdx_c[st][k] += dx;
                        }
                    }
                }
            }
        }
    }

    // Exact graph-tool 3.9 lmemb term (only at base level l=0 when independent=TRUE)
    double lmemb_dl() const {
        if (!_is_base || !_independent || _C_active <= 1) return 0.0;
        const LMath& M = lmath();
        double x = M.lbinom(_C_active, 1);
        if (x > 15.0)
            return static_cast<double>(_C_active) * x - M.lgam1(_C_active);
        return M.lbinom(2 * _C_active - 1, _C_active);
    }

    // Exact graph-tool 3.9 edges_dl term (at top level)
    double edges_dl() const {
        if (!_is_top || _E == 0) return 0.0;
        const LMath& M = lmath();
        long Bn = B();
        int64_t P = _directed ? int64_t(Bn) * Bn : int64_t(Bn) * (Bn + 1) / 2;
        double S = M.lbinom(_E - 1, _C_active - 1);
        for (int l = 0; l < _C; ++l) {
            if (_E_l[l] > 0) S += M.lbinom(P + _E_l[l] - 1, _E_l[l]);
        }
        return S;
    }

    // Full exact description length of this level
    double entropy(bool include_lmemb = true) const {
        const LMath& M = lmath();
        long Bn = B();
        // Partition DL (exact graph-tool 3.9: log(N) + lbinom(N-1, B-1) + lgamma(N+1) - sum lgamma(wr+1))
        double S = 0.0;
        if (_N > 0) {
            S += M.logn(_N) + M.lbinom(_N - 1, Bn - 1) + M.lgam1(_N);
            for (int r : _occupied) S -= M.lgam1(_wr[r]);
        }
        S += edges_dl();

        if (_is_base) {
            S += _base_const;
            if (include_lmemb) S += lmemb_dl();

            // vterm + degree sequence DL
            for (int r : _occupied) {
                if (!_independent) {
                    if (_dc) {
                        S += M.lgam1(_ep[r]);
                        if (_directed) S += M.lgam1(_em[r]);
                        S += deg_dl_group_flat(r);
                    } else {
                        long dr = _ep[r] + (_directed ? _em[r] : 0);
                        S += static_cast<double>(dr) * M.logn(_wr[r]);
                    }
                } else {
                    for (const auto& [l, nrl] : _wr_l[r].items) {
                        if (nrl <= 0) continue;
                        long epl = _ep_l[r].get(l);
                        long eml = _directed ? _em_l[r].get(l) : 0;
                        if (_dc) {
                            S += M.lgam1(epl);
                            if (_directed) S += M.lgam1(eml);
                            S += deg_dl_group_layer(l, r, nrl, epl, eml);
                        } else {
                            S += static_cast<double>(epl + eml) * M.logn(nrl);
                        }
                    }
                }
            }

            // Base-level eterm (identical for independent=FALSE and independent=TRUE!)
            for (const auto& [pk, ps] : _pairs) {
                if (ps.m_tot <= 0) continue;
                auto [r, s] = lunpack_pair(pk);
                S += base_eterm(r, s, ps);
            }
        } else {
            // Dense upper level (l >= 1)
            for (const auto& [pk, ps] : _pairs) {
                if (ps.m_tot <= 0) continue;
                auto [r, s] = lunpack_pair(pk);
                if (!_independent) {
                    for (const auto& [l, m] : ps.ml.items) {
                        if (m > 0) S += dense_eterm(r, s, m, _wr[r], _wr[s]);
                    }
                } else {
                    for (const auto& [l, m] : ps.ml.items) {
                        if (m > 0)
                            S += dense_eterm(r, s, m, _wr_l[r].get(l), _wr_l[s].get(l));
                    }
                }
            }
        }
        return S + cov_entropy(true, true);
    }

    double cov_entropy(bool include_rec = true, bool include_dl = true) const {
        if (_K == 0 || _E == 0) return 0.0;
        double S = 0.0;
        if (include_rec) {
            for (const auto& [pk, ps] : _pairs) {
                (void)pk;
                if (!_independent) {
                    if (ps.m_cov <= 0) continue;
                    double m = static_cast<double>(ps.m_cov);
                    for (int k = 0; k < _K; ++k) {
                        S += cov_pair_neg_log_P(_covs[k], m, ps.w.x(k), ps.w.x2(k));
                    }
                } else {
                    for (const auto& wi : ps.wl.items) {
                        if (wi.m_cov <= 0) continue;
                        double m = static_cast<double>(wi.m_cov);
                        for (int k = 0; k < _K; ++k) {
                            S += cov_pair_neg_log_P(_covs[k], m, wi.w.x(k), wi.w.x2(k));
                        }
                    }
                }
            }
        }
        for (int st = 0; st < _num_cov_strata; ++st) {
            for (int k = 0; k < _K; ++k) {
                const CovSpec& c = _covs[k];
                if (c.type == CovType::REAL_NORMAL) {
                    if (c.micro) {
                        if (include_dl && _is_top) {
                            S += -signed_w_log_P(static_cast<double>(_B_E_c[st]), _recsum_c[st][k],
                                                 _recx2_c[st][k], c.p0, c.p1, c.p2, c.p3, c.eps);
                        }
                        S += -positive_w_log_P(static_cast<double>(_B_E_D_c[st]), _recdx_c[st][k],
                                               c.p2, c.p3, c.eps);
                    }
                } else if (include_dl && _is_top && c.micro) {
                    if (c.type == CovType::REAL_EXPONENTIAL) {
                        S += -positive_w_log_P(static_cast<double>(_B_E_c[st]), _recsum_c[st][k],
                                               c.p0, c.p1, c.eps);
                    } else {
                        S += -geometric_w_log_P(static_cast<double>(_B_E_c[st]), _recsum_c[st][k],
                                                NAN, NAN);
                    }
                }
            }
        }
        return S;
    }

    struct PairDelta {
        uint64_t key = 0;
        long dm_tot = 0;
        long dm_cov = 0;
        LayerCountVec dml;
        WVec dw;
        LayerWVec dwl;
    };

    struct MoveBuf {
        std::vector<PairDelta> pdeltas;
        lt_hash_map<uint64_t, size_t> scratch_idx;
        long d_ep_r = 0, d_ep_s = 0, d_em_r = 0, d_em_s = 0;
        LayerCountVec d_ep_l_r, d_ep_l_s, d_em_l_r, d_em_l_s;
        std::vector<long> d_B_E_c, d_B_E_D_c;
        std::vector<std::vector<double>> d_recx2_c, d_recdx_c;
        std::vector<uint8_t> touched_strata;
    };

    // Exact O(deg(v)) incremental virtual_move(v, s)
    double virtual_move(int v, int s) const {
        int r = _b[v];
        if (r == s) return 0.0;
        return compute_move_delta(v, r, s, &_scratch_buf);
    }

    // Thread-safe read-only virtual_move(v, s) using caller-supplied MoveBuf
    double virtual_move_const(int v, int s, MoveBuf& buf) const {
        int r = _b[v];
        if (r == s) return 0.0;
        return compute_move_delta(v, r, s, &buf);
    }

    // Exact O(deg(v)) incremental move(v, s)
    double move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double dS = compute_move_delta(v, r, s, &_scratch_buf);
        apply_move(v, r, s, _scratch_buf);
        return dS;
    }

    // Apply the move most recently evaluated by virtual_move(v, s)
    void apply_last_virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s) return;
        apply_move(v, r, s, _scratch_buf);
    }

    // Direct group-level virtual_merge(r, t)
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

    // Peixoto neighbor-based ergodic proposal (using collapsed adjacency counts)
    template <class Rng>
    int sample_group(int v, double c, Rng& rng) const {
        if (_occupied.empty()) return 0;
        const auto& adj = _adj[v];
        if (adj.empty() || std::isinf(c))
            return _occupied[rng.index(_occupied.size())];
        const auto& a = adj[rng.index(adj.size())];
        int t = _b[a.nbr];
        long kt = _ep[t] + (_directed ? _em[t] : 0);
        double B_act = static_cast<double>(_occupied.size());
        double p_rand = (c * B_act) / (static_cast<double>(kt) + c * B_act);
        if (rng.unif() < p_rand)
            return _occupied[rng.index(_occupied.size())];
        double u = rng.unif() * static_cast<double>(kt);
        double acc = 0.0;
        for (const auto& [s, m] : _m_out_tot[t]) {
            acc += static_cast<double>(m);
            if (u <= acc) return s;
        }
        if (_directed) {
            for (const auto& [s, m] : _m_in_tot[t]) {
                acc += static_cast<double>(m);
                if (u <= acc) return s;
            }
        }
        return _occupied[rng.index(_occupied.size())];
    }

    double move_prob(int v, int s, double c) const {
        if (_occupied.empty()) return 0.0;
        const auto& adj = _adj[v];
        double B_act = static_cast<double>(_occupied.size());
        if (adj.empty() || std::isinf(c)) return 1.0 / B_act;
        double prob = 0.0;
        for (const auto& a : adj) {
            int t = _b[a.nbr];
            long kt = _ep[t] + (_directed ? _em[t] : 0);
            long m_ts = 0;
            auto it1 = _m_out_tot[t].find(s);
            if (it1 != _m_out_tot[t].end()) m_ts += it1->second;
            if (_directed) {
                auto it2 = _m_in_tot[t].find(s);
                if (it2 != _m_in_tot[t].end()) m_ts += it2->second;
            }
            prob += (static_cast<double>(m_ts) + c) / (static_cast<double>(kt) + c * B_act);
        }
        return prob / static_cast<double>(adj.size());
    }

    // Export 1-based contiguous block labels
    Rcpp::IntegerVector external_blocks_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        lt_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i) gmap.emplace(occ[i], static_cast<int>(i) + 1);
        Rcpp::IntegerVector out(_N);
        for (int v = 0; v < _N; ++v) out[v] = gmap[_b[v]];
        return out;
    }

private:
    int _N, _C;
    int _C_active = 0;
    bool _directed, _independent, _dc, _dist, _is_base, _is_top;
    std::vector<CovSpec> _covs;
    int _K = 0;
    int _num_cov_strata = 0;
    int _cap;
    long _E = 0;
    double _base_const = 0.0;

    std::vector<LGraphEdge> _edges;
    std::vector<std::vector<LAdjRef>> _adj;
    std::vector<long> _kout, _kin;
    std::vector<std::vector<int>> _v_layers;
    std::vector<LayerCountVec> _v_kout_l, _v_kin_l;
    std::vector<long> _E_l;

    std::vector<int> _b, _wr;
    std::vector<long> _ep, _em;
    std::vector<LayerCountVec> _wr_l, _ep_l, _em_l;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;

    lt_hash_map<uint64_t, LPairStat> _pairs;
    std::vector<lt_hash_map<int, long>> _m_out_tot, _m_in_tot;

    std::vector<lt_hash_map<long, long>> _h_out, _h_in;
    std::vector<std::vector<lt_hash_map<long, long>>> _h_out_l, _h_in_l;

    std::vector<int> _normal_idx;
    std::vector<long> _B_E_c, _B_E_D_c;
    std::vector<std::vector<double>> _recsum_c, _recx2_c, _recdx_c;

    mutable MoveBuf _scratch_buf;

    inline double base_eterm(int r, int s, const LPairStat& ps) const {
        const LMath& M = lmath();
        double val = 0.0;
        for (const auto& [l, m] : ps.ml.items) {
            if (m > 0) val -= M.lgam1(m);
        }
        if (!_directed && r == s && ps.m_tot > 0)
            val -= static_cast<double>(ps.m_tot) * L_LOG2;
        return val;
    }

    inline double base_eterm_layer_delta(int r, int s, long m0, long m1) const {
        const LMath& M = lmath();
        double d = -(M.lgam1(m1) - M.lgam1(m0));
        if (!_directed && r == s) d -= static_cast<double>(m1 - m0) * L_LOG2;
        return d;
    }

    inline double dense_eterm(int r, int s, long m, long wr_r, long wr_s) const {
        if (m <= 0) return 0.0;
        int64_t nn = (_directed || r != s) ? (int64_t(wr_r) * wr_s)
                                           : (int64_t(wr_r) * (wr_r + 1) / 2);
        return lmath().lbinom(nn + m - 1, m);
    }

    double deg_dl_group_flat(int r) const {
        if (!_dc || _wr[r] == 0) return 0.0;
        const LMath& M = lmath();
        long n = _wr[r];
        if (!_dist) {
            double s = M.lbinom(n + _ep[r] - 1, _ep[r]);
            if (_directed) s += M.lbinom(n + _em[r] - 1, _em[r]);
            return s;
        }
        const LLogQ& q = llog_q();
        double s = q(_ep[r], n) + M.lgam1(n);
        for (const auto& [dk, cnt] : _h_out[r]) s -= M.lgam1(cnt);
        if (_directed) {
            s += q(_em[r], n) + M.lgam1(n);
            for (const auto& [dk, cnt] : _h_in[r]) s -= M.lgam1(cnt);
        }
        return s;
    }

    double deg_dl_group_layer(int l, int r, long n, long epl, long eml) const {
        if (!_dc || n == 0) return 0.0;
        const LMath& M = lmath();
        if (!_dist) {
            double s = M.lbinom(n + epl - 1, epl);
            if (_directed) s += M.lbinom(n + eml - 1, eml);
            return s;
        }
        const LLogQ& q = llog_q();
        double s = q(epl, n) + M.lgam1(n);
        for (const auto& [dk, cnt] : _h_out_l[l][r]) s -= M.lgam1(cnt);
        if (_directed) {
            s += q(eml, n) + M.lgam1(n);
            for (const auto& [dk, cnt] : _h_in_l[l][r]) s -= M.lgam1(cnt);
        }
        return s;
    }

    double deg_dl_delta(const lt_hash_map<long, long>& ho,
                        const lt_hash_map<long, long>* hi,
                        long n0, long n1, long ep0, long ep1,
                        long em0, long em1, long ko_v, long ki_v, long dcnt) const {
        if (!_dc) return 0.0;
        const LMath& M = lmath();
        if (!_dist) {
            double s0 = (n0 > 0) ? (M.lbinom(n0 + ep0 - 1, ep0) +
                                    (_directed ? M.lbinom(n0 + em0 - 1, em0) : 0.0)) : 0.0;
            double s1 = (n1 > 0) ? (M.lbinom(n1 + ep1 - 1, ep1) +
                                    (_directed ? M.lbinom(n1 + em1 - 1, em1) : 0.0)) : 0.0;
            return s1 - s0;
        }
        const LLogQ& q = llog_q();
        double s0 = (n0 > 0) ? (q(ep0, n0) + M.lgam1(n0)) : 0.0;
        double s1 = (n1 > 0) ? (q(ep1, n1) + M.lgam1(n1)) : 0.0;
        auto ito = ho.find(ko_v);
        long co0 = (ito == ho.end()) ? 0 : ito->second;
        double d = (s1 - s0) - (M.lgam1(co0 + dcnt) - M.lgam1(co0));
        if (_directed && hi != nullptr) {
            double si0 = (n0 > 0) ? (q(em0, n0) + M.lgam1(n0)) : 0.0;
            double si1 = (n1 > 0) ? (q(em1, n1) + M.lgam1(n1)) : 0.0;
            auto iti = hi->find(ki_v);
            long ci0 = (iti == hi->end()) ? 0 : iti->second;
            d += (si1 - si0) - (M.lgam1(ci0 + dcnt) - M.lgam1(ci0));
        }
        return d;
    }

    double compute_move_delta(int v, int r, int s, MoveBuf* out) const {
        const LMath& M = lmath();
        MoveBuf& buf = out ? *out : _scratch_buf;
        auto& pds = buf.pdeltas;
        auto& sidx = buf.scratch_idx;
        pds.clear();
        buf.d_ep_l_r.clear();
        buf.d_ep_l_s.clear();
        buf.d_em_l_r.clear();
        buf.d_em_l_s.clear();

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
                    pds.emplace_back();
                    pds.back().key = key;
                    if (_K > 0 && !_independent) pds.back().dw.zero(_K);
                    return pds.back();
                }
            }
            auto it = sidx.find(key);
            if (it == sidx.end()) {
                it = sidx.emplace(key, pds.size()).first;
                pds.emplace_back();
                pds.back().key = key;
                if (_K > 0 && !_independent) pds.back().dw.zero(_K);
            }
            return pds[it->second];
        };

        auto add_edge_delta = [&](uint64_t key, const LGraphEdge& e, long sign) {
            PairDelta& pd = find_or_add(key);
            pd.dm_tot += sign * e.m_tot;
            if (e.single_layer >= 0) pd.dml.add(e.single_layer, sign * e.m_tot);
            else pd.dml.add_all(e.ml, sign);
            if (_K > 0) {
                double sd = static_cast<double>(sign);
                if (!_independent) {
                    pd.dm_cov += sign * e.m_cov;
                    for (int k = 0; k < _K; ++k) {
                        pd.dw.x(k) += sd * e.w.x(k);
                        pd.dw.x2(k) += sd * e.w.x2(k);
                    }
                } else if (e.single_layer >= 0) {
                    pd.dwl.add_cov(e.single_layer, sign * e.m_cov, e.w, sd, _K);
                } else {
                    for (const auto& wi : e.wl.items) {
                        pd.dwl.add_cov(wi.l, sign * wi.m_cov, wi.w, sd, _K);
                    }
                }
            }
        };

        long d_ep_r = 0, d_ep_s = 0, d_em_r = 0, d_em_s = 0;
        for (const auto& a : _adj[v]) {
            const auto& e = _edges[a.eidx];
            if (e.u == e.v) {
                if (_directed && !a.out) continue;
                add_edge_delta(lcanon_pair(r, r, _directed), e, -1);
                add_edge_delta(lcanon_pair(s, s, _directed), e, +1);
                if (_directed) {
                    d_ep_r -= e.m_tot; d_em_r -= e.m_tot;
                    d_ep_s += e.m_tot; d_em_s += e.m_tot;
                } else {
                    d_ep_r -= 2 * e.m_tot;
                    d_ep_s += 2 * e.m_tot;
                }
            } else {
                int t = _b[a.nbr];
                if (_directed) {
                    if (a.out) {
                        add_edge_delta(lpack_pair(r, t), e, -1);
                        add_edge_delta(lpack_pair(s, t), e, +1);
                        d_ep_r -= e.m_tot;
                        d_ep_s += e.m_tot;
                    } else {
                        add_edge_delta(lpack_pair(t, r), e, -1);
                        add_edge_delta(lpack_pair(t, s), e, +1);
                        d_em_r -= e.m_tot;
                        d_em_s += e.m_tot;
                    }
                } else {
                    add_edge_delta(lcanon_pair(r, t, false), e, -1);
                    add_edge_delta(lcanon_pair(s, t, false), e, +1);
                    d_ep_r -= e.m_tot;
                    d_ep_s += e.m_tot;
                }
            }
        }

        if (_independent && _is_base) {
            for (const auto& [l, k] : _v_kout_l[v].items) {
                buf.d_ep_l_r.add(l, -k);
                buf.d_ep_l_s.add(l, +k);
            }
            if (_directed) {
                for (const auto& [l, k] : _v_kin_l[v].items) {
                    buf.d_em_l_r.add(l, -k);
                    buf.d_em_l_s.add(l, +k);
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
                for (int l = 0; l < _C; ++l) {
                    if (_E_l[l] > 0)
                        dS += M.lbinom(p1 + _E_l[l] - 1, _E_l[l]) -
                              M.lbinom(p0 + _E_l[l] - 1, _E_l[l]);
                }
            }
        }
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_s1) - M.lgam1(wr_s0));

        if (_is_base) {
            if (!_independent) {
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
                    const lt_hash_map<long, long>* hi_r = (_dist && _directed) ? &_h_in[r] : nullptr;
                    const lt_hash_map<long, long>* hi_s = (_dist && _directed) ? &_h_in[s] : nullptr;
                    static const lt_hash_map<long, long> empty_h;
                    dS += deg_dl_delta(_dist ? _h_out[r] : empty_h, hi_r,
                                       wr_r0, wr_r1, ep_r0, ep_r1, em_r0, em_r1, ko_v, ki_v, -1);
                    dS += deg_dl_delta(_dist ? _h_out[s] : empty_h, hi_s,
                                       wr_s0, wr_s1, ep_s0, ep_s1, em_s0, em_s1, ko_v, ki_v, +1);
                } else {
                    long dr0 = ep_r0 + em_r0, dr1 = ep_r1 + em_r1;
                    long ds0 = ep_s0 + em_s0, ds1 = ep_s1 + em_s1;
                    dS += (wr_r1 > 0 ? dr1 * M.logn(wr_r1) : 0.0) - dr0 * M.logn(wr_r0);
                    dS += ds1 * M.logn(wr_s1) - (wr_s0 > 0 ? ds0 * M.logn(wr_s0) : 0.0);
                }
            } else {
                // independent = TRUE: update vterm + deg_dl ONLY for layers in _v_layers[v]
                static const lt_hash_map<long, long> empty_h;
                for (int l : _v_layers[v]) {
                    long nr0 = _wr_l[r].get(l), nr1 = nr0 - 1;
                    long ns0 = _wr_l[s].get(l), ns1 = ns0 + 1;
                    long ko_v = _v_kout_l[v].get(l);
                    long ki_v = _directed ? _v_kin_l[v].get(l) : 0;
                    long epr0 = _ep_l[r].get(l), epr1 = epr0 - ko_v;
                    long eps0 = _ep_l[s].get(l), eps1 = eps0 + ko_v;
                    long emr0 = _directed ? _em_l[r].get(l) : 0, emr1 = emr0 - ki_v;
                    long ems0 = _directed ? _em_l[s].get(l) : 0, ems1 = ems0 + ki_v;
                    if (_dc) {
                        dS += (M.lgam1(epr1) - M.lgam1(epr0)) + (M.lgam1(eps1) - M.lgam1(eps0));
                        if (_directed)
                            dS += (M.lgam1(emr1) - M.lgam1(emr0)) + (M.lgam1(ems1) - M.lgam1(ems0));
                        const lt_hash_map<long, long>* hi_r = (_dist && _directed) ? &_h_in_l[l][r] : nullptr;
                        const lt_hash_map<long, long>* hi_s = (_dist && _directed) ? &_h_in_l[l][s] : nullptr;
                        dS += deg_dl_delta(_dist ? _h_out_l[l][r] : empty_h, hi_r,
                                           nr0, nr1, epr0, epr1, emr0, emr1, ko_v, ki_v, -1);
                        dS += deg_dl_delta(_dist ? _h_out_l[l][s] : empty_h, hi_s,
                                           ns0, ns1, eps0, eps1, ems0, ems1, ko_v, ki_v, +1);
                    } else {
                        long dr0 = epr0 + emr0, dr1 = epr1 + emr1;
                        long ds0 = eps0 + ems0, ds1 = eps1 + ems1;
                        dS += (nr1 > 0 ? dr1 * M.logn(nr1) : 0.0) - (nr0 > 0 ? dr0 * M.logn(nr0) : 0.0);
                        dS += (ns1 > 0 ? ds1 * M.logn(ns1) : 0.0) - (ns0 > 0 ? ds0 * M.logn(ns0) : 0.0);
                    }
                }
            }

            // Base-level eterm deltas across touched pairs
            for (const auto& pd : pds) {
                auto [a, b_g] = lunpack_pair(pd.key);
                auto it = _pairs.find(pd.key);
                for (const auto& [l, dm] : pd.dml.items) {
                    long m0 = (it == _pairs.end()) ? 0 : it->second.ml.get(l);
                    long m1 = m0 + dm;
                    dS += base_eterm_layer_delta(a, b_g, m0, m1);
                }
            }
        } else {
            // Dense upper level: also touch all existing pairs incident on r or s
            auto touch_dense = [&](int g) {
                for (const auto& [h, m2] : _m_out_tot[g])
                    find_or_add(lcanon_pair(g, h, _directed));
                if (_directed) {
                    for (const auto& [h, m2] : _m_in_tot[g])
                        find_or_add(lpack_pair(h, g));
                }
            };
            touch_dense(r);
            touch_dense(s);

            if (!_independent) {
                auto wr_after = [&](int g) -> long {
                    if (g == r) return wr_r1;
                    if (g == s) return wr_s1;
                    return _wr[g];
                };
                for (const auto& pd : pds) {
                    auto [a, b_g] = lunpack_pair(pd.key);
                    auto it = _pairs.find(pd.key);
                    long wa0 = _wr[a], wb0 = _wr[b_g];
                    long wa1 = wr_after(a), wb1 = wr_after(b_g);
                    // Union of layers in existing pair and pd.dml
                    if (it != _pairs.end()) {
                        for (const auto& [l, m0] : it->second.ml.items) {
                            long m1 = m0 + pd.dml.get(l);
                            dS += dense_eterm(a, b_g, m1, wa1, wb1) -
                                  dense_eterm(a, b_g, m0, wa0, wb0);
                        }
                    }
                    for (const auto& [l, dm] : pd.dml.items) {
                        if (it != _pairs.end() && it->second.ml.get(l) != 0) continue;
                        if (dm > 0) dS += dense_eterm(a, b_g, dm, wa1, wb1);
                    }
                }
            } else {
                // independent = TRUE: wr_l[g][l] changes for l in _v_layers[v]
                auto v_has_l = [&](int l) -> bool {
                    return std::binary_search(_v_layers[v].begin(), _v_layers[v].end(), l);
                };
                auto wr_l_after = [&](int g, int l) -> long {
                    long base = _wr_l[g].get(l);
                    if (g == r && v_has_l(l)) return base - 1;
                    if (g == s && v_has_l(l)) return base + 1;
                    return base;
                };
                for (const auto& pd : pds) {
                    auto [a, b_g] = lunpack_pair(pd.key);
                    auto it = _pairs.find(pd.key);
                    if (it != _pairs.end()) {
                        for (const auto& [l, m0] : it->second.ml.items) {
                            long m1 = m0 + pd.dml.get(l);
                            long wa0 = _wr_l[a].get(l), wb0 = _wr_l[b_g].get(l);
                            long wa1 = wr_l_after(a, l), wb1 = wr_l_after(b_g, l);
                            dS += dense_eterm(a, b_g, m1, wa1, wb1) -
                                  dense_eterm(a, b_g, m0, wa0, wb0);
                        }
                    }
                    for (const auto& [l, dm] : pd.dml.items) {
                        if (it != _pairs.end() && it->second.ml.get(l) != 0) continue;
                        if (dm > 0) {
                            long wa1 = wr_l_after(a, l), wb1 = wr_l_after(b_g, l);
                            dS += dense_eterm(a, b_g, dm, wa1, wb1);
                        }
                    }
                }
            }
        }

        if (_K > 0) {
            buf.d_B_E_c.assign(_num_cov_strata, 0);
            buf.d_B_E_D_c.assign(_num_cov_strata, 0);
            if (static_cast<int>(buf.d_recx2_c.size()) != _num_cov_strata) {
                buf.d_recx2_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
                buf.d_recdx_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
            } else {
                for (int st = 0; st < _num_cov_strata; ++st) {
                    buf.d_recx2_c[st].assign(_K, 0.0);
                    buf.d_recdx_c[st].assign(_K, 0.0);
                }
            }
            buf.touched_strata.assign(_num_cov_strata, 0);

            for (const auto& pd : pds) {
                auto it = _pairs.find(pd.key);
                if (!_independent) {
                    long mc0 = (it == _pairs.end()) ? 0 : it->second.m_cov;
                    long mc1 = mc0 + pd.dm_cov;
                    if (mc0 <= 0 && mc1 <= 0) continue;
                    buf.touched_strata[0] = 1;
                    buf.d_B_E_c[0] += (mc1 > 0 ? 1 : 0) - (mc0 > 0 ? 1 : 0);
                    buf.d_B_E_D_c[0] += (mc1 > 1 ? 1 : 0) - (mc0 > 1 ? 1 : 0);
                    for (int k = 0; k < _K; ++k) {
                        double x0 = (it == _pairs.end() || mc0 <= 0) ? 0.0 : it->second.w.x(k);
                        double x20 = (it == _pairs.end() || mc0 <= 0) ? 0.0 : it->second.w.x2(k);
                        double x1 = (mc1 > 0) ? (x0 + pd.dw.x(k)) : 0.0;
                        double x21 = (mc1 > 0) ? (x20 + pd.dw.x2(k)) : 0.0;
                        dS += cov_pair_neg_log_P(_covs[k], static_cast<double>(mc1), x1, x21) -
                              cov_pair_neg_log_P(_covs[k], static_cast<double>(mc0), x0, x20);
                        if (_covs[k].type == CovType::REAL_NORMAL) {
                            buf.d_recx2_c[0][k] += (mc1 > 0 ? x1 * x1 : 0.0) - (mc0 > 0 ? x0 * x0 : 0.0);
                            double dx0 = (mc0 > 1) ? std::max(0.0, x20 - x0 * (x0 / static_cast<double>(mc0))) : 0.0;
                            double dx1 = (mc1 > 1) ? std::max(0.0, x21 - x1 * (x1 / static_cast<double>(mc1))) : 0.0;
                            buf.d_recdx_c[0][k] += dx1 - dx0;
                        }
                    }
                } else {
                    for (const auto& dwi : pd.dwl.items) {
                        const LayerWItem* old_it = (it == _pairs.end()) ? nullptr : it->second.wl.find_const(dwi.l);
                        long mc0 = old_it ? old_it->m_cov : 0;
                        long mc1 = mc0 + dwi.m_cov;
                        if (mc0 <= 0 && mc1 <= 0) continue;
                        int st = dwi.l;
                        buf.touched_strata[st] = 1;
                        buf.d_B_E_c[st] += (mc1 > 0 ? 1 : 0) - (mc0 > 0 ? 1 : 0);
                        buf.d_B_E_D_c[st] += (mc1 > 1 ? 1 : 0) - (mc0 > 1 ? 1 : 0);
                        for (int k = 0; k < _K; ++k) {
                            double x0 = (old_it && mc0 > 0) ? old_it->w.x(k) : 0.0;
                            double x20 = (old_it && mc0 > 0) ? old_it->w.x2(k) : 0.0;
                            double x1 = (mc1 > 0) ? (x0 + dwi.w.x(k)) : 0.0;
                            double x21 = (mc1 > 0) ? (x20 + dwi.w.x2(k)) : 0.0;
                            dS += cov_pair_neg_log_P(_covs[k], static_cast<double>(mc1), x1, x21) -
                                  cov_pair_neg_log_P(_covs[k], static_cast<double>(mc0), x0, x20);
                            if (_covs[k].type == CovType::REAL_NORMAL) {
                                buf.d_recx2_c[st][k] += (mc1 > 0 ? x1 * x1 : 0.0) - (mc0 > 0 ? x0 * x0 : 0.0);
                                double dx0 = (mc0 > 1) ? std::max(0.0, x20 - x0 * (x0 / static_cast<double>(mc0))) : 0.0;
                                double dx1 = (mc1 > 1) ? std::max(0.0, x21 - x1 * (x1 / static_cast<double>(mc1))) : 0.0;
                                buf.d_recdx_c[st][k] += dx1 - dx0;
                            }
                        }
                    }
                }
            }

            for (int st = 0; st < _num_cov_strata; ++st) {
                if (!buf.touched_strata[st]) continue;
                long be0 = _B_E_c[st], be1 = be0 + buf.d_B_E_c[st];
                long bed0 = _B_E_D_c[st], bed1 = bed0 + buf.d_B_E_D_c[st];
                for (int k = 0; k < _K; ++k) {
                    const CovSpec& c = _covs[k];
                    if (c.type == CovType::REAL_NORMAL) {
                        if (c.micro) {
                            if (_is_top) {
                                dS += -signed_w_log_P(static_cast<double>(be1), _recsum_c[st][k],
                                                      std::max(0.0, _recx2_c[st][k] + buf.d_recx2_c[st][k]),
                                                      c.p0, c.p1, c.p2, c.p3, c.eps) -
                                      (-signed_w_log_P(static_cast<double>(be0), _recsum_c[st][k],
                                                       _recx2_c[st][k], c.p0, c.p1, c.p2, c.p3, c.eps));
                            }
                            dS += -positive_w_log_P(static_cast<double>(bed1),
                                                    std::max(0.0, _recdx_c[st][k] + buf.d_recdx_c[st][k]),
                                                    c.p2, c.p3, c.eps) -
                                  (-positive_w_log_P(static_cast<double>(bed0), _recdx_c[st][k],
                                                     c.p2, c.p3, c.eps));
                        }
                    } else if (_is_top && c.micro && be0 != be1) {
                        if (c.type == CovType::REAL_EXPONENTIAL) {
                            dS += -positive_w_log_P(static_cast<double>(be1), _recsum_c[st][k],
                                                    c.p0, c.p1, c.eps) -
                                  (-positive_w_log_P(static_cast<double>(be0), _recsum_c[st][k],
                                                     c.p0, c.p1, c.eps));
                        } else {
                            dS += -geometric_w_log_P(static_cast<double>(be1), _recsum_c[st][k], NAN, NAN) -
                                  (-geometric_w_log_P(static_cast<double>(be0), _recsum_c[st][k], NAN, NAN));
                        }
                    }
                }
            }
        }

        buf.d_ep_r = d_ep_r;
        buf.d_ep_s = d_ep_s;
        buf.d_em_r = d_em_r;
        buf.d_em_s = d_em_s;
        return dS;
    }

    void apply_move(int v, int r, int s, const MoveBuf& buf) {
        if (_is_base && _dc && _dist) {
            if (!_independent) {
                auto ito = _h_out[r].find(_kout[v]);
                if (ito != _h_out[r].end() && --ito->second == 0) _h_out[r].erase(ito);
                _h_out[s][_kout[v]]++;
                if (_directed) {
                    auto iti = _h_in[r].find(_kin[v]);
                    if (iti != _h_in[r].end() && --iti->second == 0) _h_in[r].erase(iti);
                    _h_in[s][_kin[v]]++;
                }
            } else {
                for (int l : _v_layers[v]) {
                    long ko = _v_kout_l[v].get(l);
                    auto ito = _h_out_l[l][r].find(ko);
                    if (ito != _h_out_l[l][r].end() && --ito->second == 0)
                        _h_out_l[l][r].erase(ito);
                    _h_out_l[l][s][ko]++;
                    if (_directed) {
                        long ki = _v_kin_l[v].get(l);
                        auto iti = _h_in_l[l][r].find(ki);
                        if (iti != _h_in_l[l][r].end() && --iti->second == 0)
                            _h_in_l[l][r].erase(iti);
                        _h_in_l[l][s][ki]++;
                    }
                }
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
        if (_independent) {
            for (int l : _v_layers[v]) {
                _wr_l[r].add(l, -1);
                _wr_l[s].add(l, +1);
            }
            if (_is_base) {
                _ep_l[r].add_all(buf.d_ep_l_r, +1);
                _ep_l[s].add_all(buf.d_ep_l_s, +1);
                if (_directed) {
                    _em_l[r].add_all(buf.d_em_l_r, +1);
                    _em_l[s].add_all(buf.d_em_l_s, +1);
                }
            }
        }

        auto add_m = [](lt_hash_map<int, long>& mp, int key_g, long d) {
            if (d == 0) return;
            auto mit = mp.find(key_g);
            if (mit == mp.end()) {
                if (d > 0) mp.emplace(key_g, d);
            } else {
                mit->second += d;
                if (mit->second <= 0) mp.erase(mit);
            }
        };

        if (_K > 0) {
            for (int st = 0; st < _num_cov_strata; ++st) {
                if (!buf.touched_strata[st]) continue;
                _B_E_c[st] += buf.d_B_E_c[st];
                _B_E_D_c[st] += buf.d_B_E_D_c[st];
                for (int k : _normal_idx) {
                    _recx2_c[st][k] = std::max(0.0, _recx2_c[st][k] + buf.d_recx2_c[st][k]);
                    _recdx_c[st][k] = std::max(0.0, _recdx_c[st][k] + buf.d_recdx_c[st][k]);
                }
            }
        }

        for (const auto& pd : buf.pdeltas) {
            if (pd.dm_tot == 0 && pd.dm_cov == 0 && pd.dml.empty() && pd.dwl.empty() &&
                (_K == 0 || _independent)) continue;
            auto [a, b_g] = lunpack_pair(pd.key);
            auto it = _pairs.find(pd.key);
            if (it == _pairs.end()) {
                if (pd.dm_tot <= 0 && pd.dm_cov <= 0) continue;
                it = _pairs.emplace(pd.key, LPairStat{}).first;
            }
            LPairStat& ps = it->second;
            ps.m_tot += pd.dm_tot;
            ps.ml.add_all(pd.dml, +1);
            if (_K > 0) {
                if (!_independent) {
                    ps.m_cov += pd.dm_cov;
                    if (ps.w.k == 0) ps.w.zero(_K);
                    for (int k = 0; k < _K; ++k) {
                        ps.w.x(k) += pd.dw.x(k);
                        ps.w.x2(k) += pd.dw.x2(k);
                    }
                } else {
                    for (const auto& dwi : pd.dwl.items) {
                        ps.wl.add_cov(dwi.l, dwi.m_cov, dwi.w, +1.0, _K, true);
                    }
                }
            }
            if (ps.m_tot <= 0 && ps.m_cov <= 0) _pairs.erase(it);

            if (_directed) {
                add_m(_m_out_tot[a], b_g, pd.dm_tot);
                add_m(_m_in_tot[b_g], a, pd.dm_tot);
            } else {
                if (a == b_g) {
                    add_m(_m_out_tot[a], a, 2 * pd.dm_tot);
                } else {
                    add_m(_m_out_tot[a], b_g, pd.dm_tot);
                    add_m(_m_out_tot[b_g], a, pd.dm_tot);
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

// ---------------------------------------------------------------------------
// Helpers to build level-0 edges and parse 1-based R blocks/layers

std::vector<LGraphEdge> build_base_layered_edges(const Graph& graph,
                                                 const Rcpp::IntegerVector& layers,
                                                 const Rcpp::NumericMatrix& weights,
                                                 int K,
                                                 int& C_out) {
    int E = static_cast<int>(boost::num_edges(graph.adjacency));
    if (layers.size() != E)
        Rcpp::stop("layers must have length graph_num_edges(graph)");
    if (K > 0 && weights.nrow() != E)
        Rcpp::stop("weights must have nrow(weights) == graph_num_edges(graph)");
    int max_l = 0;
    for (int i = 0; i < E; ++i) {
        if (layers[i] == NA_INTEGER || layers[i] < 1)
            Rcpp::stop("layers must contain positive integers >= 1");
        if (layers[i] > max_l) max_l = layers[i];
    }
    C_out = std::max(1, max_l);
    std::vector<LGraphEdge> edges;
    edges.reserve(E);
    int idx = 0;
    for (auto [it, end] = boost::edges(graph.adjacency); it != end; ++it, ++idx) {
        LGraphEdge e;
        e.u = static_cast<int>(boost::source(*it, graph.adjacency));
        e.v = static_cast<int>(boost::target(*it, graph.adjacency));
        e.m_tot = 1;
        e.single_layer = layers[idx] - 1;
        e.m_cov = 1;
        if (K > 0) {
            e.w.zero(K);
            for (int k = 0; k < K; ++k) {
                double val = weights(static_cast<R_xlen_t>(idx), k);
                e.w.x(k) = val;
                e.w.x2(k) = val * val;
            }
        }
        edges.push_back(std::move(e));
    }
    return edges;
}

std::vector<LGraphEdge> build_upper_layered_edges(const LayeredLevel& lower,
                                                  lt_hash_map<int, int>& gmap,
                                                  std::vector<std::vector<int>>* next_vl = nullptr) {
    std::vector<int> occ = lower.occupied();
    std::sort(occ.begin(), occ.end());
    gmap.clear();
    gmap.reserve(occ.size());
    for (size_t i = 0; i < occ.size(); ++i)
        gmap.emplace(occ[i], static_cast<int>(i));

    int K = lower.K();
    bool indep = lower.independent();
    std::vector<LGraphEdge> next_edges;
    next_edges.reserve(lower.pairs().size());
    for (const auto& [pk, ps] : lower.pairs()) {
        if (ps.m_tot <= 0) continue;
        auto [r, s] = lunpack_pair(pk);
        LGraphEdge e;
        e.u = gmap[r];
        e.v = gmap[s];
        e.m_tot = ps.m_tot;
        e.single_layer = -1;
        e.ml = ps.ml;
        if (K > 0) {
            if (!indep) {
                e.m_cov = (ps.m_cov > 0) ? 1 : 0;
                e.w.zero(K);
                if (ps.m_cov > 0) {
                    for (int k = 0; k < K; ++k) {
                        double x = ps.w.x(k);
                        e.w.x(k) = x;
                        e.w.x2(k) = x * x;
                    }
                }
            } else {
                e.m_cov = 0;
                for (const auto& wi : ps.wl.items) {
                    if (wi.m_cov <= 0) continue;
                    WVec w_up(K);
                    for (int k = 0; k < K; ++k) {
                        double x = wi.w.x(k);
                        w_up.x(k) = x;
                        w_up.x2(k) = x * x;
                    }
                    e.wl.add_cov(wi.l, 1L, w_up, +1.0, K);
                }
            }
        }
        next_edges.push_back(std::move(e));
    }
    if (next_vl != nullptr && indep) {
        int next_N = static_cast<int>(occ.size());
        next_vl->assign(next_N, {});
        for (size_t i = 0; i < occ.size(); ++i) {
            const auto& wl = lower.wr_l(occ[i]);
            for (const auto& [lay, cnt] : wl.items) {
                if (cnt > 0) (*next_vl)[i].push_back(lay);
            }
        }
    }
    return next_edges;
}

std::vector<int> parse_l_blocks(const Rcpp::IntegerVector& blocks, int N) {
    if (blocks.size() != N)
        Rcpp::stop("blocks must have length equal to number of vertices");
    lt_hash_map<int, int> remap;
    std::vector<int> out(N);
    int next_id = 0;
    for (int i = 0; i < N; ++i) {
        int val = blocks[i];
        if (val == NA_INTEGER || val < 1)
            Rcpp::stop("blocks must contain positive integer labels");
        auto it = remap.find(val);
        if (it == remap.end()) it = remap.emplace(val, next_id++).first;
        out[i] = it->second;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Flat greedy minimization & golden-section / agglomerative search over B

void layered_greedy_sweeps(LayeredLevel& lv, LRng& rng, size_t B_min, size_t B_max,
                           int max_sweeps, double epsilon, int nthreads = 0) {
    if (lv.N() <= 1) return;
    int nw = rgt_parallel::resolve_threads(nthreads);
    std::vector<int> vs(lv.N());
    std::iota(vs.begin(), vs.end(), 0);
    for (int sw = 0; sw < max_sweeps; ++sw) {
        double tot_dS = 0.0;
        rng.shuffle(vs);
        if (nw > 1 && lv.N() >= 512 && lv.B() > 1) {
            size_t N = static_cast<size_t>(lv.N());
            size_t batch = std::min(N, std::max<size_t>(512, static_cast<size_t>(nw) * 256));
            std::vector<LayeredLevel::MoveBuf> bufs(static_cast<size_t>(nw));
            std::vector<lt_hash_set<int>> cand_sets(static_cast<size_t>(nw));
            std::vector<int> cand_s;
            for (size_t base = 0; base < N; base += batch) {
                size_t len = std::min(batch, N - base);
                cand_s.assign(len, -1);
                uint64_t b_seed = static_cast<uint64_t>(rng.index(1ULL << 30)) ^
                                  (static_cast<uint64_t>(rng.index(1ULL << 30)) << 30) ^
                                  0x9e3779b97f4a7c15ULL;
                const size_t B_snap = static_cast<size_t>(lv.B());
                rgt_parallel::parallel_for(len, nw, [&](int tid, size_t lo, size_t hi) {
                    rgt_parallel::ThreadRng trng(
                        b_seed ^ (static_cast<uint64_t>(lo + 1) * 0xbf58476d1ce4e5b9ULL));
                    auto& buf = bufs[static_cast<size_t>(tid)];
                    auto& cands = cand_sets[static_cast<size_t>(tid)];
                    for (size_t idx = lo; idx < hi; ++idx) {
                        int v = vs[base + idx];
                        int r = lv.group(v);
                        cands.clear();
                        const auto& occ = lv.occupied();
                        if (occ.size() <= 12) {
                            for (int g : occ) cands.insert(g);
                        } else {
                            for (int k = 0; k < 8; ++k) cands.insert(lv.sample_group(v, 0.5, trng));
                            for (int k = 0; k < 4; ++k) cands.insert(lv.random_occupied(trng));
                        }
                        if (B_snap < B_max && B_snap < N &&
                            lv.size(r) > 1 && lv.empty_group() >= 0) {
                            cands.insert(lv.empty_group());
                        }
                        cands.erase(r);
                        double best_dS = -1e-11;
                        int best_s = -1;
                        for (int s : cands) {
                            if (lv.size(r) == 1 && lv.size(s) > 0 && B_snap <= B_min) continue;
                            double dS = lv.virtual_move_const(v, s, buf);
                            if (dS < best_dS) { best_dS = dS; best_s = s; }
                        }
                        cand_s[idx] = best_s;
                    }
                }, 32);
                for (size_t idx = 0; idx < len; ++idx) {
                    int s = cand_s[idx];
                    if (s < 0) continue;
                    int v = vs[base + idx];
                    int r = lv.group(v);
                    if (s == r) continue;
                    size_t B_cur = static_cast<size_t>(lv.B());
                    if (lv.size(r) == 1 && lv.size(s) > 0 && B_cur <= B_min) continue;
                    if (lv.size(s) == 0 && B_cur >= B_max && lv.size(r) > 1) continue;
                    if (lv.virtual_move(v, s) < -1e-11) tot_dS += lv.move(v, s);
                }
            }
        } else {
            for (int v : vs) {
                int r = lv.group(v);
                size_t B_cur = static_cast<size_t>(lv.B());
                double best_dS = -1e-11;
                int best_s = r;

                lt_hash_set<int> cands;
                const auto& occ = lv.occupied();
                if (occ.size() <= 12) {
                    for (int g : occ) cands.insert(g);
                } else {
                    for (int k = 0; k < 8; ++k) cands.insert(lv.sample_group(v, 0.5, rng));
                    for (int k = 0; k < 4; ++k) cands.insert(lv.random_occupied(rng));
                }
                if (B_cur < B_max && B_cur < static_cast<size_t>(lv.N()) &&
                    lv.size(r) > 1 && lv.empty_group() >= 0) {
                    cands.insert(lv.empty_group());
                }
                cands.erase(r);

                for (int s : cands) {
                    if (lv.size(r) == 1 && lv.size(s) > 0 && B_cur <= B_min) continue;
                    double dS = lv.virtual_move(v, s);
                    if (dS < best_dS) {
                        best_dS = dS;
                        best_s = s;
                    }
                }
                if (best_s != r) {
                    tot_dS += lv.move(v, best_s);
                }
            }
        }
        lv.refresh_cov_globals();
        if (std::abs(tot_dS) <= epsilon) break;
    }
}

void layered_agglomerate_to(LayeredLevel& lv, LRng& rng, size_t target_B,
                            double epsilon, int nthreads = 0) {
    while (static_cast<size_t>(lv.B()) > target_B) {
        std::vector<int> occ = lv.occupied();
        size_t cur_B = occ.size();
        size_t n_merge = std::max<size_t>(1, (cur_B - target_B + 2) / 3);
        if (cur_B - target_B <= 4) n_merge = 1;

        struct MergeCand {
            double dS;
            int r, t;
        };
        std::vector<MergeCand> cands;
        cands.reserve(occ.size());
        for (int r : occ) {
            if (lv.size(r) == 0) continue;
            lt_hash_set<int> nbrs;
            if (occ.size() <= 8) {
                for (int t : occ) if (t != r) nbrs.insert(t);
            } else {
                int v0 = lv.members(r)[0];
                for (int k = 0; k < 5; ++k) {
                    int t = lv.sample_group(v0, 0.5, rng);
                    if (t != r) nbrs.insert(t);
                }
                for (int k = 0; k < 2; ++k) {
                    int t = lv.random_occupied(rng);
                    if (t != r) nbrs.insert(t);
                }
            }
            int best_t = -1;
            double best_dS = L_INF;
            for (int t : nbrs) {
                double dS = lv.virtual_merge(r, t);
                if (dS < best_dS) {
                    best_dS = dS;
                    best_t = t;
                }
            }
            if (best_t >= 0) cands.push_back({best_dS, r, best_t});
        }
        if (cands.empty()) break;
        std::sort(cands.begin(), cands.end(),
                  [](const MergeCand& a, const MergeCand& b) { return a.dS < b.dS; });

        lt_hash_set<int> touched;
        size_t merged = 0;
        for (const auto& c : cands) {
            if (merged >= n_merge || static_cast<size_t>(lv.B()) <= target_B) break;
            if (touched.count(c.r) || touched.count(c.t)) continue;
            if (lv.size(c.r) == 0 || lv.size(c.t) == 0) continue;
            lv.merge(c.r, c.t);
            touched.insert(c.r);
            touched.insert(c.t);
            ++merged;
        }
        if (merged == 0) {
            int r = lv.occupied()[0], t = lv.occupied()[1];
            lv.merge(r, t);
        }
        layered_greedy_sweeps(lv, rng, target_B, std::max(target_B, static_cast<size_t>(lv.B())),
                              2, epsilon, nthreads);
    }
}

double layered_minimize_flat(LayeredLevel& lv, LRng& rng, size_t B_min,
                             size_t B_max, bool refine, double epsilon,
                             std::vector<std::vector<int>>* checkpoints = nullptr,
                             int nthreads = 0) {
    int N = lv.N();
    if (N <= 1) {
        std::vector<int> b(N, 0);
        lv.set_partition(b);
        if (checkpoints) checkpoints->push_back(b);
        return lv.entropy();
    }
    B_min = std::max<size_t>(1, std::min<size_t>(B_min, static_cast<size_t>(N)));
    B_max = std::max<size_t>(B_min, std::min<size_t>(B_max, static_cast<size_t>(N)));

    std::vector<int> b_init(N);
    std::iota(b_init.begin(), b_init.end(), 0);
    lv.set_partition(b_init);

    if (B_max < static_cast<size_t>(N))
        layered_agglomerate_to(lv, rng, B_max, epsilon, nthreads);

    std::vector<int> best_b = lv.blocks();
    double best_S = lv.entropy();

    auto record_cp = [&]() {
        if (!checkpoints) return;
        Rcpp::IntegerVector eb = lv.external_blocks_r();
        std::vector<int> b(N);
        for (int i = 0; i < N; ++i) b[i] = eb[i] - 1;
        checkpoints->push_back(std::move(b));
    };

    size_t cur_B = static_cast<size_t>(lv.B());
    size_t next_cp_max = std::min<size_t>(B_max, 32);
    while (cur_B > B_min) {
        size_t next_B;
        if (cur_B <= 12) next_B = cur_B - 1;
        else next_B = std::max(B_min, static_cast<size_t>(std::floor(static_cast<double>(cur_B) * 0.70)));
        layered_agglomerate_to(lv, rng, next_B, epsilon, nthreads);
        if (refine) layered_greedy_sweeps(lv, rng, next_B, next_B, 3, epsilon, nthreads);
        double S = lv.entropy();
        if (S < best_S) {
            best_S = S;
            best_b = lv.blocks();
        }
        cur_B = static_cast<size_t>(lv.B());
        if (checkpoints && cur_B <= next_cp_max) {
            record_cp();
            next_cp_max = (cur_B <= 6) ? (cur_B - 1) : static_cast<size_t>(std::floor(cur_B * 0.72));
        }
    }

    lv.set_partition(best_b);
    if (refine) {
        layered_greedy_sweeps(lv, rng, B_min, B_max, 12, epsilon, nthreads);
        while (static_cast<size_t>(lv.B()) > B_min) {
            std::vector<int> occ = lv.occupied();
            if (occ.size() > 16) break;
            double best_mdS = -1e-9;
            int mr = -1, mt = -1;
            for (size_t i = 0; i < occ.size(); ++i) {
                for (size_t j = i + 1; j < occ.size(); ++j) {
                    double dS = lv.virtual_merge(occ[i], occ[j]);
                    if (dS < best_mdS) {
                        best_mdS = dS;
                        mr = occ[i];
                        mt = occ[j];
                    }
                }
            }
            if (mr < 0) break;
            lv.merge(mr, mt);
            layered_greedy_sweeps(lv, rng, B_min, B_max, 6, epsilon, nthreads);
        }
        double S_ref = lv.entropy();
        if (S_ref < best_S) {
            best_S = S_ref;
            best_b = lv.blocks();
        } else {
            lv.set_partition(best_b);
        }
    }
    if (checkpoints) record_cp();
    return best_S;
}

// ---------------------------------------------------------------------------
// Nested Layered SBM (NestedBlockState(base_state=LayeredBlockState))

class NestedLayeredSBM {
public:
    NestedLayeredSBM(const Graph& graph, const Rcpp::IntegerVector& layers,
                     bool independent, bool dc, bool dist,
                     const Rcpp::NumericMatrix& weights = Rcpp::NumericMatrix(0, 0),
                     std::vector<CovSpec> covs = {})
        : _directed(graph.directed), _independent(independent), _dc(dc),
          _dist(dist), _base_N(static_cast<int>(boost::num_vertices(graph.adjacency))),
          _base_covs(std::move(covs)),
          _up_covs(make_upper_covs(_base_covs)) {
        _base_edges = build_base_layered_edges(
            graph, layers, weights, static_cast<int>(_base_covs.size()), _C);
    }

    void build(const std::vector<std::vector<int>>& bs) {
        _levels.clear();
        int L = static_cast<int>(bs.size());
        if (L < 1) Rcpp::stop("Nested layered state requires at least 1 level");

        int cur_N = _base_N;
        std::vector<LGraphEdge> cur_edges = _base_edges;
        std::vector<std::vector<int>> cur_v_layers;

        for (int l = 0; l < L; ++l) {
            bool is_base = (l == 0);
            bool is_top = (l == L - 1);
            LayeredLevel lv(cur_N, _C, _directed, _independent,
                            is_base ? _dc : false, is_base ? _dist : false,
                            is_base, is_top,
                            is_base ? _base_covs : _up_covs);
            lv.set_edges(std::move(cur_edges),
                         (!is_base && _independent) ? &cur_v_layers : nullptr);
            lv.set_partition(bs[l]);

            if (l + 1 < L) {
                lt_hash_map<int, int> gmap;
                cur_edges = build_upper_layered_edges(lv, gmap, &cur_v_layers);
                cur_N = static_cast<int>(gmap.size());
            }
            _levels.push_back(std::move(lv));
        }
    }

    int num_levels() const { return static_cast<int>(_levels.size()); }
    LayeredLevel& level(int l) { return _levels[l]; }
    const LayeredLevel& level(int l) const { return _levels[l]; }

    std::vector<double> level_entropies() const {
        int L = static_cast<int>(_levels.size());
        std::vector<double> out(L, 0.0);
        bool has_normal = !_levels[0].normal_idx().empty();
        int S_strata = _levels[0].num_cov_strata();
        std::vector<long> L_active(S_strata, 0);
        std::vector<std::vector<double>> L_recdx(
            S_strata, std::vector<double>(_base_covs.size(), 0.0));
        if (L > 1 && has_normal) {
            for (int l = 0; l < L; ++l) {
                for (int st = 0; st < S_strata; ++st) {
                    long bed = _levels[l].B_E_D_c(st);
                    if (bed > 0) ++L_active[st];
                    for (int k : _levels[0].normal_idx()) {
                        L_recdx[st][k] += _levels[l].recdx_c(st, k) * static_cast<double>(bed);
                    }
                }
            }
        }
        for (int l = 0; l < L; ++l) {
            double s = _levels[l].entropy(l == 0);
            if (L > 1 && has_normal) {
                for (int st = 0; st < S_strata; ++st) {
                    for (int k : _levels[0].normal_idx()) {
                        s += -w_safelog(static_cast<double>(_levels[l].B_E_D_c(st)));
                        if (l == 0) {
                            const CovSpec& c = _base_covs[k];
                            s += -positive_w_log_P(static_cast<double>(L_active[st]),
                                                   L_recdx[st][k], c.p2, c.p3, c.eps);
                        }
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

    std::vector<std::vector<int>> save_bs() const {
        std::vector<std::vector<int>> bs;
        bs.reserve(_levels.size());
        for (const auto& lv : _levels) {
            Rcpp::IntegerVector eb = lv.external_blocks_r();
            std::vector<int> b(eb.size());
            for (int i = 0; i < eb.size(); ++i) b[i] = eb[i] - 1;
            bs.push_back(std::move(b));
        }
        return bs;
    }

    Rcpp::List external_bs_r() const {
        Rcpp::List out(_levels.size());
        for (size_t l = 0; l < _levels.size(); ++l)
            out[l] = _levels[l].external_blocks_r();
        return out;
    }

    // Fast rebuild of upper levels (l >= start_l + 1) from existing _levels[start_l],
    // preserving upper-level partitions without touching level 0!
    void rebuild_upper_from(int start_l) {
        int L = static_cast<int>(_levels.size());
        if (start_l + 1 >= L) return;
        std::vector<std::vector<int>> upper_bs;
        upper_bs.reserve(L - start_l - 1);
        for (int l = start_l + 1; l < L; ++l) {
            Rcpp::IntegerVector eb = _levels[l].external_blocks_r();
            std::vector<int> b(eb.size());
            for (int i = 0; i < eb.size(); ++i) b[i] = eb[i] - 1;
            upper_bs.push_back(std::move(b));
        }
        _levels.erase(_levels.begin() + start_l + 1, _levels.end());

        for (int idx = 0; idx < static_cast<int>(upper_bs.size()); ++idx) {
            int l = start_l + 1 + idx;
            const LayeredLevel& prev = _levels[l - 1];
            lt_hash_map<int, int> gmap;
            std::vector<std::vector<int>> next_vl;
            std::vector<LGraphEdge> next_edges = build_upper_layered_edges(prev, gmap, &next_vl);
            int next_N = static_cast<int>(gmap.size());
            bool is_top = (l == L - 1);
            LayeredLevel ulv(next_N, _C, _directed, _independent, false, false, false, is_top, _up_covs);
            ulv.set_edges(std::move(next_edges), _independent ? &next_vl : nullptr);
            ulv.set_partition(upper_bs[idx]);
            _levels.push_back(std::move(ulv));
        }
    }

    // Rebuild upper levels keeping current upper assignments (when B at level 0 is unchanged)
    void rebuild_upper_keep_bs() {
        rebuild_upper_from(0);
    }

    // Build optimal upper hierarchy above a given level-0 partition b0
    void optimize_from_base(const std::vector<int>& b0, LRng& rng,
                            bool refine, double epsilon, int nthreads = 0) {
        std::vector<std::vector<int>> bs;
        bs.push_back(b0);

        // Build base level first to extract coarsened graph
        LayeredLevel lv0(_base_N, _C, _directed, _independent, _dc, _dist, true, false, _base_covs);
        lv0.set_edges(_base_edges);
        lv0.set_partition(b0);

        std::vector<LGraphEdge> cur_edges;
        std::vector<std::vector<int>> cur_vl;
        lt_hash_map<int, int> gmap;
        cur_edges = build_upper_layered_edges(lv0, gmap, &cur_vl);
        int cur_N = static_cast<int>(gmap.size());

        while (cur_N > 1) {
            LayeredLevel ulv(cur_N, _C, _directed, _independent, false, false, false, true, _up_covs);
            ulv.set_edges(cur_edges, _independent ? &cur_vl : nullptr);
            layered_minimize_flat(ulv, rng, 1, static_cast<size_t>(cur_N), refine, epsilon, nullptr, nthreads);
            if (ulv.B() >= cur_N) {
                // Force collapse to B=1 if no compression at this level
                std::vector<int> b1(cur_N, 0);
                ulv.set_partition(b1);
            }
            Rcpp::IntegerVector eb = ulv.external_blocks_r();
            std::vector<int> bl(eb.size());
            for (int i = 0; i < eb.size(); ++i) bl[i] = eb[i] - 1;
            bs.push_back(std::move(bl));
            if (ulv.B() == 1) break;
            cur_edges = build_upper_layered_edges(ulv, gmap, &cur_vl);
            cur_N = static_cast<int>(gmap.size());
        }
        if (bs.size() == 1 && lv0.B() == 1) {
            bs.push_back(std::vector<int>(1, 0));
        }
        build(bs);
    }

    // Joint refinement across levels under full nested entropy
    void refine_joint(LRng& rng, size_t B_min, size_t B_max, double epsilon, int nthreads = 0) {
        double best_S = entropy();
        std::vector<std::vector<int>> best_bs = save_bs();
        std::vector<int> vs(_base_N);
        std::iota(vs.begin(), vs.end(), 0);

        for (int pass = 0; pass < 4; ++pass) {
            bool improved = false;
            LayeredLevel& lv0 = _levels[0];
            rng.shuffle(vs);
            for (int v : vs) {
                int r = lv0.group(v);
                size_t B_cur = static_cast<size_t>(lv0.B());
                lt_hash_set<int> cands;
                for (int k = 0; k < 4; ++k) cands.insert(lv0.sample_group(v, 0.5, rng));
                cands.erase(r);
                int best_t = r;
                double best_cand_S = best_S - 1e-8;
                for (int t : cands) {
                    if (lv0.size(r) == 1 && B_cur <= B_min) continue;
                    double d0 = lv0.virtual_move(v, t);
                    if (d0 > 25.0) continue;
                    int old_B0 = lv0.B();
                    lv0.apply_last_virtual_move(v, t);
                    bool changed_B0 = (lv0.B() != old_B0);
                    double cand_S = L_INF;
                    if (!changed_B0) {
                        rebuild_upper_keep_bs();
                        cand_S = entropy();
                    } else {
                        Rcpp::IntegerVector eb0 = lv0.external_blocks_r();
                        std::vector<int> b0_c(_base_N);
                        for (int i = 0; i < _base_N; ++i) b0_c[i] = eb0[i] - 1;
                        NestedLayeredSBM tmp = *this;
                        tmp.optimize_from_base(b0_c, rng, false, epsilon, nthreads);
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
                        optimize_from_base(b0_c, rng, true, epsilon, nthreads);
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
            // Try level-0 merges under full nested entropy
            while (static_cast<size_t>(_levels[0].B()) > B_min) {
                std::vector<int> occ = _levels[0].occupied();
                if (occ.size() > 30) break;
                int mr = -1;
                double mS = best_S - 1e-8;
                std::vector<std::vector<int>> m_bs;
                for (size_t i = 0; i < occ.size(); ++i) {
                    for (size_t j = i + 1; j < occ.size(); ++j) {
                        double d0 = _levels[0].virtual_merge(occ[i], occ[j]);
                        if (d0 > 25.0) continue;
                        NestedLayeredSBM tmp = *this;
                        tmp._levels[0].merge(occ[i], occ[j]);
                        Rcpp::IntegerVector eb0 = tmp._levels[0].external_blocks_r();
                        std::vector<int> b0_c(_base_N);
                        for (int k = 0; k < _base_N; ++k) b0_c[k] = eb0[k] - 1;
                        tmp.optimize_from_base(b0_c, rng, false, epsilon, nthreads);
                        double s_try = tmp.entropy();
                        if (s_try < mS) {
                            mS = s_try;
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
            // Re-optimize upper levels from best level-0
            if (best_bs[0].size() > 1 && static_cast<size_t>(_levels[0].B()) <= B_max) {
                NestedLayeredSBM tmp = *this;
                tmp.optimize_from_base(best_bs[0], rng, true, epsilon, nthreads);
                double s_up = tmp.entropy();
                if (s_up < best_S - 1e-8) {
                    best_S = s_up;
                    best_bs = tmp.save_bs();
                    build(best_bs);
                    improved = true;
                }
            }
            if (!improved) break;
        }
    }

private:
    bool _directed, _independent, _dc, _dist;
    int _base_N, _C = 1;
    std::vector<CovSpec> _base_covs, _up_covs;
    std::vector<LGraphEdge> _base_edges;
    std::vector<LayeredLevel> _levels;
};

// ---------------------------------------------------------------------------
// Finite-temperature MCMC (single-node + merge-split) for flat LayeredLevel

void layered_flat_mcmc_sweep(LayeredLevel& s, LRng& rng, double beta,
                             double c, double d, bool multiflip,
                             size_t B_min, size_t B_max,
                             double& total_dS, size_t& nattempts,
                             size_t& nmoves) {
    if (s.N() <= 1 || s.E() == 0) return;
    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    rng.shuffle(vs);

    // 1. Single-node Metropolis-Hastings sweep (exact early rejection!)
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
        double log_u = (u_rand > 0.0) ? std::log(u_rand) : -L_INF;
        // Since prev <= 1.0, log_a <= -beta * dS - log(pfwd): exact early reject!
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

    // 2. Merge-split moves (auxiliary Gibbs 2-way split / merge)
    if (multiflip && s.B() >= 1) {
        if (std::isinf(beta)) {
            // Greedy (beta = Inf) minimization: detailed balance is irrelevant,
            // so keep the Gibbs-guided splits that find improving moves.
            size_t n_ms = std::max<size_t>(1, static_cast<size_t>(s.B()));
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

                    double dS = 0.0;
                    double log_q_split = 0.0;
                    int cnt_t = 0;
                    for (int v : mem) {
                        if (rng.unif() < 0.5) {
                            dS += s.move(v, t);
                            ++cnt_t;
                        }
                    }
                    if (cnt_t == 0) {
                        dS += s.move(mem[0], t);
                        cnt_t = 1;
                    } else if (cnt_t == static_cast<int>(nr)) {
                        dS += s.move(mem[0], r);
                        --cnt_t;
                    }
                    for (int gsw = 0; gsw < 3; ++gsw) {
                        rng.shuffle(mem);
                        for (int v : mem) {
                            int cur = s.group(v);
                            int oth = (cur == r) ? t : r;
                            if (s.size(cur) <= 1) continue;
                            double d_mv = s.virtual_move(v, oth);
                            double p_oth;
                            if (std::isinf(beta)) {
                                p_oth = (d_mv < 0.0) ? 1.0 : 0.0;
                            } else {
                                double z = -beta * d_mv;
                                p_oth = (z > 35.0) ? 1.0 : (z < -35.0 ? 0.0 : 1.0 / (1.0 + std::exp(-z)));
                            }
                            bool go_oth = (rng.unif() < p_oth);
                            if (go_oth) {
                                s.apply_last_virtual_move(v, oth);
                                dS += d_mv;
                                if (gsw == 2) log_q_split += std::log(std::max(1e-300, p_oth));
                            } else {
                                if (gsw == 2) log_q_split += std::log(std::max(1e-300, 1.0 - p_oth));
                            }
                        }
                    }
                    bool accept = false;
                    if (std::isinf(beta)) {
                        accept = (dS < 0.0);
                    } else {
                        size_t B_after = static_cast<size_t>(s.B());
                        double log_pf = -std::log(static_cast<double>(B_before)) + log_q_split;
                        double log_pr = -std::log(static_cast<double>(B_after) * (B_after - 1) / 2.0);
                        double log_a = -beta * dS + log_pr - log_pf;
                        accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                    }
                    if (accept) {
                        total_dS += dS;
                        ++nmoves;
                    } else {
                        std::vector<int> mem_t = s.members(t);
                        for (int v : mem_t) s.move(v, r);
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
                    if (std::isinf(beta)) {
                        accept = (dS < 0.0);
                    } else {
                        size_t B_after = static_cast<size_t>(s.B());
                        double log_pf = -std::log(static_cast<double>(B_before) * (B_before - 1) / 2.0);
                        double log_pr = -std::log(static_cast<double>(B_after)) -
                                        static_cast<double>(mem_r.size()) * L_LOG2;
                        double log_a = -beta * dS + log_pr - log_pf;
                        accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                    }
                    if (accept) {
                        total_dS += dS;
                        ++nmoves;
                    } else {
                        for (int v : mem_r) s.move(v, r);
                    }
                }
            }
        } else {
            // Finite beta: exact MH kernel.
            // A state-independent number of merge-split attempts per sweep (as in
            // blockmodel.cpp: psplit + pmerge = 2). Using s.B() made the number of
            // pi-invariant steps depend on the current state, which breaks
            // invariance (states with more groups got more merge attempts).
            const size_t n_ms = 2;
            for (size_t step = 0; step < n_ms; ++step) {
                // Uniform random 2-way split / pair merge with exact proposal
                // probabilities (see weighted.cpp): q_split = (1/B) 2^(1 - n_r),
                // q_merge = 2 / B^2 for an unordered pair. The earlier Gibbs
                // split proposal did not account for its random launch and
                // earlier sweeps, which broke detailed balance.
                const double LN2 = std::log(2.0);
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
                    bool accept = false;
                    if (std::isinf(beta)) {
                        accept = (dS < 0.0);
                    } else {
                        double log_pf = -std::log(static_cast<double>(B_before)) +
                                        (1.0 - static_cast<double>(nr)) * LN2;
                        double log_pr = LN2 - 2.0 * std::log(static_cast<double>(B_before + 1));
                        double log_a = -beta * dS + log_pr - log_pf;
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
                    int r = s.random_occupied(rng);
                    int t = s.random_occupied(rng);
                    if (r == t) continue;
                    size_t B_before = static_cast<size_t>(s.B());
                    std::vector<int> mem_r = s.members(r);
                    size_t n_tot = mem_r.size() + static_cast<size_t>(s.size(t));
                    double dS = s.merge(r, t);
                    bool accept = false;
                    if (std::isinf(beta)) {
                        accept = (dS < 0.0);
                    } else {
                        double log_pf = LN2 - 2.0 * std::log(static_cast<double>(B_before));
                        double log_pr = -std::log(static_cast<double>(B_before - 1)) +
                                        (1.0 - static_cast<double>(n_tot)) * LN2;
                        double log_a = -beta * dS + log_pr - log_pf;
                        accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                    }
                    if (accept) {
                        total_dS += dS;
                        ++nmoves;
                    } else {
                        for (int v : mem_r) s.move(v, r);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Finite-temperature MCMC for NestedLayeredSBM

void layered_nested_mcmc_sweep(NestedLayeredSBM& nw, LRng& rng, double beta,
                               double c, double d, size_t B_min, size_t B_max,
                               double& total_dS, size_t& nattempts,
                               size_t& nmoves) {
    int L = nw.num_levels();
    double cur_S = nw.entropy();

    // Sweep each level from 0 up to L - 2
    for (int l = 0; l < L; ++l) {
        LayeredLevel& lv = nw.level(l);
        if (lv.N() <= 1) continue;
        if (l == L - 1 && lv.B() == 1) continue;
        std::vector<int> vs(lv.N());
        std::iota(vs.begin(), vs.end(), 0);
        rng.shuffle(vs);

        size_t l_bmin = (l == 0) ? B_min : 1;
        size_t l_bmax = (l == 0) ? B_max : static_cast<size_t>(lv.N());
        double cl = c * std::pow(2.0, static_cast<double>(l));

        for (int v : vs) {
            ++nattempts;
            int r = lv.group(v);
            size_t B_cur = static_cast<size_t>(lv.B());
            bool can_create = (l == 0 && B_cur < static_cast<size_t>(lv.N()) &&
                               B_cur < l_bmax && d > 0.0 && lv.empty_group() >= 0);
            double p_new = can_create ? d : 0.0;
            bool propose_new = (can_create && rng.unif() < p_new);
            int t;
            double pfwd = 0.0;
            if (propose_new) {
                if (lv.size(r) == 1) continue;
                t = lv.empty_group();
                pfwd = p_new;
            } else {
                t = lv.sample_group(v, cl, rng);
                if (t == r) continue;
                if (lv.size(r) == 1 && B_cur <= l_bmin) continue;
                pfwd = (1.0 - p_new) * lv.move_prob(v, t, cl);
            }
            if (pfwd <= 0.0) continue;

            // Fast level-l pre-check for massive rejections
            double dSl = lv.virtual_move(v, t);
            if (!std::isinf(beta) && -beta * dSl < -40.0) continue;

            bool empties = (lv.size(r) == 1);
            bool creates = (lv.size(t) == 0);

            if (!empties && !creates) {
                // B at level l unchanged: fast exact nested delta
                lv.apply_last_virtual_move(v, t);
                double prev = (1.0 - p_new) * lv.move_prob(v, r, cl);
                // Recompute upper levels (for l < L - 1) without rebuilding level 0!
                double new_S;
                if (l == L - 1 && nw.level(0).normal_idx().empty()) {
                    new_S = cur_S + dSl;
                } else {
                    if (l < L - 1) nw.rebuild_upper_from(l);
                    new_S = nw.entropy();
                }
                double dS = new_S - cur_S;
                bool accept = false;
                if (std::isinf(beta)) {
                    accept = (dS < 0.0);
                } else if (prev > 0.0) {
                    double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) {
                    cur_S = new_S;
                    total_dS += dS;
                    ++nmoves;
                } else {
                    nw.level(l).move(v, r);
                    if (l < L - 1) nw.rebuild_upper_from(l);
                }
            } else if (l == 0) {
                // Creating or emptying a level-0 block: preserve upper hierarchy via parent copy/remove
                std::vector<std::vector<int>> old_bs = nw.save_bs();
                // Map r's parent at level 1 before move
                Rcpp::IntegerVector eb_before = lv.external_blocks_r();
                int r_ext = eb_before[v] - 1;
                int parent_r = (L > 1 && r_ext < static_cast<int>(old_bs[1].size()))
                                   ? old_bs[1][r_ext] : 0;

                lv.apply_last_virtual_move(v, t);
                Rcpp::IntegerVector eb_after = lv.external_blocks_r();
                std::vector<int> new_b0(eb_after.size());
                for (int i = 0; i < eb_after.size(); ++i) new_b0[i] = eb_after[i] - 1;

                // Map old level-0 external groups to new level-0 external groups
                int new_B0 = lv.B();
                std::vector<int> new_b1(new_B0, parent_r);
                if (L > 1) {
                    for (int u = 0; u < lv.N(); ++u) {
                        if (u == v && creates) continue;
                        int old_g = eb_before[u] - 1;
                        int new_g = new_b0[u];
                        if (old_g < static_cast<int>(old_bs[1].size()))
                            new_b1[new_g] = old_bs[1][old_g];
                    }
                    // Ensure new_b1 is contiguous 0..B1-1
                    lt_hash_map<int, int> rm;
                    int nxt = 0;
                    for (int& x : new_b1) {
                        auto it = rm.find(x);
                        if (it == rm.end()) it = rm.emplace(x, nxt++).first;
                        x = it->second;
                    }
                }
                std::vector<std::vector<int>> cand_bs;
                cand_bs.push_back(new_b0);
                if (L > 1) {
                    cand_bs.push_back(new_b1);
                    int prev_B = *std::max_element(new_b1.begin(), new_b1.end()) + 1;
                    for (int ll = 2; ll < L; ++ll) {
                        if (prev_B == static_cast<int>(old_bs[ll].size())) {
                            cand_bs.push_back(old_bs[ll]);
                            prev_B = *std::max_element(old_bs[ll].begin(), old_bs[ll].end()) + 1;
                        } else {
                            cand_bs.push_back(std::vector<int>(prev_B, 0));
                            break;
                        }
                    }
                }
                if (cand_bs.back().size() > 1 &&
                    *std::max_element(cand_bs.back().begin(), cand_bs.back().end()) > 0) {
                    int topB = *std::max_element(cand_bs.back().begin(), cand_bs.back().end()) + 1;
                    cand_bs.push_back(std::vector<int>(topB, 0));
                }

                nw.build(cand_bs);
                double new_S = nw.entropy();
                double dS = new_S - cur_S;
                double prev = 0.0;
                if (empties) {
                    prev = (d > 0.0) ? d : 0.0;
                } else {
                    size_t B_after = static_cast<size_t>(nw.level(0).B());
                    bool rev_can = (B_after < static_cast<size_t>(nw.level(0).N()) &&
                                    B_after < l_bmax && d > 0.0);
                    double p_new_rev = rev_can ? d : 0.0;
                    int r_after = nw.level(0).group(v == 0 ? 1 : 0);
                    // Find block holding former neighbors of r
                    for (int u = 0; u < nw.level(0).N(); ++u) {
                        if (u != v && eb_before[u] == eb_before[v]) {
                            r_after = nw.level(0).group(u);
                            break;
                        }
                    }
                    prev = (1.0 - p_new_rev) * nw.level(0).move_prob(v, r_after, cl);
                }
                bool accept = false;
                if (std::isinf(beta)) {
                    accept = (dS < 0.0);
                } else if (prev > 0.0) {
                    double log_a = -beta * dS + std::log(prev) - std::log(pfwd);
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) {
                    cur_S = new_S;
                    total_dS += dS;
                    ++nmoves;
                    L = nw.num_levels();
                    break;  // level count or size changed
                } else {
                    nw.build(old_bs);
                    L = nw.num_levels();
                }
            }
        }
    }
}

inline Rcpp::NumericMatrix unwrap_weights(Rcpp::Nullable<Rcpp::NumericMatrix> weights) {
    return weights.isNotNull() ? Rcpp::NumericMatrix(weights.get()) : Rcpp::NumericMatrix(0, 0);
}
inline Rcpp::CharacterVector unwrap_types(Rcpp::Nullable<Rcpp::CharacterVector> types) {
    return types.isNotNull() ? Rcpp::CharacterVector(types.get()) : Rcpp::CharacterVector();
}
inline Rcpp::List unwrap_params(Rcpp::Nullable<Rcpp::List> params) {
    return params.isNotNull() ? Rcpp::List(params.get()) : Rcpp::List();
}

}  // namespace

// ---------------------------------------------------------------------------
// Rcpp exports

// [[Rcpp::export]]
double cpp_layered_entropy(SEXP pointer, Rcpp::IntegerVector blocks,
                           Rcpp::IntegerVector layers, bool independent,
                           bool degree_corrected, bool distributed,
                           bool include_lmemb, double q_max,
                           Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                           Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                           Rcpp::Nullable<Rcpp::List> params = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    int C = 1;
    auto edges = build_base_layered_edges(
        graph, layers, W, static_cast<int>(covs.size()), C);
    LayeredLevel lv(N, C, graph.directed, independent, degree_corrected,
                    distributed, true, true, covs);
    lv.set_edges(std::move(edges));
    lv.set_partition(parse_l_blocks(blocks, N));
    return lv.entropy(include_lmemb);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_layered_check_moves(SEXP pointer,
                                            Rcpp::IntegerVector blocks,
                                            Rcpp::IntegerVector layers,
                                            bool independent,
                                            bool degree_corrected,
                                            bool distributed, int nmoves,
                                            double q_max,
                                            Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                            Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                            Rcpp::Nullable<Rcpp::List> params = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    int C = 1;
    auto edges = build_base_layered_edges(
        graph, layers, W, static_cast<int>(covs.size()), C);
    LayeredLevel lv(N, C, graph.directed, independent, degree_corrected,
                    distributed, true, true, covs);
    lv.set_edges(std::move(edges));
    lv.set_partition(parse_l_blocks(blocks, N));
    LRng rng;
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
        err[i] = std::max(std::abs(dS_v - (S_new - S)), std::abs(dS_a - (S_new - S)));
        S = S_new;
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_layered_minimize(SEXP pointer, Rcpp::IntegerVector layers,
                                bool independent, bool degree_corrected,
                                bool distributed, int B_min, double B_max,
                                bool refine, double epsilon, double q_max,
                                Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                                int nthreads = 0) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    int C = 1;
    auto edges = build_base_layered_edges(
        graph, layers, W, static_cast<int>(covs.size()), C);
    LayeredLevel lv(N, C, graph.directed, independent, degree_corrected,
                    distributed, true, true, covs);
    lv.set_edges(std::move(edges));
    LRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    double S = layered_minimize_flat(lv, rng, static_cast<size_t>(B_min), bmax,
                                     refine, epsilon, nullptr, nthreads);
    return Rcpp::List::create(
        Rcpp::Named("blocks") = lv.external_blocks_r(),
        Rcpp::Named("entropy") = S,
        Rcpp::Named("B") = lv.B());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_nested_layered_entropy(SEXP pointer, Rcpp::List bs,
                                               Rcpp::IntegerVector layers,
                                               bool independent,
                                               bool degree_corrected,
                                               bool distributed, double q_max,
                                               Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                               Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                               Rcpp::Nullable<Rcpp::List> params = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    NestedLayeredSBM nw(graph, layers, independent, degree_corrected, distributed,
                        W, covs);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_l_blocks(x, cur_n);
        cur_n = *std::max_element(bl.begin(), bl.end()) + 1;
        b_vec.push_back(std::move(bl));
    }
    nw.build(b_vec);
    auto lvS = nw.level_entropies();
    return Rcpp::NumericVector(lvS.begin(), lvS.end());
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_layered_minimize(SEXP pointer, Rcpp::IntegerVector layers,
                                       bool independent, bool degree_corrected,
                                       bool distributed, int B_min, double B_max,
                                       bool refine, double epsilon, double q_max,
                                       Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                       Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                       Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                                       int nthreads = 0) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    int C = 1;
    auto edges = build_base_layered_edges(
        graph, layers, W, static_cast<int>(covs.size()), C);
    // Sweep base level without flat top-level edges_dl penalty so nested-scale B candidates are explored
    LayeredLevel lv0(N, C, graph.directed, independent, degree_corrected,
                     distributed, true, false, covs);
    lv0.set_edges(std::move(edges));
    LRng rng;
    size_t bmax = std::isinf(B_max) ? std::numeric_limits<size_t>::max()
                                    : static_cast<size_t>(B_max);
    std::vector<std::vector<int>> cps;
    layered_minimize_flat(lv0, rng, static_cast<size_t>(B_min), bmax, refine, epsilon, &cps, nthreads);

    NestedLayeredSBM best_nw(graph, layers, independent, degree_corrected, distributed,
                             W, covs);
    double best_S = L_INF;
    for (const auto& b0 : cps) {
        NestedLayeredSBM cand(graph, layers, independent, degree_corrected, distributed,
                              W, covs);
        cand.optimize_from_base(b0, rng, false, epsilon, nthreads);
        double s_cand = cand.entropy();
        if (s_cand < best_S) {
            best_S = s_cand;
            best_nw = std::move(cand);
        }
    }
    if (refine)
        best_nw.refine_joint(rng, static_cast<size_t>(B_min), bmax, epsilon, nthreads);

    auto lvS = best_nw.level_entropies();
    return Rcpp::List::create(
        Rcpp::Named("bs") = best_nw.external_bs_r(),
        Rcpp::Named("level_entropy") = lvS,
        Rcpp::Named("entropy") = std::accumulate(lvS.begin(), lvS.end(), 0.0));
}

// [[Rcpp::export]]
Rcpp::List cpp_layered_mcmc(SEXP pointer, Rcpp::IntegerVector blocks,
                            Rcpp::IntegerVector layers, bool independent,
                            bool degree_corrected, bool distributed,
                            Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                            double q_max,
                            Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                            Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                            Rcpp::Nullable<Rcpp::List> params = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    int C = 1;
    auto edges = build_base_layered_edges(
        graph, layers, W, static_cast<int>(covs.size()), C);
    LayeredLevel lv(N, C, graph.directed, independent, degree_corrected,
                    distributed, true, true, covs);
    lv.set_edges(std::move(edges));
    lv.set_partition(parse_l_blocks(blocks, N));
    LRng rng;

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
            layered_flat_mcmc_sweep(lv, rng, beta, c, d, multiflip, B_min, B_max,
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
                layered_flat_mcmc_sweep(lv, rng, cur_beta, c, d, multiflip,
                                        B_min, B_max, total_dS, nattempts, nmoves);
            }
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
        Rcpp::Named("blocks") = lv.external_blocks_r(),
        Rcpp::Named("entropy") = lv.entropy(),
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
Rcpp::List cpp_nested_layered_mcmc(SEXP pointer, Rcpp::List bs,
                                   Rcpp::IntegerVector layers, bool independent,
                                   bool degree_corrected, bool distributed,
                                   Rcpp::List mcmc_opts, Rcpp::List eq_opts,
                                   double q_max,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                   Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                   Rcpp::Nullable<Rcpp::List> params = R_NilValue) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) llog_q().init(l_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_weights(weights);
    Rcpp::CharacterVector T = unwrap_types(types);
    Rcpp::List P = unwrap_params(params);
    auto covs = parse_covs(W, T, P);
    NestedLayeredSBM nw(graph, layers, independent, degree_corrected, distributed,
                        W, covs);
    std::vector<std::vector<int>> b_vec;
    b_vec.reserve(bs.size());
    int cur_n = N;
    for (R_xlen_t l = 0; l < bs.size(); ++l) {
        Rcpp::IntegerVector x = bs[l];
        auto bl = parse_l_blocks(x, cur_n);
        cur_n = *std::max_element(bl.begin(), bl.end()) + 1;
        b_vec.push_back(std::move(bl));
    }
    nw.build(b_vec);
    LRng rng;

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
            layered_nested_mcmc_sweep(nw, rng, beta, c, d, B_min, B_max,
                                      total_dS, nattempts, nmoves);
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

    auto run_eq_loop = [&](double cur_beta, int max_it, int f_it) {
        double best_S = nw.entropy();
        double worst_S = best_S;
        int no_break = 0;
        int breaks_left = nbreaks;
        for (int step = 1; step <= max_it; ++step) {
            for (int sub = 0; sub < std::max(1, niter); ++sub) {
                layered_nested_mcmc_sweep(nw, rng, cur_beta, c, d, B_min, B_max,
                                          total_dS, nattempts, nmoves);
            }
            double S_now = nw.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_now);
            if (collect && (step % thin == 0)) {
                samples.push_back(nw.external_bs_r());
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
        Rcpp::Named("bs") = nw.external_bs_r(),
        Rcpp::Named("entropy") = nw.entropy(),
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
