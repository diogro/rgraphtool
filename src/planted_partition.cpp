// Native C++ Bayesian Planted-Partition (assortative) stochastic block model
// engine matching graph-tool 3.9's PPBlockState to machine precision.
//
// References:
//   Zhang, L. & Peixoto, T. P. (2020). Statistical inference of assortative
//   community structures. Physical Review Research 2, 043271.
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
#include "weight_cov.h"
#include "logq.h"

namespace {

using namespace rgt_cov;

struct PPAdj {
    int nbr;
    int eidx;
    int layer;
};

struct PPLayerCountVec {
    std::vector<std::pair<int, long>> items;

    void clear() { items.clear(); }
    bool empty() const { return items.empty(); }

    long get(int l) const {
        auto it = std::lower_bound(
            items.begin(), items.end(), l,
            [](const std::pair<int, long>& p, int key) { return p.first < key; });
        return (it != items.end() && it->first == l) ? it->second : 0;
    }

    void add(int l, long d) {
        if (d == 0) return;
        auto it = std::lower_bound(
            items.begin(), items.end(), l,
            [](const std::pair<int, long>& p, int key) { return p.first < key; });
        if (it == items.end() || it->first != l) {
            if (d > 0) items.emplace(it, l, d);
        } else {
            it->second += d;
            if (it->second <= 0) items.erase(it);
        }
    }
};

struct SelfStratumW {
    int id;
    WVec w;
};

inline Rcpp::NumericMatrix unwrap_pp_weights(Rcpp::Nullable<Rcpp::NumericMatrix> weights) {
    return weights.isNotNull() ? Rcpp::NumericMatrix(weights.get()) : Rcpp::NumericMatrix(0, 0);
}
inline Rcpp::CharacterVector unwrap_pp_types(Rcpp::Nullable<Rcpp::CharacterVector> types) {
    return types.isNotNull() ? Rcpp::CharacterVector(types.get()) : Rcpp::CharacterVector();
}
inline Rcpp::List unwrap_pp_params(Rcpp::Nullable<Rcpp::List> params) {
    return params.isNotNull() ? Rcpp::List(params.get()) : Rcpp::List();
}
inline Rcpp::IntegerVector unwrap_pp_layers(Rcpp::Nullable<Rcpp::IntegerVector> layers) {
    return layers.isNotNull() ? Rcpp::IntegerVector(layers.get()) : Rcpp::IntegerVector();
}

constexpr double PP_INF = std::numeric_limits<double>::infinity();
constexpr double PP_LOG2 = 0.693147180559945309417232121458176568;

template <class K, class V>
using pp_hash_map = std::unordered_map<K, V>;
template <class K>
using pp_hash_set = std::unordered_set<K>;

struct PPMath {
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

PPMath& ppmath() {
    static PPMath m;
    return m;
}
// Shared log q(n, k) table (src/logq.h).
using PPLogQ = rgt_logq::LogQ;
inline PPLogQ& pplog_q() { return rgt_logq::shared(); }

size_t pp_q_cache_max(double q_max) {
    if (!std::isfinite(q_max) || q_max <= 0.0) return 10000;
    return std::min<size_t>(static_cast<size_t>(q_max), 10000);
}

struct PPRng {
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

inline uint64_t pp_canon(int a, int b, bool directed) {
    if (!directed && a > b) std::swap(a, b);
    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(b));
}

// ---------------------------------------------------------------------------
// Native Planted-Partition state (O(B) memory, zero-allocation virtual_move!)

class PPState {
public:
    PPState(const Graph& graph, bool uniform, bool dc, bool dist,
            const Rcpp::NumericMatrix& weights = Rcpp::NumericMatrix(0, 0),
            const std::vector<CovSpec>& covs = {},
            const Rcpp::IntegerVector& layers = Rcpp::IntegerVector(),
            bool independent = false)
        : _N(static_cast<int>(boost::num_vertices(graph.adjacency))),
          _directed(graph.directed), _uniform(uniform), _dc(dc), _dist(dist),
          _covs(covs), _K(static_cast<int>(covs.size())),
          _cap(std::max(2, _N)),
          _adj(_N), _kout(_N, 0), _kin(graph.directed ? _N : 0, 0),
          _self_loops(_N, 0), _self_loops_l(_N),
          _b(_N, 0), _wr(_cap, 0), _ep(_cap, 0), _em(graph.directed ? _cap : 0, 0),
          _err(_cap, 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(_N, 0),
          _h_out(!independent && dc && dist ? _cap : 0),
          _h_in(!independent && dc && dist && graph.directed ? _cap : 0) {
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);

        long E_total = static_cast<long>(boost::num_edges(graph.adjacency));
        bool has_layers = (layers.size() > 0);
        if (has_layers) {
            if (static_cast<long>(layers.size()) != E_total)
                Rcpp::stop("layers must have length graph_num_edges(graph)");
            int max_l = 0;
            for (long i = 0; i < E_total; ++i) {
                int lv = layers[static_cast<R_xlen_t>(i)];
                if (lv == NA_INTEGER || lv < 1)
                    Rcpp::stop("layers must contain positive integers >= 1");
                if (lv > max_l) max_l = lv;
            }
            _C = std::max(1, max_l);
            _independent = independent;
        } else {
            _C = 1;
            _independent = independent && (E_total > 0);
        }
        // The degree histograms were sized from the requested `independent`;
        // resize them to the resolved mode (e.g. independent requested on an
        // edgeless graph falls back to pooled degrees).
        if (!_independent && _dc && _dist) {
            _h_out.resize(_cap);
            if (_directed) _h_in.resize(_cap);
        }
        _E_l.assign(_C, 0);
        _e_in_l.assign(_C, 0);
        _e_out_l.assign(_C, 0);
        _err_l.assign(_C, std::vector<long>(_cap, 0));
        _dm_r_l.assign(_C, 0);
        _dm_s_l.assign(_C, 0);
        _touched_layers.reserve(_C);
        _touched_l_flag.assign(_C, 0);

        if (_independent) {
            _v_layers.resize(_N);
            _v_kout_l.resize(_N);
            if (_directed) _v_kin_l.resize(_N);
            _wr_l.assign(_C, std::vector<int>(_cap, 0));
            _ep_l.assign(_C, std::vector<long>(_cap, 0));
            if (_directed) _em_l.assign(_C, std::vector<long>(_cap, 0));
            if (_dc && _dist) {
                _h_out_l.assign(_C, std::vector<pp_hash_map<long, long>>(_cap));
                if (_directed) _h_in_l.assign(_C, std::vector<pp_hash_map<long, long>>(_cap));
            }
        }

        _num_cov_strata = (_K == 0) ? 0 : (_independent ? _C : 1);
        if (_K > 0) {
            if (weights.nrow() != E_total || weights.ncol() != _K)
                Rcpp::stop("weights dimensions do not match graph edges and types");
            for (int k = 0; k < _K; ++k) {
                if (_covs[k].type == CovType::REAL_NORMAL) _normal_idx.push_back(k);
            }
            _ew_x.assign(static_cast<size_t>(E_total) * _K, 0.0);
            _ew_x2.assign(static_cast<size_t>(E_total) * _K, 0.0);
            _self_w_s.resize(_N);
            _self_w_l.resize(_N);
            _w_in_s.assign(_num_cov_strata, WVec(_K));
            _w_out_s.assign(_num_cov_strata, WVec(_K));
            if (!_uniform) _w_rr_s.assign(_num_cov_strata, std::vector<WVec>(_cap, WVec(_K)));
            _w_in_l.assign(_C, WVec(_K));
            _w_out_l.assign(_C, WVec(_K));
            _B_E_c.assign(_num_cov_strata, 0);
            _B_E_D_c.assign(_num_cov_strata, 0);
            _recsum_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
            _recx2_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
            _recdx_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));

            _w_mr_s.assign(_num_cov_strata, WVec(_K));
            _w_ms_s.assign(_num_cov_strata, WVec(_K));
            _w_mr_l.assign(_C, WVec(_K));
            _w_ms_l.assign(_C, WVec(_K));
            _last_dw_r_s.assign(_num_cov_strata, WVec(_K));
            _last_dw_s_s.assign(_num_cov_strata, WVec(_K));
            _last_dw_in_s.assign(_num_cov_strata, WVec(_K));
            _last_dw_out_s.assign(_num_cov_strata, WVec(_K));
            _last_dw_in_l.assign(_C, WVec(_K));
            _last_dw_out_l.assign(_C, WVec(_K));
            _touched_strata.reserve(_num_cov_strata);
            _touched_s_flag.assign(_num_cov_strata, 0);
            _last_d_B_E_c.assign(_num_cov_strata, 0);
            _last_d_B_E_D_c.assign(_num_cov_strata, 0);
            _last_d_recx2_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
            _last_d_recdx_c.assign(_num_cov_strata, std::vector<double>(_K, 0.0));
        }

        _E = 0;
        pp_hash_map<uint64_t, PPLayerCountVec> emult;
        int ei = 0;
        for (auto [it, end] = boost::edges(graph.adjacency); it != end; ++it, ++ei) {
            int u = static_cast<int>(boost::source(*it, graph.adjacency));
            int v = static_cast<int>(boost::target(*it, graph.adjacency));
            int l = has_layers ? (layers[static_cast<R_xlen_t>(ei)] - 1) : 0;
            int st = _independent ? l : 0;
            ++_E;
            ++_E_l[l];
            emult[pp_canon(u, v, _directed)].add(l, 1);
            if (_K > 0) {
                for (int k = 0; k < _K; ++k) {
                    double val = weights(static_cast<R_xlen_t>(ei), k);
                    _recsum_c[st][k] += val;
                    if (u == v) {
                        add_self_entry(_self_w_s[u], st, k, val, _K);
                        add_self_entry(_self_w_l[u], l, k, val, _K);
                    } else {
                        size_t slot = static_cast<size_t>(ei) * _K + k;
                        _ew_x[slot] = val;
                        _ew_x2[slot] = val * val;
                    }
                }
            }
            if (u == v) {
                _self_loops[u]++;
                _self_loops_l[u].add(l, 1);
                if (_directed) {
                    _kout[u]++;
                    _kin[u]++;
                    if (_independent) {
                        _v_kout_l[u].add(l, 1);
                        _v_kin_l[u].add(l, 1);
                    }
                } else {
                    _kout[u] += 2;
                    if (_independent) _v_kout_l[u].add(l, 2);
                }
            } else {
                _adj[u].push_back({v, ei, l});
                _adj[v].push_back({u, ei, l});
                if (_directed) {
                    _kout[u]++;
                    _kin[v]++;
                    if (_independent) {
                        _v_kout_l[u].add(l, 1);
                        _v_kin_l[v].add(l, 1);
                    }
                } else {
                    _kout[u]++;
                    _kout[v]++;
                    if (_independent) {
                        _v_kout_l[u].add(l, 1);
                        _v_kout_l[v].add(l, 1);
                    }
                }
            }
        }

        _active_layers.reserve(_C);
        for (int l = 0; l < _C; ++l) {
            if (_E_l[l] > 0) _active_layers.push_back(l);
        }
        _C_active = static_cast<int>(_active_layers.size());

        if (_independent) {
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
        }

        ppmath().ensure(static_cast<size_t>(_N + 2 * _E + 2 * _C + 64));
        if (_dc && _dist) pplog_q().init(rgt_logq::table_size(static_cast<double>(_E), 5000));

        const PPMath& M = ppmath();
        _base_const = 0.0;
        for (const auto& [pk, lcv] : emult) {
            int u = static_cast<int>(pk >> 32);
            int v = static_cast<int>(pk & 0xffffffffULL);
            for (const auto& [l, mult] : lcv.items) {
                (void)l;
                if (mult <= 0) continue;
                _base_const += M.lgam1(mult);
                if (!_directed && u == v)
                    _base_const += static_cast<double>(mult) * PP_LOG2;
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

    int N() const { return _N; }
    int C() const { return _C; }
    int C_active() const { return _C_active; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
    int K() const { return _K; }
    bool independent() const { return _independent; }
    long E_in() const { return _e_in; }
    long E_out() const { return _e_out; }
    int group(int v) const { return _b[v]; }
    int size(int r) const { return _wr[r]; }
    long err(int r) const { return _err[r]; }
    const std::vector<int>& occupied() const { return _occupied; }
    const std::vector<int>& members(int r) const { return _members[r]; }
    const std::vector<int>& blocks() const { return _b; }
    int empty_group() const { return _empty.empty() ? -1 : _empty.back(); }
    int random_occupied(PPRng& rng) const { return _occupied[rng.index(_occupied.size())]; }

    void set_partition(const std::vector<int>& b) {
        for (int r : _occupied) {
            _wr[r] = 0;
            _ep[r] = 0;
            if (_directed) _em[r] = 0;
            _err[r] = 0;
            for (int l = 0; l < _C; ++l) _err_l[l][r] = 0;
            if (_K > 0 && !_uniform) {
                for (int st = 0; st < _num_cov_strata; ++st) _w_rr_s[st][r].zero(_K);
            }
            if (_independent) {
                for (int l = 0; l < _C; ++l) {
                    _wr_l[l][r] = 0;
                    _ep_l[l][r] = 0;
                    if (_directed) _em_l[l][r] = 0;
                    if (_dc && _dist) {
                        _h_out_l[l][r].clear();
                        if (_directed) _h_in_l[l][r].clear();
                    }
                }
            }
            _members[r].clear();
            _occ_pos[r] = -1;
            if (!_independent && _dc && _dist) {
                _h_out[r].clear();
                if (_directed) _h_in[r].clear();
            }
        }
        _occupied.clear();
        _empty.clear();
        _empty_pos.assign(_cap, -1);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);
        _e_in = 0;
        _e_out = 0;
        _e_in_l.assign(_C, 0);
        _e_out_l.assign(_C, 0);
        if (_K > 0) {
            for (int st = 0; st < _num_cov_strata; ++st) {
                _w_in_s[st].zero(_K);
                _w_out_s[st].zero(_K);
            }
            for (int l = 0; l < _C; ++l) {
                _w_in_l[l].zero(_K);
                _w_out_l[l].zero(_K);
            }
        }

        for (int v = 0; v < _N; ++v) {
            int r = b[v];
            _b[v] = r;
            add_member(v, r);
            _ep[r] += _kout[v];
            if (_directed) _em[r] += _kin[v];
            if (_independent) {
                for (int l : _v_layers[v]) {
                    _wr_l[l][r]++;
                    long ko_vl = _v_kout_l[v].get(l);
                    _ep_l[l][r] += ko_vl;
                    if (_dc && _dist) _h_out_l[l][r][ko_vl]++;
                    if (_directed) {
                        long ki_vl = _v_kin_l[v].get(l);
                        _em_l[l][r] += ki_vl;
                        if (_dc && _dist) _h_in_l[l][r][ki_vl]++;
                    }
                }
            } else if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                if (_directed) _h_in[r][_kin[v]]++;
            }
            if (_self_loops[v] > 0) {
                _err[r] += _self_loops[v];
                _e_in += _self_loops[v];
                for (const auto& [l, sl_l] : _self_loops_l[v].items) {
                    _err_l[l][r] += sl_l;
                    _e_in_l[l] += sl_l;
                }
                if (_K > 0) {
                    for (const auto& sw : _self_w_s[v]) {
                        add_w(_w_in_s[sw.id], sw.w);
                        if (!_uniform) add_w(_w_rr_s[sw.id][r], sw.w);
                    }
                    for (const auto& sw : _self_w_l[v]) {
                        add_w(_w_in_l[sw.id], sw.w);
                    }
                }
            }
        }
        // Count non-self edges (each appears twice in _adj)
        long in_twice = 0, out_twice = 0;
        for (int v = 0; v < _N; ++v) {
            int rv = _b[v];
            for (const auto& a : _adj[v]) {
                int l = a.layer;
                int st = _independent ? l : 0;
                if (_b[a.nbr] == rv) {
                    if (v < a.nbr) {
                        _err[rv]++;
                        _err_l[l][rv]++;
                        _e_in_l[l]++;
                        in_twice += 2;
                        if (_K > 0) {
                            add_edge_w(_w_in_s[st], a.eidx);
                            add_edge_w(_w_in_l[l], a.eidx);
                            if (!_uniform) add_edge_w(_w_rr_s[st][rv], a.eidx);
                        }
                    }
                } else {
                    out_twice++;
                    if (v < a.nbr) {
                        _e_out_l[l]++;
                        if (_K > 0) {
                            add_edge_w(_w_out_s[st], a.eidx);
                            add_edge_w(_w_out_l[l], a.eidx);
                        }
                    }
                }
            }
        }
        _e_in += in_twice / 2;
        _e_out = out_twice / 2;
        refresh_cov_globals();
    }

    void refresh_cov_globals() {
        if (_K == 0) return;
        for (int st = 0; st < _num_cov_strata; ++st) {
            _B_E_c[st] = 0;
            _B_E_D_c[st] = 0;
            _recx2_c[st].assign(_K, 0.0);
            _recdx_c[st].assign(_K, 0.0);
            auto accum_bucket = [&](long m, const WVec& w) {
                if (m <= 0) return;
                ++_B_E_c[st];
                if (m > 1) ++_B_E_D_c[st];
                for (int k : _normal_idx) {
                    double x = w.x(k);
                    _recx2_c[st][k] += x * x;
                    if (m > 1) {
                        double dx = w.x2(k) - x * (x / static_cast<double>(m));
                        if (dx > 0.0) _recdx_c[st][k] += dx;
                    }
                }
            };
            if (_uniform) {
                accum_bucket(_independent ? _e_in_l[st] : _e_in, _w_in_s[st]);
                accum_bucket(_independent ? _e_out_l[st] : _e_out, _w_out_s[st]);
            } else {
                for (int r : _occupied) {
                    long mr = _independent ? _err_l[st][r] : _err[r];
                    accum_bucket(mr, _w_rr_s[st][r]);
                }
                accum_bucket(_independent ? _e_out_l[st] : _e_out, _w_out_s[st]);
            }
        }
    }

    double lmemb_dl() const {
        if (!_independent || _C_active <= 1) return 0.0;
        const PPMath& M = ppmath();
        double x = M.lbinom(_C_active, 1);
        if (x > 15.0)
            return static_cast<double>(_C_active) * x - M.lgam1(_C_active);
        return M.lbinom(2 * _C_active - 1, _C_active);
    }

    double cov_entropy() const {
        if (_K == 0 || _E == 0) return 0.0;
        double S = 0.0;
        for (int st = 0; st < _num_cov_strata; ++st) {
            auto add_bucket_logp = [&](long m, const WVec& w) {
                if (m <= 0) return;
                double md = static_cast<double>(m);
                for (int k = 0; k < _K; ++k)
                    S += cov_pair_neg_log_P(_covs[k], md, w.x(k), w.x2(k));
            };
            if (_uniform) {
                add_bucket_logp(_independent ? _e_in_l[st] : _e_in, _w_in_s[st]);
                add_bucket_logp(_independent ? _e_out_l[st] : _e_out, _w_out_s[st]);
            } else {
                for (int r : _occupied) {
                    long mr = _independent ? _err_l[st][r] : _err[r];
                    add_bucket_logp(mr, _w_rr_s[st][r]);
                }
                add_bucket_logp(_independent ? _e_out_l[st] : _e_out, _w_out_s[st]);
            }
            for (int k = 0; k < _K; ++k) {
                const CovSpec& c = _covs[k];
                if (c.type == CovType::REAL_NORMAL) {
                    if (c.micro) {
                        S += -signed_w_log_P(static_cast<double>(_B_E_c[st]), _recsum_c[st][k],
                                             _recx2_c[st][k], c.p0, c.p1, c.p2, c.p3, c.eps);
                        S += -positive_w_log_P(static_cast<double>(_B_E_D_c[st]), _recdx_c[st][k],
                                               c.p2, c.p3, c.eps);
                    }
                } else if (c.micro) {
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

    double entropy() const {
        const PPMath& M = ppmath();
        long Bn = B();
        double S = _base_const;

        // Partition DL
        if (_N > 0) {
            S += M.logn(_N) + M.lbinom(_N - 1, Bn - 1) + M.lgam1(_N);
            for (int r : _occupied) S -= M.lgam1(_wr[r]);
        }

        // Layer-split term + lmemb_dl
        if (_E > 0 && _C_active > 1) S += M.lbinom(_E - 1, _C_active - 1);
        S += lmemb_dl();

        // Vertex/group degree + degree DL
        if (!_independent) {
            for (int r : _occupied) {
                if (_dc) {
                    S += M.lgam1(_ep[r]);
                    if (_directed) S += M.lgam1(_em[r]);
                    S += deg_dl_group(r);
                } else {
                    long dr = _ep[r] + (_directed ? _em[r] : 0);
                    S += static_cast<double>(dr) * M.logn(_wr[r]);
                }
            }
        } else {
            for (int l : _active_layers) {
                for (int r : _occupied) {
                    long nrl = _wr_l[l][r];
                    if (nrl <= 0) continue;
                    long epl = _ep_l[l][r];
                    long eml = _directed ? _em_l[l][r] : 0;
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

        // Per-layer edges DL + assortative PP edge likelihood
        for (int l : _active_layers) {
            long el = _E_l[l];
            if (Bn > 1 && el > 0) S += M.logn(el + 1);
            long ein_l = _e_in_l[l];
            long eout_l = _e_out_l[l];
            if (_uniform) {
                S -= M.lgam1(ein_l) + M.lgam1(eout_l);
                if (_directed) {
                    S += M.logn(Bn) * static_cast<double>(ein_l);
                    S += (M.lbinom(Bn, 2) + PP_LOG2) * static_cast<double>(eout_l);
                } else {
                    S += (M.logn(Bn) - PP_LOG2) * static_cast<double>(ein_l);
                    S += M.lbinom(Bn, 2) * static_cast<double>(eout_l);
                }
            } else {
                S -= M.lgam1(eout_l);
                if (_directed)
                    S += static_cast<double>(eout_l) * (M.lbinom(Bn, 2) + PP_LOG2);
                else
                    S += static_cast<double>(eout_l) * M.lbinom(Bn, 2);
                for (int r : _occupied) {
                    long err_lr = _err_l[l][r];
                    if (_directed) S -= M.lgam1(err_lr);
                    else S -= static_cast<double>(err_lr) * PP_LOG2 + M.lgam1(err_lr);
                }
                S += M.lbinom(Bn + ein_l - 1, ein_l);
            }
        }
        return S + cov_entropy();
    }

    // Exact O(deg(v)) incremental virtual_move(v, s)
    double virtual_move(int v, int s) const {
        int r = _b[v];
        if (r == s) return 0.0;

        for (int l : _touched_layers) {
            _dm_r_l[l] = 0;
            _dm_s_l[l] = 0;
            _touched_l_flag[l] = 0;
            if (_K > 0) {
                _w_mr_l[l].zero(_K);
                _w_ms_l[l].zero(_K);
            }
        }
        _touched_layers.clear();
        if (_K > 0) {
            for (int st : _touched_strata) {
                _w_mr_s[st].zero(_K);
                _w_ms_s[st].zero(_K);
                _touched_s_flag[st] = 0;
            }
            _touched_strata.clear();
        }

        auto mark_layer = [&](int l) {
            if (!_touched_l_flag[l]) {
                _touched_l_flag[l] = 1;
                _touched_layers.push_back(l);
            }
        };
        auto mark_stratum = [&](int st) {
            if (!_touched_s_flag[st]) {
                _touched_s_flag[st] = 1;
                _touched_strata.push_back(st);
            }
        };

        long m_r = 0, m_s = 0;
        for (const auto& a : _adj[v]) {
            int bu = _b[a.nbr];
            int l = a.layer;
            int st = _independent ? l : 0;
            if (bu == r) {
                mark_layer(l);
                ++m_r;
                ++_dm_r_l[l];
                if (_K > 0) {
                    mark_stratum(st);
                    add_edge_w(_w_mr_s[st], a.eidx);
                    add_edge_w(_w_mr_l[l], a.eidx);
                }
            } else if (bu == s) {
                mark_layer(l);
                ++m_s;
                ++_dm_s_l[l];
                if (_K > 0) {
                    mark_stratum(st);
                    add_edge_w(_w_ms_s[st], a.eidx);
                    add_edge_w(_w_ms_l[l], a.eidx);
                }
            }
        }
        if (!_uniform && !_self_loops_l[v].empty()) {
            for (const auto& [l, sl_l] : _self_loops_l[v].items) {
                (void)sl_l;
                mark_layer(l);
                if (_K > 0) mark_stratum(_independent ? l : 0);
            }
        }

        long sl = _self_loops[v];
        long derr_r = -(m_r + sl);
        long derr_s = +(m_s + sl);
        long de_in = m_s - m_r;
        long de_out = -de_in;
        _last_derr_r = derr_r;
        _last_derr_s = derr_s;
        _last_de_in = de_in;
        _last_de_out = de_out;

        double dS = delta_entropy(r, s, v);
        if (_K > 0) {
            for (int l : _touched_layers) {
                for (int k = 0; k < _K; ++k) {
                    double dwin = _w_ms_l[l].x(k) - _w_mr_l[l].x(k);
                    double dwin2 = _w_ms_l[l].x2(k) - _w_mr_l[l].x2(k);
                    _last_dw_in_l[l].x(k) = dwin;
                    _last_dw_in_l[l].x2(k) = dwin2;
                    _last_dw_out_l[l].x(k) = -dwin;
                    _last_dw_out_l[l].x2(k) = -dwin2;
                }
            }
            for (int st : _touched_strata) {
                const WVec* self_st = find_self_w(_self_w_s[v], st);
                for (int k = 0; k < _K; ++k) {
                    double win_d = _w_ms_s[st].x(k) - _w_mr_s[st].x(k);
                    double win2_d = _w_ms_s[st].x2(k) - _w_mr_s[st].x2(k);
                    _last_dw_in_s[st].x(k) = win_d;
                    _last_dw_in_s[st].x2(k) = win2_d;
                    _last_dw_out_s[st].x(k) = -win_d;
                    _last_dw_out_s[st].x2(k) = -win2_d;
                    if (!_uniform) {
                        double sw_x = self_st ? self_st->x(k) : 0.0;
                        double sw_x2 = self_st ? self_st->x2(k) : 0.0;
                        _last_dw_r_s[st].x(k) = -(_w_mr_s[st].x(k) + sw_x);
                        _last_dw_r_s[st].x2(k) = -(_w_mr_s[st].x2(k) + sw_x2);
                        _last_dw_s_s[st].x(k) = +(_w_ms_s[st].x(k) + sw_x);
                        _last_dw_s_s[st].x2(k) = +(_w_ms_s[st].x2(k) + sw_x2);
                    }
                }
            }
            dS += delta_cov(r, s, v);
        }
        return dS;
    }

    double move(int v, int s) {
        int r = _b[v];
        if (r == s) return 0.0;
        double dS = virtual_move(v, s);
        apply_last_virtual_move(v, s);
        return dS;
    }

    void apply_last_virtual_move(int v, int s) {
        int r = _b[v];
        if (r == s) return;
        if (_dc && _dist) {
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
        _ep[r] -= _kout[v];
        _ep[s] += _kout[v];
        if (_directed) {
            _em[r] -= _kin[v];
            _em[s] += _kin[v];
        }
        if (_independent) {
            for (int l : _v_layers[v]) {
                _wr_l[l][r]--;
                _wr_l[l][s]++;
                long ko = _v_kout_l[v].get(l);
                _ep_l[l][r] -= ko;
                _ep_l[l][s] += ko;
                if (_directed) {
                    long ki = _v_kin_l[v].get(l);
                    _em_l[l][r] -= ki;
                    _em_l[l][s] += ki;
                }
            }
        }

        _err[r] += _last_derr_r;
        _err[s] += _last_derr_s;
        _e_in += _last_de_in;
        _e_out += _last_de_out;

        for (int l : _touched_layers) {
            long dm_r = _dm_r_l[l];
            long dm_s = _dm_s_l[l];
            long sl_l = _self_loops_l[v].get(l);
            _err_l[l][r] += -(dm_r + sl_l);
            _err_l[l][s] += +(dm_s + sl_l);
            long dein_l = dm_s - dm_r;
            _e_in_l[l] += dein_l;
            _e_out_l[l] += -dein_l;
            if (_K > 0) {
                apply_bucket_delta(_w_in_l[l], _last_dw_in_l[l], _e_in_l[l]);
                apply_bucket_delta(_w_out_l[l], _last_dw_out_l[l], _e_out_l[l]);
            }
        }

        if (_K > 0) {
            for (int st : _touched_strata) {
                long min_new = _independent ? _e_in_l[st] : _e_in;
                long mout_new = _independent ? _e_out_l[st] : _e_out;
                apply_bucket_delta(_w_in_s[st], _last_dw_in_s[st], min_new);
                apply_bucket_delta(_w_out_s[st], _last_dw_out_s[st], mout_new);
                if (!_uniform) {
                    long mr_new = _independent ? _err_l[st][r] : _err[r];
                    long ms_new = _independent ? _err_l[st][s] : _err[s];
                    apply_bucket_delta(_w_rr_s[st][r], _last_dw_r_s[st], mr_new);
                    apply_bucket_delta(_w_rr_s[st][s], _last_dw_s_s[st], ms_new);
                }
                _B_E_c[st] += _last_d_B_E_c[st];
                _B_E_D_c[st] += _last_d_B_E_D_c[st];
                for (int k : _normal_idx) {
                    _recx2_c[st][k] = std::max(0.0, _recx2_c[st][k] + _last_d_recx2_c[st][k]);
                    _recdx_c[st][k] = std::max(0.0, _recdx_c[st][k] + _last_d_recdx_c[st][k]);
                }
            }
        }
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

    // Assortative neighbor-based ergodic proposal (strongly favors neighbor groups!)
    int sample_group(int v, double c, PPRng& rng) const {
        if (_occupied.empty()) return 0;
        const auto& adj = _adj[v];
        if (adj.empty() || std::isinf(c))
            return _occupied[rng.index(_occupied.size())];
        double B_act = static_cast<double>(_occupied.size());
        double deg = static_cast<double>(adj.size());
        double p_rand = (c * B_act) / (deg + c * B_act);
        if (rng.unif() < p_rand)
            return _occupied[rng.index(_occupied.size())];
        return _b[adj[rng.index(adj.size())].nbr];
    }

    double move_prob(int v, int s, double c) const {
        if (_occupied.empty()) return 0.0;
        const auto& adj = _adj[v];
        double B_act = static_cast<double>(_occupied.size());
        if (adj.empty() || std::isinf(c)) return 1.0 / B_act;
        long cnt = 0;
        for (const auto& a : adj) if (_b[a.nbr] == s) ++cnt;
        double deg = static_cast<double>(adj.size());
        return (static_cast<double>(cnt) + c) / (deg + c * B_act);
    }

    Rcpp::NumericVector w_in_mean_r() const {
        Rcpp::NumericVector out(_K);
        for (int k = 0; k < _K; ++k)
            out[k] = (_e_in > 0) ? (_w_in_s_tot(true, k) / static_cast<double>(_e_in)) : NA_REAL;
        return out;
    }

    Rcpp::NumericVector w_out_mean_r() const {
        Rcpp::NumericVector out(_K);
        for (int k = 0; k < _K; ++k)
            out[k] = (_e_out > 0) ? (_w_in_s_tot(false, k) / static_cast<double>(_e_out)) : NA_REAL;
        return out;
    }

    Rcpp::NumericVector e_in_l_r() const {
        Rcpp::NumericVector out(_C);
        for (int l = 0; l < _C; ++l) out[l] = static_cast<double>(_e_in_l[l]);
        return out;
    }

    Rcpp::NumericVector e_out_l_r() const {
        Rcpp::NumericVector out(_C);
        for (int l = 0; l < _C; ++l) out[l] = static_cast<double>(_e_out_l[l]);
        return out;
    }

    Rcpp::NumericVector e_tot_l_r() const {
        Rcpp::NumericVector out(_C);
        for (int l = 0; l < _C; ++l) out[l] = static_cast<double>(_E_l[l]);
        return out;
    }

    Rcpp::NumericMatrix external_err_l_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        Rcpp::NumericMatrix out(static_cast<int>(occ.size()), _C);
        for (int l = 0; l < _C; ++l) {
            for (size_t i = 0; i < occ.size(); ++i)
                out(static_cast<R_xlen_t>(i), l) = static_cast<double>(_err_l[l][occ[i]]);
        }
        return out;
    }

    Rcpp::NumericMatrix w_in_mean_l_r() const {
        Rcpp::NumericMatrix out(_C, _K);
        for (int l = 0; l < _C; ++l) {
            for (int k = 0; k < _K; ++k)
                out(l, k) = (_e_in_l[l] > 0) ? (_w_in_l[l].x(k) / static_cast<double>(_e_in_l[l])) : NA_REAL;
        }
        return out;
    }

    Rcpp::NumericMatrix w_out_mean_l_r() const {
        Rcpp::NumericMatrix out(_C, _K);
        for (int l = 0; l < _C; ++l) {
            for (int k = 0; k < _K; ++k)
                out(l, k) = (_e_out_l[l] > 0) ? (_w_out_l[l].x(k) / static_cast<double>(_e_out_l[l])) : NA_REAL;
        }
        return out;
    }

    Rcpp::IntegerVector external_blocks_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        pp_hash_map<int, int> gmap;
        for (size_t i = 0; i < occ.size(); ++i) gmap.emplace(occ[i], static_cast<int>(i) + 1);
        Rcpp::IntegerVector out(_N);
        for (int v = 0; v < _N; ++v) out[v] = gmap[_b[v]];
        return out;
    }

    Rcpp::NumericVector external_err_r() const {
        std::vector<int> occ = _occupied;
        std::sort(occ.begin(), occ.end());
        Rcpp::NumericVector out(occ.size());
        for (size_t i = 0; i < occ.size(); ++i) out[i] = static_cast<double>(_err[occ[i]]);
        return out;
    }

private:
    int _N, _C = 1, _C_active = 0;
    bool _directed, _uniform, _dc, _dist, _independent = false;
    std::vector<CovSpec> _covs;
    int _K = 0;
    int _num_cov_strata = 0;
    std::vector<int> _normal_idx;
    int _cap;
    long _E = 0, _e_in = 0, _e_out = 0;
    double _base_const = 0.0;

    std::vector<std::vector<PPAdj>> _adj;
    std::vector<long> _kout, _kin, _self_loops;
    std::vector<PPLayerCountVec> _self_loops_l;
    std::vector<long> _E_l, _e_in_l, _e_out_l;
    std::vector<int> _active_layers;
    std::vector<std::vector<int>> _v_layers;
    std::vector<PPLayerCountVec> _v_kout_l, _v_kin_l;

    std::vector<double> _ew_x, _ew_x2;
    std::vector<std::vector<SelfStratumW>> _self_w_s, _self_w_l;
    std::vector<int> _b, _wr;
    std::vector<std::vector<int>> _wr_l;
    std::vector<long> _ep, _em, _err;
    std::vector<std::vector<long>> _ep_l, _em_l, _err_l;
    std::vector<WVec> _w_in_s, _w_out_s, _w_in_l, _w_out_l;
    std::vector<std::vector<WVec>> _w_rr_s;
    std::vector<long> _B_E_c, _B_E_D_c;
    std::vector<std::vector<double>> _recsum_c, _recx2_c, _recdx_c;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;
    std::vector<pp_hash_map<long, long>> _h_out, _h_in;
    std::vector<std::vector<pp_hash_map<long, long>>> _h_out_l, _h_in_l;

    mutable long _last_derr_r = 0, _last_derr_s = 0, _last_de_in = 0, _last_de_out = 0;
    mutable std::vector<long> _dm_r_l, _dm_s_l;
    mutable std::vector<int> _touched_layers, _touched_strata;
    mutable std::vector<uint8_t> _touched_l_flag, _touched_s_flag;
    mutable std::vector<WVec> _w_mr_s, _w_ms_s, _w_mr_l, _w_ms_l;
    mutable std::vector<WVec> _last_dw_r_s, _last_dw_s_s, _last_dw_in_s, _last_dw_out_s;
    mutable std::vector<WVec> _last_dw_in_l, _last_dw_out_l;
    mutable std::vector<long> _last_d_B_E_c, _last_d_B_E_D_c;
    mutable std::vector<std::vector<double>> _last_d_recx2_c, _last_d_recdx_c;

    double _w_in_s_tot(bool in_bucket, int k) const {
        double s = 0.0;
        for (int l = 0; l < _C; ++l) s += in_bucket ? _w_in_l[l].x(k) : _w_out_l[l].x(k);
        return s;
    }

    static void add_self_entry(std::vector<SelfStratumW>& vec, int id, int k, double val, int K) {
        for (auto& e : vec) {
            if (e.id == id) {
                e.w.x(k) += val;
                e.w.x2(k) += val * val;
                return;
            }
        }
        SelfStratumW nw{id, WVec(K)};
        nw.w.x(k) = val;
        nw.w.x2(k) = val * val;
        vec.push_back(std::move(nw));
    }

    static const WVec* find_self_w(const std::vector<SelfStratumW>& vec, int id) {
        for (const auto& e : vec) {
            if (e.id == id) return &e.w;
        }
        return nullptr;
    }

    inline void add_w(WVec& dst, const WVec& src) const {
        for (int k = 0; k < _K; ++k) {
            dst.x(k) += src.x(k);
            dst.x2(k) += src.x2(k);
        }
    }

    inline void add_edge_w(WVec& dst, int eidx) const {
        size_t base = static_cast<size_t>(eidx) * _K;
        for (int k = 0; k < _K; ++k) {
            dst.x(k) += _ew_x[base + k];
            dst.x2(k) += _ew_x2[base + k];
        }
    }

    inline void apply_bucket_delta(WVec& dst, const WVec& dw, long m_new) const {
        if (m_new <= 0) {
            dst.zero(_K);
        } else {
            for (int k = 0; k < _K; ++k) {
                dst.x(k) += dw.x(k);
                dst.x2(k) += dw.x2(k);
            }
        }
    }

    double delta_cov(int r, int s, int v) const {
        (void)v;
        double dS = 0.0;
        for (int st : _touched_strata) {
            _last_d_B_E_c[st] = 0;
            _last_d_B_E_D_c[st] = 0;
            for (int k = 0; k < _K; ++k) {
                _last_d_recx2_c[st][k] = 0.0;
                _last_d_recdx_c[st][k] = 0.0;
            }
            auto eval_bucket = [&](long m0, const WVec& w0, long dm, const WVec& dw) {
                long m1 = m0 + dm;
                if (m0 <= 0 && m1 <= 0) return;
                _last_d_B_E_c[st] += (m1 > 0 ? 1 : 0) - (m0 > 0 ? 1 : 0);
                _last_d_B_E_D_c[st] += (m1 > 1 ? 1 : 0) - (m0 > 1 ? 1 : 0);
                for (int k = 0; k < _K; ++k) {
                    double x0 = (m0 > 0) ? w0.x(k) : 0.0;
                    double x20 = (m0 > 0) ? w0.x2(k) : 0.0;
                    double x1 = (m1 > 0) ? (x0 + dw.x(k)) : 0.0;
                    double x21 = (m1 > 0) ? (x20 + dw.x2(k)) : 0.0;
                    dS += cov_pair_neg_log_P(_covs[k], static_cast<double>(m1), x1, x21) -
                          cov_pair_neg_log_P(_covs[k], static_cast<double>(m0), x0, x20);
                    if (_covs[k].type == CovType::REAL_NORMAL) {
                        _last_d_recx2_c[st][k] += (m1 > 0 ? x1 * x1 : 0.0) - (m0 > 0 ? x0 * x0 : 0.0);
                        double dx0 = (m0 > 1) ? std::max(0.0, x20 - x0 * (x0 / static_cast<double>(m0))) : 0.0;
                        double dx1 = (m1 > 1) ? std::max(0.0, x21 - x1 * (x1 / static_cast<double>(m1))) : 0.0;
                        _last_d_recdx_c[st][k] += dx1 - dx0;
                    }
                }
            };
            if (_uniform) {
                long ein0 = _independent ? _e_in_l[st] : _e_in;
                long eout0 = _independent ? _e_out_l[st] : _e_out;
                long dein = _independent ? (_dm_s_l[st] - _dm_r_l[st]) : _last_de_in;
                eval_bucket(ein0, _w_in_s[st], dein, _last_dw_in_s[st]);
                eval_bucket(eout0, _w_out_s[st], -dein, _last_dw_out_s[st]);
            } else {
                long mr0 = _independent ? _err_l[st][r] : _err[r];
                long ms0 = _independent ? _err_l[st][s] : _err[s];
                long eout0 = _independent ? _e_out_l[st] : _e_out;
                long sl_st = _independent ? _self_loops_l[v].get(st) : _self_loops[v];
                long dm_r_st = _independent ? _dm_r_l[st] : (-_last_derr_r - sl_st);
                long dm_s_st = _independent ? _dm_s_l[st] : (_last_derr_s - sl_st);
                long dr_st = -(dm_r_st + sl_st);
                long ds_st = +(dm_s_st + sl_st);
                long dout_st = dm_r_st - dm_s_st;
                eval_bucket(mr0, _w_rr_s[st][r], dr_st, _last_dw_r_s[st]);
                eval_bucket(ms0, _w_rr_s[st][s], ds_st, _last_dw_s_s[st]);
                eval_bucket(eout0, _w_out_s[st], dout_st, _last_dw_out_s[st]);
            }
            long be0 = _B_E_c[st], be1 = be0 + _last_d_B_E_c[st];
            long bed0 = _B_E_D_c[st], bed1 = bed0 + _last_d_B_E_D_c[st];
            for (int k = 0; k < _K; ++k) {
                const CovSpec& c = _covs[k];
                if (c.type == CovType::REAL_NORMAL) {
                    if (c.micro) {
                        dS += -signed_w_log_P(static_cast<double>(be1), _recsum_c[st][k],
                                              std::max(0.0, _recx2_c[st][k] + _last_d_recx2_c[st][k]),
                                              c.p0, c.p1, c.p2, c.p3, c.eps) -
                              (-signed_w_log_P(static_cast<double>(be0), _recsum_c[st][k],
                                               _recx2_c[st][k], c.p0, c.p1, c.p2, c.p3, c.eps));
                        dS += -positive_w_log_P(static_cast<double>(bed1),
                                                std::max(0.0, _recdx_c[st][k] + _last_d_recdx_c[st][k]),
                                                c.p2, c.p3, c.eps) -
                              (-positive_w_log_P(static_cast<double>(bed0), _recdx_c[st][k],
                                                 c.p2, c.p3, c.eps));
                    }
                } else if (c.micro && be0 != be1) {
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
        return dS;
    }

    double deg_dl_group(int r) const {
        if (!_dc || _wr[r] == 0) return 0.0;
        const PPMath& M = ppmath();
        long n = _wr[r];
        if (!_dist) {
            double s = M.lbinom(n + _ep[r] - 1, _ep[r]);
            if (_directed) s += M.lbinom(n + _em[r] - 1, _em[r]);
            return s;
        }
        const PPLogQ& q = pplog_q();
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
        const PPMath& M = ppmath();
        if (!_dist) {
            double s = M.lbinom(n + epl - 1, epl);
            if (_directed) s += M.lbinom(n + eml - 1, eml);
            return s;
        }
        const PPLogQ& q = pplog_q();
        double s = q(epl, n) + M.lgam1(n);
        for (const auto& [dk, cnt] : _h_out_l[l][r]) s -= M.lgam1(cnt);
        if (_directed) {
            s += q(eml, n) + M.lgam1(n);
            for (const auto& [dk, cnt] : _h_in_l[l][r]) s -= M.lgam1(cnt);
        }
        return s;
    }

    double deg_dl_delta_maps(const pp_hash_map<long, long>* ho,
                             const pp_hash_map<long, long>* hi,
                             long n0, long n1, long ep0, long ep1,
                             long em0, long em1, long ko_v, long ki_v, long dcnt) const {
        if (!_dc) return 0.0;
        const PPMath& M = ppmath();
        if (!_dist) {
            double s0 = (n0 > 0) ? (M.lbinom(n0 + ep0 - 1, ep0) +
                                    (_directed ? M.lbinom(n0 + em0 - 1, em0) : 0.0)) : 0.0;
            double s1 = (n1 > 0) ? (M.lbinom(n1 + ep1 - 1, ep1) +
                                    (_directed ? M.lbinom(n1 + em1 - 1, em1) : 0.0)) : 0.0;
            return s1 - s0;
        }
        const PPLogQ& q = pplog_q();
        double s0 = (n0 > 0) ? (q(ep0, n0) + M.lgam1(n0)) : 0.0;
        double s1 = (n1 > 0) ? (q(ep1, n1) + M.lgam1(n1)) : 0.0;
        auto ito = ho->find(ko_v);
        long co0 = (ito == ho->end()) ? 0 : ito->second;
        double d = (s1 - s0) - (M.lgam1(co0 + dcnt) - M.lgam1(co0));
        if (_directed) {
            double si0 = (n0 > 0) ? (q(em0, n0) + M.lgam1(n0)) : 0.0;
            double si1 = (n1 > 0) ? (q(em1, n1) + M.lgam1(n1)) : 0.0;
            auto iti = hi->find(ki_v);
            long ci0 = (iti == hi->end()) ? 0 : iti->second;
            d += (si1 - si0) - (M.lgam1(ci0 + dcnt) - M.lgam1(ci0));
        }
        return d;
    }

    double delta_entropy(int r, int s, int v) const {
        const PPMath& M = ppmath();
        long wr_r0 = _wr[r], wr_r1 = wr_r0 - 1;
        long wr_s0 = _wr[s], wr_s1 = wr_s0 + 1;
        long B0 = B();
        long B1 = B0 - (wr_r1 == 0 ? 1 : 0) + (wr_s0 == 0 ? 1 : 0);

        double dS = 0.0;
        if (B0 != B1) dS += M.lbinom(_N - 1, B1 - 1) - M.lbinom(_N - 1, B0 - 1);
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_s1) - M.lgam1(wr_s0));

        if (!_independent) {
            long ko = _kout[v], ki = _directed ? _kin[v] : 0;
            long ep_r0 = _ep[r], ep_r1 = ep_r0 - ko;
            long ep_s0 = _ep[s], ep_s1 = ep_s0 + ko;
            long em_r0 = _directed ? _em[r] : 0, em_r1 = em_r0 - ki;
            long em_s0 = _directed ? _em[s] : 0, em_s1 = em_s0 + ki;
            if (_dc) {
                dS += (M.lgam1(ep_r1) - M.lgam1(ep_r0)) + (M.lgam1(ep_s1) - M.lgam1(ep_s0));
                if (_directed)
                    dS += (M.lgam1(em_r1) - M.lgam1(em_r0)) + (M.lgam1(em_s1) - M.lgam1(em_s0));
                const pp_hash_map<long, long>* hor = _dist ? &_h_out[r] : nullptr;
                const pp_hash_map<long, long>* hos = _dist ? &_h_out[s] : nullptr;
                const pp_hash_map<long, long>* hir = (_dist && _directed) ? &_h_in[r] : nullptr;
                const pp_hash_map<long, long>* his = (_dist && _directed) ? &_h_in[s] : nullptr;
                dS += deg_dl_delta_maps(hor, hir, wr_r0, wr_r1, ep_r0, ep_r1, em_r0, em_r1, ko, ki, -1);
                dS += deg_dl_delta_maps(hos, his, wr_s0, wr_s1, ep_s0, ep_s1, em_s0, em_s1, ko, ki, +1);
            } else {
                long dr0 = ep_r0 + em_r0, dr1 = ep_r1 + em_r1;
                long ds0 = ep_s0 + em_s0, ds1 = ep_s1 + em_s1;
                dS += (wr_r1 > 0 ? dr1 * M.logn(wr_r1) : 0.0) - dr0 * M.logn(wr_r0);
                dS += ds1 * M.logn(wr_s1) - (wr_s0 > 0 ? ds0 * M.logn(wr_s0) : 0.0);
            }
        } else {
            for (int l : _v_layers[v]) {
                long nr0 = _wr_l[l][r], nr1 = nr0 - 1;
                long ns0 = _wr_l[l][s], ns1 = ns0 + 1;
                long ko = _v_kout_l[v].get(l);
                long ki = _directed ? _v_kin_l[v].get(l) : 0;
                long epr0 = _ep_l[l][r], epr1 = epr0 - ko;
                long eps0 = _ep_l[l][s], eps1 = eps0 + ko;
                long emr0 = _directed ? _em_l[l][r] : 0, emr1 = emr0 - ki;
                long ems0 = _directed ? _em_l[l][s] : 0, ems1 = ems0 + ki;
                if (_dc) {
                    dS += (M.lgam1(epr1) - M.lgam1(epr0)) + (M.lgam1(eps1) - M.lgam1(eps0));
                    if (_directed)
                        dS += (M.lgam1(emr1) - M.lgam1(emr0)) + (M.lgam1(ems1) - M.lgam1(ems0));
                    const pp_hash_map<long, long>* hor = _dist ? &_h_out_l[l][r] : nullptr;
                    const pp_hash_map<long, long>* hos = _dist ? &_h_out_l[l][s] : nullptr;
                    const pp_hash_map<long, long>* hir = (_dist && _directed) ? &_h_in_l[l][r] : nullptr;
                    const pp_hash_map<long, long>* his = (_dist && _directed) ? &_h_in_l[l][s] : nullptr;
                    dS += deg_dl_delta_maps(hor, hir, nr0, nr1, epr0, epr1, emr0, emr1, ko, ki, -1);
                    dS += deg_dl_delta_maps(hos, his, ns0, ns1, eps0, eps1, ems0, ems1, ko, ki, +1);
                } else {
                    long dr0 = epr0 + emr0, dr1 = epr1 + emr1;
                    long ds0 = eps0 + ems0, ds1 = eps1 + ems1;
                    dS += (nr1 > 0 ? dr1 * M.logn(nr1) : 0.0) - (nr0 > 0 ? dr0 * M.logn(nr0) : 0.0);
                    dS += (ns1 > 0 ? ds1 * M.logn(ns1) : 0.0) - (ns0 > 0 ? ds0 * M.logn(ns0) : 0.0);
                }
            }
        }

        auto eval_layer_pp = [&](int l) {
            long el = _E_l[l];
            if (el == 0) return;
            long dm_r = _touched_l_flag[l] ? _dm_r_l[l] : 0;
            long dm_s = _touched_l_flag[l] ? _dm_s_l[l] : 0;
            long sl_l = _self_loops_l[v].get(l);
            long derr_rl = -(dm_r + sl_l);
            long derr_sl = +(dm_s + sl_l);
            long dein_l = dm_s - dm_r;
            long deout_l = -dein_l;
            long ein0 = _e_in_l[l], ein1 = ein0 + dein_l;
            long eout0 = _e_out_l[l], eout1 = eout0 + deout_l;
            if (B0 != B1) {
                dS += ((B1 > 1) ? M.logn(el + 1) : 0.0) -
                      ((B0 > 1) ? M.logn(el + 1) : 0.0);
            }
            if (_uniform) {
                dS -= (M.lgam1(ein1) - M.lgam1(ein0)) + (M.lgam1(eout1) - M.lgam1(eout0));
                if (_directed) {
                    dS += M.logn(B1) * ein1 - M.logn(B0) * ein0;
                    dS += (M.lbinom(B1, 2) + PP_LOG2) * eout1 -
                          (M.lbinom(B0, 2) + PP_LOG2) * eout0;
                } else {
                    dS += (M.logn(B1) - PP_LOG2) * ein1 - (M.logn(B0) - PP_LOG2) * ein0;
                    dS += M.lbinom(B1, 2) * eout1 - M.lbinom(B0, 2) * eout0;
                }
            } else {
                dS -= M.lgam1(eout1) - M.lgam1(eout0);
                if (_directed) {
                    dS += eout1 * (M.lbinom(B1, 2) + PP_LOG2) -
                          eout0 * (M.lbinom(B0, 2) + PP_LOG2);
                    dS -= (M.lgam1(_err_l[l][r] + derr_rl) - M.lgam1(_err_l[l][r])) +
                          (M.lgam1(_err_l[l][s] + derr_sl) - M.lgam1(_err_l[l][s]));
                } else {
                    dS += eout1 * M.lbinom(B1, 2) - eout0 * M.lbinom(B0, 2);
                    dS -= (derr_rl * PP_LOG2 + M.lgam1(_err_l[l][r] + derr_rl) - M.lgam1(_err_l[l][r])) +
                          (derr_sl * PP_LOG2 + M.lgam1(_err_l[l][s] + derr_sl) - M.lgam1(_err_l[l][s]));
                }
                dS += M.lbinom(B1 + ein1 - 1, ein1) - M.lbinom(B0 + ein0 - 1, ein0);
            }
        };
        if (B0 != B1) {
            for (int l : _active_layers) eval_layer_pp(l);
        } else {
            for (int l : _touched_layers) eval_layer_pp(l);
        }
        return dS;
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

std::vector<int> parse_pp_blocks(const Rcpp::IntegerVector& blocks, int N) {
    if (blocks.size() != N)
        Rcpp::stop("blocks must have length equal to number of vertices");
    pp_hash_map<int, int> remap;
    std::vector<int> out(N);
    int nxt = 0;
    for (int i = 0; i < N; ++i) {
        int val = blocks[i];
        if (val == NA_INTEGER || val < 1)
            Rcpp::stop("blocks must contain positive integer labels");
        auto it = remap.find(val);
        if (it == remap.end()) it = remap.emplace(val, nxt++).first;
        out[i] = it->second;
    }
    return out;
}

void pp_greedy_sweeps(PPState& st, PPRng& rng, size_t B_min, size_t B_max,
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
            pp_hash_set<int> cands;
            const auto& occ = st.occupied();
            if (occ.size() <= 16) {
                for (int g : occ) cands.insert(g);
            } else {
                for (int k = 0; k < 10; ++k) cands.insert(st.sample_group(v, 0.5, rng));
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
                if (dS < best_dS) {
                    best_dS = dS;
                    best_s = s;
                }
            }
            if (best_s != r) tot_dS += st.move(v, best_s);
        }
        if (std::abs(tot_dS) <= epsilon) break;
    }
}

void pp_mcmc_sweep_once(PPState& s, PPRng& rng, double beta, double c, double d,
                        bool multiflip, size_t B_min, size_t B_max,
                        double& total_dS, size_t& nattempts, size_t& nmoves) {
    if (s.N() <= 1 || s.E() == 0) return;
    std::vector<int> vs(s.N());
    std::iota(vs.begin(), vs.end(), 0);
    rng.shuffle(vs);

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
        double log_u = (u_rand > 0.0) ? std::log(u_rand) : -PP_INF;
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
                for (int gsw = 0; gsw < 8; ++gsw) {
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
                            if (gsw == 7) log_q += std::log(std::max(1e-300, p_oth));
                        } else {
                            if (gsw == 7) log_q += std::log(std::max(1e-300, 1.0 - p_oth));
                        }
                    }
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
                                    static_cast<double>(mem_r.size()) * PP_LOG2;
                    double log_a = -beta * dS + log_pr - log_pf;
                    accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
                }
                if (accept) { total_dS += dS; ++nmoves; }
                else { for (int v : mem_r) s.move(v, r); }
            }
        }
    }
}

}  // namespace

// [[Rcpp::export]]
double cpp_pp_entropy(SEXP pointer, Rcpp::IntegerVector blocks, bool uniform,
                      bool degree_corrected, bool distributed, double q_max,
                      Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                      Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                      Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                      Rcpp::Nullable<Rcpp::IntegerVector> layers = R_NilValue,
                      bool independent = false) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_pp_weights(weights);
    Rcpp::CharacterVector T = unwrap_pp_types(types);
    Rcpp::List P = unwrap_pp_params(params);
    Rcpp::IntegerVector L = unwrap_pp_layers(layers);
    auto covs = parse_covs(W, T, P);
    PPState st(graph, uniform, degree_corrected, distributed, W, covs, L, independent);
    st.set_partition(parse_pp_blocks(blocks, N));
    return st.entropy();
}

// [[Rcpp::export]]
Rcpp::List cpp_pp_stats(SEXP pointer, Rcpp::IntegerVector blocks, bool uniform,
                        bool degree_corrected, bool distributed, double q_max,
                        Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                        Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerVector> layers = R_NilValue,
                        bool independent = false) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_pp_weights(weights);
    Rcpp::CharacterVector T = unwrap_pp_types(types);
    Rcpp::List P = unwrap_pp_params(params);
    Rcpp::IntegerVector L = unwrap_pp_layers(layers);
    auto covs = parse_covs(W, T, P);
    PPState st(graph, uniform, degree_corrected, distributed, W, covs, L, independent);
    st.set_partition(parse_pp_blocks(blocks, N));
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("E_in") = static_cast<double>(st.E_in()),
        Rcpp::Named("E_out") = static_cast<double>(st.E_out()),
        Rcpp::Named("E_total") = static_cast<double>(st.E()),
        Rcpp::Named("B") = st.B(),
        Rcpp::Named("err") = st.external_err_r());
    if (L.size() > 0) {
        out["E_in_layer"] = st.e_in_l_r();
        out["E_out_layer"] = st.e_out_l_r();
        out["E_total_layer"] = st.e_tot_l_r();
        out["err_layer"] = st.external_err_l_r();
    }
    if (st.K() > 0) {
        out["weight_in_mean"] = st.w_in_mean_r();
        out["weight_out_mean"] = st.w_out_mean_r();
        if (L.size() > 0) {
            out["weight_in_mean_layer"] = st.w_in_mean_l_r();
            out["weight_out_mean_layer"] = st.w_out_mean_l_r();
        }
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_pp_check_moves(SEXP pointer, Rcpp::IntegerVector blocks,
                                       bool uniform, bool degree_corrected,
                                       bool distributed, int nmoves, double q_max,
                                       Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                                       Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                                       Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                                       Rcpp::Nullable<Rcpp::IntegerVector> layers = R_NilValue,
                                       bool independent = false) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_pp_weights(weights);
    Rcpp::CharacterVector T = unwrap_pp_types(types);
    Rcpp::List P = unwrap_pp_params(params);
    Rcpp::IntegerVector L = unwrap_pp_layers(layers);
    auto covs = parse_covs(W, T, P);
    PPState st(graph, uniform, degree_corrected, distributed, W, covs, L, independent);
    st.set_partition(parse_pp_blocks(blocks, N));
    PPRng rng;
    Rcpp::NumericVector err(nmoves);
    if (N <= 1) return err;
    double S = st.entropy();
    for (int i = 0; i < nmoves; ++i) {
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
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_pp_minimize(SEXP pointer, bool uniform, bool degree_corrected,
                           bool distributed, int B_min, double B_max,
                           bool refine, double epsilon, double q_max,
                           Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                           Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                           Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerVector> layers = R_NilValue,
                           bool independent = false) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_pp_weights(weights);
    Rcpp::CharacterVector T = unwrap_pp_types(types);
    Rcpp::List P = unwrap_pp_params(params);
    Rcpp::IntegerVector L = unwrap_pp_layers(layers);
    auto covs = parse_covs(W, T, P);
    PPState st(graph, uniform, degree_corrected, distributed, W, covs, L, independent);
    PPRng rng;
    if (N <= 1) {
        std::vector<int> b(N, 0);
        st.set_partition(b);
        return Rcpp::List::create(
            Rcpp::Named("blocks") = st.external_blocks_r(),
            Rcpp::Named("entropy") = st.entropy(),
            Rcpp::Named("B") = st.B());
    }
    size_t bmin = std::max<size_t>(1, std::min<size_t>(static_cast<size_t>(B_min), static_cast<size_t>(N)));
    size_t bmax = std::isinf(B_max) ? static_cast<size_t>(N)
                                    : std::max(bmin, std::min<size_t>(static_cast<size_t>(B_max), static_cast<size_t>(N)));

    std::vector<int> best_b(N, 0);
    double best_S = PP_INF;

    // Multi-start warm thermal mixing (beta = 1 -> 2 -> 4) + greedy quench at beta = Inf
    std::vector<size_t> init_Bs = {bmin};
    for (size_t cand_b : {size_t(4), size_t(12)}) {
        size_t clamped = std::min<size_t>(bmax, std::max<size_t>(bmin, std::min<size_t>(N, cand_b)));
        if (std::find(init_Bs.begin(), init_Bs.end(), clamped) == init_Bs.end())
            init_Bs.push_back(clamped);
    }

    for (size_t b0_k : init_Bs) {
        std::vector<int> b_init(N);
        for (int v = 0; v < N; ++v) b_init[v] = static_cast<int>(v % b0_k);
        st.set_partition(b_init);
        double tot_dS = 0.0;
        size_t na = 0, nm = 0;

        int warm_sw = refine ? 40 : 15;
        for (int sw = 0; sw < warm_sw; ++sw) {
            double b_cur = (sw < warm_sw / 2) ? 1.0 : ((sw < 3 * warm_sw / 4) ? 2.0 : 4.0);
            pp_mcmc_sweep_once(st, rng, b_cur, 0.5, 0.05, true, bmin, bmax, tot_dS, na, nm);
        }

        double cur_S = st.entropy();
        int no_imp = 0;
        int max_sw = refine ? 60 : 20;
        for (int sw = 0; sw < max_sw; ++sw) {
            pp_mcmc_sweep_once(st, rng, PP_INF, 0.5, 0.05, true, bmin, bmax, tot_dS, na, nm);
            double new_S = st.entropy();
            if (new_S < cur_S - epsilon) {
                cur_S = new_S;
                no_imp = 0;
            } else if (++no_imp >= 6) {
                break;
            }
        }
        if (refine) {
            pp_greedy_sweeps(st, rng, bmin, bmax, 10, epsilon);
            cur_S = st.entropy();
        }
        if (cur_S < best_S) {
            best_S = cur_S;
            best_b = st.blocks();
        }
    }
    st.set_partition(best_b);
    best_S = st.entropy();

    return Rcpp::List::create(
        Rcpp::Named("blocks") = st.external_blocks_r(),
        Rcpp::Named("entropy") = best_S,
        Rcpp::Named("B") = st.B());
}

// [[Rcpp::export]]
Rcpp::List cpp_pp_mcmc(SEXP pointer, Rcpp::IntegerVector blocks, bool uniform,
                       bool degree_corrected, bool distributed,
                       Rcpp::List mcmc_opts, Rcpp::List eq_opts, double q_max,
                       Rcpp::Nullable<Rcpp::NumericMatrix> weights = R_NilValue,
                       Rcpp::Nullable<Rcpp::CharacterVector> types = R_NilValue,
                       Rcpp::Nullable<Rcpp::List> params = R_NilValue,
                       Rcpp::Nullable<Rcpp::IntegerVector> layers = R_NilValue,
                       bool independent = false) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    Rcpp::NumericMatrix W = unwrap_pp_weights(weights);
    Rcpp::CharacterVector T = unwrap_pp_types(types);
    Rcpp::List P = unwrap_pp_params(params);
    Rcpp::IntegerVector L = unwrap_pp_layers(layers);
    auto covs = parse_covs(W, T, P);
    PPState st(graph, uniform, degree_corrected, distributed, W, covs, L, independent);
    st.set_partition(parse_pp_blocks(blocks, N));
    PPRng rng;

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
            pp_mcmc_sweep_once(st, rng, beta, c, d, multiflip, B_min, B_max,
                               total_dS, nattempts, nmoves);
            if ((it + 1) % 25 == 0) st.refresh_cov_globals();
            ++steps_done;
            Rcpp::checkUserInterrupt();
        }
        return Rcpp::List::create(
            Rcpp::Named("blocks") = st.external_blocks_r(),
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
                pp_mcmc_sweep_once(st, rng, cur_beta, c, d, multiflip,
                                   B_min, B_max, total_dS, nattempts, nmoves);
            }
            double S_now = st.entropy();
            ++steps_done;
            if (keep_hist) hist_S.push_back(S_now);
            if (collect && (step % thin == 0)) {
                samples.push_back(st.external_blocks_r());
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
        Rcpp::Named("blocks") = st.external_blocks_r(),
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
