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

namespace {

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

struct PPLogQ {
    size_t max_n = 0;
    std::vector<double> table;
    mutable pp_hash_map<uint64_t, double> memo;

    void init(size_t n) {
        if (n <= max_n) return;
        n = std::min<size_t>(n, 2500);
        if (n <= max_n) return;
        size_t stride = n + 1;
        std::vector<double> nt(stride * stride, 0.0);
        for (size_t i = 1; i <= n; ++i) {
            nt[i * stride + 1] = 0.0;
            for (size_t j = 2; j <= i; ++j) {
                double a = nt[i * stride + (j - 1)];
                double b = nt[(i - j) * stride + std::min(i - j, j)];
                double mx = std::max(a, b);
                nt[i * stride + j] = mx + std::log1p(std::exp(-std::abs(a - b)));
            }
        }
        table.swap(nt);
        max_n = n;
    }

    double operator()(long n, long k) const {
        if (n <= 0 || k <= 0) return 0.0;
        if (k > n) k = n;
        if (static_cast<size_t>(n) <= max_n)
            return table[static_cast<size_t>(n) * (max_n + 1) + static_cast<size_t>(k)];
        uint64_t key = (static_cast<uint64_t>(n) << 32) | static_cast<uint64_t>(k);
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;
        double dn = static_cast<double>(n);
        double c = M_PI / std::sqrt(6.0);
        double res = 2.0 * c * std::sqrt(dn) - std::log(4.0 * std::sqrt(3.0) * dn);
        memo.emplace(key, res);
        return res;
    }
};

PPLogQ& pplog_q() {
    static PPLogQ q;
    return q;
}

size_t pp_q_cache_max(double q_max) {
    if (!std::isfinite(q_max) || q_max <= 0.0) return 2000;
    return std::min<size_t>(static_cast<size_t>(q_max), 2000);
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
    PPState(const Graph& graph, bool uniform, bool dc, bool dist)
        : _N(static_cast<int>(boost::num_vertices(graph.adjacency))),
          _directed(graph.directed), _uniform(uniform), _dc(dc), _dist(dist),
          _cap(std::max(2, _N)),
          _adj(_N), _kout(_N, 0), _kin(graph.directed ? _N : 0, 0),
          _self_loops(_N, 0),
          _b(_N, 0), _wr(_cap, 0), _ep(_cap, 0), _em(graph.directed ? _cap : 0, 0),
          _err(_cap, 0),
          _occ_pos(_cap, -1), _empty_pos(_cap, -1),
          _members(_cap), _pos(_N, 0),
          _h_out(dc && dist ? _cap : 0),
          _h_in(dc && dist && graph.directed ? _cap : 0) {
        _empty.reserve(_cap);
        for (int r = _cap - 1; r >= 0; --r) push_empty(r);

        _E = 0;
        pp_hash_map<uint64_t, long> emult;
        for (auto [it, end] = boost::edges(graph.adjacency); it != end; ++it) {
            int u = static_cast<int>(boost::source(*it, graph.adjacency));
            int v = static_cast<int>(boost::target(*it, graph.adjacency));
            ++_E;
            emult[pp_canon(u, v, _directed)]++;
            if (u == v) {
                _self_loops[u]++;
                if (_directed) {
                    _kout[u]++;
                    _kin[u]++;
                } else {
                    _kout[u] += 2;
                }
            } else {
                _adj[u].push_back(v);
                _adj[v].push_back(u);
                if (_directed) {
                    _kout[u]++;
                    _kin[v]++;
                } else {
                    _kout[u]++;
                    _kout[v]++;
                }
            }
        }

        ppmath().ensure(static_cast<size_t>(_N + 2 * _E + 64));
        if (_dc && _dist) pplog_q().init(std::min<size_t>(static_cast<size_t>(2 * _E + 10), 2000));

        const PPMath& M = ppmath();
        _base_const = 0.0;
        for (const auto& [pk, mult] : emult) {
            int u = static_cast<int>(pk >> 32);
            int v = static_cast<int>(pk & 0xffffffffULL);
            _base_const += M.lgam1(mult);
            if (!_directed && u == v)
                _base_const += static_cast<double>(mult) * PP_LOG2;
        }
        if (_dc) {
            for (int v = 0; v < _N; ++v) {
                _base_const -= M.lgam1(_kout[v]);
                if (_directed) _base_const -= M.lgam1(_kin[v]);
            }
        }
    }

    int N() const { return _N; }
    long E() const { return _E; }
    int B() const { return static_cast<int>(_occupied.size()); }
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
            _members[r].clear();
            _occ_pos[r] = -1;
            if (_dc && _dist) {
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

        for (int v = 0; v < _N; ++v) {
            int r = b[v];
            _b[v] = r;
            add_member(v, r);
            _ep[r] += _kout[v];
            if (_directed) _em[r] += _kin[v];
            if (_self_loops[v] > 0) {
                _err[r] += _self_loops[v];
                _e_in += _self_loops[v];
            }
            if (_dc && _dist) {
                _h_out[r][_kout[v]]++;
                if (_directed) _h_in[r][_kin[v]]++;
            }
        }
        // Count non-self edges (each appears twice in _adj)
        long in_twice = 0, out_twice = 0;
        for (int v = 0; v < _N; ++v) {
            int rv = _b[v];
            for (int u : _adj[v]) {
                if (_b[u] == rv) {
                    if (v < u) {
                        _err[rv]++;
                        in_twice += 2;
                    }
                } else {
                    out_twice++;
                }
            }
        }
        _e_in += in_twice / 2;
        _e_out = out_twice / 2;
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

        // Edges DL
        if (Bn > 1 && _E > 0) S += M.logn(_E + 1);

        // Vertex/group degree + degree DL
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

        // Assortative edge likelihood
        if (_uniform) {
            S -= M.lgam1(_e_in) + M.lgam1(_e_out);
            if (_directed) {
                S += M.logn(Bn) * static_cast<double>(_e_in);
                S += (M.lbinom(Bn, 2) + PP_LOG2) * static_cast<double>(_e_out);
            } else {
                S += (M.logn(Bn) - PP_LOG2) * static_cast<double>(_e_in);
                S += M.lbinom(Bn, 2) * static_cast<double>(_e_out);
            }
        } else {
            S -= M.lgam1(_e_out);
            if (_directed)
                S += static_cast<double>(_e_out) * (M.lbinom(Bn, 2) + PP_LOG2);
            else
                S += static_cast<double>(_e_out) * M.lbinom(Bn, 2);
            for (int r : _occupied) {
                if (_directed) S -= M.lgam1(_err[r]);
                else S -= static_cast<double>(_err[r]) * PP_LOG2 + M.lgam1(_err[r]);
            }
            S += M.lbinom(Bn + _e_in - 1, _e_in);
        }
        return S;
    }

    // Exact O(deg(v)) zero-allocation virtual_move(v, s)
    double virtual_move(int v, int s) const {
        int r = _b[v];
        if (r == s) return 0.0;
        long m_r = 0, m_s = 0;
        for (int u : _adj[v]) {
            int bu = _b[u];
            if (bu == r) ++m_r;
            else if (bu == s) ++m_s;
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
        return delta_entropy(r, s, 1, _kout[v], _directed ? _kin[v] : 0,
                             derr_r, derr_s, de_in, de_out, true, v);
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
        _ep[r] -= _kout[v];
        _ep[s] += _kout[v];
        if (_directed) {
            _em[r] -= _kin[v];
            _em[s] += _kin[v];
        }
        _err[r] += _last_derr_r;
        _err[s] += _last_derr_s;
        _e_in += _last_de_in;
        _e_out += _last_de_out;
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
        return _b[adj[rng.index(adj.size())]];
    }

    double move_prob(int v, int s, double c) const {
        if (_occupied.empty()) return 0.0;
        const auto& adj = _adj[v];
        double B_act = static_cast<double>(_occupied.size());
        if (adj.empty() || std::isinf(c)) return 1.0 / B_act;
        long cnt = 0;
        for (int u : adj) if (_b[u] == s) ++cnt;
        double deg = static_cast<double>(adj.size());
        return (static_cast<double>(cnt) + c) / (deg + c * B_act);
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
    int _N;
    bool _directed, _uniform, _dc, _dist;
    int _cap;
    long _E = 0, _e_in = 0, _e_out = 0;
    double _base_const = 0.0;

    std::vector<std::vector<int>> _adj;
    std::vector<long> _kout, _kin, _self_loops;
    std::vector<int> _b, _wr;
    std::vector<long> _ep, _em, _err;
    std::vector<int> _occ_pos, _occupied, _empty, _empty_pos;
    std::vector<std::vector<int>> _members;
    std::vector<int> _pos;
    std::vector<pp_hash_map<long, long>> _h_out, _h_in;

    mutable long _last_derr_r = 0, _last_derr_s = 0, _last_de_in = 0, _last_de_out = 0;

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

    double deg_dl_delta_group(int r, long n0, long n1, long ep0, long ep1,
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

    double delta_entropy(int r, int s, long n_mv, long ko, long ki,
                         long derr_r, long derr_s, long de_in, long de_out,
                         bool single_node, int v) const {
        const PPMath& M = ppmath();
        long wr_r0 = _wr[r], wr_r1 = wr_r0 - n_mv;
        long wr_s0 = _wr[s], wr_s1 = wr_s0 + n_mv;
        long B0 = B();
        long B1 = B0 - (wr_r1 == 0 ? 1 : 0) + (wr_s0 == 0 ? 1 : 0);

        double dS = 0.0;
        if (B0 != B1) {
            dS += M.lbinom(_N - 1, B1 - 1) - M.lbinom(_N - 1, B0 - 1);
            if (_E > 0) {
                dS += ((B1 > 1) ? M.logn(_E + 1) : 0.0) -
                      ((B0 > 1) ? M.logn(_E + 1) : 0.0);
            }
        }
        dS -= (M.lgam1(wr_r1) - M.lgam1(wr_r0)) + (M.lgam1(wr_s1) - M.lgam1(wr_s0));

        long ep_r0 = _ep[r], ep_r1 = ep_r0 - ko;
        long ep_s0 = _ep[s], ep_s1 = ep_s0 + ko;
        long em_r0 = _directed ? _em[r] : 0, em_r1 = em_r0 - ki;
        long em_s0 = _directed ? _em[s] : 0, em_s1 = em_s0 + ki;

        if (_dc) {
            dS += (M.lgam1(ep_r1) - M.lgam1(ep_r0)) + (M.lgam1(ep_s1) - M.lgam1(ep_s0));
            if (_directed)
                dS += (M.lgam1(em_r1) - M.lgam1(em_r0)) + (M.lgam1(em_s1) - M.lgam1(em_s0));
            if (single_node) {
                dS += deg_dl_delta_group(r, wr_r0, wr_r1, ep_r0, ep_r1, em_r0, em_r1, ko, ki, -1);
                dS += deg_dl_delta_group(s, wr_s0, wr_s1, ep_s0, ep_s1, em_s0, em_s1, ko, ki, +1);
            }
            (void)v;
        } else {
            long dr0 = ep_r0 + em_r0, dr1 = ep_r1 + em_r1;
            long ds0 = ep_s0 + em_s0, ds1 = ep_s1 + em_s1;
            dS += (wr_r1 > 0 ? dr1 * M.logn(wr_r1) : 0.0) - dr0 * M.logn(wr_r0);
            dS += ds1 * M.logn(wr_s1) - (wr_s0 > 0 ? ds0 * M.logn(wr_s0) : 0.0);
        }

        long ein0 = _e_in, ein1 = ein0 + de_in;
        long eout0 = _e_out, eout1 = eout0 + de_out;

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
                dS -= (M.lgam1(_err[r] + derr_r) - M.lgam1(_err[r])) +
                      (M.lgam1(_err[s] + derr_s) - M.lgam1(_err[s]));
            } else {
                dS += eout1 * M.lbinom(B1, 2) - eout0 * M.lbinom(B0, 2);
                dS -= (derr_r * PP_LOG2 + M.lgam1(_err[r] + derr_r) - M.lgam1(_err[r])) +
                      (derr_s * PP_LOG2 + M.lgam1(_err[s] + derr_s) - M.lgam1(_err[s]));
            }
            dS += M.lbinom(B1 + ein1 - 1, ein1) - M.lbinom(B0 + ein0 - 1, ein0);
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
                      bool degree_corrected, bool distributed, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    PPState st(graph, uniform, degree_corrected, distributed);
    st.set_partition(parse_pp_blocks(blocks, N));
    return st.entropy();
}

// [[Rcpp::export]]
Rcpp::List cpp_pp_stats(SEXP pointer, Rcpp::IntegerVector blocks, bool uniform,
                        bool degree_corrected, bool distributed, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    PPState st(graph, uniform, degree_corrected, distributed);
    st.set_partition(parse_pp_blocks(blocks, N));
    return Rcpp::List::create(
        Rcpp::Named("E_in") = static_cast<double>(st.E_in()),
        Rcpp::Named("E_out") = static_cast<double>(st.E_out()),
        Rcpp::Named("E_total") = static_cast<double>(st.E()),
        Rcpp::Named("B") = st.B(),
        Rcpp::Named("err") = st.external_err_r());
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_pp_check_moves(SEXP pointer, Rcpp::IntegerVector blocks,
                                       bool uniform, bool degree_corrected,
                                       bool distributed, int nmoves, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    PPState st(graph, uniform, degree_corrected, distributed);
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
                           bool refine, double epsilon, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    PPState st(graph, uniform, degree_corrected, distributed);
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
                       Rcpp::List mcmc_opts, Rcpp::List eq_opts, double q_max) {
    Graph& graph = get_graph(pointer);
    int N = static_cast<int>(boost::num_vertices(graph.adjacency));
    if (degree_corrected && distributed) pplog_q().init(pp_q_cache_max(q_max));
    PPState st(graph, uniform, degree_corrected, distributed);
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
