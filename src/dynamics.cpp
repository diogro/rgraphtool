// rgraphtool: native C++ network reconstruction from dynamics and behavior
// (Peixoto 2019 PRL, Peixoto 2024, Peixoto 2025 PRX), plus forward simulators
// for Ising, epidemics (SI/SIS/SIR/SIRS), Gaussian, and Potts dynamics.
//
// Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
// Copyright (C) 2026 Diogo Melo (R port)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>
#include "dyn_sbm.h"
#include "graph_types.h"

namespace {

using namespace rgt_dyn;

constexpr double HALF_LOG_2PI = 0.9189385332046727417803297364056176398614;

enum DynModel {
    DYN_PSEUDO_ISING = 0,
    DYN_ISING_GLAUBER = 1,
    DYN_PSEUDO_CISING = 2,
    DYN_CISING_GLAUBER = 3,
    DYN_EPIDEMICS = 4,
    DYN_PSEUDO_NORMAL = 5,
    DYN_NORMAL_GLAUBER = 6,
    DYN_LINEAR_NORMAL = 7,
    DYN_PSEUDO_POTTS = 8,
    DYN_POTTS_GLAUBER = 9
};

inline bool is_time_series_model(int m) {
    return m == DYN_ISING_GLAUBER || m == DYN_CISING_GLAUBER ||
           m == DYN_EPIDEMICS || m == DYN_NORMAL_GLAUBER ||
           m == DYN_LINEAR_NORMAL || m == DYN_POTTS_GLAUBER;
}

inline double log_2cosh(double x) {
    double ax = std::abs(x);
    return ax + std::log1p(std::exp(-2.0 * ax));
}

inline double log_1p2cosh(double x) {
    double ax = std::abs(x);
    return ax + std::log1p(std::exp(-ax) + std::exp(-2.0 * ax));
}

inline double log_2sinh_over_x(double x) {
    double ax = std::abs(x);
    if (ax < 1e-5) {
        double x2 = ax * ax;
        return LOG2 + x2 / 6.0 - (x2 * x2) / 180.0;
    }
    return ax + std::log(-std::expm1(-2.0 * ax)) - std::log(ax);
}

inline double log_sinh_half(double z) {
    // log(sinh(z / 2)) for z > 0
    double hz = 0.5 * z;
    if (hz < 1e-5) return std::log(hz) + (hz * hz) / 6.0;
    if (hz > 20.0) return hz - LOG2;
    return std::log(std::sinh(hz));
}

inline double log_1mexp_neg_half(double z) {
    // log(1 - exp(-z / 2)) for z > 0
    double hz = 0.5 * z;
    if (hz < 1e-5) return std::log(hz) - 0.5 * hz + (hz * hz) / 24.0;
    return std::log(-std::expm1(-hz));
}

// Discretized Laplace NLL for a single value v with rate lam > 0 and bin width delta > 0
inline double laplace_bin_nll(double v, double lam, double delta) {
    if (lam <= 0.0) return 0.0;
    if (v == 0.0) return -log_1mexp_neg_half(lam * delta);
    return lam * std::abs(v) - log_sinh_half(lam * delta);
}

// Quantized MDL distribution prior matching Peixoto (2025 PRX) & graph-tool 3.9
// (including the 128-nat floating-point range encoding constant when E > 0)
double quantized_dist_nll(const std::vector<double>& vals, double delta, double lam,
                          bool uniform, bool exclude_zero) {
    size_t E = vals.size();
    if (E == 0) return 0.0;
    gt_hash_map<int64_t, int64_t> cnt;
    int64_t k_min = std::numeric_limits<int64_t>::max();
    int64_t k_max = std::numeric_limits<int64_t>::min();
    for (double v : vals) {
        int64_t q = static_cast<int64_t>(std::llround(v / delta));
        cnt[q]++;
        if (q < k_min) k_min = q;
        if (q > k_max) k_max = q;
    }
    int64_t K = static_cast<int64_t>(cnt.size());
    int64_t Ei = static_cast<int64_t>(E);

    double S = 0.0;
    if (uniform) {
        S += static_cast<double>(Ei) * std::log(static_cast<double>(K));
    } else {
        S += std::log(static_cast<double>(Ei)) +
             lbinom_d(static_cast<double>(Ei - 1), static_cast<double>(K - 1)) +
             lgam(static_cast<double>(Ei + 1));
        for (const auto& [q, m] : cnt) {
            S -= lgam(static_cast<double>(m + 1));
        }
    }

    int64_t span = k_max - k_min;
    bool has_zero_between = exclude_zero && (k_min < 0 && k_max > 0);
    if (K > 2) {
        int64_t avail = span - 1 - (has_zero_between ? 1 : 0);
        if (avail < K - 2) avail = K - 2;
        S += lbinom_d(static_cast<double>(avail), static_cast<double>(K - 2));
    }

    double z1 = static_cast<double>(k_min) * delta;
    double zK = static_cast<double>(k_max) * delta;
    if (lam > 0.0) {
        double z_shift = exclude_zero ? delta : 0.0;
        S += 128.0 + lam * (std::abs(z1) + std::abs(zK) - z_shift) -
             2.0 * log_sinh_half(lam * delta);
        if (K >= 2) S -= LOG2;
    }
    return S;
}

struct PriorOpts {
    bool sbm = true;
    bool latent_edges = true;
    bool density = false;
    double aE = 1.0;
    bool xdist = false;
    bool xdist_uniform = false;
    double xl1 = 1.0;
    double xdelta = 1e-8;
    std::vector<bool> tdist;
    std::vector<bool> tdist_uniform;
    std::vector<double> tl1;
    std::vector<double> tdelta;
    bool normal_prior = false;
    double normal_mu = 0.0;
    double normal_sigma = 1.0;
};

class DynamicsEngine {
public:
    DynamicsEngine(int model, int N, bool directed, bool self_loops,
                   int q_t, bool has_zero, bool exposed, bool positive, double pslack,
                   double xmin, double xmax,
                   const std::vector<double>& tmin, const std::vector<double>& tmax,
                   bool dc, bool dist, size_t q_max, int sbm_levels)
        : _model(model), _N(N), _directed(directed), _self_loops(self_loops),
          _q_t(q_t), _has_zero(has_zero), _exposed(exposed),
          _positive(positive), _pslack(pslack),
          _xmin(xmin), _xmax(xmax), _tmin(tmin), _tmax(tmax),
          _prior(N, directed, dc, dist, q_max, sbm_levels),
          _adj_in(N), _adj_out(directed ? N : 0),
          _abs_deg(N, 0.0),
          _theta(N, std::vector<double>(q_t, 0.0)),
          _node_nll(N, 0.0) {
        if (directed) _M_dyads = self_loops ? int64_t(N) * N : int64_t(N) * (N - 1);
        else _M_dyads = self_loops ? int64_t(N) * (N + 1) / 2 : int64_t(N) * (N - 1) / 2;
    }

    int N() const { return _N; }
    int model() const { return _model; }
    int q_t() const { return _q_t; }
    bool directed() const { return _directed; }
    bool self_loops() const { return _self_loops; }
    int64_t M_dyads() const { return _M_dyads; }
    DynBlockPrior& prior() { return _prior; }
    const DynBlockPrior& prior() const { return _prior; }
    const gt_hash_map<uint64_t, double>& weights() const { return _x; }
    const std::vector<std::vector<double>>& theta() const { return _theta; }

    uint64_t key(int u, int v) const {
        if (!_directed && u > v) std::swap(u, v);
        return pair_key(u, v);
    }

    bool valid_dyad(int u, int v) const {
        if (u < 0 || u >= _N || v < 0 || v >= _N) return false;
        if (!_self_loops && u == v) return false;
        return true;
    }

    double get_x(int u, int v) const {
        auto it = _x.find(key(u, v));
        return it == _x.end() ? 0.0 : it->second;
    }

    bool has_edge(int u, int v) const {
        return _x.find(key(u, v)) != _x.end();
    }

    // Load data from list of N x M_c matrices (supports single matrix or multiple cascades)
    void set_data(const Rcpp::List& s_list, const Rcpp::NumericMatrix& f_mat) {
        _f_q = f_mat.nrow();
        _f.assign(_f_q * _f_q, 0.0);
        for (int r = 0; r < _f_q; ++r)
            for (int c = 0; c < _f_q; ++c)
                _f[r * _f_q + c] = f_mat(r, c);

        bool ts = is_time_series_model(_model);
        _M = 0;
        for (R_xlen_t c = 0; c < s_list.size(); ++c) {
            Rcpp::NumericMatrix sm = s_list[c];
            if (sm.nrow() != _N) Rcpp::stop("Data matrix nrow must equal N.");
            int mc = sm.ncol();
            if (ts) _M += std::max(0, mc - 1);
            else _M += mc;
        }
        _s_prev.assign(static_cast<size_t>(_N) * _M, 0.0);
        _s_next.assign(static_cast<size_t>(_N) * _M, 0.0);
        if (_model == DYN_PSEUDO_POTTS || _model == DYN_POTTS_GLAUBER) {
            _si_prev.assign(static_cast<size_t>(_N) * _M, 0);
            _si_next.assign(static_cast<size_t>(_N) * _M, 0);
        }

        int pos = 0;
        for (R_xlen_t c = 0; c < s_list.size(); ++c) {
            Rcpp::NumericMatrix sm = s_list[c];
            int mc = sm.ncol();
            if (ts) {
                for (int t = 0; t + 1 < mc; ++t) {
                    for (int i = 0; i < _N; ++i) {
                        double sp = sm(i, t), sn = sm(i, t + 1);
                        _s_prev[static_cast<size_t>(i) * _M + pos] = sp;
                        _s_next[static_cast<size_t>(i) * _M + pos] = sn;
                        if (!_si_prev.empty()) {
                            _si_prev[static_cast<size_t>(i) * _M + pos] = static_cast<int>(sp);
                            _si_next[static_cast<size_t>(i) * _M + pos] = static_cast<int>(sn);
                        }
                    }
                    ++pos;
                }
            } else {
                for (int t = 0; t < mc; ++t) {
                    for (int i = 0; i < _N; ++i) {
                        double sv = sm(i, t);
                        _s_prev[static_cast<size_t>(i) * _M + pos] = sv;
                        _s_next[static_cast<size_t>(i) * _M + pos] = sv;
                        if (!_si_prev.empty()) {
                            _si_prev[static_cast<size_t>(i) * _M + pos] = static_cast<int>(sv);
                            _si_next[static_cast<size_t>(i) * _M + pos] = static_cast<int>(sv);
                        }
                    }
                    ++pos;
                }
            }
        }
    }

    void init_state(const std::vector<std::tuple<int, int, double>>& init_edges,
                    const std::vector<std::vector<double>>& theta_init,
                    const std::vector<std::vector<int>>& bs) {
        _x.clear();
        for (int i = 0; i < _N; ++i) {
            _adj_in[i].clear();
            if (_directed) _adj_out[i].clear();
            _abs_deg[i] = 0.0;
        }
        std::vector<std::tuple<int, int, long>> sbm_edges;
        sbm_edges.reserve(init_edges.size());
        for (const auto& [u0, v0, w] : init_edges) {
            if (!valid_dyad(u0, v0)) continue;
            int u = u0, v = v0;
            if (!_directed && u > v) std::swap(u, v);
            uint64_t k = pair_key(u, v);
            if (_x.find(k) != _x.end()) continue;
            _x[k] = w;
            sbm_edges.emplace_back(u, v, 1L);
            if (_directed) {
                _adj_in[v][u] = w;
                _adj_out[u][v] = w;
                _abs_deg[v] += std::abs(w);
            } else {
                _adj_in[u][v] = w;
                _abs_deg[u] += std::abs(w);
                if (u != v) {
                    _adj_in[v][u] = w;
                    _abs_deg[v] += std::abs(w);
                }
            }
        }
        _theta = theta_init;
        _prior.init(sbm_edges, bs);
        rebuild_active();
        recompute_all_fields_and_nll();
    }

    void rebuild_active() {
        _active.clear();
        _active_pos.clear();
        for (const auto& [k, w] : _x) {
            _active_pos[k] = _active.size();
            _active.push_back(k);
        }
    }

    // Recompute cached fields and per-vertex NLL from scratch
    void recompute_all_fields_and_nll() {
        if (_model == DYN_PSEUDO_POTTS || _model == DYN_POTTS_GLAUBER) {
            _h_potts.assign(static_cast<size_t>(_N) * _f_q * _M, 0.0);
        } else {
            _h.assign(static_cast<size_t>(_N) * _M, 0.0);
        }
        _total_nll = 0.0;
        for (int i = 0; i < _N; ++i) {
            recompute_vertex_field(i);
            _node_nll[i] = compute_vertex_nll(i);
            _total_nll += _node_nll[i];
        }
    }

    double data_nll() const { return _total_nll; }

    double geom_edge_prior(long E, bool density, double aE) const {
        if (density) {
            return static_cast<double>(E) * std::log(aE) - lgam(static_cast<double>(E) + 1.0) + aE;
        }
        if (_M_dyads <= 0) return 0.0;
        double Md = static_cast<double>(_M_dyads);
        return std::log(Md + 1.0) + static_cast<double>(E) * std::log((Md + 1.0) / Md);
    }

    double x_prior(const PriorOpts& opts) const {
        if (opts.xdist) {
            std::vector<double> vals;
            vals.reserve(_x.size());
            for (const auto& [k, w] : _x) vals.push_back(w);
            return quantized_dist_nll(vals, opts.xdelta, opts.xl1, opts.xdist_uniform, true);
        }
        double S = 0.0;
        if (opts.normal_prior) {
            double inv_s2 = 1.0 / (opts.normal_sigma * opts.normal_sigma);
            double c = std::log(opts.normal_sigma) + HALF_LOG_2PI;
            int64_t E = static_cast<int64_t>(_x.size());
            int64_t n_zero = _M_dyads - E;
            // Active edges
            for (const auto& [k, w] : _x) {
                double d = w - opts.normal_mu;
                S += 0.5 * d * d * inv_s2 + c;
            }
            if (n_zero > 0) {
                double d0 = -opts.normal_mu;
                S += static_cast<double>(n_zero) * (0.5 * d0 * d0 * inv_s2 + c);
            }
        }
        if (opts.xl1 > 0.0) {
            if (opts.sbm) {
                for (const auto& [k, w] : _x) {
                    S += opts.xl1 * (std::abs(w) - 0.5 * opts.xdelta) -
                         log_sinh_half(opts.xl1 * opts.xdelta);
                }
            } else {
                int64_t E = static_cast<int64_t>(_x.size());
                int64_t n_zero = _M_dyads - E;
                if (n_zero > 0) S += static_cast<double>(n_zero) * laplace_bin_nll(0.0, opts.xl1, opts.xdelta);
                for (const auto& [k, w] : _x) S += laplace_bin_nll(w, opts.xl1, opts.xdelta);
            }
        }
        return S;
    }

    double t_prior(const PriorOpts& opts) const {
        double S = 0.0;
        for (int l = 0; l < _q_t; ++l) {
            bool td = (l < static_cast<int>(opts.tdist.size())) ? opts.tdist[l] : false;
            bool tu = (l < static_cast<int>(opts.tdist_uniform.size())) ? opts.tdist_uniform[l] : false;
            double lam = (l < static_cast<int>(opts.tl1.size())) ? opts.tl1[l] : 0.0;
            double del = (l < static_cast<int>(opts.tdelta.size())) ? opts.tdelta[l] : 1e-8;
            if (td) {
                std::vector<double> vals(_N);
                for (int i = 0; i < _N; ++i) vals[i] = _theta[i][l];
                S += quantized_dist_nll(vals, del, lam, tu, false);
            } else if (lam > 0.0) {
                for (int i = 0; i < _N; ++i) S += laplace_bin_nll(_theta[i][l], lam, del);
            }
        }
        return S;
    }

    double total_entropy(const PriorOpts& opts) const {
        double S = geom_edge_prior(static_cast<long>(_x.size()), opts.density, opts.aE);
        if (opts.latent_edges) S += _total_nll;
        S += x_prior(opts);
        S += t_prior(opts);
        if (opts.sbm) S += _prior.entropy();
        return S;
    }

    // Exact incremental delta for changing edge (u, v) from current state to (present_new, w_new)
    double edge_move_dS(int u0, int v0, bool present_new, double w_new,
                        const PriorOpts& opts) {
        if (!valid_dyad(u0, v0)) return INF;
        if (present_new && (w_new < _xmin || w_new > _xmax)) return INF;
        int u = u0, v = v0;
        if (!_directed && u > v) std::swap(u, v);
        uint64_t k = pair_key(u, v);
        auto it = _x.find(k);
        bool present_old = (it != _x.end());
        double w_old = present_old ? it->second : 0.0;
        if (!present_new) w_new = 0.0;
        if (present_old == present_new && w_old == w_new) return 0.0;

        double dw = w_new - w_old;
        double dS_nll = 0.0;
        if (opts.latent_edges && dw != 0.0) {
            if (_directed) {
                // Edge u -> v only affects target v's field!
                double d_abs = std::abs(w_new) - std::abs(w_old);
                dS_nll += eval_vertex_nll_with_edge_delta(v, u, dw, d_abs) - _node_nll[v];
            } else {
                double d_abs = std::abs(w_new) - std::abs(w_old);
                dS_nll += eval_vertex_nll_with_edge_delta(u, v, dw, d_abs) - _node_nll[u];
                if (u != v) {
                    dS_nll += eval_vertex_nll_with_edge_delta(v, u, dw, d_abs) - _node_nll[v];
                }
            }
        }

        long E0 = static_cast<long>(_x.size());
        long dE = (present_new ? 1L : 0L) - (present_old ? 1L : 0L);
        double dS_geom = (dE != 0)
            ? (geom_edge_prior(E0 + dE, opts.density, opts.aE) - geom_edge_prior(E0, opts.density, opts.aE))
            : 0.0;

        double dS_sbm = (opts.sbm && dE != 0) ? _prior.edge_delta(u, v, dE, false) : 0.0;

        double dS_xp = 0.0;
        if (opts.xdist) {
            double old_xp = x_prior(opts);
            // Temporarily evaluate new x_prior
            std::vector<double> vals;
            vals.reserve(_x.size() + 1);
            for (const auto& [ek, ew] : _x) {
                if (ek == k) {
                    if (present_new) vals.push_back(w_new);
                } else {
                    vals.push_back(ew);
                }
            }
            if (!present_old && present_new) vals.push_back(w_new);
            dS_xp = quantized_dist_nll(vals, opts.xdelta, opts.xl1, opts.xdist_uniform, true) - old_xp;
        } else {
            if (opts.normal_prior) {
                double inv_s2 = 1.0 / (opts.normal_sigma * opts.normal_sigma);
                double d0 = (present_old ? w_old : 0.0) - opts.normal_mu;
                double d1 = (present_new ? w_new : 0.0) - opts.normal_mu;
                dS_xp += 0.5 * (d1 * d1 - d0 * d0) * inv_s2;
            }
            if (opts.xl1 > 0.0) {
                if (opts.sbm) {
                    double old_l = present_old
                        ? (opts.xl1 * (std::abs(w_old) - 0.5 * opts.xdelta) - log_sinh_half(opts.xl1 * opts.xdelta))
                        : 0.0;
                    double new_l = present_new
                        ? (opts.xl1 * (std::abs(w_new) - 0.5 * opts.xdelta) - log_sinh_half(opts.xl1 * opts.xdelta))
                        : 0.0;
                    dS_xp += (new_l - old_l);
                } else {
                    double old_l = laplace_bin_nll(present_old ? w_old : 0.0, opts.xl1, opts.xdelta);
                    double new_l = laplace_bin_nll(present_new ? w_new : 0.0, opts.xl1, opts.xdelta);
                    dS_xp += (new_l - old_l);
                }
            }
        }
        return dS_nll + dS_geom + dS_sbm + dS_xp;
    }

    void apply_edge_move(int u0, int v0, bool present_new, double w_new, const PriorOpts& opts) {
        int u = u0, v = v0;
        if (!_directed && u > v) std::swap(u, v);
        uint64_t k = pair_key(u, v);
        auto it = _x.find(k);
        bool present_old = (it != _x.end());
        double w_old = present_old ? it->second : 0.0;
        if (!present_new) w_new = 0.0;
        if (present_old == present_new && w_old == w_new) return;

        double dw = w_new - w_old;
        double d_abs = std::abs(w_new) - std::abs(w_old);
        long dE = (present_new ? 1L : 0L) - (present_old ? 1L : 0L);

        if (present_new) {
            _x[k] = w_new;
            if (!present_old) active_insert(k);
            if (_directed) {
                _adj_in[v][u] = w_new;
                _adj_out[u][v] = w_new;
                _abs_deg[v] += d_abs;
            } else {
                _adj_in[u][v] = w_new;
                _abs_deg[u] += d_abs;
                if (u != v) {
                    _adj_in[v][u] = w_new;
                    _abs_deg[v] += d_abs;
                }
            }
        } else {
            _x.erase(it);
            active_erase(k);
            if (_directed) {
                _adj_in[v].erase(u);
                _adj_out[u].erase(v);
                _abs_deg[v] += d_abs;
            } else {
                _adj_in[u].erase(v);
                _abs_deg[u] += d_abs;
                if (u != v) {
                    _adj_in[v].erase(u);
                    _abs_deg[v] += d_abs;
                }
            }
        }

        if (dw != 0.0) {
            if (_directed) {
                apply_vertex_field_delta(v, u, dw);
                double nll_new = compute_vertex_nll(v);
                _total_nll += (nll_new - _node_nll[v]);
                _node_nll[v] = nll_new;
            } else {
                apply_vertex_field_delta(u, v, dw);
                double nll_u = compute_vertex_nll(u);
                _total_nll += (nll_u - _node_nll[u]);
                _node_nll[u] = nll_u;
                if (u != v) {
                    apply_vertex_field_delta(v, u, dw);
                    double nll_v = compute_vertex_nll(v);
                    _total_nll += (nll_v - _node_nll[v]);
                    _node_nll[v] = nll_v;
                }
            }
        }

        if (dE != 0) _prior.edge_delta(u, v, dE, true);
    }

    // Exact incremental delta for updating node i's parameter theta[i][l] -> t_new
    double theta_move_dS(int i, int l, double t_new, const PriorOpts& opts) {
        if (i < 0 || i >= _N || l < 0 || l >= _q_t) return INF;
        if (t_new < _tmin[l] || t_new > _tmax[l]) return INF;
        double t_old = _theta[i][l];
        if (t_new == t_old) return 0.0;

        double dS_nll = 0.0;
        if (opts.latent_edges) {
            _theta[i][l] = t_new;
            double nll_new = compute_vertex_nll(i);
            _theta[i][l] = t_old;
            dS_nll = nll_new - _node_nll[i];
        }

        double dS_tp = 0.0;
        bool td = (l < static_cast<int>(opts.tdist.size())) ? opts.tdist[l] : false;
        bool tu = (l < static_cast<int>(opts.tdist_uniform.size())) ? opts.tdist_uniform[l] : false;
        double lam = (l < static_cast<int>(opts.tl1.size())) ? opts.tl1[l] : 0.0;
        double del = (l < static_cast<int>(opts.tdelta.size())) ? opts.tdelta[l] : 1e-8;
        if (td) {
            std::vector<double> v0(_N), v1(_N);
            for (int j = 0; j < _N; ++j) { v0[j] = _theta[j][l]; v1[j] = _theta[j][l]; }
            v1[i] = t_new;
            dS_tp = quantized_dist_nll(v1, del, lam, tu, false) -
                    quantized_dist_nll(v0, del, lam, tu, false);
        } else if (lam > 0.0) {
            dS_tp = laplace_bin_nll(t_new, lam, del) - laplace_bin_nll(t_old, lam, del);
        }
        return dS_nll + dS_tp;
    }

    void apply_theta_move(int i, int l, double t_new) {
        if (_theta[i][l] == t_new) return;
        _theta[i][l] = t_new;
        double nll_new = compute_vertex_nll(i);
        _total_nll += (nll_new - _node_nll[i]);
        _node_nll[i] = nll_new;
    }

    // 1D golden-section minimization on [lo, hi] for edge weight (u, v) or node theta(i, l)
    double optimize_edge_weight(int u, int v, const PriorOpts& opts, double& best_dS) {
        double w0 = get_x(u, v);
        // Bracket search around w0 (or around 0 if adding new edge)
        double lo = std::max(_xmin, w0 - 2.0);
        double hi = std::min(_xmax, w0 + 2.0);
        if (_model == DYN_EPIDEMICS) {
            lo = std::max(_xmin, -6.0);
            hi = std::min(_xmax, -1e-4);
        } else if (std::isinf(lo) || std::isinf(hi)) {
            lo = std::max(_xmin, -5.0);
            hi = std::min(_xmax, 5.0);
        }
        // Quick coarse scan over 9 points to bracket basin cleanly
        double best_w = (w0 != 0.0) ? w0 : 0.5 * (lo + hi);
        best_dS = edge_move_dS(u, v, true, best_w, opts);
        int n_grid = 9;
        for (int g = 0; g <= n_grid; ++g) {
            double wg = lo + (hi - lo) * static_cast<double>(g) / static_cast<double>(n_grid);
            if (wg == 0.0) continue;
            double d = edge_move_dS(u, v, true, wg, opts);
            if (d < best_dS) { best_dS = d; best_w = wg; }
        }
        double step = (hi - lo) / static_cast<double>(n_grid);
        double a = std::max(_xmin, best_w - step);
        double b = std::min(_xmax, best_w + step);
        constexpr double gr = 0.6180339887498948482;
        double c = b - gr * (b - a);
        double d = a + gr * (b - a);
        double fc = edge_move_dS(u, v, true, c, opts);
        double fd = edge_move_dS(u, v, true, d, opts);
        for (int it = 0; it < 18; ++it) {
            if (fc < fd) {
                b = d; d = c; fd = fc;
                c = b - gr * (b - a);
                fc = edge_move_dS(u, v, true, c, opts);
            } else {
                a = c; c = d; fc = fd;
                d = a + gr * (b - a);
                fd = edge_move_dS(u, v, true, d, opts);
            }
        }
        double mid = 0.5 * (a + b);
        double fmid = edge_move_dS(u, v, true, mid, opts);
        if (fmid < best_dS) { best_dS = fmid; best_w = mid; }
        return best_w;
    }

    double optimize_theta_val(int i, int l, const PriorOpts& opts, double& best_dS) {
        double t0 = _theta[i][l];
        double lo = std::max(_tmin[l], t0 - 2.0);
        double hi = std::min(_tmax[l], t0 + 2.0);
        if (_model == DYN_EPIDEMICS) {
            lo = std::max(_tmin[l], -12.0);
            hi = std::min(_tmax[l], -1e-4);
        }
        double best_t = t0;
        best_dS = 0.0;
        int n_grid = 8;
        for (int g = 0; g <= n_grid; ++g) {
            double tg = lo + (hi - lo) * static_cast<double>(g) / static_cast<double>(n_grid);
            double d = theta_move_dS(i, l, tg, opts);
            if (d < best_dS) { best_dS = d; best_t = tg; }
        }
        double step = (hi - lo) / static_cast<double>(n_grid);
        double a = std::max(_tmin[l], best_t - step);
        double b = std::min(_tmax[l], best_t + step);
        constexpr double gr = 0.6180339887498948482;
        double c = b - gr * (b - a);
        double d = a + gr * (b - a);
        double fc = theta_move_dS(i, l, c, opts);
        double fd = theta_move_dS(i, l, d, opts);
        for (int it = 0; it < 16; ++it) {
            if (fc < fd) {
                b = d; d = c; fd = fc;
                c = b - gr * (b - a);
                fc = theta_move_dS(i, l, c, opts);
            } else {
                a = c; c = d; fc = fd;
                d = a + gr * (b - a);
                fd = theta_move_dS(i, l, d, opts);
            }
        }
        double mid = 0.5 * (a + b);
        double fmid = theta_move_dS(i, l, mid, opts);
        if (fmid < best_dS) { best_dS = fmid; best_t = mid; }
        return best_t;
    }

    // Build top-k candidate dyads using fast empirical correlation screening
    std::vector<uint64_t> build_candidates(int k_per_node) const {
        std::vector<uint64_t> cands;
        if (_N <= 1) return cands;
        int64_t cap_all = _M_dyads;
        if (cap_all <= static_cast<int64_t>(_N) * std::max(k_per_node, 20) || _N <= 45) {
            // Enumerate all valid dyads directly for small/medium graphs
            cands.reserve(static_cast<size_t>(cap_all));
            for (int u = 0; u < _N; ++u) {
                int v_start = _directed ? 0 : u;
                for (int v = v_start; v < _N; ++v) {
                    if (!_self_loops && u == v) continue;
                    cands.push_back(key(u, v));
                }
            }
            return cands;
        }

        // For larger N: screen top k_per_node partners per vertex via normalized correlation
        std::vector<double> mean_p(_N, 0.0), sd_p(_N, 0.0);
        std::vector<double> mean_n(_N, 0.0), sd_n(_N, 0.0);
        for (int i = 0; i < _N; ++i) {
            double sp = 0.0, sn = 0.0;
            for (int m = 0; m < _M; ++m) {
                sp += _s_prev[static_cast<size_t>(i) * _M + m];
                sn += _s_next[static_cast<size_t>(i) * _M + m];
            }
            mean_p[i] = sp / std::max(1, _M);
            mean_n[i] = sn / std::max(1, _M);
            double vp = 0.0, vn = 0.0;
            for (int m = 0; m < _M; ++m) {
                double dp = _s_prev[static_cast<size_t>(i) * _M + m] - mean_p[i];
                double dn = _s_next[static_cast<size_t>(i) * _M + m] - mean_n[i];
                vp += dp * dp;
                vn += dn * dn;
            }
            sd_p[i] = std::sqrt(vp) + 1e-12;
            sd_n[i] = std::sqrt(vn) + 1e-12;
        }

        gt_hash_map<uint64_t, int> seen;
        std::vector<std::pair<double, int>> scores;
        scores.reserve(_N);
        int k_keep = std::min(_N - 1, std::max(4, k_per_node));
        for (int v = 0; v < _N; ++v) {
            scores.clear();
            for (int u = 0; u < _N; ++u) {
                if (!_self_loops && u == v) continue;
                double dot = 0.0;
                const double* pu = &_s_prev[static_cast<size_t>(u) * _M];
                const double* pv = &_s_next[static_cast<size_t>(v) * _M];
                double mu = mean_p[u], mv = mean_n[v];
                for (int m = 0; m < _M; ++m) dot += (pu[m] - mu) * (pv[m] - mv);
                double corr = std::abs(dot) / (sd_p[u] * sd_n[v]);
                scores.emplace_back(corr, u);
            }
            if (static_cast<int>(scores.size()) > k_keep) {
                std::nth_element(scores.begin(), scores.begin() + k_keep, scores.end(),
                                 [](const auto& a, const auto& b) { return a.first > b.first; });
                scores.resize(k_keep);
            }
            for (const auto& [sc, u] : scores) {
                uint64_t k = key(u, v);
                if (seen.emplace(k, 1).second) cands.push_back(k);
            }
        }
        return cands;
    }

    // Greedy MAP / coordinate-descent sweep (beta = Inf)
    void greedy_sweep(Rng& rng, const std::vector<uint64_t>& cands,
                      const PriorOpts& opts, bool update_sbm,
                      double& total_dS, size_t& nattempts, size_t& nmoves) {
        // 1. Update node parameters theta[i][l]
        std::vector<int> vord(_N);
        std::iota(vord.begin(), vord.end(), 0);
        rng.shuffle(vord);
        for (int i : vord) {
            for (int l = 0; l < _q_t; ++l) {
                ++nattempts;
                double best_dS = 0.0;
                double best_t = optimize_theta_val(i, l, opts, best_dS);
                if (best_dS < -1e-9) {
                    apply_theta_move(i, l, best_t);
                    total_dS += best_dS;
                    ++nmoves;
                }
            }
        }

        // 2. Sweep candidate dyads (add / refine / delete)
        std::vector<uint64_t> ord = cands;
        // Also ensure all currently active edges are visited
        gt_hash_map<uint64_t, int> in_ord;
        for (uint64_t k : ord) in_ord[k] = 1;
        for (uint64_t k : _active) {
            if (in_ord.find(k) == in_ord.end()) ord.push_back(k);
        }
        rng.shuffle(ord);

        for (uint64_t k : ord) {
            auto [u, v] = unpack_key(k);
            ++nattempts;
            bool cur_pres = has_edge(u, v);
            double dS_del = cur_pres ? edge_move_dS(u, v, false, 0.0, opts) : 0.0;
            double dS_opt = INF;
            double w_opt = optimize_edge_weight(u, v, opts, dS_opt);
            if (cur_pres) {
                if (dS_del < dS_opt && dS_del < -1e-9) {
                    apply_edge_move(u, v, false, 0.0, opts);
                    total_dS += dS_del;
                    ++nmoves;
                } else if (dS_opt < -1e-9) {
                    apply_edge_move(u, v, true, w_opt, opts);
                    total_dS += dS_opt;
                    ++nmoves;
                }
            } else {
                if (dS_opt < -1e-9) {
                    apply_edge_move(u, v, true, w_opt, opts);
                    total_dS += dS_opt;
                    ++nmoves;
                }
            }
        }

        // 3. SBM partition sweep (if sbm = TRUE)
        if (opts.sbm && update_sbm && _N > 1) {
            _prior.mcmc_sweep(rng, INF, 1.0, 0.01, total_dS, nattempts, nmoves);
        }
    }

    // Finite-temperature posterior MCMC sweep (beta = 1)
    void mcmc_sweep(Rng& rng, double beta, const std::vector<uint64_t>& cands,
                    const PriorOpts& opts, double w_step, double t_step,
                    bool update_sbm,
                    double& total_dS, size_t& nattempts, size_t& nmoves) {
        // 1. Node parameters MH sweep
        for (int i = 0; i < _N; ++i) {
            for (int l = 0; l < _q_t; ++l) {
                ++nattempts;
                double t0 = _theta[i][l];
                double t1 = t0 + (rng.unif() * 2.0 - 1.0) * t_step;
                if (t1 < _tmin[l] || t1 > _tmax[l]) continue;
                double dS = theta_move_dS(i, l, t1, opts);
                if (!std::isfinite(dS)) continue;
                if (dS <= 0.0 || rng.unif() < std::exp(-beta * dS)) {
                    apply_theta_move(i, l, t1);
                    total_dS += dS;
                    ++nmoves;
                }
            }
        }

        // 2. Weight update sweep on active edges
        std::vector<uint64_t> act = _active;
        rng.shuffle(act);
        for (uint64_t k : act) {
            auto [u, v] = unpack_key(k);
            ++nattempts;
            double w0 = get_x(u, v);
            double w1 = w0 + (rng.unif() * 2.0 - 1.0) * w_step;
            if (w1 == 0.0 || w1 < _xmin || w1 > _xmax) continue;
            double dS = edge_move_dS(u, v, true, w1, opts);
            if (!std::isfinite(dS)) continue;
            if (dS <= 0.0 || rng.unif() < std::exp(-beta * dS)) {
                apply_edge_move(u, v, true, w1, opts);
                total_dS += dS;
                ++nmoves;
            }
        }

        // 3. Reversible add/delete edge moves over candidate dyads (symmetric fixed list!)
        // Propose weight w_new ~ Uniform([lo_p, hi_p]) on add, reverse is delete with prob 1/(hi_p - lo_p)
        double lo_p = (_model == DYN_EPIDEMICS) ? std::max(_xmin, -5.0) : std::max(_xmin, -3.0);
        double hi_p = (_model == DYN_EPIDEMICS) ? std::min(_xmax, -1e-4) : std::min(_xmax, 3.0);
        double span_p = std::max(1e-6, hi_p - lo_p);
        double log_span = std::log(span_p);

        std::vector<uint64_t> ord = cands;
        rng.shuffle(ord);
        for (uint64_t k : ord) {
            auto [u, v] = unpack_key(k);
            ++nattempts;
            bool cur_pres = has_edge(u, v);
            if (!cur_pres) {
                double w1 = lo_p + rng.unif() * span_p;
                if (w1 == 0.0) continue;
                double dS = edge_move_dS(u, v, true, w1, opts);
                if (!std::isfinite(dS)) continue;
                // P_fwd(0 -> w1) = 1 / span_p; P_rev(w1 -> 0) = 1 -> log(P_rev / P_fwd) = +log(span_p)
                double log_a = -beta * dS + log_span;
                if (log_a >= 0.0 || rng.unif() < std::exp(log_a)) {
                    apply_edge_move(u, v, true, w1, opts);
                    total_dS += dS;
                    ++nmoves;
                }
            } else {
                double w0 = get_x(u, v);
                if (w0 < lo_p || w0 > hi_p) continue;  // maintain exact reversibility
                double dS = edge_move_dS(u, v, false, 0.0, opts);
                if (!std::isfinite(dS)) continue;
                // P_fwd(w0 -> 0) = 1; P_rev(0 -> w0) = 1 / span_p -> log(P_rev / P_fwd) = -log(span_p)
                double log_a = -beta * dS - log_span;
                if (log_a >= 0.0 || rng.unif() < std::exp(log_a)) {
                    apply_edge_move(u, v, false, 0.0, opts);
                    total_dS += dS;
                    ++nmoves;
                }
            }
        }

        // 4. SBM partition sweep
        if (opts.sbm && update_sbm && _N > 1) {
            _prior.mcmc_sweep(rng, beta, 1.0, 0.01, total_dS, nattempts, nmoves);
        }
    }

private:
    int _model, _N;
    bool _directed, _self_loops;
    int _q_t;
    bool _has_zero, _exposed, _positive;
    double _pslack;
    double _xmin, _xmax;
    std::vector<double> _tmin, _tmax;
    int64_t _M_dyads = 0;
    int _M = 0;

    DynBlockPrior _prior;
    gt_hash_map<uint64_t, double> _x;
    std::vector<gt_hash_map<int, double>> _adj_in, _adj_out;
    std::vector<double> _abs_deg;
    std::vector<std::vector<double>> _theta;
    std::vector<double> _node_nll;
    double _total_nll = 0.0;

    // Flattened data [i * _M + m]
    std::vector<double> _s_prev, _s_next;
    std::vector<int> _si_prev, _si_next;
    int _f_q = 0;
    std::vector<double> _f;

    // Cached local fields:
    // - Ising / CIsing / Epidemics / Normal / Linear: _h[i * _M + m] = sum_j W_{ij} s_prev[j, m]
    //   (for Epidemics: sum_{j->i, s_prev[j,m]==1} x_{ji})
    // - Potts: _h_potts[(i * _f_q + r) * _M + m] = sum_j W_{ij} f[r, s_prev[j, m]]
    std::vector<double> _h, _h_potts;

    std::vector<uint64_t> _active;
    gt_hash_map<uint64_t, size_t> _active_pos;

    void active_insert(uint64_t k) {
        if (_active_pos.find(k) != _active_pos.end()) return;
        _active_pos[k] = _active.size();
        _active.push_back(k);
    }

    void active_erase(uint64_t k) {
        auto it = _active_pos.find(k);
        if (it == _active_pos.end()) return;
        size_t pos = it->second;
        uint64_t last = _active.back();
        _active[pos] = last;
        _active_pos[last] = pos;
        _active.pop_back();
        _active_pos.erase(it);
    }

    void recompute_vertex_field(int i) {
        if (_model == DYN_PSEUDO_POTTS || _model == DYN_POTTS_GLAUBER) {
            double* hi = &_h_potts[static_cast<size_t>(i) * _f_q * _M];
            std::fill(hi, hi + static_cast<size_t>(_f_q) * _M, 0.0);
            for (const auto& [j, w] : _adj_in[i]) {
                const int* sj = &_si_prev[static_cast<size_t>(j) * _M];
                for (int r = 0; r < _f_q; ++r) {
                    double* hir = hi + static_cast<size_t>(r) * _M;
                    const double* fr = &_f[r * _f_q];
                    for (int m = 0; m < _M; ++m) hir[m] += w * fr[sj[m]];
                }
            }
            return;
        }
        double* hi = &_h[static_cast<size_t>(i) * _M];
        std::fill(hi, hi + _M, 0.0);
        if (_model == DYN_EPIDEMICS) {
            for (const auto& [j, w] : _adj_in[i]) {
                const double* sj = &_s_prev[static_cast<size_t>(j) * _M];
                for (int m = 0; m < _M; ++m) {
                    if (sj[m] == 1.0) hi[m] += w;
                }
            }
        } else {
            for (const auto& [j, w] : _adj_in[i]) {
                const double* sj = &_s_prev[static_cast<size_t>(j) * _M];
                for (int m = 0; m < _M; ++m) hi[m] += w * sj[m];
            }
        }
    }

    void apply_vertex_field_delta(int i, int j, double dw) {
        if (_model == DYN_PSEUDO_POTTS || _model == DYN_POTTS_GLAUBER) {
            double* hi = &_h_potts[static_cast<size_t>(i) * _f_q * _M];
            const int* sj = &_si_prev[static_cast<size_t>(j) * _M];
            for (int r = 0; r < _f_q; ++r) {
                double* hir = hi + static_cast<size_t>(r) * _M;
                const double* fr = &_f[r * _f_q];
                for (int m = 0; m < _M; ++m) hir[m] += dw * fr[sj[m]];
            }
            return;
        }
        double* hi = &_h[static_cast<size_t>(i) * _M];
        const double* sj = &_s_prev[static_cast<size_t>(j) * _M];
        if (_model == DYN_EPIDEMICS) {
            for (int m = 0; m < _M; ++m) {
                if (sj[m] == 1.0) hi[m] += dw;
            }
        } else {
            for (int m = 0; m < _M; ++m) hi[m] += dw * sj[m];
        }
    }

    double eval_vertex_nll_with_edge_delta(int i, int j, double dw, double d_abs) const {
        if (_model == DYN_PSEUDO_POTTS || _model == DYN_POTTS_GLAUBER) {
            const double* hi = &_h_potts[static_cast<size_t>(i) * _f_q * _M];
            const int* sj = &_si_prev[static_cast<size_t>(j) * _M];
            const int* si = &_si_next[static_cast<size_t>(i) * _M];
            const auto& thi = _theta[i];
            double nll = 0.0;
            std::vector<double> lg(_f_q);
            for (int m = 0; m < _M; ++m) {
                double mx = -INF;
                int sjm = sj[m];
                for (int r = 0; r < _f_q; ++r) {
                    double v = hi[static_cast<size_t>(r) * _M + m] + dw * _f[r * _f_q + sjm] + thi[r];
                    lg[r] = v;
                    if (v > mx) mx = v;
                }
                double sm = 0.0;
                for (int r = 0; r < _f_q; ++r) sm += std::exp(lg[r] - mx);
                nll += -(lg[si[m]] - (mx + std::log(sm)));
            }
            return nll;
        }

        const double* hi = &_h[static_cast<size_t>(i) * _M];
        const double* sj = &_s_prev[static_cast<size_t>(j) * _M];
        const double* sp = &_s_prev[static_cast<size_t>(i) * _M];
        const double* sn = &_s_next[static_cast<size_t>(i) * _M];

        if (_model == DYN_PSEUDO_ISING || _model == DYN_ISING_GLAUBER) {
            double th = _theta[i][0];
            double nll = 0.0;
            if (_has_zero) {
                for (int m = 0; m < _M; ++m) {
                    double fld = hi[m] + dw * sj[m] + th;
                    nll += -sn[m] * fld + log_1p2cosh(fld);
                }
            } else {
                for (int m = 0; m < _M; ++m) {
                    double fld = hi[m] + dw * sj[m] + th;
                    nll += -sn[m] * fld + log_2cosh(fld);
                }
            }
            return nll;
        }

        if (_model == DYN_PSEUDO_CISING || _model == DYN_CISING_GLAUBER) {
            double th = _theta[i][0];
            double nll = 0.0;
            for (int m = 0; m < _M; ++m) {
                double fld = hi[m] + dw * sj[m] + th;
                nll += -sn[m] * fld + log_2sinh_over_x(fld);
            }
            return nll;
        }

        if (_model == DYN_EPIDEMICS) {
            double th = _theta[i][0];
            double log_1m_r = log1mexp_d(th);  // log(1 - exp(theta_i))
            double target_inf = _exposed ? -1.0 : 1.0;
            double nll = 0.0;
            for (int m = 0; m < _M; ++m) {
                if (sp[m] == 0.0) {
                    double xsum = hi[m] + ((sj[m] == 1.0) ? dw : 0.0);
                    double log_stay = log_1m_r + xsum;  // <= 0
                    if (sn[m] == 0.0) {
                        nll -= log_stay;
                    } else if (sn[m] == target_inf) {
                        nll -= log1mexp_d(log_stay);
                    }
                }
            }
            return nll;
        }

        if (_model == DYN_PSEUDO_NORMAL || _model == DYN_NORMAL_GLAUBER) {
            double th0 = _theta[i][0];
            double mu = (_model == DYN_PSEUDO_NORMAL) ? _theta[i][1] : 0.0;
            if (_model == DYN_PSEUDO_NORMAL && _positive) {
                double k_abs = _abs_deg[i] + d_abs;
                if (k_abs > 0.0) th0 = std::min(-0.5 * std::log(k_abs) - _pslack, th0);
            }
            double s2 = std::exp(2.0 * th0);
            double inv_s2 = std::exp(-2.0 * th0);
            double c = th0 + HALF_LOG_2PI;
            double nll = 0.0;
            for (int m = 0; m < _M; ++m) {
                double ws = hi[m] + dw * sj[m];
                double mc = mu - s2 * ws;
                double diff = sn[m] - mc;
                nll += 0.5 * diff * diff * inv_s2 + c;
            }
            return nll;
        }

        if (_model == DYN_LINEAR_NORMAL) {
            double th0 = _theta[i][0];
            double inv_s2 = std::exp(-2.0 * th0);
            double c = th0 + HALF_LOG_2PI;
            double nll = 0.0;
            for (int m = 0; m < _M; ++m) {
                double ws = hi[m] + dw * sj[m];
                double diff = sn[m] - sp[m] - ws;
                nll += 0.5 * diff * diff * inv_s2 + c;
            }
            return nll;
        }
        return 0.0;
    }

    double compute_vertex_nll(int i) const {
        return eval_vertex_nll_with_edge_delta(i, 0, 0.0, 0.0);
    }
};

PriorOpts parse_prior_opts(Rcpp::List opts_list, int q_t) {
    PriorOpts p;
    if (opts_list.containsElementNamed("sbm")) p.sbm = Rcpp::as<bool>(opts_list["sbm"]);
    if (opts_list.containsElementNamed("latent_edges")) p.latent_edges = Rcpp::as<bool>(opts_list["latent_edges"]);
    if (opts_list.containsElementNamed("density")) p.density = Rcpp::as<bool>(opts_list["density"]);
    if (opts_list.containsElementNamed("aE")) p.aE = Rcpp::as<double>(opts_list["aE"]);
    if (opts_list.containsElementNamed("xdist")) p.xdist = Rcpp::as<bool>(opts_list["xdist"]);
    if (opts_list.containsElementNamed("xdist_uniform")) p.xdist_uniform = Rcpp::as<bool>(opts_list["xdist_uniform"]);
    if (opts_list.containsElementNamed("xl1")) p.xl1 = Rcpp::as<double>(opts_list["xl1"]);
    if (opts_list.containsElementNamed("xdelta")) p.xdelta = Rcpp::as<double>(opts_list["xdelta"]);
    if (opts_list.containsElementNamed("normal")) p.normal_prior = Rcpp::as<bool>(opts_list["normal"]);
    if (opts_list.containsElementNamed("mu")) p.normal_mu = Rcpp::as<double>(opts_list["mu"]);
    if (opts_list.containsElementNamed("sigma")) p.normal_sigma = Rcpp::as<double>(opts_list["sigma"]);

    auto read_bool_vec = [&](const char* nm, bool def) {
        std::vector<bool> v(q_t, def);
        if (opts_list.containsElementNamed(nm)) {
            Rcpp::LogicalVector lv = opts_list[nm];
            for (int l = 0; l < q_t; ++l) v[l] = (lv[l % lv.size()] != 0);
        }
        return v;
    };
    auto read_dbl_vec = [&](const char* nm, double def) {
        std::vector<double> v(q_t, def);
        if (opts_list.containsElementNamed(nm)) {
            Rcpp::NumericVector nv = opts_list[nm];
            for (int l = 0; l < q_t; ++l) v[l] = nv[l % nv.size()];
        }
        return v;
    };
    p.tdist = read_bool_vec("tdist", false);
    p.tdist_uniform = read_bool_vec("tdist_uniform", false);
    p.tl1 = read_dbl_vec("tl1", 1.0);
    p.tdelta = read_dbl_vec("tdelta", 1e-8);
    return p;
}

DynamicsEngine build_dynamics_engine(
    int model, int N, bool directed, bool self_loops,
    Rcpp::List s_list, Rcpp::NumericMatrix f_mat,
    bool has_zero, bool exposed, bool positive, double pslack,
    double xmin, double xmax, Rcpp::NumericVector tmin, Rcpp::NumericVector tmax,
    Rcpp::IntegerMatrix init_edges, Rcpp::NumericVector init_x,
    Rcpp::NumericMatrix init_theta, Rcpp::List bs_list,
    bool dc, bool dist, double q_max) {
    int q_t = init_theta.ncol();
    int sbm_L = static_cast<int>(bs_list.size());
    size_t qm = std::isinf(q_max) ? 5000 : static_cast<size_t>(q_max);
    std::vector<double> tmin_v(tmin.begin(), tmin.end());
    std::vector<double> tmax_v(tmax.begin(), tmax.end());
    DynamicsEngine eng(model, N, directed, self_loops, q_t, has_zero, exposed,
                       positive, pslack, xmin, xmax, tmin_v, tmax_v,
                       dc, dist, qm, sbm_L);
    eng.set_data(s_list, f_mat);

    std::vector<std::tuple<int, int, double>> edges;
    edges.reserve(init_edges.nrow());
    for (int i = 0; i < init_edges.nrow(); ++i) {
        edges.emplace_back(init_edges(i, 0) - 1, init_edges(i, 1) - 1, init_x[i]);
    }
    std::vector<std::vector<double>> th(N, std::vector<double>(q_t, 0.0));
    for (int i = 0; i < N; ++i)
        for (int l = 0; l < q_t; ++l) th[i][l] = init_theta(i, l);

    std::vector<std::vector<int>> bs;
    bs.reserve(sbm_L);
    for (int l = 0; l < sbm_L; ++l) {
        Rcpp::IntegerVector x = bs_list[l];
        std::vector<int> v(x.size());
        for (R_xlen_t i = 0; i < x.size(); ++i) v[i] = x[i] - 1;
        bs.push_back(std::move(v));
    }
    eng.init_state(edges, th, bs);
    return eng;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_dynamics_entropy(
    int model, int N, bool directed, bool self_loops,
    Rcpp::List s_list, Rcpp::NumericMatrix f_mat,
    bool has_zero, bool exposed, bool positive, double pslack,
    double xmin, double xmax, Rcpp::NumericVector tmin, Rcpp::NumericVector tmax,
    Rcpp::IntegerMatrix init_edges, Rcpp::NumericVector init_x,
    Rcpp::NumericMatrix init_theta, Rcpp::List bs_list,
    bool dc, bool dist, Rcpp::List prior_opts, double q_max) {
    DynamicsEngine eng = build_dynamics_engine(
        model, N, directed, self_loops, s_list, f_mat, has_zero, exposed, positive,
        pslack, xmin, xmax, tmin, tmax, init_edges, init_x, init_theta, bs_list,
        dc, dist, q_max);
    PriorOpts opts = parse_prior_opts(prior_opts, init_theta.ncol());
    double nll = opts.latent_edges ? eng.data_nll() : 0.0;
    double g_pr = eng.geom_edge_prior(static_cast<long>(eng.weights().size()), opts.density, opts.aE);
    double x_pr = eng.x_prior(opts);
    double t_pr = eng.t_prior(opts);
    double sbm_s = opts.sbm ? eng.prior().entropy() : 0.0;
    return Rcpp::List::create(
        Rcpp::Named("total") = nll + g_pr + x_pr + t_pr + sbm_s,
        Rcpp::Named("nll") = nll,
        Rcpp::Named("geom_prior") = g_pr,
        Rcpp::Named("x_prior") = x_pr,
        Rcpp::Named("t_prior") = t_pr,
        Rcpp::Named("sbm") = sbm_s);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_dynamics_check_moves(
    int model, int N, bool directed, bool self_loops,
    Rcpp::List s_list, Rcpp::NumericMatrix f_mat,
    bool has_zero, bool exposed, bool positive, double pslack,
    double xmin, double xmax, Rcpp::NumericVector tmin, Rcpp::NumericVector tmax,
    Rcpp::IntegerMatrix init_edges, Rcpp::NumericVector init_x,
    Rcpp::NumericMatrix init_theta, Rcpp::List bs_list,
    bool dc, bool dist, Rcpp::List prior_opts, int nmoves, double q_max) {
    DynamicsEngine eng = build_dynamics_engine(
        model, N, directed, self_loops, s_list, f_mat, has_zero, exposed, positive,
        pslack, xmin, xmax, tmin, tmax, init_edges, init_x, init_theta, bs_list,
        dc, dist, q_max);
    PriorOpts opts = parse_prior_opts(prior_opts, init_theta.ncol());
    Rng rng;
    Rcpp::NumericVector errs(nmoves);
    double S = eng.total_entropy(opts);
    for (int it = 0; it < nmoves; ++it) {
        double r = rng.unif();
        if (r < 0.5) {
            int u = static_cast<int>(rng.index(static_cast<size_t>(N)));
            int v = static_cast<int>(rng.index(static_cast<size_t>(N)));
            if (!self_loops && u == v) v = (v + 1) % N;
            bool cur = eng.has_edge(u, v);
            bool new_p = cur ? (rng.unif() < 0.5) : true;
            double new_w = 0.0;
            if (new_p) {
                if (model == DYN_EPIDEMICS) new_w = -0.1 - 1.5 * rng.unif();
                else new_w = (rng.unif() * 2.0 - 1.0) * 1.2;
                if (new_w == 0.0) new_w = 0.3;
            }
            double dS = eng.edge_move_dS(u, v, new_p, new_w, opts);
            if (!std::isfinite(dS)) continue;
            eng.apply_edge_move(u, v, new_p, new_w, opts);
            double S_new = eng.total_entropy(opts);
            errs[it] = std::abs(dS - (S_new - S));
            S = S_new;
        } else if (r < 0.8) {
            int i = static_cast<int>(rng.index(static_cast<size_t>(N)));
            int l = static_cast<int>(rng.index(static_cast<size_t>(eng.q_t())));
            double t_new = (model == DYN_EPIDEMICS) ? (-0.2 - 2.0 * rng.unif()) : (rng.unif() * 1.0 - 0.5);
            double dS = eng.theta_move_dS(i, l, t_new, opts);
            if (!std::isfinite(dS)) continue;
            eng.apply_theta_move(i, l, t_new);
            double S_new = eng.total_entropy(opts);
            errs[it] = std::abs(dS - (S_new - S));
            S = S_new;
        } else if (opts.sbm) {
            int v = static_cast<int>(rng.index(static_cast<size_t>(N)));
            int s = eng.prior().base().random_occupied(rng);
            int cur_g = eng.prior().base().group(v);
            if (s == cur_g) continue;
            if (eng.prior().base().size(cur_g) == 1 && !eng.prior().same_parent(0, cur_g, s)) continue;
            double dS = eng.prior().vertex_move_dS(v, s);
            eng.prior().vertex_move(v, s);
            double S_new = eng.total_entropy(opts);
            errs[it] = std::abs(dS - (S_new - S));
            S = S_new;
        }
    }
    return errs;
}

// [[Rcpp::export]]
Rcpp::List cpp_dynamics_reconstruct(
    int model, int N, bool directed, bool self_loops,
    Rcpp::List s_list, Rcpp::NumericMatrix f_mat,
    bool has_zero, bool exposed, bool positive, double pslack,
    double xmin, double xmax, Rcpp::NumericVector tmin, Rcpp::NumericVector tmax,
    Rcpp::IntegerMatrix init_edges, Rcpp::NumericVector init_x,
    Rcpp::NumericMatrix init_theta, Rcpp::List bs_list,
    bool dc, bool dist, Rcpp::List prior_opts,
    int niter, int nsamples, int burnin, int thin, double beta,
    int k_cand, double w_step, double t_step, double tol, double q_max) {
    DynamicsEngine eng = build_dynamics_engine(
        model, N, directed, self_loops, s_list, f_mat, has_zero, exposed, positive,
        pslack, xmin, xmax, tmin, tmax, init_edges, init_x, init_theta, bs_list,
        dc, dist, q_max);
    PriorOpts opts = parse_prior_opts(prior_opts, init_theta.ncol());
    Rng rng;
    std::vector<uint64_t> cands = eng.build_candidates(k_cand);

    double total_dS = 0.0;
    size_t nattempts = 0, nmoves = 0;
    std::vector<double> entropies;

    // Phase 1: Greedy MAP optimization sweeps (when niter > 0)
    double prev_S = eng.total_entropy(opts);
    for (int it = 0; it < niter; ++it) {
        double step_dS = 0.0;
        eng.greedy_sweep(rng, cands, opts, true, step_dS, nattempts, nmoves);
        total_dS += step_dS;
        double cur_S = eng.total_entropy(opts);
        entropies.push_back(cur_S);
        if ((it & 3) == 0) Rcpp::checkUserInterrupt();
        if (it >= 2 && std::abs(prev_S - cur_S) < tol && step_dS >= -tol) break;
        prev_S = cur_S;
    }

    // Phase 2: Optional posterior MCMC sampling (when nsamples > 0)
    gt_hash_map<uint64_t, long> edge_cnt;
    gt_hash_map<uint64_t, double> edge_wsum, edge_w2sum;
    std::vector<std::vector<double>> th_sum(N, std::vector<double>(eng.q_t(), 0.0));
    std::vector<std::vector<double>> th_2sum(N, std::vector<double>(eng.q_t(), 0.0));

    if (nsamples > 0) {
        double beta_samp = std::isinf(beta) ? 1.0 : beta;
        for (int b = 0; b < burnin; ++b) {
            eng.mcmc_sweep(rng, beta_samp, cands, opts, w_step, t_step, true, total_dS, nattempts, nmoves);
            if ((b & 7) == 0) Rcpp::checkUserInterrupt();
        }
        int thin_n = std::max(1, thin);
        for (int s = 0; s < nsamples; ++s) {
            for (int t = 0; t < thin_n; ++t) {
                eng.mcmc_sweep(rng, beta_samp, cands, opts, w_step, t_step, true, total_dS, nattempts, nmoves);
            }
            for (const auto& [k, w] : eng.weights()) {
                edge_cnt[k] += 1;
                edge_wsum[k] += w;
                edge_w2sum[k] += w * w;
            }
            const auto& th = eng.theta();
            for (int i = 0; i < N; ++i) {
                for (int l = 0; l < eng.q_t(); ++l) {
                    th_sum[i][l] += th[i][l];
                    th_2sum[i][l] += th[i][l] * th[i][l];
                }
            }
            entropies.push_back(eng.total_entropy(opts));
            if ((s & 7) == 0) Rcpp::checkUserInterrupt();
        }
    }

    // Export final MAP/MCMC edges
    const auto& fw = eng.weights();
    Rcpp::IntegerMatrix fin_edges(fw.size(), 2);
    Rcpp::NumericVector fin_x(fw.size());
    size_t fi = 0;
    for (const auto& [k, w] : fw) {
        auto [u, v] = unpack_key(k);
        fin_edges(fi, 0) = u + 1;
        fin_edges(fi, 1) = v + 1;
        fin_x[fi] = w;
        ++fi;
    }

    // Export final theta
    Rcpp::NumericMatrix fin_theta(N, eng.q_t());
    for (int i = 0; i < N; ++i)
        for (int l = 0; l < eng.q_t(); ++l) fin_theta(i, l) = eng.theta()[i][l];

    // Export marginal edges (if nsamples > 0, else MAP edges with prob = 1)
    size_t nm = (nsamples > 0) ? edge_cnt.size() : fw.size();
    Rcpp::IntegerVector m_u(nm), m_v(nm);
    Rcpp::NumericVector m_prob(nm), m_w(nm), m_wsd(nm);
    if (nsamples > 0) {
        size_t idx = 0;
        for (const auto& [k, cnt] : edge_cnt) {
            auto [u, v] = unpack_key(k);
            m_u[idx] = u + 1;
            m_v[idx] = v + 1;
            m_prob[idx] = static_cast<double>(cnt) / static_cast<double>(nsamples);
            double mw = edge_wsum[k] / static_cast<double>(cnt);
            double var = std::max(0.0, edge_w2sum[k] / static_cast<double>(cnt) - mw * mw);
            m_w[idx] = mw;
            m_wsd[idx] = std::sqrt(var);
            ++idx;
        }
    } else {
        size_t idx = 0;
        for (const auto& [k, w] : fw) {
            auto [u, v] = unpack_key(k);
            m_u[idx] = u + 1;
            m_v[idx] = v + 1;
            m_prob[idx] = 1.0;
            m_w[idx] = w;
            m_wsd[idx] = 0.0;
            ++idx;
        }
    }

    Rcpp::NumericMatrix th_mean(N, eng.q_t()), th_sd(N, eng.q_t());
    if (nsamples > 0) {
        for (int i = 0; i < N; ++i) {
            for (int l = 0; l < eng.q_t(); ++l) {
                double m = th_sum[i][l] / static_cast<double>(nsamples);
                double v = std::max(0.0, th_2sum[i][l] / static_cast<double>(nsamples) - m * m);
                th_mean(i, l) = m;
                th_sd(i, l) = std::sqrt(v);
            }
        }
    } else {
        th_mean = Rcpp::clone(fin_theta);
    }

    return Rcpp::List::create(
        Rcpp::Named("from") = m_u,
        Rcpp::Named("to") = m_v,
        Rcpp::Named("prob") = m_prob,
        Rcpp::Named("weight") = m_w,
        Rcpp::Named("weight_sd") = m_wsd,
        Rcpp::Named("final_edges") = fin_edges,
        Rcpp::Named("final_x") = fin_x,
        Rcpp::Named("final_theta") = fin_theta,
        Rcpp::Named("theta_mean") = th_mean,
        Rcpp::Named("theta_sd") = th_sd,
        Rcpp::Named("final_bs") = eng.prior().export_bs(),
        Rcpp::Named("entropies") = entropies,
        Rcpp::Named("entropy") = eng.total_entropy(opts),
        Rcpp::Named("nll") = eng.data_nll(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves));
}

// ---------------------------------------------------------------------------
// Native C++ forward dynamical simulators

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_sim_ising(
    SEXP pointer, Rcpp::NumericVector weights, Rcpp::NumericVector theta,
    int n_samples, int burnin, int thin, bool async_eq, bool has_zero,
    Rcpp::IntegerVector init_state) {
    Graph& graph = get_graph(pointer);
    const Adjacency& g = graph.adjacency;
    int N = static_cast<int>(boost::num_vertices(g));
    bool directed = graph.directed;

    std::vector<std::vector<std::pair<int, double>>> adj_in(N);
    auto er = boost::edges(g);
    int ei = 0;
    for (auto it = er.first; it != er.second; ++it, ++ei) {
        int u = static_cast<int>(boost::source(*it, g));
        int v = static_cast<int>(boost::target(*it, g));
        double w = weights[ei];
        adj_in[v].emplace_back(u, w);
        if (!directed && u != v) adj_in[u].emplace_back(v, w);
    }

    std::vector<int> s(N);
    for (int i = 0; i < N; ++i) s[i] = init_state[i];

    auto draw_spin = [&](double h) -> int {
        if (!has_zero) {
            double p_pos = 1.0 / (1.0 + std::exp(-2.0 * h));
            return (unif_rand() < p_pos) ? 1 : -1;
        }
        double l_neg = -h, l_zero = 0.0, l_pos = h;
        double mx = std::max({l_neg, l_zero, l_pos});
        double w_neg = std::exp(l_neg - mx), w_zero = std::exp(l_zero - mx), w_pos = std::exp(l_pos - mx);
        double u = unif_rand() * (w_neg + w_zero + w_pos);
        if (u < w_neg) return -1;
        if (u < w_neg + w_zero) return 0;
        return 1;
    };

    std::vector<int> s_next(N);
    auto do_step = [&]() {
        if (async_eq) {
            for (int k = 0; k < N; ++k) {
                int i = static_cast<int>(R_unif_index(static_cast<double>(N)));
                double h = theta[i];
                for (const auto& [j, w] : adj_in[i]) h += w * static_cast<double>(s[j]);
                s[i] = draw_spin(h);
            }
        } else {
            for (int i = 0; i < N; ++i) {
                double h = theta[i];
                for (const auto& [j, w] : adj_in[i]) h += w * static_cast<double>(s[j]);
                s_next[i] = draw_spin(h);
            }
            s.swap(s_next);
        }
    };

    for (int b = 0; b < burnin; ++b) do_step();
    Rcpp::IntegerMatrix out(N, n_samples);
    int thin_n = std::max(1, thin);
    for (int m = 0; m < n_samples; ++m) {
        if (m > 0) {
            for (int t = 0; t < thin_n; ++t) do_step();
        }
        for (int i = 0; i < N; ++i) out(i, m) = s[i];
        if ((m & 63) == 0) Rcpp::checkUserInterrupt();
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_sim_epidemics(
    SEXP pointer, Rcpp::NumericVector log1m_t, Rcpp::NumericVector spont,
    Rcpp::NumericVector gamma_rec, Rcpp::NumericVector mu_resusc,
    int n_steps, bool exposed, double expose_rate, Rcpp::IntegerVector init_state) {
    Graph& graph = get_graph(pointer);
    const Adjacency& g = graph.adjacency;
    int N = static_cast<int>(boost::num_vertices(g));
    bool directed = graph.directed;

    std::vector<std::vector<std::pair<int, double>>> adj_in(N);
    auto er = boost::edges(g);
    int ei = 0;
    for (auto it = er.first; it != er.second; ++it, ++ei) {
        int u = static_cast<int>(boost::source(*it, g));
        int v = static_cast<int>(boost::target(*it, g));
        double x = log1m_t[ei];  // log(1 - t_e) <= 0
        adj_in[v].emplace_back(u, x);
        if (!directed && u != v) adj_in[u].emplace_back(v, x);
    }

    std::vector<int> s(N), sn(N);
    for (int i = 0; i < N; ++i) s[i] = init_state[i];

    Rcpp::IntegerMatrix out(N, n_steps);
    for (int t = 0; t < n_steps; ++t) {
        for (int i = 0; i < N; ++i) out(i, t) = s[i];
        if (t + 1 == n_steps) break;
        for (int i = 0; i < N; ++i) {
            if (s[i] == 0) {
                double log_stay = std::log1p(-spont[i]);
                for (const auto& [j, x] : adj_in[i]) {
                    if (s[j] == 1) log_stay += x;
                }
                double p_inf = -std::expm1(log_stay);
                sn[i] = (unif_rand() < p_inf) ? (exposed ? -1 : 1) : 0;
            } else if (s[i] == -1) {
                sn[i] = (unif_rand() < expose_rate) ? 1 : -1;
            } else if (s[i] == 1) {
                if (gamma_rec[i] > 0.0 && unif_rand() < gamma_rec[i]) {
                    sn[i] = (mu_resusc[i] < 0.0) ? 0 : 2;  // -1 signals SIS (back to 0), >=0 signals SIR/SIRS (to 2)
                } else {
                    sn[i] = 1;
                }
            } else {
                sn[i] = (mu_resusc[i] > 0.0 && unif_rand() < mu_resusc[i]) ? 0 : 2;
            }
        }
        s.swap(sn);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_sim_normal_dyn(
    SEXP pointer, Rcpp::NumericVector weights, Rcpp::NumericVector sigma,
    Rcpp::NumericVector mu, int n_samples, int burnin, int thin, int mode) {
    // mode: 0 = equilibrium Gibbs (PseudoNormal), 1 = Glauber sync, 2 = LinearNormal
    Graph& graph = get_graph(pointer);
    const Adjacency& g = graph.adjacency;
    int N = static_cast<int>(boost::num_vertices(g));
    bool directed = graph.directed;

    std::vector<std::vector<std::pair<int, double>>> adj_in(N);
    auto er = boost::edges(g);
    int ei = 0;
    for (auto it = er.first; it != er.second; ++it, ++ei) {
        int u = static_cast<int>(boost::source(*it, g));
        int v = static_cast<int>(boost::target(*it, g));
        double w = weights[ei];
        adj_in[v].emplace_back(u, w);
        if (!directed && u != v) adj_in[u].emplace_back(v, w);
    }

    std::vector<double> s(N), sn(N);
    for (int i = 0; i < N; ++i) s[i] = mu[i] + sigma[i] * norm_rand();

    auto do_step = [&]() {
        if (mode == 0) {
            for (int k = 0; k < N; ++k) {
                int i = static_cast<int>(R_unif_index(static_cast<double>(N)));
                double ws = 0.0;
                for (const auto& [j, w] : adj_in[i]) ws += w * s[j];
                double mc = mu[i] - sigma[i] * sigma[i] * ws;
                s[i] = mc + sigma[i] * norm_rand();
            }
        } else if (mode == 1) {
            for (int i = 0; i < N; ++i) {
                double ws = 0.0;
                for (const auto& [j, w] : adj_in[i]) ws += w * s[j];
                double mc = mu[i] - sigma[i] * sigma[i] * ws;
                sn[i] = mc + sigma[i] * norm_rand();
            }
            s.swap(sn);
        } else {
            for (int i = 0; i < N; ++i) {
                double ws = 0.0;
                for (const auto& [j, w] : adj_in[i]) ws += w * s[j];
                sn[i] = s[i] + ws + sigma[i] * norm_rand();
            }
            s.swap(sn);
        }
    };

    for (int b = 0; b < burnin; ++b) do_step();
    Rcpp::NumericMatrix out(N, n_samples);
    int thin_n = std::max(1, thin);
    for (int m = 0; m < n_samples; ++m) {
        if (m > 0) {
            for (int t = 0; t < thin_n; ++t) do_step();
        }
        for (int i = 0; i < N; ++i) out(i, m) = s[i];
        if ((m & 63) == 0) Rcpp::checkUserInterrupt();
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_sim_potts(
    SEXP pointer, Rcpp::NumericVector weights, Rcpp::NumericMatrix f_mat,
    Rcpp::NumericMatrix theta, int n_samples, int burnin, int thin, bool async_eq,
    Rcpp::IntegerVector init_state) {
    Graph& graph = get_graph(pointer);
    const Adjacency& g = graph.adjacency;
    int N = static_cast<int>(boost::num_vertices(g));
    int q = f_mat.nrow();
    bool directed = graph.directed;

    std::vector<std::vector<std::pair<int, double>>> adj_in(N);
    auto er = boost::edges(g);
    int ei = 0;
    for (auto it = er.first; it != er.second; ++it, ++ei) {
        int u = static_cast<int>(boost::source(*it, g));
        int v = static_cast<int>(boost::target(*it, g));
        double w = weights[ei];
        adj_in[v].emplace_back(u, w);
        if (!directed && u != v) adj_in[u].emplace_back(v, w);
    }

    std::vector<int> s(N), sn(N);
    for (int i = 0; i < N; ++i) s[i] = init_state[i];
    std::vector<double> probs(q);

    auto draw_cat = [&](int i) -> int {
        double mx = -INF;
        for (int r = 0; r < q; ++r) {
            double h = theta(i, r);
            for (const auto& [j, w] : adj_in[i]) h += w * f_mat(r, s[j]);
            probs[r] = h;
            if (h > mx) mx = h;
        }
        double tot = 0.0;
        for (int r = 0; r < q; ++r) { probs[r] = std::exp(probs[r] - mx); tot += probs[r]; }
        double u = unif_rand() * tot;
        double acc = 0.0;
        for (int r = 0; r < q; ++r) {
            acc += probs[r];
            if (u <= acc) return r;
        }
        return q - 1;
    };

    auto do_step = [&]() {
        if (async_eq) {
            for (int k = 0; k < N; ++k) {
                int i = static_cast<int>(R_unif_index(static_cast<double>(N)));
                s[i] = draw_cat(i);
            }
        } else {
            for (int i = 0; i < N; ++i) sn[i] = draw_cat(i);
            s.swap(sn);
        }
    };

    for (int b = 0; b < burnin; ++b) do_step();
    Rcpp::IntegerMatrix out(N, n_samples);
    int thin_n = std::max(1, thin);
    for (int m = 0; m < n_samples; ++m) {
        if (m > 0) {
            for (int t = 0; t < thin_n; ++t) do_step();
        }
        for (int i = 0; i < N; ++i) out(i, m) = s[i];
        if ((m & 63) == 0) Rcpp::checkUserInterrupt();
    }
    return out;
}
