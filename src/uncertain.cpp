// rgraphtool: native C++ uncertain network reconstruction, measured networks
// (homogeneous and heterogeneous errors), extraneous edge probability estimates,
// latent Poisson multigraphs (EM and Bayesian SBM), and conditional edge
// probabilities (Peixoto 2018 PRX, Peixoto 2020 PRE).
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

enum MeasModel {
    MODEL_MEASURED = 0,
    MODEL_MIXED = 1,
    MODEL_UNCERTAIN = 2,
    MODEL_LATENT_MULTI = 3
};

class UncertainEngine {
public:
    UncertainEngine(int model, int N, bool directed, bool self_loops,
                    bool dc, bool dist, size_t q_max, int levels)
        : _model(model), _N(N), _directed(directed), _self_loops(self_loops),
          _prior(N, directed, dc, dist, q_max, levels) {
        if (directed) _M_dyads = self_loops ? int64_t(N) * N : int64_t(N) * (N - 1);
        else _M_dyads = self_loops ? int64_t(N) * (N + 1) / 2 : int64_t(N) * (N - 1) / 2;
    }

    int N() const { return _N; }
    bool directed() const { return _directed; }
    bool self_loops() const { return _self_loops; }
    int model() const { return _model; }
    DynBlockPrior& prior() { return _prior; }
    const DynBlockPrior& prior() const { return _prior; }

    uint64_t key(int u, int v) const {
        if (!_directed && u > v) std::swap(u, v);
        return pair_key(u, v);
    }

    bool valid_dyad(int u, int v) const {
        if (u < 0 || u >= _N || v < 0 || v >= _N) return false;
        if (!_self_loops && u == v) return false;
        return true;
    }

    // Uniformly sample a valid dyad (u, v) from all _M_dyads with exact prob 1/_M_dyads
    std::pair<int, int> sample_uniform_dyad(Rng& rng) const {
        if (_directed) {
            int u = static_cast<int>(rng.index(static_cast<size_t>(_N)));
            if (_self_loops) {
                int v = static_cast<int>(rng.index(static_cast<size_t>(_N)));
                return {u, v};
            }
            int v = static_cast<int>(rng.index(static_cast<size_t>(_N - 1)));
            if (v >= u) ++v;
            return {u, v};
        } else {
            if (_self_loops) {
                // Exact bijection from [0..N] x [0..N-1] (size N(N+1)=2M) to 0<=u<=v<N
                int a = static_cast<int>(rng.index(static_cast<size_t>(_N + 1)));
                int b = static_cast<int>(rng.index(static_cast<size_t>(_N)));
                if (a <= b) return {a, b};
                return {b, a - 1};
            }
            int u = static_cast<int>(rng.index(static_cast<size_t>(_N)));
            int v = static_cast<int>(rng.index(static_cast<size_t>(_N - 1)));
            if (v >= u) ++v;
            if (u > v) std::swap(u, v);
            return {u, v};
        }
    }

    void setup_measured(const std::vector<std::tuple<int, int, long, long>>& obs,
                        long n_def, long x_def, double alpha, double beta,
                        double mu, double nu, double lp, double lq) {
        _n_def = n_def; _x_def = x_def;
        _alpha = alpha; _beta = beta; _mu = mu; _nu = nu;
        _lp = lp; _lq = lq;
        _obs_nx.clear();
        _obs_keys.clear();
        _N_meas = 0; _X_meas = 0; _lbinom_const = 0.0;
        for (const auto& [u0, v0, ni, xi] : obs) {
            if (!valid_dyad(u0, v0)) continue;
            uint64_t k = key(u0, v0);
            if (_obs_nx.find(k) == _obs_nx.end()) _obs_keys.push_back(k);
            _obs_nx[k] = {ni, xi};
        }
        for (const auto& [k, nx] : _obs_nx) {
            _N_meas += nx.first;
            _X_meas += nx.second;
            _lbinom_const += lbinom_d(static_cast<double>(nx.first), static_cast<double>(nx.second));
        }
        int64_t n_unobs = _M_dyads - static_cast<int64_t>(_obs_nx.size());
        _N_meas += n_unobs * n_def;
        _X_meas += n_unobs * x_def;
        _lbinom_const += static_cast<double>(n_unobs) *
                         lbinom_d(static_cast<double>(n_def), static_cast<double>(x_def));
        if (_model == MODEL_MIXED) sync_mixed_q();
    }

    void sync_mixed_q() {
        auto tf = [&](double n, double x) -> std::pair<double, double> {
            double ra = lbeta_d(n - x + _alpha, x + _beta) - lbeta_d(_alpha, _beta);
            double rb = lbeta_d(x + _mu, n - x + _nu) - lbeta_d(_mu, _nu);
            return {ra, rb};
        };
        _obs_q.clear();
        double rb_sum = 0.0;
        for (const auto& [k, nx] : _obs_nx) {
            auto [ra, rb] = tf(static_cast<double>(nx.first), static_cast<double>(nx.second));
            _obs_q[k] = ra - rb;
            rb_sum += rb;
        }
        auto [dra, drb] = tf(static_cast<double>(_n_def), static_cast<double>(_x_def));
        _q_def = dra - drb;
        int64_t n_unobs = _M_dyads - static_cast<int64_t>(_obs_nx.size());
        _S_const = static_cast<double>(n_unobs) * drb + rb_sum;
        recompute_q_sum();
    }

    void setup_uncertain(const std::vector<std::tuple<int, int, double>>& obs_prob,
                         double q_default) {
        _obs_q.clear();
        _obs_keys.clear();
        for (const auto& [u0, v0, qv] : obs_prob) {
            if (!valid_dyad(u0, v0)) continue;
            uint64_t k = key(u0, v0);
            if (_obs_q.find(k) == _obs_q.end()) _obs_keys.push_back(k);
            _obs_q[k] = qv;
        }
        double sum_q = 0.0;
        for (const auto& [k, qv] : _obs_q) sum_q += qv;
        int64_t n_unobs = _M_dyads - static_cast<int64_t>(_obs_q.size());
        double p = (sum_q + static_cast<double>(n_unobs) * q_default) / static_cast<double>(_M_dyads);
        double lp_odds = (p > 0.0 && p < 1.0) ? (std::log(p) - std::log1p(-p)) : 0.0;
        double s_c = 0.0;
        for (auto& [k, qv] : _obs_q) {
            if (qv < 1.0) s_c += std::log1p(-qv);
            if (qv <= 0.0) qv = -INF;
            else if (qv >= 1.0) qv = INF;
            else qv = std::log(qv) - std::log1p(-qv) - lp_odds;
        }
        if (q_default > 0.0 && q_default < 1.0) {
            _q_def = std::log(q_default) - std::log1p(-q_default) - lp_odds;
            s_c += static_cast<double>(n_unobs) * std::log1p(-q_default);
        } else if (q_default >= 1.0) {
            _q_def = INF;
        } else {
            _q_def = -INF;
        }
        if (p < 1.0) s_c -= static_cast<double>(_M_dyads) * std::log1p(-p);
        _S_const = s_c;
        recompute_q_sum();
    }

    void setup_latent_multi(const std::vector<std::pair<int, int>>& obs_edges) {
        _obs_q.clear();
        _obs_keys.clear();
        for (const auto& [u0, v0] : obs_edges) {
            if (!valid_dyad(u0, v0)) continue;
            uint64_t k = key(u0, v0);
            if (_obs_q.find(k) == _obs_q.end()) {
                _obs_q[k] = INF;
                _obs_keys.push_back(k);
            }
        }
        _q_def = -INF;
        _S_const = 0.0;
        _q_sum = 0.0;
    }

    void init_latent(const std::vector<std::tuple<int, int, long>>& latent_edges,
                     const std::vector<std::vector<int>>& bs) {
        _prior.init(latent_edges, bs);
        recompute_meas_stats();
        rebuild_active();
    }

    void rebuild_active() {
        _active.clear();
        _active_pos.clear();
        for (const auto& [k, w] : _prior.base().edge_weights()) {
            if (w > 0) {
                _active_pos[k] = _active.size();
                _active.push_back(k);
            }
        }
    }

    void recompute_meas_stats() {
        if (_model == MODEL_MEASURED) {
            _M_pos = 0; _T_pos = 0;
            for (const auto& [k, w] : _prior.base().edge_weights()) {
                if (w <= 0) continue;
                auto [u, v] = unpack_key(k);
                if (!valid_dyad(u, v)) continue;
                auto it = _obs_nx.find(k);
                if (it != _obs_nx.end()) {
                    _M_pos += it->second.first;
                    _T_pos += it->second.second;
                } else {
                    _M_pos += _n_def;
                    _T_pos += _x_def;
                }
            }
        } else if (_model == MODEL_MIXED || _model == MODEL_UNCERTAIN) {
            recompute_q_sum();
        }
    }

    void recompute_q_sum() {
        _q_sum = 0.0;
        _n_inf_q = 0;
        _n_ninf_q = 0;
        for (const auto& [k, w] : _prior.base().edge_weights()) {
            if (w <= 0) continue;
            auto [u, v] = unpack_key(k);
            if (!valid_dyad(u, v)) continue;
            double qv = get_q_val(k);
            if (std::isinf(qv)) {
                if (qv > 0) ++_n_inf_q;
                else ++_n_ninf_q;
            } else {
                _q_sum += qv;
            }
        }
    }

    double get_q_val(uint64_t k) const {
        auto it = _obs_q.find(k);
        return it == _obs_q.end() ? _q_def : it->second;
    }

    std::pair<long, long> get_nx_val(uint64_t k) const {
        auto it = _obs_nx.find(k);
        return it == _obs_nx.end() ? std::make_pair(_n_def, _x_def) : it->second;
    }

    double meas_entropy(bool density = false, double aE = 1.0) const {
        double S = 0.0;
        if (_model == MODEL_MEASURED) {
            S = meas_entropy_from(_M_pos, _T_pos);
        } else if (_model == MODEL_MIXED || _model == MODEL_UNCERTAIN) {
            if (_n_ninf_q > 0) return INF;
            S = -(_S_const + _q_sum);
        } else if (_model == MODEL_LATENT_MULTI) {
            S = 0.0;
        }
        if (density) {
            long E = _prior.base().E();
            S += static_cast<double>(E) * std::log(aE) - lgam(static_cast<double>(E) + 1.0) + aE;
        }
        return S;
    }

    double total_entropy(bool sbm = true, bool density = false, double aE = 1.0) const {
        double S = meas_entropy(density, aE);
        if (sbm) S += _prior.entropy();
        return S;
    }

    std::pair<double, double> p_posterior() const {
        return {static_cast<double>(_M_pos - _T_pos) + _alpha,
                static_cast<double>(_T_pos) + _beta};
    }
    std::pair<double, double> q_posterior() const {
        return {static_cast<double>(_X_meas - _T_pos) + _mu,
                static_cast<double>(_N_meas - _X_meas - (_M_pos - _T_pos)) + _nu};
    }

    double edge_move_dS(int u, int v, long dm, bool density = false, double aE = 1.0) {
        if (!valid_dyad(u, v)) return INF;
        long w0 = _prior.base().get_w(u, v);
        long w1 = w0 + dm;
        if (w1 < 0) return INF;
        if (dm == 0) return 0.0;

        uint64_t k = key(u, v);
        double dS_meas = 0.0;
        if (_model == MODEL_LATENT_MULTI) {
            bool is_obs = (_obs_q.find(k) != _obs_q.end());
            if (u == v && _self_loops) {
                if (is_obs && w1 == 0) return INF;
            } else {
                if (is_obs && w1 == 0) return INF;
                if (!is_obs && w1 > 0) return INF;
            }
        } else if (_model == MODEL_MEASURED) {
            if ((w0 == 0 && w1 > 0) || (w0 > 0 && w1 == 0)) {
                auto [ni, xi] = get_nx_val(k);
                int sgn = (w1 > 0) ? 1 : -1;
                dS_meas = meas_entropy_from(_M_pos + sgn * ni, _T_pos + sgn * xi) -
                          meas_entropy_from(_M_pos, _T_pos);
            }
        } else if (_model == MODEL_MIXED || _model == MODEL_UNCERTAIN) {
            if ((w0 == 0 && w1 > 0) || (w0 > 0 && w1 == 0)) {
                double qv = get_q_val(k);
                if (w1 > 0) {
                    if (std::isinf(qv)) return (qv < 0) ? INF : 0.0;
                    dS_meas = -qv;
                } else {
                    if (std::isinf(qv)) return (qv > 0) ? INF : 0.0;
                    dS_meas = +qv;
                }
            }
        }
        if (std::isinf(dS_meas)) return dS_meas;

        double dS_dens = 0.0;
        if (density) {
            long E0 = _prior.base().E(), E1 = E0 + dm;
            dS_dens = (static_cast<double>(E1) * std::log(aE) - lgam(static_cast<double>(E1) + 1.0)) -
                      (static_cast<double>(E0) * std::log(aE) - lgam(static_cast<double>(E0) + 1.0));
        }

        double dS_sbm = _prior.edge_delta(u, v, dm, false);
        return dS_meas + dS_dens + dS_sbm;
    }

    double apply_edge_move(int u, int v, long dm, bool density = false, double aE = 1.0) {
        long w0 = _prior.base().get_w(u, v);
        long w1 = w0 + dm;
        uint64_t k = key(u, v);
        double dS_meas = 0.0;
        if (_model == MODEL_MEASURED) {
            if ((w0 == 0 && w1 > 0) || (w0 > 0 && w1 == 0)) {
                auto [ni, xi] = get_nx_val(k);
                int sgn = (w1 > 0) ? 1 : -1;
                double s_old = meas_entropy_from(_M_pos, _T_pos);
                _M_pos += sgn * ni;
                _T_pos += sgn * xi;
                dS_meas = meas_entropy_from(_M_pos, _T_pos) - s_old;
            }
        } else if (_model == MODEL_MIXED || _model == MODEL_UNCERTAIN) {
            if ((w0 == 0 && w1 > 0) || (w0 > 0 && w1 == 0)) {
                double qv = get_q_val(k);
                int sgn = (w1 > 0) ? 1 : -1;
                if (std::isinf(qv)) {
                    if (qv > 0) _n_inf_q += sgn;
                    else _n_ninf_q += sgn;
                } else {
                    _q_sum += static_cast<double>(sgn) * qv;
                    dS_meas = -static_cast<double>(sgn) * qv;
                }
            }
        }
        double dS_dens = 0.0;
        if (density) {
            long E0 = _prior.base().E(), E1 = E0 + dm;
            dS_dens = (static_cast<double>(E1) * std::log(aE) - lgam(static_cast<double>(E1) + 1.0)) -
                      (static_cast<double>(E0) * std::log(aE) - lgam(static_cast<double>(E0) + 1.0));
        }
        double dS_sbm = _prior.edge_delta(u, v, dm, true);
        if (w0 == 0 && w1 > 0) active_insert(k);
        else if (w0 > 0 && w1 == 0) active_erase(k);
        return dS_meas + dS_dens + dS_sbm;
    }

    double get_edge_lprob(int u, int v, double epsilon = 1e-8) {
        if (!valid_dyad(u, v)) return -INF;
        long w0 = _prior.base().get_w(u, v);
        if (w0 > 0) apply_edge_move(u, v, -w0);

        double lse_pos = -INF;
        double lse_all = 0.0;
        double prev_term = 0.0;
        for (long m = 1; m <= 1024; ++m) {
            double dS = edge_move_dS(u, v, m);
            if (std::isinf(dS) && dS > 0) {
                if (m == 1) {
                    if (w0 > 0) apply_edge_move(u, v, w0);
                    return -INF;
                }
                break;
            }
            double lt = -dS;
            double new_pos = log_sum_exp(lse_pos, lt);
            lse_all = log_sum_exp(lse_all, lt);
            bool conv = (m > 2 && lt < prev_term && std::exp(lt - new_pos) < epsilon);
            lse_pos = new_pos;
            prev_term = lt;
            if (conv) break;
        }
        if (w0 > 0) apply_edge_move(u, v, w0);
        return lse_pos - lse_all;
    }

    // Exact detailed-balance edge MCMC sweep:
    // - LATENT_MULTI: visits observed edges (min_w=1) and self-loops (min_w=0 if enabled).
    // - MEASURED/MIXED/UNCERTAIN: mixture of:
    //     Move A (p=0.5): pick fixed observed dyad (u,v) in _obs_keys, symmetric +/-1 reflected at 0.
    //     Move B (p=0.5): reversible add/delete pair:
    //       - with prob 0.5, pick uniform dyad from all M dyads, propose w -> w + 1
    //         (reverse picks uniform active edge from A', propose w' -> w' - 1)
    //       - with prob 0.5, pick uniform active edge from A, propose w -> w - 1
    //         (reverse picks uniform dyad from M, propose w' -> w' + 1)
    void edge_mcmc_sweep(Rng& rng, double beta, int n_steps,
                         double& total_dS, size_t& nattempts, size_t& nmoves) {
        if (_N <= 1 || _M_dyads == 0) return;

        if (_model == MODEL_LATENT_MULTI) {
            int steps = std::max(1, n_steps);
            size_t n_obs = _obs_keys.size();
            size_t n_tot = n_obs + (_self_loops ? static_cast<size_t>(_N) : 0);
            if (n_tot == 0) return;
            for (int st = 0; st < steps; ++st) {
                size_t pick = rng.index(n_tot);
                int u, v;
                long min_w;
                if (pick < n_obs) {
                    auto pr = unpack_key(_obs_keys[pick]);
                    u = pr.first; v = pr.second;
                    min_w = 1L;
                } else {
                    u = v = static_cast<int>(pick - n_obs);
                    min_w = (_obs_q.find(key(u, u)) != _obs_q.end()) ? 1L : 0L;
                }
                ++nattempts;
                long w0 = _prior.base().get_w(u, v);
                long w1;
                double log_ratio;
                if (w0 == min_w) {
                    w1 = min_w + 1;
                    log_ratio = -LOG2;
                } else if (rng.unif() < 0.5) {
                    w1 = w0 - 1;
                    log_ratio = (w1 == min_w) ? LOG2 : 0.0;
                } else {
                    w1 = w0 + 1;
                    log_ratio = 0.0;
                }
                long dm = w1 - w0;
                double dS = edge_move_dS(u, v, dm);
                if (std::isinf(dS) && dS > 0) continue;
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS + log_ratio));
                if (accept) {
                    total_dS += apply_edge_move(u, v, dm);
                    ++nmoves;
                }
            }
            return;
        }

        int steps = std::max(1, n_steps);
        double log_M = std::log(static_cast<double>(_M_dyads));
        for (int st = 0; st < steps; ++st) {
            if (!_obs_keys.empty() && rng.unif() < 0.5) {
                // Move A: observed dyad (state-independent selection)
                auto [u, v] = unpack_key(_obs_keys[rng.index(_obs_keys.size())]);
                if (!valid_dyad(u, v)) continue;
                ++nattempts;
                long w0 = _prior.base().get_w(u, v);
                long w1;
                double log_ratio;
                if (w0 == 0) {
                    w1 = 1;
                    log_ratio = -LOG2;
                } else if (rng.unif() < 0.5) {
                    w1 = w0 - 1;
                    log_ratio = (w1 == 0) ? LOG2 : 0.0;
                } else {
                    w1 = w0 + 1;
                    log_ratio = 0.0;
                }
                long dm = w1 - w0;
                double dS = edge_move_dS(u, v, dm);
                if (std::isinf(dS) && dS > 0) continue;
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS + log_ratio));
                if (accept) {
                    total_dS += apply_edge_move(u, v, dm);
                    ++nmoves;
                }
            } else {
                // Move B: exact reversible add (uniform dyad) / delete (uniform active edge)
                ++nattempts;
                if (rng.unif() < 0.5) {
                    // Add move: w0 -> w0 + 1
                    auto [u, v] = sample_uniform_dyad(rng);
                    long w0 = _prior.base().get_w(u, v);
                    size_t A_after = _active.size() + (w0 == 0 ? 1 : 0);
                    double log_ratio = log_M - std::log(static_cast<double>(A_after));
                    double dS = edge_move_dS(u, v, +1);
                    if (std::isinf(dS) && dS > 0) continue;
                    bool accept = std::isinf(beta)
                        ? (dS < 0.0)
                        : (rng.unif() < std::exp(-beta * dS + log_ratio));
                    if (accept) {
                        total_dS += apply_edge_move(u, v, +1);
                        ++nmoves;
                    }
                } else {
                    // Delete move: pick uniform active edge, w0 -> w0 - 1
                    if (_active.empty()) continue;
                    size_t A_before = _active.size();
                    uint64_t k = _active[rng.index(A_before)];
                    auto [u, v] = unpack_key(k);
                    double log_ratio = std::log(static_cast<double>(A_before)) - log_M;
                    double dS = edge_move_dS(u, v, -1);
                    if (std::isinf(dS) && dS > 0) continue;
                    bool accept = std::isinf(beta)
                        ? (dS < 0.0)
                        : (rng.unif() < std::exp(-beta * dS + log_ratio));
                    if (accept) {
                        total_dS += apply_edge_move(u, v, -1);
                        ++nmoves;
                    }
                }
            }
        }
    }

    void h_mcmc_step(Rng& rng, double beta_temp, double hstep,
                     double& total_dS, size_t& nattempts, size_t& nmoves) {
        if (_model != MODEL_MIXED) return;
        ++nattempts;
        double hs[4] = {_alpha, _beta, _mu, _nu};
        int j = static_cast<int>(rng.index(4));
        double old_h = hs[j];
        double lo_f = std::max(1e-4, old_h - hstep);
        double hi_f = old_h + hstep;
        double pf = 1.0 / (hi_f - lo_f);
        double new_h = lo_f + rng.unif() * (hi_f - lo_f);
        double lo_b = std::max(1e-4, new_h - hstep);
        double hi_b = new_h + hstep;
        double pb = 1.0 / (hi_b - lo_b);

        double Sb = meas_entropy();
        hs[j] = new_h;
        _alpha = hs[0]; _beta = hs[1]; _mu = hs[2]; _nu = hs[3];
        sync_mixed_q();
        double Sa = meas_entropy();
        double dS = Sa - Sb;

        bool accept = false;
        if (std::isinf(beta_temp)) accept = (dS < 0.0);
        else {
            double log_a = -beta_temp * dS + std::log(pb) - std::log(pf);
            accept = (log_a >= 0.0 || rng.unif() < std::exp(log_a));
        }
        if (accept) {
            total_dS += dS;
            ++nmoves;
        } else {
            hs[j] = old_h;
            _alpha = hs[0]; _beta = hs[1]; _mu = hs[2]; _nu = hs[3];
            sync_mixed_q();
        }
    }

    double alpha() const { return _alpha; }
    double beta_hp() const { return _beta; }
    double mu() const { return _mu; }
    double nu() const { return _nu; }

private:
    int _model;
    int _N;
    bool _directed, _self_loops;
    int64_t _M_dyads = 0;
    DynBlockPrior _prior;

    gt_hash_map<uint64_t, std::pair<long, long>> _obs_nx;
    std::vector<uint64_t> _obs_keys;
    long _n_def = 1, _x_def = 0;
    int64_t _N_meas = 0, _X_meas = 0, _M_pos = 0, _T_pos = 0;
    double _lbinom_const = 0.0;
    double _alpha = 1.0, _beta = 1.0, _mu = 1.0, _nu = 1.0;
    double _lp = std::numeric_limits<double>::quiet_NaN();
    double _lq = std::numeric_limits<double>::quiet_NaN();

    gt_hash_map<uint64_t, double> _obs_q;
    double _q_def = -INF;
    double _S_const = 0.0;
    double _q_sum = 0.0;
    int64_t _n_inf_q = 0, _n_ninf_q = 0;

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

    double meas_entropy_from(int64_t M_pos, int64_t T_pos) const {
        double log_p_term;
        if (std::isnan(_lp)) {
            log_p_term = lbeta_d(static_cast<double>(M_pos - T_pos) + _alpha,
                                 static_cast<double>(T_pos) + _beta) -
                         lbeta_d(_alpha, _beta);
        } else {
            log_p_term = static_cast<double>(M_pos - T_pos) * _lp +
                         static_cast<double>(T_pos) * log1mexp_d(_lp);
        }
        double log_q_term;
        int64_t fp_pos = _X_meas - T_pos;
        int64_t fp_neg = _N_meas - _X_meas - (M_pos - T_pos);
        if (std::isnan(_lq)) {
            log_q_term = lbeta_d(static_cast<double>(fp_pos) + _mu,
                                 static_cast<double>(fp_neg) + _nu) -
                         lbeta_d(_mu, _nu);
        } else {
            log_q_term = static_cast<double>(fp_pos) * _lq +
                         static_cast<double>(fp_neg) * log1mexp_d(_lq);
        }
        return -(_lbinom_const + log_p_term + log_q_term);
    }
};

UncertainEngine build_uncertain_engine(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist, double q_max) {
    int L = static_cast<int>(bs_list.size());
    size_t qm = std::isinf(q_max) ? 5000 : static_cast<size_t>(q_max);
    UncertainEngine eng(model, N, directed, self_loops, dc, dist, qm, L);

    if (model == MODEL_MEASURED || model == MODEL_MIXED) {
        std::vector<std::tuple<int, int, long, long>> obs;
        obs.reserve(obs_edges.nrow());
        for (int i = 0; i < obs_edges.nrow(); ++i) {
            obs.emplace_back(obs_edges(i, 0) - 1, obs_edges(i, 1) - 1,
                             static_cast<long>(obs_v1[i]), static_cast<long>(obs_v2[i]));
        }
        eng.setup_measured(obs, static_cast<long>(def1), static_cast<long>(def2),
                           hparams[0], hparams[1], hparams[2], hparams[3], lp, lq);
    } else if (model == MODEL_UNCERTAIN) {
        std::vector<std::tuple<int, int, double>> obs;
        obs.reserve(obs_edges.nrow());
        for (int i = 0; i < obs_edges.nrow(); ++i) {
            obs.emplace_back(obs_edges(i, 0) - 1, obs_edges(i, 1) - 1, obs_v1[i]);
        }
        eng.setup_uncertain(obs, def1);
    } else {
        std::vector<std::pair<int, int>> obs;
        obs.reserve(obs_edges.nrow());
        for (int i = 0; i < obs_edges.nrow(); ++i) {
            obs.emplace_back(obs_edges(i, 0) - 1, obs_edges(i, 1) - 1);
        }
        eng.setup_latent_multi(obs);
    }

    std::vector<std::tuple<int, int, long>> lat;
    lat.reserve(latent_edges.nrow());
    for (int i = 0; i < latent_edges.nrow(); ++i) {
        lat.emplace_back(latent_edges(i, 0) - 1, latent_edges(i, 1) - 1,
                         static_cast<long>(latent_w[i]));
    }
    std::vector<std::vector<int>> bs;
    bs.reserve(L);
    for (int l = 0; l < L; ++l) {
        Rcpp::IntegerVector x = bs_list[l];
        std::vector<int> v(x.size());
        for (R_xlen_t i = 0; i < x.size(); ++i) v[i] = x[i] - 1;
        bs.push_back(std::move(v));
    }
    eng.init_latent(lat, bs);
    return eng;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_uncertain_entropy(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist, bool sbm, bool density, double aE, double q_max) {
    UncertainEngine eng = build_uncertain_engine(
        model, N, directed, self_loops, obs_edges, obs_v1, obs_v2, def1, def2,
        hparams, lp, lq, latent_edges, latent_w, bs_list, dc, dist, q_max);
    double mS = eng.meas_entropy(density, aE);
    double sS = sbm ? eng.prior().entropy() : 0.0;
    auto pp = eng.p_posterior();
    auto qp = eng.q_posterior();
    return Rcpp::List::create(
        Rcpp::Named("total") = mS + sS,
        Rcpp::Named("measurement") = mS,
        Rcpp::Named("sbm") = sS,
        Rcpp::Named("p_posterior") = Rcpp::NumericVector::create(pp.first, pp.second),
        Rcpp::Named("q_posterior") = Rcpp::NumericVector::create(qp.first, qp.second));
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_uncertain_edge_probs(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist,
    Rcpp::IntegerMatrix query_edges, double epsilon, double q_max) {
    UncertainEngine eng = build_uncertain_engine(
        model, N, directed, self_loops, obs_edges, obs_v1, obs_v2, def1, def2,
        hparams, lp, lq, latent_edges, latent_w, bs_list, dc, dist, q_max);
    int m = query_edges.nrow();
    Rcpp::NumericVector lprobs(m);
    for (int i = 0; i < m; ++i) {
        lprobs[i] = eng.get_edge_lprob(query_edges(i, 0) - 1, query_edges(i, 1) - 1, epsilon);
    }
    return lprobs;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_uncertain_edge_dS(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist,
    Rcpp::IntegerMatrix query_edges, Rcpp::IntegerVector dm_vec, double q_max) {
    UncertainEngine eng = build_uncertain_engine(
        model, N, directed, self_loops, obs_edges, obs_v1, obs_v2, def1, def2,
        hparams, lp, lq, latent_edges, latent_w, bs_list, dc, dist, q_max);
    int m = query_edges.nrow();
    Rcpp::NumericVector out(m);
    for (int i = 0; i < m; ++i) {
        out[i] = eng.edge_move_dS(query_edges(i, 0) - 1, query_edges(i, 1) - 1,
                                  static_cast<long>(dm_vec[i]));
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_uncertain_check_moves(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist, int nmoves, double q_max) {
    UncertainEngine eng = build_uncertain_engine(
        model, N, directed, self_loops, obs_edges, obs_v1, obs_v2, def1, def2,
        hparams, lp, lq, latent_edges, latent_w, bs_list, dc, dist, q_max);
    Rng rng;
    Rcpp::NumericVector err(nmoves);
    double S = eng.total_entropy();
    for (int i = 0; i < nmoves; ++i) {
        if (rng.unif() < 0.6) {
            int u = static_cast<int>(rng.index(static_cast<size_t>(N)));
            int v = static_cast<int>(rng.index(static_cast<size_t>(N)));
            if (!self_loops && u == v) v = (v + 1) % N;
            long w0 = eng.prior().base().get_w(u, v);
            long dm = (w0 > 0 && rng.unif() < 0.5) ? -1L : 1L;
            if (model == MODEL_LATENT_MULTI) {
                if (obs_edges.nrow() == 0) continue;
                int ei = static_cast<int>(rng.index(static_cast<size_t>(obs_edges.nrow())));
                u = obs_edges(ei, 0) - 1; v = obs_edges(ei, 1) - 1;
                w0 = eng.prior().base().get_w(u, v);
                dm = (w0 > 1 && rng.unif() < 0.5) ? -1L : 1L;
            }
            double dS = eng.edge_move_dS(u, v, dm);
            if (std::isinf(dS)) continue;
            eng.apply_edge_move(u, v, dm);
            double S_new = eng.total_entropy();
            err[i] = std::abs(dS - (S_new - S));
            S = S_new;
        } else {
            int v = static_cast<int>(rng.index(static_cast<size_t>(N)));
            int s = eng.prior().base().random_occupied(rng);
            int r = eng.prior().base().group(v);
            if (s == r) continue;
            if (eng.prior().base().size(r) == 1 && !eng.prior().same_parent(0, r, s)) continue;
            double dS = eng.prior().vertex_move_dS(v, s);
            eng.prior().vertex_move(v, s);
            double S_new = eng.total_entropy();
            err[i] = std::abs(dS - (S_new - S));
            S = S_new;
        }
    }
    return err;
}

// [[Rcpp::export]]
Rcpp::List cpp_uncertain_reconstruct(
    int model, int N, bool directed, bool self_loops,
    Rcpp::IntegerMatrix obs_edges, Rcpp::NumericVector obs_v1, Rcpp::NumericVector obs_v2,
    double def1, double def2, Rcpp::NumericVector hparams, double lp, double lq,
    Rcpp::IntegerMatrix latent_edges, Rcpp::IntegerVector latent_w,
    Rcpp::List bs_list, bool dc, bool dist,
    int nsamples, int burnin, int thin, double beta,
    double pedges, double ph, double hstep, double c_prop, double d_prop,
    bool collect_partitions, double q_max) {
    UncertainEngine eng = build_uncertain_engine(
        model, N, directed, self_loops, obs_edges, obs_v1, obs_v2, def1, def2,
        hparams, lp, lq, latent_edges, latent_w, bs_list, dc, dist, q_max);
    Rng rng;
    double total_dS = 0.0;
    size_t nattempts = 0, nmoves = 0;

    int edge_steps = std::max(N * 2, static_cast<int>(obs_edges.nrow()) * 2);
    auto do_sweep = [&]() {
        if (model == MODEL_MIXED && rng.unif() < ph) {
            eng.h_mcmc_step(rng, beta, hstep, total_dS, nattempts, nmoves);
        } else if (rng.unif() < pedges) {
            eng.edge_mcmc_sweep(rng, beta, edge_steps, total_dS, nattempts, nmoves);
        } else {
            eng.prior().mcmc_sweep(rng, beta, c_prop, d_prop, total_dS, nattempts, nmoves);
        }
    };

    for (int it = 0; it < burnin; ++it) {
        do_sweep();
        if ((it & 15) == 0) Rcpp::checkUserInterrupt();
    }

    gt_hash_map<uint64_t, long> edge_present_cnt;
    gt_hash_map<uint64_t, long> edge_weight_sum;
    std::vector<double> entropies;
    std::vector<double> p_means, q_means;
    std::vector<Rcpp::List> part_samples;
    entropies.reserve(nsamples);
    if (model == MODEL_MEASURED) { p_means.reserve(nsamples); q_means.reserve(nsamples); }
    if (collect_partitions) part_samples.reserve(nsamples);

    int thin_n = std::max(1, thin);
    for (int s = 0; s < nsamples; ++s) {
        for (int t = 0; t < thin_n; ++t) do_sweep();
        for (const auto& [k, w] : eng.prior().base().edge_weights()) {
            if (w > 0) {
                edge_present_cnt[k] += 1;
                edge_weight_sum[k] += w;
            }
        }
        entropies.push_back(eng.total_entropy());
        if (model == MODEL_MEASURED) {
            auto pp = eng.p_posterior();
            auto qp = eng.q_posterior();
            p_means.push_back(pp.first / (pp.first + pp.second));
            q_means.push_back(qp.first / (qp.first + qp.second));
        }
        if (collect_partitions) part_samples.push_back(eng.prior().export_bs());
        if ((s & 15) == 0) Rcpp::checkUserInterrupt();
    }

    size_t nm = edge_present_cnt.size();
    Rcpp::IntegerVector m_u(nm), m_v(nm);
    Rcpp::NumericVector m_prob(nm), m_mw(nm);
    size_t idx = 0;
    for (const auto& [k, cnt] : edge_present_cnt) {
        auto [u, v] = unpack_key(k);
        m_u[idx] = u + 1;
        m_v[idx] = v + 1;
        m_prob[idx] = static_cast<double>(cnt) / static_cast<double>(nsamples);
        m_mw[idx] = static_cast<double>(edge_weight_sum[k]) / static_cast<double>(nsamples);
        ++idx;
    }

    const auto& fw = eng.prior().base().edge_weights();
    Rcpp::IntegerMatrix fin_edges(fw.size(), 2);
    Rcpp::IntegerVector fin_w(fw.size());
    size_t fi = 0;
    for (const auto& [k, w] : fw) {
        auto [u, v] = unpack_key(k);
        fin_edges(fi, 0) = u + 1;
        fin_edges(fi, 1) = v + 1;
        fin_w[fi] = static_cast<int>(w);
        ++fi;
    }

    auto pp = eng.p_posterior();
    auto qp = eng.q_posterior();
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("from") = m_u,
        Rcpp::Named("to") = m_v,
        Rcpp::Named("prob") = m_prob,
        Rcpp::Named("mean_weight") = m_mw,
        Rcpp::Named("final_edges") = fin_edges,
        Rcpp::Named("final_weights") = fin_w,
        Rcpp::Named("final_bs") = eng.prior().export_bs(),
        Rcpp::Named("final_hparams") = Rcpp::NumericVector::create(
            eng.alpha(), eng.beta_hp(), eng.mu(), eng.nu()),
        Rcpp::Named("p_posterior") = Rcpp::NumericVector::create(pp.first, pp.second),
        Rcpp::Named("q_posterior") = Rcpp::NumericVector::create(qp.first, qp.second),
        Rcpp::Named("p_samples") = p_means,
        Rcpp::Named("q_samples") = q_means,
        Rcpp::Named("entropies") = entropies,
        Rcpp::Named("entropy") = eng.total_entropy(),
        Rcpp::Named("dS") = total_dS,
        Rcpp::Named("nattempts") = static_cast<double>(nattempts),
        Rcpp::Named("nmoves") = static_cast<double>(nmoves));
    if (collect_partitions) {
        Rcpp::List pl(part_samples.size());
        for (size_t i = 0; i < part_samples.size(); ++i) pl[i] = part_samples[i];
        out["partitions"] = pl;
    }
    return out;
}

// Fast O(V + E) Expectation-Maximization algorithm for latent Poisson multigraph
// inference from an erased simple graph (Peixoto 2020 PRE).
// [[Rcpp::export]]
Rcpp::List cpp_latent_multigraph_em(SEXP pointer, double epsilon, int max_niter) {
    Graph& graph = get_graph(pointer);
    const Adjacency& g = graph.adjacency;
    int N = static_cast<int>(boost::num_vertices(g));
    bool directed = graph.directed;

    gt_hash_map<uint64_t, int> seen;
    std::vector<std::pair<int, int>> edges;
    std::vector<double> k_out(N, 0.0), k_in(N, 0.0);
    auto er = boost::edges(g);
    for (auto it = er.first; it != er.second; ++it) {
        int u = static_cast<int>(boost::source(*it, g));
        int v = static_cast<int>(boost::target(*it, g));
        if (u == v) continue;
        if (!directed && u > v) std::swap(u, v);
        uint64_t k = pair_key(u, v);
        if (seen.emplace(k, 1).second) {
            edges.emplace_back(u, v);
            k_out[u] += 1.0;
            k_in[v] += 1.0;
            if (!directed) { k_out[v] += 1.0; k_in[u] += 1.0; }
        }
    }

    double sum_ko = 0.0, sum_ki = 0.0;
    for (int i = 0; i < N; ++i) { sum_ko += k_out[i]; sum_ki += k_in[i]; }
    std::vector<double> th_o(N, 0.0), th_i(N, 0.0);
    if (sum_ko > 0.0) {
        double inv_so = 1.0 / std::sqrt(sum_ko);
        double inv_si = (directed && sum_ki > 0.0) ? (1.0 / std::sqrt(sum_ki)) : inv_so;
        for (int i = 0; i < N; ++i) {
            th_o[i] = k_out[i] * inv_so;
            th_i[i] = directed ? (k_in[i] * inv_si) : th_o[i];
        }
    }

    int iter = 0;
    int limit = (max_niter > 0) ? max_niter : 100000;
    std::vector<double> d_o(N, 0.0), d_i(N, 0.0);
    double max_diff = 0.0;

    if (sum_ko > 0.0) {
        for (iter = 1; iter <= limit; ++iter) {
            std::fill(d_o.begin(), d_o.end(), 0.0);
            if (directed) std::fill(d_i.begin(), d_i.end(), 0.0);

            for (const auto& [u, v] : edges) {
                double lam = th_o[u] * th_i[v];
                double w = (lam > 1e-12) ? (lam / (-std::expm1(-lam))) : 1.0;
                d_o[u] += w;
                if (directed) d_i[v] += w;
                else d_o[v] += w;
            }
            for (int u = 0; u < N; ++u) {
                double lam = th_o[u] * th_i[u];
                d_o[u] += lam;
                if (directed) d_i[u] += lam;
            }
            double tot_o = 0.0, tot_i = 0.0;
            for (int u = 0; u < N; ++u) {
                tot_o += d_o[u];
                if (directed) tot_i += d_i[u];
            }
            double inv_o = (tot_o > 0.0) ? (1.0 / std::sqrt(tot_o)) : 0.0;
            double inv_i = directed ? ((tot_i > 0.0) ? (1.0 / std::sqrt(tot_i)) : 0.0) : inv_o;
            max_diff = 0.0;
            for (int u = 0; u < N; ++u) {
                double no = d_o[u] * inv_o;
                double ni = directed ? (d_i[u] * inv_i) : no;
                max_diff = std::max(max_diff, std::abs(no - th_o[u]));
                if (directed) max_diff = std::max(max_diff, std::abs(ni - th_i[u]));
                th_o[u] = no;
                th_i[u] = ni;
            }
            if (max_diff < epsilon) break;
            if ((iter & 255) == 0) Rcpp::checkUserInterrupt();
        }
    }

    size_t total_out = edges.size() + static_cast<size_t>(N);
    Rcpp::IntegerMatrix out_edges(total_out, 2);
    Rcpp::NumericVector out_w(total_out);
    size_t pos = 0;
    for (const auto& [u, v] : edges) {
        double lam = th_o[u] * th_i[v];
        out_edges(pos, 0) = u + 1;
        out_edges(pos, 1) = v + 1;
        out_w[pos] = (lam > 1e-12) ? (lam / (-std::expm1(-lam))) : 1.0;
        ++pos;
    }
    for (int u = 0; u < N; ++u) {
        out_edges(pos, 0) = u + 1;
        out_edges(pos, 1) = u + 1;
        out_w[pos] = th_o[u] * th_i[u] * (directed ? 1.0 : 0.5);
        ++pos;
    }

    return Rcpp::List::create(
        Rcpp::Named("edges") = out_edges,
        Rcpp::Named("weights") = out_w,
        Rcpp::Named("theta_out") = th_o,
        Rcpp::Named("theta_in") = th_i,
        Rcpp::Named("niter") = iter,
        Rcpp::Named("delta") = max_diff);
}
