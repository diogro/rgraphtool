// rgraphtool: native C++ LatentClosureBlockState (Peixoto 2022 PRX)
// Disentangling homophily / community structure (flat or nested SBM substrate)
// and L >= 1 generations of zero-inflated triadic closures.
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

namespace {

using namespace rgt_dyn;

struct ClosureDelta {
    int l;
    int k;
    int64_t dM;
    int64_t dE;
};

class ClosureEngine {
public:
    ClosureEngine(int N, int L, bool dc, bool dist, size_t q_max, int sbm_levels)
        : _N(N), _L(L), _prior(N, /*directed=*/false, dc, dist, q_max, sbm_levels),
          _obs_adj(N),
          _M_cnt(L + 1, std::vector<int64_t>(N, 0)),
          _E_cnt(L + 1, std::vector<int64_t>(N, 0)),
          _NM(L + 1, 0), _NE(L + 1, 0) {}

    int N() const { return _N; }
    int L() const { return _L; }
    int num_edges() const { return static_cast<int>(_u.size()); }
    DynBlockPrior& prior() { return _prior; }
    const DynBlockPrior& prior() const { return _prior; }
    const std::vector<int>& gen() const { return _gen; }
    const std::vector<std::vector<int>>& ex() const { return _ex; }
    const std::vector<int64_t>& M_layer(int l) const { return _M_cnt[l]; }
    const std::vector<int64_t>& E_layer(int l) const { return _E_cnt[l]; }

    static uint64_t ukey(int a, int b) {
        if (a > b) std::swap(a, b);
        return pair_key(a, b);
    }

    void init(const std::vector<std::pair<int, int>>& obs_edges,
              const std::vector<int>& gen_init,
              const std::vector<long>& w0_init,
              const std::vector<long>& loops_init,
              const std::vector<std::vector<int>>& ex_init,
              const std::vector<std::vector<int>>& bs) {
        int m = static_cast<int>(obs_edges.size());
        _u.resize(m);
        _v.resize(m);
        _gen.resize(m);
        _ex.resize(m);
        _edge_id.clear();
        for (int i = 0; i < _N; ++i) _obs_adj[i].clear();

        for (int e = 0; e < m; ++e) {
            int a = obs_edges[e].first, b = obs_edges[e].second;
            if (a == b) Rcpp::stop("LatentClosureBlockState expects simple undirected edges without self-loops.");
            if (a > b) std::swap(a, b);
            uint64_t k = pair_key(a, b);
            if (_edge_id.find(k) != _edge_id.end())
                Rcpp::stop("Duplicate observed edge in LatentClosureBlockState.");
            _edge_id[k] = e;
            _u[e] = a;
            _v[e] = b;
            _obs_adj[a].emplace_back(b, e);
            _obs_adj[b].emplace_back(a, e);
            _gen[e] = gen_init[e];
            _ex[e] = ex_init[e];
            std::sort(_ex[e].begin(), _ex[e].end());
            _ex[e].erase(std::unique(_ex[e].begin(), _ex[e].end()), _ex[e].end());
        }
        for (int i = 0; i < _N; ++i) {
            std::sort(_obs_adj[i].begin(), _obs_adj[i].end());
        }

        // Precompute triangles (common neighbors) for every observed edge
        _common.assign(m, {});
        for (int e = 0; e < m; ++e) {
            const auto& au = _obs_adj[_u[e]];
            const auto& av = _obs_adj[_v[e]];
            size_t i = 0, j = 0;
            while (i < au.size() && j < av.size()) {
                if (au[i].first == av[j].first) {
                    _common[e].emplace_back(au[i].first, au[i].second, av[j].second);
                    ++i; ++j;
                } else if (au[i].first < av[j].first) {
                    ++i;
                } else {
                    ++j;
                }
            }
        }

        // Initialize substrate SBM with layer-0 edges + self-loops
        std::vector<std::tuple<int, int, long>> sub_edges;
        sub_edges.reserve(m + _N);
        for (int e = 0; e < m; ++e) {
            if (_gen[e] == 0) {
                long w = std::max(1L, w0_init[e]);
                sub_edges.emplace_back(_u[e], _v[e], w);
            }
        }
        for (int i = 0; i < _N; ++i) {
            if (i < static_cast<int>(loops_init.size()) && loops_init[i] > 0) {
                sub_edges.emplace_back(i, i, loops_init[i]);
            }
        }
        _prior.init(sub_edges, bs);

        // Compute M_cnt, E_cnt and validate initial state
        recompute_counts_and_validate();
    }

    double layer_entropy(int l) const {
        if (l < 1 || l > _L) return 0.0;
        const MathTab& MT = mtab();
        int64_t nm = _NM[l], ne = _NE[l];
        double S = MT.logn(nm + 1) + MT.lbinom(nm, ne);
        for (int k = 0; k < _N; ++k) {
            int64_t ek = _E_cnt[l][k];
            if (ek > 0) {
                int64_t mk = _M_cnt[l][k];
                S += MT.logn(mk) + MT.lbinom(mk, ek);
            }
        }
        return S;
    }

    double closure_entropy() const {
        double S = 0.0;
        for (int l = 1; l <= _L; ++l) S += layer_entropy(l);
        return S;
    }

    double total_entropy() const {
        return _prior.entropy() + closure_entropy();
    }

    // Candidate egos C_l(u, v) for observed edge e at generation l >= 1
    void candidate_egos(int e, int l, std::vector<int>& out) const {
        out.clear();
        if (l < 1 || l > _L) return;
        for (const auto& [k, e_uk, e_vk] : _common[e]) {
            if (std::max(_gen[e_uk], _gen[e_vk]) == l - 1) out.push_back(k);
        }
    }

    // Check if assigning gen[e] = l1 breaks any existing ego using edge e as a leg
    bool layer_valid_for_legs(int e, int l1) const {
        int u = _u[e], v = _v[e];
        for (const auto& [w, e_uw, e_vw] : _common[e]) {
            // Check if u is an ego closing (v, w)
            if (_gen[e_vw] >= 1 && has_ego(e_vw, u)) {
                if (std::max(l1, _gen[e_uw]) != _gen[e_vw] - 1) return false;
            }
            // Check if v is an ego closing (u, w)
            if (_gen[e_uw] >= 1 && has_ego(e_uw, v)) {
                if (std::max(l1, _gen[e_vw]) != _gen[e_uw] - 1) return false;
            }
        }
        return true;
    }

    // Collect exact sparse (l, k, dM, dE) deltas for changing edge e from
    // (_gen[e], _ex[e]) to (l1, new_ex).
    void collect_closure_deltas(int e, int l1, const std::vector<int>& new_ex,
                                std::vector<ClosureDelta>& deltas) const {
        deltas.clear();
        int l0 = _gen[e];
        auto add = [&](int l, int k, int64_t dM, int64_t dE) {
            if (l < 1 || l > _L || (dM == 0 && dE == 0)) return;
            for (auto& d : deltas) {
                if (d.l == l && d.k == k) { d.dM += dM; d.dE += dE; return; }
            }
            deltas.push_back({l, k, dM, dE});
        };

        if (l0 >= 1) {
            for (int k : _ex[e]) add(l0, k, 0, -1);
        }
        if (l1 >= 1) {
            for (int k : new_ex) add(l1, k, 0, +1);
        }

        if (l0 != l1) {
            int u = _u[e], v = _v[e];
            // Ego u: pairs {v, w} for w in N_obs(u) \ {v}
            for (const auto& [w, e_uw] : _obs_adj[u]) {
                if (w == v) continue;
                int g_uw = _gen[e_uw];
                auto it = _edge_id.find(ukey(v, w));
                int g_vw = (it == _edge_id.end()) ? (_L + 1) : _gen[it->second];
                int lo = std::max(l0, g_uw) + 1;
                int ln = std::max(l1, g_uw) + 1;
                if (lo <= _L && g_vw >= lo) add(lo, u, -1, 0);
                if (ln <= _L && g_vw >= ln) add(ln, u, +1, 0);
            }
            // Ego v: pairs {u, w} for w in N_obs(v) \ {u}
            for (const auto& [w, e_vw] : _obs_adj[v]) {
                if (w == u) continue;
                int g_vw = _gen[e_vw];
                auto it = _edge_id.find(ukey(u, w));
                int g_uw = (it == _edge_id.end()) ? (_L + 1) : _gen[it->second];
                int lo = std::max(l0, g_vw) + 1;
                int ln = std::max(l1, g_vw) + 1;
                if (lo <= _L && g_uw >= lo) add(lo, v, -1, 0);
                if (ln <= _L && g_uw >= ln) add(ln, v, +1, 0);
            }
            // Common neighbors k: pair {u, v}
            for (const auto& [k, e_uk, e_vk] : _common[e]) {
                int ls = std::max(_gen[e_uk], _gen[e_vk]) + 1;
                if (ls <= _L) {
                    int old_open = (l0 >= ls) ? 1 : 0;
                    int new_open = (l1 >= ls) ? 1 : 0;
                    if (new_open != old_open) add(ls, k, new_open - old_open, 0);
                }
            }
        }
    }

    double eval_closure_deltas(const std::vector<ClosureDelta>& deltas) const {
        const MathTab& MT = mtab();
        std::vector<int64_t> dNM(_L + 1, 0), dNE(_L + 1, 0);
        std::vector<bool> touched(_L + 1, false);
        double dS = 0.0;
        for (const auto& d : deltas) {
            if (d.dM == 0 && d.dE == 0) continue;
            int l = d.l, k = d.k;
            int64_t m0 = _M_cnt[l][k], m1 = m0 + d.dM;
            int64_t e0 = _E_cnt[l][k], e1 = e0 + d.dE;
            if (m1 < 0 || e1 < 0 || e1 > m1) return INF;
            touched[l] = true;
            if (m0 == 0 && m1 > 0) dNM[l]++;
            else if (m0 > 0 && m1 == 0) dNM[l]--;
            if (e0 == 0 && e1 > 0) dNE[l]++;
            else if (e0 > 0 && e1 == 0) dNE[l]--;
            double t0 = (e0 > 0) ? (MT.logn(m0) + MT.lbinom(m0, e0)) : 0.0;
            double t1 = (e1 > 0) ? (MT.logn(m1) + MT.lbinom(m1, e1)) : 0.0;
            dS += (t1 - t0);
        }
        for (int l = 1; l <= _L; ++l) {
            if (!touched[l]) continue;
            if (dNM[l] != 0 || dNE[l] != 0) {
                int64_t nm0 = _NM[l], ne0 = _NE[l];
                int64_t nm1 = nm0 + dNM[l], ne1 = ne0 + dNE[l];
                dS += (MT.logn(nm1 + 1) + MT.lbinom(nm1, ne1)) -
                      (MT.logn(nm0 + 1) + MT.lbinom(nm0, ne0));
            }
        }
        return dS;
    }

    void apply_closure_deltas(const std::vector<ClosureDelta>& deltas) {
        for (const auto& d : deltas) {
            if (d.dM == 0 && d.dE == 0) continue;
            int l = d.l, k = d.k;
            int64_t m0 = _M_cnt[l][k], m1 = m0 + d.dM;
            int64_t e0 = _E_cnt[l][k], e1 = e0 + d.dE;
            if (m0 == 0 && m1 > 0) _NM[l]++;
            else if (m0 > 0 && m1 == 0) _NM[l]--;
            if (e0 == 0 && e1 > 0) _NE[l]++;
            else if (e0 > 0 && e1 == 0) _NE[l]--;
            _M_cnt[l][k] = m1;
            _E_cnt[l][k] = e1;
        }
    }

    // Evaluate (and optionally apply) moving edge e to (l1, w1_sub, new_ex)
    double edge_state_move(int e, int l1, long w1_sub, const std::vector<int>& new_ex,
                           bool apply, std::vector<ClosureDelta>& scratch) {
        int l0 = _gen[e];
        int u = _u[e], v = _v[e];
        long w0 = (l0 == 0) ? _prior.base().get_w(u, v) : 0L;
        long w1 = (l1 == 0) ? w1_sub : 0L;
        long dm = w1 - w0;

        collect_closure_deltas(e, l1, new_ex, scratch);
        double dS_cl = eval_closure_deltas(scratch);
        if (std::isinf(dS_cl)) return INF;
        double dS_sbm = (dm != 0) ? _prior.edge_delta(u, v, dm, apply) : 0.0;
        if (apply) {
            apply_closure_deltas(scratch);
            _gen[e] = l1;
            _ex[e] = new_ex;
        }
        return dS_cl + dS_sbm;
    }

    // Exact detailed-balance MCMC sweep over edges (layer switches, ego toggles,
    // substrate edge multiplicities) and substrate self-loops.
    void edge_mcmc_sweep(Rng& rng, double beta, bool sample_loops,
                         double& total_dS, size_t& nattempts, size_t& nmoves) {
        int m = static_cast<int>(_u.size());
        if (m == 0) return;

        std::vector<int> eord(m);
        std::iota(eord.begin(), eord.end(), 0);
        rng.shuffle(eord);

        std::vector<int> valid_layers, cands0, cands1, new_ex;
        std::vector<ClosureDelta> deltas;

        for (int e : eord) {
            ++nattempts;
            int l0 = _gen[e];
            int u = _u[e], v = _v[e];
            double move_pick = rng.unif();

            if (move_pick < 0.5) {
                // Sub-move A: Layer switch (only between (0, w=1) and (l>=1, |ex|=1))
                if (l0 == 0 && _prior.base().get_w(u, v) != 1) continue;
                if (l0 >= 1 && _ex[e].size() != 1) continue;

                valid_layers.clear();
                for (int l = 0; l <= _L; ++l) {
                    if (!layer_valid_for_legs(e, l)) continue;
                    if (l >= 1) {
                        candidate_egos(e, l, cands1);
                        if (cands1.empty()) continue;
                    }
                    valid_layers.push_back(l);
                }
                if (valid_layers.size() <= 1) continue;

                // Sample l1 uniformly from valid_layers \ {l0}
                size_t pick = rng.index(valid_layers.size() - 1);
                int l1 = valid_layers[pick];
                if (l1 == l0) l1 = valid_layers.back();

                size_t n_c0 = 1, n_c1 = 1;
                if (l0 >= 1) {
                    candidate_egos(e, l0, cands0);
                    n_c0 = cands0.size();
                }
                new_ex.clear();
                if (l1 >= 1) {
                    candidate_egos(e, l1, cands1);
                    n_c1 = cands1.size();
                    new_ex.push_back(cands1[rng.index(n_c1)]);
                }
                // P_fwd = (1 / (|L_valid| - 1)) * (1 / n_c1)
                // P_rev = (1 / (|L_valid| - 1)) * (1 / n_c0)
                // => log(P_rev / P_fwd) = log(n_c1) - log(n_c0)
                double log_ratio = std::log(static_cast<double>(n_c1)) -
                                   std::log(static_cast<double>(n_c0));
                long w1 = (l1 == 0) ? 1L : 0L;
                double dS = edge_state_move(e, l1, w1, new_ex, false, deltas);
                if (std::isinf(dS) && dS > 0) continue;
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS + log_ratio));
                if (accept) {
                    total_dS += edge_state_move(e, l1, w1, new_ex, true, deltas);
                    ++nmoves;
                }
            } else if (l0 == 0) {
                // Sub-move B1: Substrate multiplicity +/-1 (reflected at w = 1)
                long w0 = _prior.base().get_w(u, v);
                long w1;
                double log_ratio;
                if (w0 == 1) {
                    w1 = 2;
                    log_ratio = -LOG2;
                } else if (rng.unif() < 0.5) {
                    w1 = w0 - 1;
                    log_ratio = (w1 == 1) ? LOG2 : 0.0;
                } else {
                    w1 = w0 + 1;
                    log_ratio = 0.0;
                }
                double dS = _prior.edge_delta(u, v, w1 - w0, false);
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS + log_ratio));
                if (accept) {
                    total_dS += _prior.edge_delta(u, v, w1 - w0, true);
                    ++nmoves;
                }
            } else {
                // Sub-move B2: Toggle one ego k in C_{l0}(u, v) (symmetric!)
                candidate_egos(e, l0, cands0);
                if (cands0.size() <= 1) continue;
                int k = cands0[rng.index(cands0.size())];
                new_ex = _ex[e];
                auto it = std::lower_bound(new_ex.begin(), new_ex.end(), k);
                if (it != new_ex.end() && *it == k) {
                    if (new_ex.size() <= 1) continue;  // ego set must stay non-empty
                    new_ex.erase(it);
                } else {
                    new_ex.insert(it, k);
                }
                double dS = edge_state_move(e, l0, 0, new_ex, false, deltas);
                if (std::isinf(dS) && dS > 0) continue;
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS));
                if (accept) {
                    total_dS += edge_state_move(e, l0, 0, new_ex, true, deltas);
                    ++nmoves;
                }
            }
        }

        // Substrate self-loop moves (symmetric +/-1 reflected at 0)
        if (sample_loops) {
            for (int u = 0; u < _N; ++u) {
                ++nattempts;
                long w0 = _prior.base().get_w(u, u);
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
                double dS = _prior.edge_delta(u, u, w1 - w0, false);
                bool accept = std::isinf(beta)
                    ? (dS < 0.0)
                    : (rng.unif() < std::exp(-beta * dS + log_ratio));
                if (accept) {
                    total_dS += _prior.edge_delta(u, u, w1 - w0, true);
                    ++nmoves;
                }
            }
        }
    }

private:
    int _N, _L;
    DynBlockPrior _prior;
    std::vector<int> _u, _v, _gen;
    std::vector<std::vector<int>> _ex;
    gt_hash_map<uint64_t, int> _edge_id;
    std::vector<std::vector<std::pair<int, int>>> _obs_adj;
    std::vector<std::vector<std::tuple<int, int, int>>> _common;

    std::vector<std::vector<int64_t>> _M_cnt, _E_cnt;
    std::vector<int64_t> _NM, _NE;

    bool has_ego(int e, int k) const {
        const auto& v = _ex[e];
        return std::binary_search(v.begin(), v.end(), k);
    }

    void recompute_counts_and_validate() {
        for (int l = 0; l <= _L; ++l) {
            std::fill(_M_cnt[l].begin(), _M_cnt[l].end(), 0);
            std::fill(_E_cnt[l].begin(), _E_cnt[l].end(), 0);
            _NM[l] = 0;
            _NE[l] = 0;
        }
        int m = static_cast<int>(_u.size());
        // Validate every closure edge's ego set and populate E_cnt
        for (int e = 0; e < m; ++e) {
            int l = _gen[e];
            if (l < 0 || l > _L) Rcpp::stop("Edge generation out of range [0, L].");
            if (l == 0) {
                if (!_ex[e].empty()) Rcpp::stop("Layer-0 substrate edge cannot have closure egos.");
                continue;
            }
            if (_ex[e].empty()) Rcpp::stop("Closure edge at generation l >= 1 must have at least one ego.");
            std::vector<int> cands;
            candidate_egos(e, l, cands);
            for (int k : _ex[e]) {
                if (!std::binary_search(cands.begin(), cands.end(), k))
                    Rcpp::stop("Invalid initial state: ego does not form a valid open triad at generation l.");
                _E_cnt[l][k]++;
            }
        }
        // Populate M_cnt
        for (int k = 0; k < _N; ++k) {
            const auto& adj = _obs_adj[k];
            size_t deg = adj.size();
            for (size_t i = 0; i < deg; ++i) {
                int a = adj[i].first, g_ka = _gen[adj[i].second];
                for (size_t j = i + 1; j < deg; ++j) {
                    int b = adj[j].first, g_kb = _gen[adj[j].second];
                    int ls = std::max(g_ka, g_kb) + 1;
                    if (ls > _L) continue;
                    auto it = _edge_id.find(ukey(a, b));
                    int g_ab = (it == _edge_id.end()) ? (_L + 1) : _gen[it->second];
                    if (g_ab >= ls) _M_cnt[ls][k]++;
                }
            }
        }
        int64_t max_M = 0;
        for (int l = 1; l <= _L; ++l) {
            for (int k = 0; k < _N; ++k) {
                if (_M_cnt[l][k] > max_M) max_M = _M_cnt[l][k];
                if (_M_cnt[l][k] > 0) _NM[l]++;
                if (_E_cnt[l][k] > 0) _NE[l]++;
                if (_E_cnt[l][k] > _M_cnt[l][k])
                    Rcpp::stop("Invalid initial state: E_k > M_k.");
            }
        }
        mtab().ensure(static_cast<size_t>(std::max<int64_t>(max_M + _N + 64, 4096)));
    }
};

ClosureEngine build_closure_engine(
    int N, int L, Rcpp::IntegerMatrix obs_edges,
    Rcpp::IntegerVector gen_init, Rcpp::IntegerVector w0_init,
    Rcpp::IntegerVector loops_init, Rcpp::List ex_init,
    Rcpp::List bs_list, bool dc, bool dist, double q_max) {
    int sbm_L = static_cast<int>(bs_list.size());
    size_t qm = std::isinf(q_max) ? 5000 : static_cast<size_t>(q_max);
    ClosureEngine eng(N, L, dc, dist, qm, sbm_L);

    int m = obs_edges.nrow();
    std::vector<std::pair<int, int>> edges(m);
    std::vector<int> gen(m);
    std::vector<long> w0(m);
    std::vector<std::vector<int>> ex(m);
    for (int i = 0; i < m; ++i) {
        edges[i] = {obs_edges(i, 0) - 1, obs_edges(i, 1) - 1};
        gen[i] = gen_init[i];
        w0[i] = static_cast<long>(w0_init[i]);
        Rcpp::IntegerVector ev = ex_init[i];
        ex[i].resize(ev.size());
        for (R_xlen_t j = 0; j < ev.size(); ++j) ex[i][j] = ev[j] - 1;
    }
    std::vector<long> loops(loops_init.size());
    for (R_xlen_t i = 0; i < loops_init.size(); ++i) loops[i] = static_cast<long>(loops_init[i]);

    std::vector<std::vector<int>> bs;
    bs.reserve(sbm_L);
    for (int l = 0; l < sbm_L; ++l) {
        Rcpp::IntegerVector x = bs_list[l];
        std::vector<int> v(x.size());
        for (R_xlen_t i = 0; i < x.size(); ++i) v[i] = x[i] - 1;
        bs.push_back(std::move(v));
    }

    eng.init(edges, gen, w0, loops, ex, bs);
    return eng;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_closure_entropy(
    int N, int L, Rcpp::IntegerMatrix obs_edges,
    Rcpp::IntegerVector gen_init, Rcpp::IntegerVector w0_init,
    Rcpp::IntegerVector loops_init, Rcpp::List ex_init,
    Rcpp::List bs_list, bool dc, bool dist, double q_max) {
    ClosureEngine eng = build_closure_engine(
        N, L, obs_edges, gen_init, w0_init, loops_init, ex_init, bs_list, dc, dist, q_max);
    double sS = eng.prior().entropy();
    Rcpp::NumericVector cl_S(L);
    Rcpp::List M_list(L), E_list(L);
    double tot_cl = 0.0;
    for (int l = 1; l <= L; ++l) {
        cl_S[l - 1] = eng.layer_entropy(l);
        tot_cl += cl_S[l - 1];
        const auto& Mv = eng.M_layer(l);
        const auto& Ev = eng.E_layer(l);
        Rcpp::NumericVector rm(N), re(N);
        for (int k = 0; k < N; ++k) { rm[k] = static_cast<double>(Mv[k]); re[k] = static_cast<double>(Ev[k]); }
        M_list[l - 1] = rm;
        E_list[l - 1] = re;
    }
    return Rcpp::List::create(
        Rcpp::Named("total") = sS + tot_cl,
        Rcpp::Named("sbm") = sS,
        Rcpp::Named("closure") = tot_cl,
        Rcpp::Named("layer_entropies") = cl_S,
        Rcpp::Named("M") = M_list,
        Rcpp::Named("E") = E_list);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_closure_check_moves(
    int N, int L, Rcpp::IntegerMatrix obs_edges,
    Rcpp::IntegerVector gen_init, Rcpp::IntegerVector w0_init,
    Rcpp::IntegerVector loops_init, Rcpp::List ex_init,
    Rcpp::List bs_list, bool dc, bool dist, int nmoves, double q_max) {
    ClosureEngine eng = build_closure_engine(
        N, L, obs_edges, gen_init, w0_init, loops_init, ex_init, bs_list, dc, dist, q_max);
    Rng rng;
    Rcpp::NumericVector errs(nmoves);
    double S = eng.total_entropy();
    int m = eng.num_edges();
    std::vector<int> cands, new_ex;
    std::vector<ClosureDelta> deltas;

    for (int it = 0; it < nmoves; ++it) {
        if (m == 0) break;
        if (rng.unif() < 0.7) {
            int e = static_cast<int>(rng.index(static_cast<size_t>(m)));
            int l0 = eng.gen()[e];
            std::vector<int> valid_l;
            for (int l = 0; l <= L; ++l) {
                if (!eng.layer_valid_for_legs(e, l)) continue;
                if (l >= 1) {
                    eng.candidate_egos(e, l, cands);
                    if (cands.empty()) continue;
                }
                valid_l.push_back(l);
            }
            if (valid_l.empty()) continue;
            int l1 = valid_l[rng.index(valid_l.size())];
            new_ex.clear();
            long w1 = 0;
            if (l1 == 0) {
                w1 = (l0 == 0 && rng.unif() < 0.5) ? 2L : 1L;
            } else {
                eng.candidate_egos(e, l1, cands);
                // Pick random non-empty subset (1 or 2 egos)
                int k1 = cands[rng.index(cands.size())];
                new_ex.push_back(k1);
                if (cands.size() > 1 && rng.unif() < 0.4) {
                    int k2 = cands[rng.index(cands.size())];
                    if (k2 != k1) {
                        new_ex.push_back(k2);
                        std::sort(new_ex.begin(), new_ex.end());
                    }
                }
            }
            double dS = eng.edge_state_move(e, l1, w1, new_ex, false, deltas);
            if (std::isinf(dS)) continue;
            eng.edge_state_move(e, l1, w1, new_ex, true, deltas);
            double S_new = eng.total_entropy();
            errs[it] = std::abs(dS - (S_new - S));
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
            errs[it] = std::abs(dS - (S_new - S));
            S = S_new;
        }
    }
    return errs;
}

// [[Rcpp::export]]
Rcpp::List cpp_closure_infer(
    int N, int L, Rcpp::IntegerMatrix obs_edges,
    Rcpp::IntegerVector gen_init, Rcpp::IntegerVector w0_init,
    Rcpp::IntegerVector loops_init, Rcpp::List ex_init,
    Rcpp::List bs_list, bool dc, bool dist,
    int nsamples, int burnin, int thin, double beta,
    double pedges, double c_prop, double d_prop, bool sample_loops,
    bool collect_partitions, double q_max) {
    ClosureEngine eng = build_closure_engine(
        N, L, obs_edges, gen_init, w0_init, loops_init, ex_init, bs_list, dc, dist, q_max);
    Rng rng;
    double total_dS = 0.0;
    size_t nattempts = 0, nmoves = 0;

    auto do_sweep = [&]() {
        if (rng.unif() < pedges) {
            eng.edge_mcmc_sweep(rng, beta, sample_loops, total_dS, nattempts, nmoves);
        } else {
            eng.prior().mcmc_sweep(rng, beta, c_prop, d_prop, total_dS, nattempts, nmoves);
        }
    };

    for (int it = 0; it < burnin; ++it) {
        do_sweep();
        if ((it & 15) == 0) Rcpp::checkUserInterrupt();
    }

    int m = eng.num_edges();
    std::vector<std::vector<long>> gen_cnt(L + 1, std::vector<long>(m, 0));
    std::vector<std::vector<double>> t_sum(L, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> m_sum(L, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> c_sum(L, std::vector<double>(N, 0.0));
    std::vector<double> entropies;
    std::vector<Rcpp::List> part_samples;
    entropies.reserve(nsamples);
    if (collect_partitions) part_samples.reserve(nsamples);

    int thin_n = std::max(1, thin);
    for (int s = 0; s < nsamples; ++s) {
        for (int t = 0; t < thin_n; ++t) do_sweep();
        const auto& gvec = eng.gen();
        for (int e = 0; e < m; ++e) gen_cnt[gvec[e]][e]++;
        for (int l = 1; l <= L; ++l) {
            const auto& Mv = eng.M_layer(l);
            const auto& Ev = eng.E_layer(l);
            for (int k = 0; k < N; ++k) {
                t_sum[l - 1][k] += static_cast<double>(Ev[k]);
                m_sum[l - 1][k] += static_cast<double>(Mv[k]);
                if (Mv[k] > 0) c_sum[l - 1][k] += static_cast<double>(Ev[k]) / static_cast<double>(Mv[k]);
            }
        }
        entropies.push_back(eng.total_entropy());
        if (collect_partitions) part_samples.push_back(eng.prior().export_bs());
        if ((s & 15) == 0) Rcpp::checkUserInterrupt();
    }

    double inv_n = 1.0 / static_cast<double>(std::max(1, nsamples));
    Rcpp::NumericMatrix layer_probs(m, L + 1);
    Rcpp::NumericVector sub_prob(m), cl_prob(m);
    for (int e = 0; e < m; ++e) {
        double p0 = static_cast<double>(gen_cnt[0][e]) * inv_n;
        layer_probs(e, 0) = p0;
        sub_prob[e] = p0;
        double pc = 0.0;
        for (int l = 1; l <= L; ++l) {
            double pl = static_cast<double>(gen_cnt[l][e]) * inv_n;
            layer_probs(e, l) = pl;
            pc += pl;
        }
        cl_prob[e] = pc;
    }

    Rcpp::List vstats(L);
    for (int l = 0; l < L; ++l) {
        Rcpp::NumericVector tv(N), mv(N), cv(N);
        for (int k = 0; k < N; ++k) {
            tv[k] = t_sum[l][k] * inv_n;
            mv[k] = m_sum[l][k] * inv_n;
            cv[k] = c_sum[l][k] * inv_n;
        }
        vstats[l] = Rcpp::DataFrame::create(
            Rcpp::Named("vertex") = Rcpp::seq(1, N),
            Rcpp::Named("t") = tv,
            Rcpp::Named("m") = mv,
            Rcpp::Named("c") = cv);
    }

    // Export final state
    Rcpp::IntegerVector fin_gen(m), fin_w0(m);
    Rcpp::List fin_ex(m);
    for (int e = 0; e < m; ++e) {
        fin_gen[e] = eng.gen()[e];
        fin_w0[e] = (eng.gen()[e] == 0)
            ? static_cast<int>(eng.prior().base().get_w(obs_edges(e, 0) - 1, obs_edges(e, 1) - 1))
            : 0;
        const auto& ev = eng.ex()[e];
        Rcpp::IntegerVector rv(ev.size());
        for (size_t j = 0; j < ev.size(); ++j) rv[j] = ev[j] + 1;
        fin_ex[e] = rv;
    }
    Rcpp::IntegerVector fin_loops(N);
    for (int u = 0; u < N; ++u) fin_loops[u] = static_cast<int>(eng.prior().base().get_w(u, u));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("substrate_prob") = sub_prob,
        Rcpp::Named("closure_prob") = cl_prob,
        Rcpp::Named("layer_probs") = layer_probs,
        Rcpp::Named("vertex_stats") = vstats,
        Rcpp::Named("final_gen") = fin_gen,
        Rcpp::Named("final_w0") = fin_w0,
        Rcpp::Named("final_loops") = fin_loops,
        Rcpp::Named("final_ex") = fin_ex,
        Rcpp::Named("final_bs") = eng.prior().export_bs(),
        Rcpp::Named("entropies") = entropies,
        Rcpp::Named("entropy") = eng.total_entropy(),
        Rcpp::Named("sbm_entropy") = eng.prior().entropy(),
        Rcpp::Named("closure_entropy") = eng.closure_entropy(),
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
