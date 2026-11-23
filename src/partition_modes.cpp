// graph-tool-R: partition label alignment, random label model
// (PartitionModeState), vertex marginals, consensus partitions, and posterior
// mode clustering (ModeClusterState), following Peixoto (2021) PRX 11, 021003.
//
// Copyright (C) 2021-2026 Tiago de Paula Peixoto (algorithms)
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
#include <numeric>
#include <vector>
#include "vendor/hash_map_wrap.hh"

namespace {

inline double lgam(double x) { return std::lgamma(x); }

inline double lbinom(double n, double k) {
    if (k == 0 || n == k)
        return 0;
    if (n < k || k < 0)
        return 0;
    return lgam(n + 1) - lgam(k + 1) - lgam(n - k + 1);
}

struct PMRng {
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

// Exact O(n^2 * m) Hungarian algorithm maximizing sum_i W[i][match[i]].
// Returns column index (0-based) matched to each row i, or -1 if unmatched.
std::vector<int> max_weight_matching(const std::vector<std::vector<long>>& W) {
    int nr = static_cast<int>(W.size());
    if (nr == 0)
        return {};
    int nc = static_cast<int>(W[0].size());
    if (nc == 0)
        return std::vector<int>(nr, -1);

    bool transposed = (nr > nc);
    int n = transposed ? nc : nr;
    int m = transposed ? nr : nc;
    long wmax = 0;
    for (int i = 0; i < nr; ++i)
        for (int j = 0; j < nc; ++j)
            if (W[i][j] > wmax)
                wmax = W[i][j];

    std::vector<long> u(n + 1, 0), v(m + 1, 0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    constexpr long INF = std::numeric_limits<long>::max() / 4;

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<long> minv(m + 1, INF);
        std::vector<char> used(m + 1, 0);
        do {
            used[j0] = 1;
            int i0 = p[j0];
            int j1 = 0;
            long delta = INF;
            for (int j = 1; j <= m; ++j) {
                if (!used[j]) {
                    long wij = transposed ? W[j - 1][i0 - 1] : W[i0 - 1][j - 1];
                    long cur = (wmax - wij) - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    std::vector<int> match(nr, -1);
    for (int j = 1; j <= m; ++j) {
        if (p[j] != 0) {
            if (transposed)
                match[j - 1] = p[j] - 1;
            else
                match[p[j] - 1] = j - 1;
        }
    }
    return match;
}

// Relabel a single vector of non-negative labels (with -1 for absent nodes) to
// align with reference vector y (also 0-based with -1 for absent).
// Returns pair (aligned_x, mapping from old x label -> new aligned label).
std::pair<std::vector<int>, gt_hash_map<int, int>>
align_pair(const std::vector<int>& x, const std::vector<int>& y) {
    size_t n = std::min(x.size(), y.size());
    std::vector<int> ux, uy;
    gt_hash_map<int, int> ix, iy;
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] < 0)
            continue;
        if (ix.find(x[i]) == ix.end()) {
            ix.emplace(x[i], static_cast<int>(ux.size()));
            ux.push_back(x[i]);
        }
    }
    int max_y = -1;
    for (size_t i = 0; i < y.size(); ++i) {
        if (y[i] < 0)
            continue;
        if (y[i] > max_y)
            max_y = y[i];
        if (iy.find(y[i]) == iy.end()) {
            iy.emplace(y[i], static_cast<int>(uy.size()));
            uy.push_back(y[i]);
        }
    }
    std::vector<std::vector<long>> W(ux.size(), std::vector<long>(uy.size(), 0));
    for (size_t i = 0; i < n; ++i) {
        if (x[i] >= 0 && y[i] >= 0)
            W[ix[x[i]]][iy[y[i]]]++;
    }
    std::vector<int> match = max_weight_matching(W);
    gt_hash_set<int> used_labels;
    gt_hash_map<int, int> remap;
    for (size_t r = 0; r < ux.size(); ++r) {
        int c = match[r];
        if (c >= 0 && W[r][c] > 0) {
            int lbl = uy[c];
            remap[ux[r]] = lbl;
            used_labels.insert(lbl);
        }
    }
    int next_lbl = 0;
    for (size_t r = 0; r < ux.size(); ++r) {
        if (remap.find(ux[r]) != remap.end())
            continue;
        while (used_labels.find(next_lbl) != used_labels.end() ||
               iy.find(next_lbl) != iy.end())
            ++next_lbl;
        remap[ux[r]] = next_lbl;
        used_labels.insert(next_lbl);
    }
    std::vector<int> out(x.size(), -1);
    for (size_t i = 0; i < x.size(); ++i)
        if (x[i] >= 0)
            out[i] = remap[x[i]];
    return {out, remap};
}

// Order labels in decreasing order of frequency (0, 1, 2, ...).
std::vector<int> order_labels_vec(const std::vector<int>& x, gt_hash_map<int, int>* remap_out = nullptr) {
    gt_hash_map<int, long> cnt;
    for (int v : x)
        if (v >= 0)
            cnt[v]++;
    std::vector<std::pair<long, int>> ord;
    ord.reserve(cnt.size());
    for (const auto& [k, c] : cnt)
        ord.emplace_back(-c, k);
    std::sort(ord.begin(), ord.end());
    gt_hash_map<int, int> remap;
    for (size_t i = 0; i < ord.size(); ++i)
        remap[ord[i].second] = static_cast<int>(i);
    std::vector<int> out(x.size(), -1);
    for (size_t i = 0; i < x.size(); ++i)
        if (x[i] >= 0)
            out[i] = remap[x[i]];
    if (remap_out)
        *remap_out = std::move(remap);
    return out;
}

typedef std::vector<std::vector<int>> HierPart;

HierPart order_nested_labels(const HierPart& x) {
    HierPart out(x.size());
    gt_hash_map<int, int> prev_remap;
    for (size_t l = 0; l < x.size(); ++l) {
        std::vector<int> cur;
        if (l == 0) {
            cur = x[0];
        } else {
            int max_idx = -1;
            for (const auto& [old_g, new_g] : prev_remap)
                if (new_g > max_idx)
                    max_idx = new_g;
            cur.assign(max_idx + 1, -1);
            for (const auto& [old_g, new_g] : prev_remap)
                if (old_g >= 0 && static_cast<size_t>(old_g) < x[l].size())
                    cur[new_g] = x[l][old_g];
        }
        gt_hash_map<int, int> remap;
        out[l] = order_labels_vec(cur, &remap);
        prev_remap = std::move(remap);
    }
    return out;
}

HierPart align_nested_pair(const HierPart& x, const HierPart& y) {
    size_t L = x.size();
    HierPart out(L);
    gt_hash_map<int, int> prev_remap;
    for (size_t l = 0; l < L; ++l) {
        std::vector<int> cur;
        if (l == 0) {
            cur = x[0];
        } else {
            int max_idx = -1;
            for (const auto& [old_g, new_g] : prev_remap)
                if (new_g > max_idx)
                    max_idx = new_g;
            cur.assign(max_idx + 1, -1);
            for (const auto& [old_g, new_g] : prev_remap)
                if (old_g >= 0 && static_cast<size_t>(old_g) < x[l].size())
                    cur[new_g] = x[l][old_g];
        }
        if (l < y.size()) {
            auto [aligned, remap] = align_pair(cur, y[l]);
            out[l] = std::move(aligned);
            prev_remap = std::move(remap);
        } else {
            gt_hash_map<int, int> remap;
            out[l] = order_labels_vec(cur, &remap);
            prev_remap = std::move(remap);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// PartitionMode: Peixoto (2021) Random Label Model for flat or nested
// partition ensembles.
//
// Each level l stores:
//   - aligned partition vectors (0-based labels, -1 for absent nodes)
//   - per-node marginal counts c[i][r] and total appearances M_i[i]
//   - total count per label label_count[r] and number of active labels B

class PartitionMode {
public:
    explicit PartitionMode(bool nested = false) : _nested(nested) {}

    size_t M() const { return _parts.size(); }
    size_t levels() const { return _levels.size(); }
    bool nested() const { return _nested; }

    size_t B(size_t l = 0) const {
        return l < _levels.size() ? _levels[l].B : 0;
    }

    size_t active_nodes(size_t l = 0) const {
        return l < _levels.size() ? _levels[l].N_active : 0;
    }

    const std::vector<HierPart>& partitions() const { return _parts; }

    // Aligns raw hierarchy `b` against current mode counts without mutating.
    HierPart align_partition(const HierPart& b) const {
        if (_parts.empty())
            return b;
        HierPart out(b.size());
        gt_hash_map<int, int> prev_remap;
        for (size_t l = 0; l < b.size(); ++l) {
            std::vector<int> cur;
            if (l == 0) {
                cur = b[0];
            } else {
                int max_idx = -1;
                for (const auto& [old_g, new_g] : prev_remap)
                    if (new_g > max_idx)
                        max_idx = new_g;
                cur.assign(max_idx + 1, -1);
                for (const auto& [old_g, new_g] : prev_remap)
                    if (old_g >= 0 && static_cast<size_t>(old_g) < b[l].size())
                        cur[new_g] = b[l][old_g];
            }
            if (l >= _levels.size() || _levels[l].B == 0) {
                // Keep labels contiguous starting from 0
                gt_hash_map<int, int> remap;
                out[l].assign(cur.size(), -1);
                for (size_t i = 0; i < cur.size(); ++i) {
                    if (cur[i] < 0)
                        continue;
                    auto it = remap.find(cur[i]);
                    if (it == remap.end())
                        it = remap.emplace(cur[i], static_cast<int>(remap.size())).first;
                    out[l][i] = it->second;
                }
                prev_remap = std::move(remap);
                continue;
            }
            const Level& lv = _levels[l];
            std::vector<int> ux;
            gt_hash_map<int, int> ix;
            for (int v : cur) {
                if (v >= 0 && ix.find(v) == ix.end()) {
                    ix.emplace(v, static_cast<int>(ux.size()));
                    ux.push_back(v);
                }
            }
            std::vector<int> active_s;
            gt_hash_map<int, int> is;
            for (const auto& [s, cnt] : lv.label_count) {
                if (cnt > 0) {
                    is.emplace(s, static_cast<int>(active_s.size()));
                    active_s.push_back(s);
                }
            }
            std::vector<std::vector<long>> W(ux.size(), std::vector<long>(active_s.size(), 0));
            size_t n = std::min(cur.size(), lv.counts.size());
            for (size_t i = 0; i < n; ++i) {
                if (cur[i] < 0)
                    continue;
                int row = ix[cur[i]];
                for (const auto& [s, c] : lv.counts[i]) {
                    auto it = is.find(s);
                    if (it != is.end())
                        W[row][it->second] += c;
                }
            }
            std::vector<int> match = max_weight_matching(W);
            gt_hash_map<int, int> remap;
            gt_hash_set<int> used;
            for (size_t r = 0; r < ux.size(); ++r) {
                int col = match[r];
                if (col >= 0 && W[r][col] > 0) {
                    int lbl = active_s[col];
                    remap[ux[r]] = lbl;
                    used.insert(lbl);
                }
            }
            int next_lbl = 0;
            for (size_t r = 0; r < ux.size(); ++r) {
                if (remap.find(ux[r]) != remap.end())
                    continue;
                while (used.find(next_lbl) != used.end() ||
                       (lv.label_count.find(next_lbl) != lv.label_count.end() &&
                        lv.label_count.at(next_lbl) > 0))
                    ++next_lbl;
                remap[ux[r]] = next_lbl;
                used.insert(next_lbl);
            }
            out[l].assign(cur.size(), -1);
            for (size_t i = 0; i < cur.size(); ++i)
                if (cur[i] >= 0)
                    out[l][i] = remap[cur[i]];
            prev_remap = std::move(remap);
        }
        return out;
    }

    size_t add_partition(const HierPart& b, bool relabel = true) {
        HierPart aligned = relabel ? align_partition(b) : b;
        size_t idx = _parts.size();
        _parts.push_back(aligned);
        update_counts(aligned, +1);
        return idx;
    }

    void remove_partition(size_t idx) {
        update_counts(_parts[idx], -1);
        if (idx + 1 < _parts.size())
            _parts[idx] = std::move(_parts.back());
        _parts.pop_back();
    }

    // Exact entropy change of adding `b` (with optional relabelling).
    double virtual_add(const HierPart& b, bool relabel = true) const {
        HierPart aligned = relabel ? align_partition(b) : b;
        double dS = 0.0;
        size_t L = std::max(_levels.size(), aligned.size());
        for (size_t l = 0; l < L; ++l) {
            if (l >= aligned.size())
                continue;
            const std::vector<int>& bl = aligned[l];
            if (l >= _levels.size() || _levels[l].N_active == 0) {
                // Brand new level
                gt_hash_set<int> uniq;
                size_t n_act = 0;
                for (int r : bl) {
                    if (r >= 0) {
                        uniq.insert(r);
                        ++n_act;
                    }
                }
                if (n_act > 0) {
                    double B1 = static_cast<double>(uniq.size());
                    // Each active node has M_i = 1, c_{i,r} = 1 -> lgamma(1 + B) - lgamma(B) - lgamma(2) = log(B)
                    dS += static_cast<double>(n_act) * std::log(B1) +
                          std::log(static_cast<double>(n_act));
                }
                continue;
            }
            const Level& lv = _levels[l];
            // Check if B or N_active changes
            gt_hash_set<int> new_labels;
            size_t added_nodes = 0;
            for (size_t i = 0; i < bl.size(); ++i) {
                int r = bl[i];
                if (r < 0)
                    continue;
                auto it = lv.label_count.find(r);
                if (it == lv.label_count.end() || it->second == 0)
                    new_labels.insert(r);
                if (i >= lv.node_total.size() || lv.node_total[i] == 0)
                    ++added_nodes;
            }
            size_t B0 = lv.B;
            size_t B1 = B0 + new_labels.size();
            size_t N0 = lv.N_active;
            size_t N1 = N0 + added_nodes;

            if (B1 == B0) {
                // Fast O(N_l) delta: only nodes present in bl change!
                if (N1 != N0)
                    dS += std::log(static_cast<double>(N1)) -
                          (N0 > 0 ? std::log(static_cast<double>(N0)) : 0.0);
                double Bd = static_cast<double>(B0);
                for (size_t i = 0; i < bl.size(); ++i) {
                    int r = bl[i];
                    if (r < 0)
                        continue;
                    long m0 = (i < lv.node_total.size()) ? lv.node_total[i] : 0;
                    long c0 = 0;
                    if (i < lv.counts.size()) {
                        auto it = lv.counts[i].find(r);
                        if (it != lv.counts[i].end())
                            c0 = it->second;
                    }
                    // (lgam(m0 + 1 + B) - lgam(m0 + B)) - (lgam(c0 + 2) - lgam(c0 + 1))
                    // = log(m0 + B) - log(c0 + 1)
                    dS += std::log(static_cast<double>(m0) + Bd) -
                          std::log(static_cast<double>(c0 + 1));
                }
            } else {
                // B changed: recompute level entropy difference directly
                double S0 = level_entropy(l);
                double S1 = (N1 > 0) ? std::log(static_cast<double>(N1)) : 0.0;
                double lgB1 = lgam(static_cast<double>(B1));
                size_t n_all = std::max(lv.node_total.size(), bl.size());
                for (size_t i = 0; i < n_all; ++i) {
                    long m0 = (i < lv.node_total.size()) ? lv.node_total[i] : 0;
                    int r_add = (i < bl.size()) ? bl[i] : -1;
                    long m1 = m0 + (r_add >= 0 ? 1 : 0);
                    if (m1 == 0)
                        continue;
                    S1 += lgam(static_cast<double>(m1 + B1)) - lgB1;
                    if (i < lv.counts.size()) {
                        for (const auto& [r, c] : lv.counts[i]) {
                            long c1 = c + (r == r_add ? 1 : 0);
                            S1 -= lgam(static_cast<double>(c1 + 1));
                        }
                        if (r_add >= 0 && lv.counts[i].find(r_add) == lv.counts[i].end())
                            S1 -= lgam(2.0);
                    } else if (r_add >= 0) {
                        S1 -= lgam(2.0);
                    }
                }
                dS += S1 - S0;
            }
        }
        return dS;
    }

    double replace_partitions(PMRng& rng) {
        size_t m = _parts.size();
        if (m <= 1)
            return 0.0;
        double S0 = entropy();
        std::vector<size_t> order(m);
        std::iota(order.begin(), order.end(), 0);
        rng.shuffle(order);
        for (size_t idx : order) {
            HierPart cur = _parts[idx];
            update_counts(cur, -1);
            HierPart rel = align_partition(cur);
            _parts[idx] = std::move(rel);
            update_counts(_parts[idx], +1);
        }
        return entropy() - S0;
    }

    void converge_alignment(PMRng& rng, size_t max_iter = 100) {
        for (size_t it = 0; it < max_iter; ++it) {
            double delta = replace_partitions(rng);
            if (std::abs(delta) <= 1e-6)
                break;
        }
    }

    // Relabels this mode to maximize overlap with `other` mode.
    void align_to(const PartitionMode& other) {
        if (_parts.empty() || other._parts.empty())
            return;
        size_t L = _levels.size();
        gt_hash_map<int, int> prev_remap;
        for (size_t l = 0; l < L; ++l) {
            if (l >= other._levels.size())
                break;
            // First permute node indices at level l > 0 according to prev_remap
            if (l > 0) {
                int max_idx = -1;
                for (const auto& [old_g, new_g] : prev_remap)
                    if (new_g > max_idx)
                        max_idx = new_g;
                for (auto& p : _parts) {
                    if (l >= p.size())
                        continue;
                    std::vector<int> cur(max_idx + 1, -1);
                    for (const auto& [old_g, new_g] : prev_remap)
                        if (old_g >= 0 && static_cast<size_t>(old_g) < p[l].size())
                            cur[new_g] = p[l][old_g];
                    p[l] = std::move(cur);
                }
                rebuild_level(l);
            }
            const Level& lv = _levels[l];
            const Level& ov = other._levels[l];
            std::vector<int> ux, uy;
            gt_hash_map<int, int> ix, iy;
            for (const auto& [r, c] : lv.label_count) {
                if (c > 0) {
                    ix.emplace(r, static_cast<int>(ux.size()));
                    ux.push_back(r);
                }
            }
            for (const auto& [s, c] : ov.label_count) {
                if (c > 0) {
                    iy.emplace(s, static_cast<int>(uy.size()));
                    uy.push_back(s);
                }
            }
            std::vector<std::vector<long>> W(ux.size(), std::vector<long>(uy.size(), 0));
            size_t n = std::min(lv.counts.size(), ov.counts.size());
            for (size_t i = 0; i < n; ++i) {
                for (const auto& [r, cr] : lv.counts[i]) {
                    auto itr = ix.find(r);
                    if (itr == ix.end())
                        continue;
                    for (const auto& [s, cs] : ov.counts[i]) {
                        auto its = iy.find(s);
                        if (its != iy.end())
                            W[itr->second][its->second] += cr * cs;
                    }
                }
            }
            std::vector<int> match = max_weight_matching(W);
            gt_hash_map<int, int> remap;
            gt_hash_set<int> used;
            for (size_t r = 0; r < ux.size(); ++r) {
                int col = match[r];
                if (col >= 0 && W[r][col] > 0) {
                    int lbl = uy[col];
                    remap[ux[r]] = lbl;
                    used.insert(lbl);
                }
            }
            int next_lbl = 0;
            for (size_t r = 0; r < ux.size(); ++r) {
                if (remap.find(ux[r]) != remap.end())
                    continue;
                while (used.find(next_lbl) != used.end() ||
                       (ov.label_count.find(next_lbl) != ov.label_count.end() &&
                        ov.label_count.at(next_lbl) > 0))
                    ++next_lbl;
                remap[ux[r]] = next_lbl;
                used.insert(next_lbl);
            }
            for (auto& p : _parts) {
                if (l >= p.size())
                    continue;
                for (int& v : p[l])
                    if (v >= 0)
                        v = remap[v];
            }
            rebuild_level(l);
            prev_remap = std::move(remap);
        }
    }

    double level_entropy(size_t l) const {
        if (l >= _levels.size() || _levels[l].N_active == 0)
            return 0.0;
        const Level& lv = _levels[l];
        double S = std::log(static_cast<double>(lv.N_active));
        double B = static_cast<double>(lv.B);
        double lgB = lgam(B);
        for (size_t i = 0; i < lv.node_total.size(); ++i) {
            long mi = lv.node_total[i];
            if (mi == 0)
                continue;
            S += lgam(static_cast<double>(mi) + B) - lgB;
            for (const auto& [r, c] : lv.counts[i])
                S -= lgam(static_cast<double>(c + 1));
        }
        return S;
    }

    double entropy() const {
        double S = 0.0;
        for (size_t l = 0; l < _levels.size(); ++l)
            S += level_entropy(l);
        return S;
    }

    double posterior_entropy(bool MLE, size_t l = 0) const {
        if (l >= _levels.size() || _levels[l].N_active == 0)
            return 0.0;
        const Level& lv = _levels[l];
        double B = static_cast<double>(lv.B);
        double H = 0.0;
        for (size_t i = 0; i < lv.node_total.size(); ++i) {
            long mi = lv.node_total[i];
            if (mi == 0)
                continue;
            if (MLE) {
                double denom = static_cast<double>(mi);
                for (const auto& [r, c] : lv.counts[i]) {
                    if (c > 0) {
                        double p = static_cast<double>(c) / denom;
                        H -= p * std::log(p);
                    }
                }
            } else {
                double denom = static_cast<double>(mi) + B;
                size_t seen = 0;
                for (const auto& [r, c] : lv.counts[i]) {
                    if (c > 0) {
                        double p = static_cast<double>(c + 1) / denom;
                        H -= p * std::log(p);
                        ++seen;
                    }
                }
                if (lv.B > seen) {
                    double p0 = 1.0 / denom;
                    H -= static_cast<double>(lv.B - seen) * p0 * std::log(p0);
                }
            }
        }
        return H;
    }

    double posterior_cdev(bool MLE, size_t l = 0) const {
        if (l >= _levels.size() || _levels[l].N_active == 0)
            return 0.0;
        const Level& lv = _levels[l];
        double B = static_cast<double>(lv.B);
        double sum_max = 0.0;
        for (size_t i = 0; i < lv.node_total.size(); ++i) {
            long mi = lv.node_total[i];
            if (mi == 0)
                continue;
            long best = 0;
            for (const auto& [r, c] : lv.counts[i])
                if (c > best)
                    best = c;
            if (MLE)
                sum_max += static_cast<double>(best) / static_cast<double>(mi);
            else
                sum_max += static_cast<double>(best + 1) / (static_cast<double>(mi) + B);
        }
        return 1.0 - sum_max / static_cast<double>(lv.N_active);
    }

    // Export contiguous 1-based R summary of this mode.
    Rcpp::List to_r() const {
        size_t L = _levels.size();
        // Map active labels at each level to contiguous 0..B_l-1 ordered by decreasing total count.
        std::vector<std::vector<int>> active_labels(L);
        std::vector<gt_hash_map<int, int>> label_to_col(L);
        for (size_t l = 0; l < L; ++l) {
            std::vector<std::pair<long, int>> ord;
            for (const auto& [r, c] : _levels[l].label_count)
                if (c > 0)
                    ord.emplace_back(-c, r);
            std::sort(ord.begin(), ord.end());
            active_labels[l].resize(ord.size());
            for (size_t j = 0; j < ord.size(); ++j) {
                active_labels[l][j] = ord[j].second;
                label_to_col[l][ord[j].second] = static_cast<int>(j);
            }
        }

        Rcpp::List marg_list(L);
        Rcpp::List cons_list(L);
        Rcpp::NumericVector level_ent(L);
        Rcpp::IntegerVector B_vec(L);

        for (size_t l = 0; l < L; ++l) {
            const Level& lv = _levels[l];
            level_ent[l] = level_entropy(l);
            int Bl = static_cast<int>(active_labels[l].size());
            B_vec[l] = Bl;
            // Active node indices in order (for l=0: 0..N-1; for l>0: active_labels[l-1])
            std::vector<int> node_ids;
            if (l == 0) {
                node_ids.resize(lv.node_total.size());
                std::iota(node_ids.begin(), node_ids.end(), 0);
            } else {
                node_ids = active_labels[l - 1];
            }
            int Nl = static_cast<int>(node_ids.size());
            Rcpp::IntegerMatrix mat(Nl, Bl);
            Rcpp::IntegerVector cons(Nl);
            for (int row = 0; row < Nl; ++row) {
                int i = node_ids[row];
                long best_c = -1;
                int best_col = 0;
                if (i >= 0 && static_cast<size_t>(i) < lv.counts.size()) {
                    for (const auto& [r, c] : lv.counts[i]) {
                        auto it = label_to_col[l].find(r);
                        if (it != label_to_col[l].end()) {
                            int col = it->second;
                            mat(row, col) = static_cast<int>(c);
                            if (c > best_c || (c == best_c && col < best_col)) {
                                best_c = c;
                                best_col = col;
                            }
                        }
                    }
                }
                cons[row] = (best_c > 0) ? (best_col + 1) : NA_INTEGER;
            }
            marg_list[l] = mat;
            cons_list[l] = cons;
        }

        // Export aligned partitions (1-based contiguous via label_to_col)
        Rcpp::List parts_r(_parts.size());
        for (size_t m = 0; m < _parts.size(); ++m) {
            if (!_nested) {
                const auto& b0 = _parts[m][0];
                Rcpp::IntegerVector v(b0.size());
                for (size_t i = 0; i < b0.size(); ++i)
                    v[i] = (b0[i] >= 0) ? (label_to_col[0][b0[i]] + 1) : NA_INTEGER;
                parts_r[m] = v;
            } else {
                Rcpp::List hl(L);
                for (size_t l = 0; l < L; ++l) {
                    std::vector<int> node_ids;
                    if (l == 0) {
                        node_ids.resize(_levels[0].node_total.size());
                        std::iota(node_ids.begin(), node_ids.end(), 0);
                    } else {
                        node_ids = active_labels[l - 1];
                    }
                    Rcpp::IntegerVector v(node_ids.size());
                    for (size_t row = 0; row < node_ids.size(); ++row) {
                        int i = node_ids[row];
                        if (l < _parts[m].size() && i >= 0 &&
                            static_cast<size_t>(i) < _parts[m][l].size() &&
                            _parts[m][l][i] >= 0) {
                            v[row] = label_to_col[l][_parts[m][l][i]] + 1;
                        } else {
                            v[row] = NA_INTEGER;
                        }
                    }
                    hl[l] = v;
                }
                parts_r[m] = hl;
            }
        }

        return Rcpp::List::create(
            Rcpp::Named("nested") = _nested,
            Rcpp::Named("M") = static_cast<int>(_parts.size()),
            Rcpp::Named("B") = _nested ? SEXP(B_vec) : SEXP(Rcpp::wrap(B_vec.size() > 0 ? B_vec[0] : 0)),
            Rcpp::Named("entropy") = entropy(),
            Rcpp::Named("level_entropy") = level_ent,
            Rcpp::Named("posterior_entropy") = posterior_entropy(true, 0),
            Rcpp::Named("posterior_entropy_mean") = posterior_entropy(false, 0),
            Rcpp::Named("posterior_cdev") = posterior_cdev(true, 0),
            Rcpp::Named("posterior_cdev_mean") = posterior_cdev(false, 0),
            Rcpp::Named("marginals") = _nested ? SEXP(marg_list) : SEXP(marg_list[0]),
            Rcpp::Named("consensus") = _nested ? SEXP(cons_list) : SEXP(cons_list[0]),
            Rcpp::Named("partitions") = parts_r);
    }

private:
    struct Level {
        std::vector<gt_hash_map<int, long>> counts;
        std::vector<long> node_total;
        gt_hash_map<int, long> label_count;
        size_t B = 0;
        size_t N_active = 0;
    };

    bool _nested;
    std::vector<HierPart> _parts;
    std::vector<Level> _levels;

    void update_counts(const HierPart& p, long delta) {
        if (delta > 0 && p.size() > _levels.size())
            _levels.resize(p.size());
        for (size_t l = 0; l < p.size(); ++l) {
            Level& lv = _levels[l];
            if (delta > 0 && p[l].size() > lv.counts.size()) {
                lv.counts.resize(p[l].size());
                lv.node_total.resize(p[l].size(), 0);
            }
            for (size_t i = 0; i < p[l].size(); ++i) {
                int r = p[l][i];
                if (r < 0)
                    continue;
                long old_nt = lv.node_total[i];
                lv.node_total[i] += delta;
                if (old_nt == 0 && lv.node_total[i] > 0)
                    ++lv.N_active;
                else if (old_nt > 0 && lv.node_total[i] == 0)
                    --lv.N_active;

                auto it = lv.counts[i].find(r);
                if (it == lv.counts[i].end()) {
                    if (delta > 0)
                        lv.counts[i].emplace(r, delta);
                } else {
                    it->second += delta;
                    if (it->second == 0)
                        lv.counts[i].erase(it);
                }

                long& lc = lv.label_count[r];
                long old_lc = lc;
                lc += delta;
                if (old_lc == 0 && lc > 0)
                    ++lv.B;
                else if (old_lc > 0 && lc == 0)
                    --lv.B;
            }
        }
    }

    void rebuild_level(size_t l) {
        Level& lv = _levels[l];
        lv.counts.clear();
        lv.node_total.clear();
        lv.label_count.clear();
        lv.B = 0;
        lv.N_active = 0;
        for (const auto& p : _parts) {
            if (l >= p.size())
                continue;
            if (p[l].size() > lv.counts.size()) {
                lv.counts.resize(p[l].size());
                lv.node_total.resize(p[l].size(), 0);
            }
            for (size_t i = 0; i < p[l].size(); ++i) {
                int r = p[l][i];
                if (r < 0)
                    continue;
                if (lv.node_total[i]++ == 0)
                    ++lv.N_active;
                lv.counts[i][r]++;
                if (lv.label_count[r]++ == 0)
                    ++lv.B;
            }
        }
    }
};

// Parse R partition or nested hierarchy into 0-based HierPart (-1 for NA).
HierPart parse_partition(SEXP item, bool nested) {
    auto parse_vec = [](Rcpp::IntegerVector x) {
        std::vector<int> v(x.size(), -1);
        for (R_xlen_t i = 0; i < x.size(); ++i)
            if (x[i] != NA_INTEGER)
                v[i] = x[i] - 1;
        return v;
    };
    if (!nested) {
        return {parse_vec(Rcpp::as<Rcpp::IntegerVector>(item))};
    }
    Rcpp::List hl = Rcpp::as<Rcpp::List>(item);
    HierPart out;
    out.reserve(hl.size());
    for (R_xlen_t l = 0; l < hl.size(); ++l)
        out.push_back(parse_vec(Rcpp::as<Rcpp::IntegerVector>(hl[l])));
    return out;
}

double cluster_prior(size_t M, const std::vector<size_t>& sizes) {
    if (M == 0)
        return 0.0;
    size_t K = 0;
    double S = 0.0;
    for (size_t sz : sizes) {
        if (sz > 0) {
            ++K;
            S -= lgam(static_cast<double>(sz + 1));
        }
    }
    S += lbinom(static_cast<double>(M - 1), static_cast<double>(K - 1)) +
         lgam(static_cast<double>(M + 1)) + std::log(static_cast<double>(M));
    return S;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_align_partition_labels(Rcpp::IntegerVector x, Rcpp::IntegerVector y) {
    std::vector<int> vx(x.size(), -1), vy(y.size(), -1);
    for (R_xlen_t i = 0; i < x.size(); ++i)
        if (x[i] != NA_INTEGER) vx[i] = x[i] - 1;
    for (R_xlen_t i = 0; i < y.size(); ++i)
        if (y[i] != NA_INTEGER) vy[i] = y[i] - 1;
    auto [aligned, remap] = align_pair(vx, vy);
    Rcpp::IntegerVector out(aligned.size());
    for (size_t i = 0; i < aligned.size(); ++i)
        out[i] = (aligned[i] >= 0) ? (aligned[i] + 1) : NA_INTEGER;
    return out;
}

// [[Rcpp::export]]
double cpp_partition_overlap(Rcpp::IntegerVector x, Rcpp::IntegerVector y, bool norm) {
    if (x.size() != y.size())
        Rcpp::stop("x and y must have the same length");
    size_t N = static_cast<size_t>(x.size());
    if (N == 0)
        return norm ? 1.0 : 0.0;
    std::vector<int> vx(N, -1), vy(N, -1);
    for (size_t i = 0; i < N; ++i) {
        if (x[i] != NA_INTEGER) vx[i] = x[i] - 1;
        if (y[i] != NA_INTEGER) vy[i] = y[i] - 1;
    }
    auto [aligned, remap] = align_pair(vx, vy);
    size_t m = 0;
    for (size_t i = 0; i < N; ++i)
        if (aligned[i] >= 0 && aligned[i] == vy[i])
            ++m;
    return norm ? static_cast<double>(m) / static_cast<double>(N) : static_cast<double>(m);
}

// [[Rcpp::export]]
Rcpp::List cpp_align_nested_partition_labels(Rcpp::List x, Rcpp::List y) {
    HierPart hx = parse_partition(x, true);
    HierPart hy = parse_partition(y, true);
    HierPart aligned = align_nested_pair(hx, hy);
    Rcpp::List out(aligned.size());
    for (size_t l = 0; l < aligned.size(); ++l) {
        Rcpp::IntegerVector v(aligned[l].size());
        for (size_t i = 0; i < aligned[l].size(); ++i)
            v[i] = (aligned[l][i] >= 0) ? (aligned[l][i] + 1) : NA_INTEGER;
        out[l] = v;
    }
    return out;
}

// [[Rcpp::export]]
double cpp_nested_partition_overlap(Rcpp::List x, Rcpp::List y, bool norm) {
    HierPart hx = parse_partition(x, true);
    HierPart hy = order_nested_labels(parse_partition(y, true));
    HierPart ax = align_nested_pair(hx, hy);
    size_t L = std::min(ax.size(), hy.size());
    size_t m = 0, N = 0;
    for (size_t l = 0; l < L; ++l) {
        const auto& xl = ax[l];
        const auto& yl = hy[l];
        size_t n_min = std::min(xl.size(), yl.size());
        for (size_t i = 0; i < n_min; ++i)
            if (xl[i] >= 0 && xl[i] == yl[i])
                ++m;
        size_t nx = 0, ny = 0;
        for (int v : xl) if (v >= 0) ++nx;
        for (int v : yl) if (v >= 0) ++ny;
        N += std::max(nx, ny);
    }
    if (norm) {
        if (N <= L)
            return 1.0;
        return 1.0 - static_cast<double>(N - m) / static_cast<double>(N - L);
    }
    return static_cast<double>(m);
}

// [[Rcpp::export]]
Rcpp::List cpp_partition_mode(Rcpp::List partitions, bool relabel, bool nested, bool converge) {
    PartitionMode mode(nested);
    for (R_xlen_t i = 0; i < partitions.size(); ++i)
        mode.add_partition(parse_partition(partitions[i], nested), relabel);
    if (converge && mode.M() > 1) {
        PMRng rng;
        mode.converge_alignment(rng);
    }
    return mode.to_r();
}

// [[Rcpp::export]]
Rcpp::List cpp_mode_cluster(Rcpp::List partitions, bool nested, bool relabel,
                            Rcpp::Nullable<Rcpp::IntegerVector> initial_b,
                            int niter, bool converge) {
    size_t M = static_cast<size_t>(partitions.size());
    if (M == 0)
        Rcpp::stop("partitions must be non-empty");
    std::vector<HierPart> raw(M);
    for (size_t i = 0; i < M; ++i)
        raw[i] = parse_partition(partitions[i], nested);

    PMRng rng;
    std::vector<int> assign(M, 0);
    if (initial_b.isNotNull()) {
        Rcpp::IntegerVector ib(initial_b);
        if (static_cast<size_t>(ib.size()) != M)
            Rcpp::stop("initial_b must have one entry per partition");
        gt_hash_map<int, int> map;
        for (size_t i = 0; i < M; ++i) {
            auto it = map.find(ib[i]);
            if (it == map.end())
                it = map.emplace(ib[i], static_cast<int>(map.size())).first;
            assign[i] = it->second;
        }
    }

    auto build_modes = [&](const std::vector<int>& b, bool do_relabel, bool do_conv) {
        int K = 0;
        for (int c : b)
            if (c + 1 > K)
                K = c + 1;
        std::vector<PartitionMode> modes;
        modes.reserve(K);
        for (int k = 0; k < K; ++k)
            modes.emplace_back(nested);
        for (size_t i = 0; i < M; ++i)
            modes[b[i]].add_partition(raw[i], do_relabel);
        if (do_conv)
            for (auto& md : modes)
                if (md.M() > 1)
                    md.converge_alignment(rng, 20);
        return modes;
    };

    auto eval_total_entropy = [&](const std::vector<PartitionMode>& modes) {
        std::vector<size_t> sizes;
        sizes.reserve(modes.size());
        double S = 0.0;
        for (const auto& md : modes) {
            if (md.M() > 0) {
                sizes.push_back(md.M());
                S += md.entropy();
            }
        }
        return S + cluster_prior(M, sizes);
    };

    // If niter == 0, evaluate exact entropy/modes for initial_b without optimizing.
    if (niter <= 0) {
        std::vector<PartitionMode> modes = build_modes(assign, relabel, converge);
        double S = eval_total_entropy(modes);
        Rcpp::List mode_list(modes.size());
        Rcpp::NumericVector weights(modes.size());
        double H_w = 0.0;
        for (size_t k = 0; k < modes.size(); ++k) {
            mode_list[k] = modes[k].to_r();
            double w = static_cast<double>(modes[k].M()) / static_cast<double>(M);
            weights[k] = w;
            if (w > 0)
                H_w -= w * std::log(w);
        }
        Rcpp::IntegerVector cl_out(M);
        for (size_t i = 0; i < M; ++i)
            cl_out[i] = assign[i] + 1;
        return Rcpp::List::create(
            Rcpp::Named("entropy") = S,
            Rcpp::Named("clusters") = cl_out,
            Rcpp::Named("K") = static_cast<int>(modes.size()),
            Rcpp::Named("Be") = std::exp(H_w),
            Rcpp::Named("weights") = weights,
            Rcpp::Named("modes") = mode_list);
    }

    // Optimize ModeClusterState via greedy single-partition moves + cluster splits & merges.
    // Maintain a pool of PartitionMode objects and track partition index inside its mode.
    struct ClusterPool {
        bool nested;
        size_t M;
        const std::vector<HierPart>& raw;
        std::vector<PartitionMode> modes;
        std::vector<std::vector<size_t>> members;  // partition indices inside each mode
        std::vector<int> cl;
        std::vector<size_t> pos;                   // index inside modes[cl[i]]

        ClusterPool(bool nst, const std::vector<HierPart>& r, const std::vector<int>& init_b)
            : nested(nst), M(r.size()), raw(r), cl(init_b), pos(r.size(), 0) {
            int K = 0;
            for (int c : init_b)
                if (c + 1 > K) K = c + 1;
            modes.reserve(K + 4);
            members.resize(K);
            for (int k = 0; k < K; ++k)
                modes.emplace_back(nested);
            for (size_t i = 0; i < M; ++i) {
                int k = cl[i];
                pos[i] = modes[k].add_partition(raw[i], true);
                members[k].push_back(i);
            }
        }

        double total_entropy() const {
            std::vector<size_t> sizes;
            double S = 0.0;
            for (const auto& md : modes) {
                if (md.M() > 0) {
                    sizes.push_back(md.M());
                    S += md.entropy();
                }
            }
            return S + cluster_prior(M, sizes);
        }

        void move_part(size_t i, int dst) {
            int src = cl[i];
            if (src == dst)
                return;
            if (dst >= static_cast<int>(modes.size())) {
                modes.emplace_back(nested);
                members.emplace_back();
            }
            size_t p_idx = pos[i];
            modes[src].remove_partition(p_idx);
            // Fix swap-and-pop bookkeeping for members[src]
            size_t last_i = members[src].back();
            if (p_idx + 1 < members[src].size()) {
                members[src][p_idx] = last_i;
                pos[last_i] = p_idx;
            }
            members[src].pop_back();

            cl[i] = dst;
            pos[i] = modes[dst].add_partition(raw[i], true);
            members[dst].push_back(i);
        }

        // Greedy single-partition sweep allowing creation of a new cluster.
        bool greedy_sweep(PMRng& rng) {
            std::vector<size_t> order(M);
            std::iota(order.begin(), order.end(), 0);
            rng.shuffle(order);
            bool changed = false;

            std::vector<size_t> sizes;
            auto prior_delta = [&](size_t s_src, size_t s_dst, size_t K) {
                size_t K1 = K + (s_dst == 0 ? 1 : 0) - (s_src == 1 ? 1 : 0);
                double d = lbinom(static_cast<double>(M - 1), static_cast<double>(K1 - 1)) -
                           lbinom(static_cast<double>(M - 1), static_cast<double>(K - 1));
                d += lgam(static_cast<double>(s_src + 1)) - lgam(static_cast<double>(s_src));
                d += lgam(static_cast<double>(s_dst + 1)) - lgam(static_cast<double>(s_dst + 2));
                return d;
            };

            size_t K_act = 0;
            for (const auto& md : modes)
                if (md.M() > 0) ++K_act;

            for (size_t i : order) {
                int src = cl[i];
                size_t s_src = modes[src].M();
                double S_src_before = modes[src].entropy();

                // Temporarily remove i from src
                size_t p_idx = pos[i];
                modes[src].remove_partition(p_idx);
                size_t last_i = members[src].back();
                if (p_idx + 1 < members[src].size()) {
                    members[src][p_idx] = last_i;
                    pos[last_i] = p_idx;
                }
                members[src].pop_back();
                double dS_remove = modes[src].entropy() - S_src_before;

                int best_dst = src;
                double best_dS = 0.0;
                int empty_slot = -1;

                for (int k = 0; k < static_cast<int>(modes.size()); ++k) {
                    if (k == src)
                        continue;
                    if (modes[k].M() == 0) {
                        if (empty_slot < 0)
                            empty_slot = k;
                        continue;
                    }
                    double dS = dS_remove + modes[k].virtual_add(raw[i], true) +
                                prior_delta(s_src, modes[k].M(), K_act);
                    if (dS < best_dS - 1e-9) {
                        best_dS = dS;
                        best_dst = k;
                    }
                }
                // Try creating a new cluster (only if src wasn't singleton)
                if (s_src > 1) {
                    PartitionMode tmp(nested);
                    double dS_new = dS_remove + tmp.virtual_add(raw[i], false) +
                                    prior_delta(s_src, 0, K_act);
                    if (dS_new < best_dS - 1e-9) {
                        best_dS = dS_new;
                        best_dst = (empty_slot >= 0) ? empty_slot : static_cast<int>(modes.size());
                    }
                }

                if (best_dst >= static_cast<int>(modes.size())) {
                    modes.emplace_back(nested);
                    members.emplace_back();
                }
                bool src_emptied = (modes[src].M() == 0);
                bool dst_was_empty = (modes[best_dst].M() == 0);
                cl[i] = best_dst;
                pos[i] = modes[best_dst].add_partition(raw[i], true);
                members[best_dst].push_back(i);
                if (best_dst != src) {
                    if (src_emptied) --K_act;
                    if (dst_was_empty) ++K_act;
                    changed = true;
                }
            }
            return changed;
        }

        // Try splitting each cluster with >= 2 partitions using 2-seed split + greedy sweeps.
        bool try_splits(PMRng& rng) {
            bool any = false;
            std::vector<int> act;
            for (int k = 0; k < static_cast<int>(modes.size()); ++k)
                if (modes[k].M() >= 2)
                    act.push_back(k);
            for (int k : act) {
                if (modes[k].M() < 2)
                    continue;
                double S_before = total_entropy();
                std::vector<size_t> mem = members[k];
                if (mem.size() < 2)
                    continue;
                // Pick two distant seeds inside cluster k
                size_t s1 = mem[rng.index(mem.size())];
                size_t s2 = s1;
                double worst_ov = 2.0;
                // Sample up to 32 candidates to find the most dissimilar seed
                size_t ncheck = std::min<size_t>(mem.size(), 32);
                for (size_t t = 0; t < ncheck; ++t) {
                    size_t cand = (ncheck == mem.size()) ? mem[t] : mem[rng.index(mem.size())];
                    if (cand == s1)
                        continue;
                    auto [al, rm] = align_pair(raw[cand][0], raw[s1][0]);
                    size_t match_cnt = 0;
                    for (size_t i = 0; i < al.size(); ++i)
                        if (al[i] >= 0 && al[i] == raw[s1][0][i])
                            ++match_cnt;
                    double ov = static_cast<double>(match_cnt) / std::max<size_t>(1, al.size());
                    if (ov < worst_ov) {
                        worst_ov = ov;
                        s2 = cand;
                    }
                }
                if (s2 == s1 || worst_ov >= 1.0 - 1e-12)
                    continue;

                // Create temporary split into k and new_k
                int new_k = -1;
                for (int j = 0; j < static_cast<int>(modes.size()); ++j) {
                    if (modes[j].M() == 0) {
                        new_k = j;
                        break;
                    }
                }
                if (new_k < 0) {
                    new_k = static_cast<int>(modes.size());
                    modes.emplace_back(nested);
                    members.emplace_back();
                }
                move_part(s2, new_k);
                // Seed initial 2-way split by maximum overlap with s1 vs s2
                auto ov_to = [&](size_t a, size_t b_idx) {
                    auto [al, rm] = align_pair(raw[a][0], raw[b_idx][0]);
                    size_t match_cnt = 0;
                    for (size_t j = 0; j < al.size(); ++j)
                        if (al[j] >= 0 && al[j] == raw[b_idx][0][j])
                            ++match_cnt;
                    return match_cnt;
                };
                for (size_t i : mem) {
                    if (i == s1 || i == s2)
                        continue;
                    if (ov_to(i, s2) > ov_to(i, s1))
                        move_part(i, new_k);
                }
                // Refine split with local sweeps between k and new_k
                for (int sw = 0; sw < 5; ++sw) {
                    bool moved = false;
                    for (size_t i : mem) {
                        int cur = cl[i];
                        int oth = (cur == k) ? new_k : k;
                        if (modes[cur].M() <= 1)
                            continue;
                        double sc_before = modes[cur].entropy() + modes[oth].entropy();
                        move_part(i, oth);
                        double sc_after = modes[cur].entropy() + modes[oth].entropy();
                        if (sc_after < sc_before - 1e-9) {
                            moved = true;
                        } else {
                            move_part(i, cur);
                        }
                    }
                    if (!moved)
                        break;
                }
                double S_after = total_entropy();
                if (S_after < S_before - 1e-8) {
                    any = true;
                } else {
                    // Undo split: move everything in new_k back to k
                    std::vector<size_t> back = members[new_k];
                    for (size_t i : back)
                        move_part(i, k);
                }
            }
            return any;
        }

        // Try merging pairs of active clusters.
        bool try_merges() {
            bool any = false;
            while (true) {
                std::vector<int> act;
                for (int k = 0; k < static_cast<int>(modes.size()); ++k)
                    if (modes[k].M() > 0)
                        act.push_back(k);
                if (act.size() <= 1)
                    break;
                double S_cur = total_entropy();
                int best_a = -1, best_b = -1;
                double best_S = S_cur - 1e-8;

                std::vector<size_t> base_sizes;
                double base_mode_sum = 0.0;
                for (int k : act) {
                    base_sizes.push_back(modes[k].M());
                    base_mode_sum += modes[k].entropy();
                }

                for (size_t ia = 0; ia < act.size(); ++ia) {
                    for (size_t ib = ia + 1; ib < act.size(); ++ib) {
                        int ka = act[ia], kb = act[ib];
                        PartitionMode merged(nested);
                        for (size_t i : members[ka])
                            merged.add_partition(raw[i], true);
                        for (size_t i : members[kb])
                            merged.add_partition(raw[i], true);
                        std::vector<size_t> new_sizes;
                        new_sizes.reserve(act.size() - 1);
                        for (size_t ic = 0; ic < act.size(); ++ic) {
                            if (ic == ia)
                                new_sizes.push_back(base_sizes[ia] + base_sizes[ib]);
                            else if (ic != ib)
                                new_sizes.push_back(base_sizes[ic]);
                        }
                        double S_cand = base_mode_sum - modes[ka].entropy() -
                                        modes[kb].entropy() + merged.entropy() +
                                        cluster_prior(M, new_sizes);
                        if (S_cand < best_S) {
                            best_S = S_cand;
                            best_a = ka;
                            best_b = kb;
                        }
                    }
                }
                if (best_a < 0)
                    break;
                std::vector<size_t> to_move = members[best_b];
                for (size_t i : to_move)
                    move_part(i, best_a);
                any = true;
            }
            return any;
        }
    };

    ClusterPool pool(nested, raw, assign);
    for (int it = 0; it < niter; ++it) {
        bool c1 = pool.try_splits(rng);
        bool c2 = pool.greedy_sweep(rng);
        bool c3 = pool.try_merges();
        Rcpp::checkUserInterrupt();
        if (!c1 && !c2 && !c3)
            break;
    }

    // Collect nonempty clusters, converge alignments, sort by decreasing size,
    // and align secondary modes to the dominant mode.
    std::vector<int> act;
    for (int k = 0; k < static_cast<int>(pool.modes.size()); ++k)
        if (pool.modes[k].M() > 0)
            act.push_back(k);
    std::sort(act.begin(), act.end(), [&](int a, int b) {
        return pool.modes[a].M() > pool.modes[b].M();
    });

    if (converge) {
        for (int k : act)
            if (pool.modes[k].M() > 1)
                pool.modes[k].converge_alignment(rng);
    }
    for (size_t idx = 1; idx < act.size(); ++idx)
        pool.modes[act[idx]].align_to(pool.modes[act[0]]);

    gt_hash_map<int, int> old_to_new;
    Rcpp::List mode_list(act.size());
    Rcpp::NumericVector weights(act.size());
    std::vector<size_t> final_sizes(act.size());
    double total_S = 0.0;
    double H_w = 0.0;
    for (size_t idx = 0; idx < act.size(); ++idx) {
        int k = act[idx];
        old_to_new[k] = static_cast<int>(idx);
        final_sizes[idx] = pool.modes[k].M();
        total_S += pool.modes[k].entropy();
        mode_list[idx] = pool.modes[k].to_r();
        double w = static_cast<double>(final_sizes[idx]) / static_cast<double>(M);
        weights[idx] = w;
        if (w > 0)
            H_w -= w * std::log(w);
    }
    total_S += cluster_prior(M, final_sizes);

    Rcpp::IntegerVector cl_out(M);
    for (size_t i = 0; i < M; ++i)
        cl_out[i] = old_to_new[pool.cl[i]] + 1;

    return Rcpp::List::create(
        Rcpp::Named("entropy") = total_S,
        Rcpp::Named("clusters") = cl_out,
        Rcpp::Named("K") = static_cast<int>(act.size()),
        Rcpp::Named("Be") = std::exp(H_w),
        Rcpp::Named("weights") = weights,
        Rcpp::Named("modes") = mode_list);
}
