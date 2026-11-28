// rgraphtool: shared log restricted-partition table q(n, k) for degree
// description lengths (graph-tool's q cache and log_q_approx).
//
// Copyright (C) 2026 Diogo Melo
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// Every SBM family uses this single implementation: log q(n, k), the log of
// the number of partitions of n into at most k parts, is tabulated exactly up
// to n_max and approximated beyond it with graph-tool's log_q_approx (the
// Szekeres saddle-point formula, with the small-k limit for k < n^(1/4)).
// Separate per-model copies used to cap the exact table at 2000 and fall back
// to an incorrect approximation, giving degree DLs hundreds of nats off.

#ifndef RGRAPHTOOL_LOGQ_H
#define RGRAPHTOOL_LOGQ_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

namespace rgt_logq {

class LogQ {
public:
    // Grow the exact table to cover n <= n_max (never shrinks). Not thread
    // safe: call from the main thread, before parallel regions.
    void ensure(size_t n_max) {
        if (n_max <= _n_max && !_table.empty())
            return;
        _n_max = n_max;
        _table.assign(n_max * (n_max + 1) / 2 + 1, -inf());
        for (size_t n = 1; n <= n_max; ++n) {
            for (size_t k = 1; k <= n; ++k) {
                double a = (k == 1) ? -inf() : get(n, k - 1);
                double b = get_raw(n - k, k);
                double hi = std::max(a, b), lo = std::min(a, b);
                at(n, k) = (hi == -inf()) ? -inf() : hi + std::log1p(std::exp(lo - hi));
            }
        }
    }
    void init(size_t n_max) { ensure(n_max); }
    size_t n_max() const { return _n_max; }

    // log q(n, k); q(0, k) = 1. Thread safe for concurrent reads.
    double operator()(long n, long k) const {
        if (k > n)
            k = n;
        if (n == 0 && k == 0)
            return 0;
        if (n < 0 || k <= 0)
            return -inf();
        if (static_cast<size_t>(n) <= _n_max)
            return get(n, k);
        thread_local std::unordered_map<uint64_t, double> approx_cache;
        uint64_t key = (static_cast<uint64_t>(n) << 32) | (static_cast<uint64_t>(k) & 0xffffffffULL);
        auto it = approx_cache.find(key);
        if (it != approx_cache.end())
            return it->second;
        double val = approx(static_cast<double>(n), static_cast<double>(k));
        if (approx_cache.size() >= 65536)
            approx_cache.clear();
        approx_cache.emplace(key, val);
        return val;
    }

    // graph-tool's log_q_approx.
    static double approx(double n, double k) {
        if (k < std::pow(n, 0.25))
            return lbinom(n - 1, k - 1) - std::lgamma(k + 1);
        double u = k / std::sqrt(n);
        double v = u, delta = 1;
        while (delta > 1e-8) {
            double n_v = u * std::sqrt(spence(std::exp(-v)));
            delta = std::abs(n_v - v);
            v = n_v;
        }
        double lf = std::log(v) - std::log1p(-std::exp(-v) * (1 + u * u / 2)) / 2 -
                    std::log(2.0) * 1.5 - std::log(u) - std::log(M_PI);
        double g = 2 * v / u - u * std::log1p(-std::exp(-v));
        return lf - std::log(n) + std::sqrt(n) * g;
    }

private:
    size_t _n_max = 0;
    std::vector<double> _table;

    static double inf() { return std::numeric_limits<double>::infinity(); }
    static double lbinom(double n, double k) {
        if (k == 0 || n == k || n < k || k < 0)
            return 0;
        return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
    }
    // spence(x) = Li2(1 - x) for 0 <= x <= 1 (scipy/cephes convention)
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

// The one table shared by every model (an inline function's static is a
// single object across translation units).
inline LogQ& shared() {
    static LogQ cache;
    return cache;
}

// Exact-table size for a graph with E edges: graph-tool's
// min(max(2E, 100), 10000), further capped by the requested q_max.
inline size_t table_size(double E, size_t q_cap) {
    size_t n = static_cast<size_t>(std::min(std::max(2.0 * E, 100.0), 10000.0));
    return std::min(n, q_cap);
}

}  // namespace rgt_logq

#endif
