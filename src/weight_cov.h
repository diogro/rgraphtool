// Shared edge covariate likelihood and specification helpers for weighted
// and combined layered-weighted SBMs (Peixoto 2018, Phys. Rev. E 97, 012306).
// LGPL-3.0-or-later.

#ifndef RGRAPHTOOL_WEIGHT_COV_H
#define RGRAPHTOOL_WEIGHT_COV_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace rgt_cov {

constexpr double W_INF = std::numeric_limits<double>::infinity();

enum class CovType : int {
    REAL_EXPONENTIAL = 0,
    REAL_NORMAL = 1,
    DISCRETE_GEOMETRIC = 2,
    DISCRETE_POISSON = 3,
    DISCRETE_BINOMIAL = 4
};

struct CovSpec {
    CovType type;
    bool micro;              // true = microcanonical (NaN hyperparams), false = Bayesian
    double p0 = NAN;         // alpha (exp/geo/poi), m0 (normal), N (binomial)
    double p1 = NAN;         // beta  (exp/geo/poi), k0 (normal), alpha (binomial)
    double p2 = NAN;         // v0 (normal), beta (binomial)
    double p3 = NAN;         // nu0 (normal)
    double eps = 1e-8;       // 0.1 * min_{w != 0} |w|
    double const_term = 0.0; // constant edge-level factor (poisson lgamma(w+1), binomial -lbinom(N,w))
};

// Compact covariate accumulator per block pair, layer entry, or delta entry.
// Inline storage for up to 4 covariates avoids heap allocations in hot loops.
struct WVec {
    static constexpr int INLINE_K = 4;
    int k = 0;
    double s_x[INLINE_K] = {0, 0, 0, 0};
    double s_x2[INLINE_K] = {0, 0, 0, 0};
    std::vector<double> d_x, d_x2;

    WVec() = default;
    explicit WVec(int K) : k(K) {
        if (K > INLINE_K) {
            d_x.assign(K, 0.0);
            d_x2.assign(K, 0.0);
        }
    }
    void zero(int K) {
        k = K;
        if (K <= INLINE_K) {
            for (int i = 0; i < K; ++i) { s_x[i] = 0.0; s_x2[i] = 0.0; }
        } else {
            d_x.assign(K, 0.0);
            d_x2.assign(K, 0.0);
        }
    }
    inline double& x(int i) { return (k <= INLINE_K) ? s_x[i] : d_x[i]; }
    inline double x(int i) const { return (k <= INLINE_K) ? s_x[i] : d_x[i]; }
    inline double& x2(int i) { return (k <= INLINE_K) ? s_x2[i] : d_x2[i]; }
    inline double x2(int i) const { return (k <= INLINE_K) ? s_x2[i] : d_x2[i]; }
};

inline double w_safe_lgamma(double x) {
    if (x <= 0.0) return 0.0;
    return std::lgamma(x);
}

inline double w_lbinom_d(double n, double k) {
    if (k <= 0.0 || k >= n) return 0.0;
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

inline double w_lbeta(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

inline double w_safelog(double x) {
    return (x > 0.0) ? std::log(x) : 0.0;
}

// Exact graph-tool 3.9 weighted log-probabilities (from weighted/entropy.hh)
inline double positive_w_log_P(double N, double x, double alpha, double beta,
                               double eps) {
    if (N <= 0.0) return 0.0;
    if (std::isnan(alpha) && std::isnan(beta)) {
        if (x < eps || N <= 1.0) return 0.0;
        return w_safe_lgamma(N) - (N - 1.0) * std::log(x);
    }
    return (w_safe_lgamma(N + alpha) - w_safe_lgamma(alpha) +
            alpha * std::log(beta) - (alpha + N) * std::log(beta + x));
}

inline double signed_w_log_P(double N, double x, double x2, double m0,
                             double k0, double v0, double nu0, double eps) {
    if (N <= 0.0) return 0.0;
    if (std::isnan(m0) && std::isnan(k0)) {
        double smu1 = x * (x / N);
        if (N < 2.0 || smu1 >= x2 || (x2 - smu1) < eps * eps)
            return 0.0;
        return (w_safe_lgamma((N - 1.0) / 2.0) + 0.5 * std::log(N) -
                ((static_cast<int>(N) - 3) / 2.0) * std::log(x2 - smu1) -
                ((N - 1.0) / 2.0) * std::log(M_PI));
    }
    double v = std::max(0.0, x2 - x * (x / N));
    double k_n = k0 + N;
    double nu_n = nu0 + N;
    double diff = m0 - x / N;
    double nu_v_n = v0 * nu0 + v + ((N * k0) / k_n) * diff * diff;
    return (w_safe_lgamma(nu_n / 2.0) - w_safe_lgamma(nu0 / 2.0) +
            0.5 * (std::log(k0) - std::log(k_n)) +
            0.5 * nu0 * std::log(nu0 * v0) - 0.5 * nu_n * std::log(nu_v_n) -
            0.5 * N * std::log(M_PI));
}

inline double geometric_w_log_P(double N, double x, double alpha, double beta) {
    if (N <= 0.0) return 0.0;
    if (std::isnan(alpha) && std::isnan(beta))
        return -w_lbinom_d((N - 1.0) + x, x);
    return w_lbeta(N + alpha, x + beta) - w_lbeta(alpha, beta);
}

inline double poisson_w_log_P(double N, double x, double alpha, double beta) {
    if (N <= 0.0) return 0.0;
    if (std::isnan(alpha) && std::isnan(beta))
        return w_safe_lgamma(x + 1.0) - x * std::log(N);
    return (w_safe_lgamma(x + alpha) - (x + alpha) * std::log(N + beta) -
            w_safe_lgamma(alpha) + alpha * std::log(beta));
}

inline double binomial_w_log_P(double N, double x, double trials,
                               double alpha, double beta) {
    if (N <= 0.0) return 0.0;
    if (std::isnan(alpha) && std::isnan(beta))
        return -w_lbinom_d(N * trials, x);
    return w_lbeta(x + alpha, N * trials - x + beta) - w_lbeta(alpha, beta);
}

inline double cov_pair_neg_log_P(const CovSpec& c, double m, double x, double x2) {
    if (m <= 0.0) return 0.0;
    switch (c.type) {
    case CovType::REAL_EXPONENTIAL:
        return -positive_w_log_P(m, x, c.p0, c.p1, c.eps);
    case CovType::REAL_NORMAL:
        return -signed_w_log_P(m, x, x2, c.p0, c.p1, c.p2, c.p3, c.eps);
    case CovType::DISCRETE_GEOMETRIC:
        return -geometric_w_log_P(m, x, c.p0, c.p1);
    case CovType::DISCRETE_POISSON:
        return -poisson_w_log_P(m, x, c.p0, c.p1);
    case CovType::DISCRETE_BINOMIAL:
        return -binomial_w_log_P(m, x, c.p0, c.p1, c.p2);
    }
    return 0.0;
}

inline std::vector<CovSpec> parse_covs(const Rcpp::NumericMatrix& weights,
                                       const Rcpp::CharacterVector& types,
                                       const Rcpp::List& params) {
    int K = static_cast<int>(types.size());
    if (K == 0) return {};
    int E = weights.nrow();
    if (weights.ncol() != K || params.size() != K)
        Rcpp::stop("weights, types, and params must agree on the number of covariates");
    std::vector<CovSpec> out(K);
    for (int k = 0; k < K; ++k) {
        std::string t = Rcpp::as<std::string>(types[k]);
        CovSpec& c = out[k];
        if (t == "real-exponential") c.type = CovType::REAL_EXPONENTIAL;
        else if (t == "real-normal") c.type = CovType::REAL_NORMAL;
        else if (t == "discrete-geometric") c.type = CovType::DISCRETE_GEOMETRIC;
        else if (t == "discrete-poisson") c.type = CovType::DISCRETE_POISSON;
        else if (t == "discrete-binomial") c.type = CovType::DISCRETE_BINOMIAL;
        else Rcpp::stop("Unknown covariate type: '%s'", t.c_str());

        double min_nz = W_INF;
        for (int e = 0; e < E; ++e) {
            double w = weights(e, k);
            if (!std::isfinite(w)) Rcpp::stop("Edge weights must be finite");
            if (w != 0.0) min_nz = std::min(min_nz, std::abs(w));
        }
        c.eps = std::isfinite(min_nz) ? (min_nz / 10.0) : 1e-8;

        Rcpp::NumericVector p = params[k];
        if (c.type == CovType::REAL_EXPONENTIAL ||
            c.type == CovType::DISCRETE_GEOMETRIC ||
            c.type == CovType::DISCRETE_POISSON) {
            c.p0 = (p.size() > 0) ? p[0] : NAN;
            c.p1 = (p.size() > 1) ? p[1] : NAN;
            c.micro = (std::isnan(c.p0) && std::isnan(c.p1));
        } else if (c.type == CovType::DISCRETE_BINOMIAL) {
            c.p0 = (p.size() > 0) ? p[0] : 1.0;
            c.p1 = (p.size() > 1) ? p[1] : NAN;
            c.p2 = (p.size() > 2) ? p[2] : NAN;
            c.micro = (std::isnan(c.p1) && std::isnan(c.p2));
        } else if (c.type == CovType::REAL_NORMAL) {
            c.p0 = (p.size() > 0) ? p[0] : NAN;
            c.p1 = (p.size() > 1) ? p[1] : NAN;
            c.p2 = (p.size() > 2) ? p[2] : NAN;
            c.p3 = (p.size() > 3) ? p[3] : NAN;
            c.micro = (std::isnan(c.p0) && std::isnan(c.p1));
        }

        c.const_term = 0.0;
        if (c.type == CovType::DISCRETE_POISSON) {
            for (int e = 0; e < E; ++e)
                c.const_term += w_safe_lgamma(weights(e, k) + 1.0);
        } else if (c.type == CovType::DISCRETE_BINOMIAL) {
            for (int e = 0; e < E; ++e)
                c.const_term -= w_lbinom_d(c.p0, weights(e, k));
        }
    }
    return out;
}

inline std::vector<CovSpec> make_upper_covs(const std::vector<CovSpec>& base_covs) {
    std::vector<CovSpec> up;
    up.reserve(base_covs.size());
    for (const auto& c : base_covs) {
        CovSpec u;
        u.micro = true;
        u.eps = c.eps;
        u.const_term = 0.0;
        if (c.type == CovType::DISCRETE_GEOMETRIC ||
            c.type == CovType::DISCRETE_POISSON ||
            c.type == CovType::DISCRETE_BINOMIAL) {
            u.type = CovType::DISCRETE_GEOMETRIC;
        } else {
            u.type = c.type;
        }
        up.push_back(u);
    }
    return up;
}

}  // namespace rgt_cov

#endif  // RGRAPHTOOL_WEIGHT_COV_H
