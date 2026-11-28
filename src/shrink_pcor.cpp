// Shrinkage partial correlations (Schafer & Strimmer 2005) without forming
// N x N matrices unless requested.
//
// Let A (N x M) hold the standardized variables scaled so that R = A A^T is
// the sample correlation matrix. The shrunk correlation is
//   R* = (1 - lambda) R + lambda I.
// With c = lambda / (1 - lambda) and C = c I_M + A^T A = L^T L (Cholesky,
// upper L), Woodbury gives R*^{-1} proportional to I - U U^T with U = A L^{-1},
// so the partial correlations are
//   pcor_ij = v_i . v_j,   v_i = u_i / sqrt(1 - |u_i|^2)   (i != j).
// The Schafer-Strimmer shrinkage intensity only needs A^T A (M x M) and
// per-observation sums, so it costs O(N M^2) time and O(N M + M^2) memory.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include <algorithm>
#include <cmath>
#include <vector>

#include "parallel.h"

namespace {

struct Standardized {
    std::vector<double> A;  // N x M, column-major
    int N = 0;
    int M = 0;
};

Standardized standardize(const Rcpp::NumericMatrix& x, bool center) {
    Standardized s;
    s.N = x.nrow();
    s.M = x.ncol();
    const int N = s.N, M = s.M;
    s.A.assign(static_cast<size_t>(N) * M, 0.0);
    for (int i = 0; i < N; ++i) {
        double mu = 0.0;
        if (center) {
            for (int k = 0; k < M; ++k) mu += x(i, k);
            mu /= M;
        }
        double ss = 0.0;
        for (int k = 0; k < M; ++k) {
            const double d = x(i, k) - mu;
            ss += d * d;
        }
        if (!(ss > 0.0) || !std::isfinite(ss)) {
            Rcpp::stop("variable %d has zero or non-finite variance", i + 1);
        }
        // A A^T = correlation: a_i = (x_i - mu_i) / sqrt(ss_i)
        const double inv = 1.0 / std::sqrt(ss);
        for (int k = 0; k < M; ++k) s.A[i + static_cast<size_t>(k) * N] = (x(i, k) - mu) * inv;
    }
    return s;
}

// G = A^T A (M x M, full symmetric)
std::vector<double> gram_mm(const Standardized& s) {
    const int N = s.N, M = s.M;
    std::vector<double> G(static_cast<size_t>(M) * M, 0.0);
    const double one = 1.0, zero = 0.0;
    F77_CALL(dsyrk)("U", "T", &M, &N, &one, s.A.data(), &N, &zero, G.data(), &M FCONE FCONE);
    for (int j = 0; j < M; ++j)
        for (int i = j + 1; i < M; ++i) G[i + static_cast<size_t>(j) * M] = G[j + static_cast<size_t>(i) * M];
    return G;
}

double estimate_lambda(const Standardized& s, const std::vector<double>& G, bool center) {
    const int N = s.N, M = s.M;
    if (N < 2) return 0.0;
    const double n = static_cast<double>(M);
    const double df = center ? n - 1.0 : n;
    // Standardized data xs = A * sqrt(df) (unit variance with df denominator)
    double s4diag = 0.0, sum_sk2 = 0.0;
    for (int k = 0; k < M; ++k) {
        double sk = 0.0;
        for (int i = 0; i < N; ++i) {
            const double a = s.A[i + static_cast<size_t>(k) * N];
            const double xs2 = a * a * df;
            sk += xs2;
            s4diag += xs2 * xs2;
        }
        sum_sk2 += sk * sk;
    }
    const double s4off = sum_sk2 - s4diag;  // sum_{i!=j} sum_k w_kij^2
    double g2 = 0.0;
    for (double g : G) g2 += g * g;          // sum_{ij} r_ij^2
    const double r2off = g2 - N;             // sum_{i!=j} r_ij^2
    if (!(r2off > 0.0)) return 1.0;
    double var_sum;
    if (center) {
        // w_bar = (n - 1)/n r ; Var(r_ij) = n/(n-1)^3 sum_k (w - w_bar)^2
        const double wbar2 = ((n - 1.0) / n) * ((n - 1.0) / n) * r2off;
        var_sum = n / ((n - 1.0) * (n - 1.0) * (n - 1.0)) * (s4off - n * wbar2);
    } else {
        // known zero mean: w_bar = r ; Var(r_ij) = 1/(n (n-1)) sum_k (w - w_bar)^2
        var_sum = (s4off - n * r2off) / (n * (n - 1.0));
    }
    double lam = var_sum / r2off;
    if (!std::isfinite(lam)) lam = 1.0;
    return std::min(1.0, std::max(0.0, lam));
}

// Returns Vt (M x N, column-major: gene i is column i) with pcor_ij = Vt_i . Vt_j
std::vector<double> pcor_factor(const Standardized& s, const std::vector<double>& G, double lam) {
    const int N = s.N, M = s.M;
    std::vector<double> C = G;
    const double c = lam / (1.0 - lam);
    for (int k = 0; k < M; ++k) C[k + static_cast<size_t>(k) * M] += c;
    int info = 0;
    F77_CALL(dpotrf)("U", &M, C.data(), &M, &info FCONE);
    if (info != 0) Rcpp::stop("Cholesky factorization failed (info = %d)", info);
    std::vector<double> U = s.A;  // N x M; solve U L = A
    const double one = 1.0;
    F77_CALL(dtrsm)("R", "U", "N", "N", &N, &M, &one, C.data(), &M, U.data(), &N
                    FCONE FCONE FCONE FCONE);
    std::vector<double> Vt(static_cast<size_t>(M) * N);
    for (int i = 0; i < N; ++i) {
        double r = 0.0;
        for (int k = 0; k < M; ++k) {
            const double u = U[i + static_cast<size_t>(k) * N];
            r += u * u;
        }
        const double denom = std::sqrt(std::max(1.0 - r, 1e-300));
        for (int k = 0; k < M; ++k) Vt[k + static_cast<size_t>(i) * M] = U[i + static_cast<size_t>(k) * N] / denom;
    }
    return Vt;
}

// Direct path for lambda == 0 (requires an invertible correlation matrix).
Rcpp::NumericMatrix pcor_direct(const Standardized& s) {
    const int N = s.N, M = s.M;
    std::vector<double> R(static_cast<size_t>(N) * N, 0.0);
    const double one = 1.0, zero = 0.0;
    F77_CALL(dsyrk)("U", "N", &N, &M, &one, s.A.data(), &N, &zero, R.data(), &N FCONE FCONE);
    int info = 0;
    F77_CALL(dpotrf)("U", &N, R.data(), &N, &info FCONE);
    if (info != 0) {
        Rcpp::stop("lambda = 0 but the correlation matrix is singular (N >= M?); use lambda > 0");
    }
    F77_CALL(dpotri)("U", &N, R.data(), &N, &info FCONE);
    if (info != 0) Rcpp::stop("matrix inversion failed (info = %d)", info);
    Rcpp::NumericMatrix P(N, N);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i <= j; ++i) {
            const double pij = R[i + static_cast<size_t>(j) * N];
            const double pii = R[i + static_cast<size_t>(i) * N];
            const double pjj = R[j + static_cast<size_t>(j) * N];
            const double v = (i == j) ? 0.0 : -pij / std::sqrt(pii * pjj);
            P(i, j) = v;
            P(j, i) = v;
        }
    }
    return P;
}

struct Cand {
    int i, j;
    double v;
};

// Dot product with four independent accumulators (vectorizes without
// -ffast-math).
inline double dot4(const double* a, const double* b, int M) {
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int k = 0;
    for (; k + 3 < M; k += 4) {
        s0 += a[k] * b[k];
        s1 += a[k + 1] * b[k + 1];
        s2 += a[k + 2] * b[k + 2];
        s3 += a[k + 3] * b[k + 3];
    }
    for (; k < M; ++k) s0 += a[k] * b[k];
    return (s0 + s1) + (s2 + s3);
}

}  // namespace

// [[Rcpp::export(.shrink_pcor_cpp)]]
Rcpp::List shrink_pcor_cpp(const Rcpp::NumericMatrix& x, double lambda, bool center,
                           const Rcpp::IntegerMatrix& pairs, int top_k, bool want_matrix,
                           bool lambda_only, int nthreads) {
    const int N = x.nrow(), M = x.ncol();
    if (M < 3) Rcpp::stop("x must have at least 3 observations (columns)");
    Standardized s = standardize(x, center);
    std::vector<double> G = gram_mm(s);
    double lam = lambda;
    if (Rcpp::NumericVector::is_na(lam)) lam = estimate_lambda(s, G, center);
    if (lam < 0.0 || lam > 1.0) Rcpp::stop("lambda must be in [0, 1]");

    Rcpp::List out = Rcpp::List::create(Rcpp::Named("lambda") = lam, Rcpp::Named("n_obs") = M);
    if (lambda_only) return out;

    const bool all_zero = lam >= 1.0;
    const bool direct = lam <= 0.0;
    Rcpp::NumericMatrix Pdir;
    std::vector<double> Vt;
    if (direct) {
        Pdir = pcor_direct(s);
    } else if (!all_zero) {
        Vt = pcor_factor(s, G, lam);
    }
    auto pc = [&](int i, int j) -> double {
        if (all_zero || i == j) return 0.0;
        if (direct) return Pdir(i, j);
        return dot4(Vt.data() + static_cast<size_t>(i) * M,
                    Vt.data() + static_cast<size_t>(j) * M, M);
    };

    if (want_matrix) {
        if (direct) {
            out["matrix"] = Pdir;
        } else {
            Rcpp::NumericMatrix P(N, N);
            if (!all_zero) {
                std::vector<double> S(static_cast<size_t>(N) * N, 0.0);
                const double one = 1.0, zero = 0.0;
                F77_CALL(dsyrk)("U", "T", &N, &M, &one, Vt.data(), &M, &zero, S.data(), &N
                                FCONE FCONE);
                for (int j = 0; j < N; ++j)
                    for (int i = 0; i < j; ++i) {
                        const double v = S[i + static_cast<size_t>(j) * N];
                        P(i, j) = v;
                        P(j, i) = v;
                    }
            }
            out["matrix"] = P;
        }
    }

    if (pairs.nrow() > 0) {
        if (pairs.ncol() != 2) Rcpp::stop("pairs must be a two-column matrix");
        Rcpp::NumericVector pv(pairs.nrow());
        for (int r = 0; r < pairs.nrow(); ++r) {
            const int i = pairs(r, 0) - 1, j = pairs(r, 1) - 1;
            if (i < 0 || j < 0 || i >= N || j >= N) Rcpp::stop("pairs index out of range");
            pv[r] = pc(i, j);
        }
        out["pairs"] = pv;
    }

    if (top_k > 0) {
        const long long n_pairs = static_cast<long long>(N) * (N - 1) / 2;
        const int k = static_cast<int>(std::min<long long>(top_k, n_pairs));
        std::vector<Cand> cand;
        cand.reserve(static_cast<size_t>(2) * k + 16);
        double thr = 0.0;
        auto keep_top = [&]() {
            if (static_cast<int>(cand.size()) <= k) return;
            std::nth_element(cand.begin(), cand.begin() + (k - 1), cand.end(),
                             [](const Cand& a, const Cand& b) { return std::fabs(a.v) > std::fabs(b.v); });
            cand.resize(k);
            thr = std::fabs(cand[k - 1].v);
        };
        const int blk = 256;
        const int n_blocks = (N + blk - 1) / blk;
        if (direct || all_zero) {
            for (int i = 0; i < N; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    const double v = pc(i, j);
                    if (std::fabs(v) > thr || static_cast<int>(cand.size()) < k) cand.push_back({i, j, v});
                }
                if (static_cast<int>(cand.size()) > 2 * k + 1024) keep_top();
            }
        } else {
            // Row blocks are scanned in parallel; each worker keeps its own
            // candidate buffer and running threshold (no R API on workers).
            const int nw = rgt_parallel::resolve_threads(nthreads);
            std::vector<std::vector<Cand>> local(std::max(1, nw));
            std::vector<double> local_thr(std::max(1, nw), 0.0);
            auto keep_local = [&](int tid) {
                auto& cv = local[tid];
                if (static_cast<int>(cv.size()) <= k) return;
                std::nth_element(cv.begin(), cv.begin() + (k - 1), cv.end(),
                                 [](const Cand& a, const Cand& b) { return std::fabs(a.v) > std::fabs(b.v); });
                cv.resize(k);
                local_thr[tid] = std::fabs(cv[k - 1].v);
            };
            rgt_parallel::parallel_for(
                static_cast<size_t>(n_blocks), nw,
                [&](int tid, size_t lo, size_t hi) {
                    // Plain C++ dot products, never BLAS: BLAS must not be
                    // called from pool workers (a multithreaded OpenBLAS would
                    // oversubscribe cores, and some builds are not reentrant).
                    // One kernel for every thread count keeps results
                    // bit-identical across nthreads.
                    auto& cv = local[tid];
                    for (size_t bi = lo; bi < hi; ++bi) {
                        const int b0 = static_cast<int>(bi) * blk;
                        const int nb = std::min(blk, N - b0);
                        for (int ii = 0; ii < nb; ++ii) {
                            const int i = b0 + ii;
                            const double* a = Vt.data() + static_cast<size_t>(i) * M;
                            for (int j = i + 1; j < N; ++j) {
                                const double v = dot4(a, Vt.data() + static_cast<size_t>(j) * M, M);
                                if (std::fabs(v) > local_thr[tid] || static_cast<int>(cv.size()) < k)
                                    cv.push_back({i, j, v});
                            }
                            if (static_cast<int>(cv.size()) > 2 * k + 1024) keep_local(tid);
                        }
                        keep_local(tid);
                    }
                },
                1);
            for (auto& cv : local) {
                cand.insert(cand.end(), cv.begin(), cv.end());
                std::vector<Cand>().swap(cv);
            }
        }
        keep_top();
        std::sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b) {
            const double fa = std::fabs(a.v), fb = std::fabs(b.v);
            if (fa != fb) return fa > fb;
            if (a.i != b.i) return a.i < b.i;
            return a.j < b.j;
        });
        Rcpp::IntegerVector from(cand.size()), to(cand.size());
        Rcpp::NumericVector val(cand.size());
        for (size_t r = 0; r < cand.size(); ++r) {
            from[r] = cand[r].i + 1;
            to[r] = cand[r].j + 1;
            val[r] = cand[r].v;
        }
        out["top"] = Rcpp::DataFrame::create(Rcpp::Named("from") = from, Rcpp::Named("to") = to,
                                             Rcpp::Named("pcor") = val);
    }
    return out;
}
