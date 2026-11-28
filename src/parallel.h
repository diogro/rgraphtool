// rgraphtool: portable C++23 thread pool and thread-safe PRNG for parallel
// large-network inference. Works identically on macOS (Apple clang without
// libomp), Linux, and Windows with no external threading library dependency.
//
// Copyright (C) 2026 Diogo Melo
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.

#ifndef RGRAPHTOOL_PARALLEL_HPP
#define RGRAPHTOOL_PARALLEL_HPP

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace rgt_parallel {

// ---------------------------------------------------------------------------
// Fast deterministic thread-local PRNG (SplitMix64 seeder + xoshiro256**).
// Workers never touch R's RNG or R API; seeds are drawn on the main thread
// from R's RNG before spawning parallel regions so set.seed() is reproducible.

struct ThreadRng {
    uint64_t s[4];

    explicit ThreadRng(uint64_t seed = 0x9e3779b97f4a7c15ULL) {
        auto sm64 = [&seed]() {
            uint64_t z = (seed += 0x9e3779b97f4a7c15ULL);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        };
        s[0] = sm64();
        s[1] = sm64();
        s[2] = sm64();
        s[3] = sm64();
    }

    static inline uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    inline uint64_t next_u64() {
        uint64_t result = rotl(s[1] * 5, 7) * 9;
        uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    // High and low 64-bit halves of a * b, without relying on a 128-bit type.
    static inline uint64_t mul_hi_lo(uint64_t a, uint64_t b, uint64_t& lo) {
#if defined(__SIZEOF_INT128__)
        __extension__ typedef unsigned __int128 u128;
        u128 p = static_cast<u128>(a) * static_cast<u128>(b);
        lo = static_cast<uint64_t>(p);
        return static_cast<uint64_t>(p >> 64);
#else
        uint64_t a_lo = a & 0xffffffffULL, a_hi = a >> 32;
        uint64_t b_lo = b & 0xffffffffULL, b_hi = b >> 32;
        uint64_t ll = a_lo * b_lo, lh = a_lo * b_hi, hl = a_hi * b_lo, hh = a_hi * b_hi;
        uint64_t mid = (ll >> 32) + (lh & 0xffffffffULL) + (hl & 0xffffffffULL);
        lo = (mid << 32) | (ll & 0xffffffffULL);
        return hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
#endif
    }

    inline size_t index(size_t n) {
        if (n <= 1) return 0;
        // Lemire's multiply-shift range reduction, with the rejection step
        // that makes it exactly unbiased (taken with probability < n / 2^64).
        uint64_t n64 = static_cast<uint64_t>(n), lo;
        uint64_t hi = mul_hi_lo(next_u64(), n64, lo);
        if (lo < n64) {
            uint64_t threshold = (0 - n64) % n64;
            while (lo < threshold) hi = mul_hi_lo(next_u64(), n64, lo);
        }
        return static_cast<size_t>(hi);
    }

    inline double unif() {
        return static_cast<double>(next_u64() >> 11) * (1.0 / 9007199254740992.0);
    }

    template <class T>
    void shuffle(std::vector<T>& x) {
        for (size_t i = x.size(); i > 1; --i)
            std::swap(x[i - 1], x[index(i)]);
    }
};

// ---------------------------------------------------------------------------
// Global thread count configuration (defaults to 1 thread unless user opts in).

inline int hardware_threads() {
    unsigned int hc = std::thread::hardware_concurrency();
    return hc > 0 ? static_cast<int>(hc) : 1;
}

inline std::atomic<int>& global_threads_ref() {
    static std::atomic<int> nthr{1};
    return nthr;
}

inline int get_num_threads() {
    return global_threads_ref().load(std::memory_order_relaxed);
}

inline int set_num_threads(int n) {
    if (n < 1) n = 1;
    if (n > 256) n = 256;
    global_threads_ref().store(n, std::memory_order_relaxed);
    return n;
}

inline int resolve_threads(int nthreads) {
    if (nthreads <= 0) nthreads = get_num_threads();
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 256) nthreads = 256;
    return nthreads;
}

// True on a thread while it runs a parallel_for body (including the calling
// thread, which works as worker 0).
inline bool& in_parallel_region() {
    thread_local bool flag = false;
    return flag;
}

inline long current_pid() {
#ifdef _WIN32
    return static_cast<long>(_getpid());
#else
    return static_cast<long>(getpid());
#endif
}

// ---------------------------------------------------------------------------
// Lightweight persistent C++23 thread pool (calling thread participates as
// worker 0; workers are spawned lazily only when nthreads > 1).
//
// One pool per process. A forked child (parallel::mclapply) inherits a copy of
// the parent's pool object, but none of its worker threads: dispatching to it
// would wait forever for workers that do not exist. instance() therefore
// starts a fresh pool when the process id changes. The stale copy is leaked,
// never destroyed, since joining threads that do not exist is undefined.

class ThreadPool {
public:
    static ThreadPool& instance() {
        struct Holder {
            std::atomic<ThreadPool*> pool{nullptr};
            std::atomic<long> owner{-1};
            std::mutex mtx;
            // Joins the workers when the library is unloaded, so no thread is
            // left running code that dlclose() unmaps.
            ~Holder() {
                ThreadPool* p = pool.load();
                if (p && owner.load() == current_pid()) delete p;
            }
        };
        static Holder h;
        const long pid = current_pid();
        // owner is published after pool: seeing our pid guarantees seeing our pool
        if (h.owner.load(std::memory_order_acquire) == pid) {
            return *h.pool.load(std::memory_order_acquire);
        }
        std::lock_guard<std::mutex> lk(h.mtx);
        if (h.owner.load(std::memory_order_acquire) != pid) {
            h.pool.store(new ThreadPool(), std::memory_order_release);
            h.owner.store(pid, std::memory_order_release);
        }
        return *h.pool.load(std::memory_order_acquire);
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(_mtx);
            _stop = true;
            ++_gen;
        }
        _cv_work.notify_all();
        for (auto& w : _workers) {
            if (w.joinable()) w.join();
        }
    }

    template <class Fn>
    void parallel_for(size_t n, int nthreads, Fn&& fn, size_t min_chunk = 1) {
        int nw = resolve_threads(nthreads);
        // A parallel_for called from inside a parallel region runs serially on
        // the calling thread: the pool is busy, and its run lock is not
        // recursive, so dispatching would deadlock.
        if (in_parallel_region()) nw = 1;
        if (nw <= 1 || n <= min_chunk) {
            if (n > 0) fn(0, size_t(0), n);
            return;
        }
        size_t max_tasks = (n + min_chunk - 1) / min_chunk;
        if (static_cast<size_t>(nw) > max_tasks) {
            nw = static_cast<int>(max_tasks);
        }
        if (nw <= 1) {
            fn(0, size_t(0), n);
            return;
        }

        // Serialize concurrent parallel_for entries from different threads.
        std::unique_lock<std::mutex> run_lk(_run_mtx);
        ensure_workers(nw - 1);

        size_t chunk = std::max(min_chunk, n / (static_cast<size_t>(nw) * 8));
        std::atomic<size_t> cursor{0};
        std::exception_ptr exc = nullptr;
        std::mutex exc_mtx;

        auto worker_body = [&](int tid) {
            bool& in_region = in_parallel_region();
            bool was_in = in_region;
            in_region = true;
            try {
                while (true) {
                    size_t lo = cursor.fetch_add(chunk, std::memory_order_relaxed);
                    if (lo >= n) break;
                    size_t hi = std::min(lo + chunk, n);
                    fn(tid, lo, hi);
                }
            } catch (...) {
                std::lock_guard<std::mutex> elk(exc_mtx);
                if (!exc) exc = std::current_exception();
            }
            in_region = was_in;
        };

        {
            std::lock_guard<std::mutex> lk(_mtx);
            _active_workers = nw - 1;
            _remaining = nw - 1;
            _task = [&](int tid) { worker_body(tid); };
            ++_gen;
        }
        _cv_work.notify_all();

        // Main thread participates immediately as worker 0
        worker_body(0);

        {
            std::unique_lock<std::mutex> lk(_mtx);
            _cv_done.wait(lk, [&]() { return _remaining == 0; });
            _task = nullptr;
        }

        if (exc) std::rethrow_exception(exc);
    }

private:
    ThreadPool() = default;
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    void ensure_workers(int needed) {
        int cur = static_cast<int>(_workers.size());
        if (needed <= cur) return;
        for (int id = cur + 1; id <= needed; ++id) {
            _workers.emplace_back([this, id]() { worker_loop(id); });
        }
    }

    void worker_loop(int tid) {
        uint64_t seen_gen = 0;
        while (true) {
            std::function<void(int)> local_task;
            bool participate = false;
            {
                std::unique_lock<std::mutex> lk(_mtx);
                _cv_work.wait(lk, [&]() { return _stop || _gen != seen_gen; });
                if (_stop) return;
                seen_gen = _gen;
                if (tid <= _active_workers) {
                    participate = true;
                    local_task = _task;
                }
            }
            if (participate) {
                local_task(tid);
                local_task = nullptr;
                std::lock_guard<std::mutex> lk(_mtx);
                if (--_remaining == 0) {
                    _cv_done.notify_one();
                }
            }
        }
    }

    std::mutex _run_mtx;
    std::mutex _mtx;
    std::condition_variable _cv_work;
    std::condition_variable _cv_done;
    std::vector<std::thread> _workers;
    std::function<void(int)> _task;
    uint64_t _gen = 0;
    int _active_workers = 0;
    int _remaining = 0;
    bool _stop = false;
};

template <class Fn>
inline void parallel_for(size_t n, int nthreads, Fn&& fn, size_t min_chunk = 1) {
    ThreadPool::instance().parallel_for(n, nthreads, std::forward<Fn>(fn), min_chunk);
}

}  // namespace rgt_parallel

#endif  // RGRAPHTOOL_PARALLEL_HPP
