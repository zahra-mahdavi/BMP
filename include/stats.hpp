#pragma once
#include <vector>
#include <cmath>
#include <cstddef>

struct Stats {
    double mean = 0.0;
    double stddev = 0.0;
    double se = 0.0;
    std::size_t n = 0;
};

inline Stats compute_stats(const std::vector<double>& xs) {
    Stats s; 
    s.n = xs.size();
    if (s.n == 0) return s;
    long double sum = 0.0L, sumsq = 0.0L;
    for (double v : xs) { sum += v; sumsq += v * v; }
    s.mean = static_cast<double>(sum / s.n);
    long double var = 0.0L;
    if (s.n > 1) {
        var = (sumsq - sum * sum / s.n) / (s.n - 1); // unbiased
    }
    s.stddev = var > 0.0L ? std::sqrt(static_cast<double>(var)) : 0.0;
    // PDF: sigma(avg) = sigma / sqrt(M - 1)
    s.se = (s.n > 1) ? s.stddev / std::sqrt(static_cast<double>(s.n - 1)) : 0.0;
    return s;
}
