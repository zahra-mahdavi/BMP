#pragma once
#include <chrono>
#include <functional>

struct TimedResult {
    double seconds = 0.0;
};

template<class F>
inline TimedResult time_call(F&& f) {
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();
    f();
    auto t1 = Clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    return TimedResult{dt.count()};
}
