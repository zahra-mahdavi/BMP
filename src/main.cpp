// src/main.cpp
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <cassert>
#include <algorithm>

#include "bitpack.hpp"
#include "traversal.hpp"
#include "contractor.hpp"

using Clock = std::chrono::steady_clock;

static std::vector<int> random_bits(size_t n, std::mt19937_64& rng){
    std::vector<int> b(n);
    std::uniform_int_distribution<int> bit01(0, 1);
    for(size_t i=0;i<n;++i) b[i] = bit01(rng);
    return b;
}

int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Defaults
    size_t nvars      = 16;   // number of matrices in chain
    size_t left_rows  = 256;  // rows of leftmost matrix
    size_t right_cols = 128;  // size of R 

    
    for(int i=1;i<argc;i++){
        std::string s = argv[i];
        auto eat = [&](const std::string& key, size_t& dst){
            if(s.rfind(key,0)==0){
                dst = static_cast<size_t>(std::stoull(s.substr(key.size())));
                return true;
            }
            return false;
        };
        if(eat("--nvars=", nvars))      continue;
        if(eat("--left=", left_rows))   continue;
        if(eat("--right=", right_cols)) continue;
    }

    // RNG 
    std::mt19937_64 rng(42);

    // Build a random row-switch chain 
    std::vector<PackedMatrix> M0, M1;
    gen_random_rowswitch_chain(nvars, left_rows, right_cols, M0, M1);

    
    PackedVector R(right_cols);
    {
        std::uniform_int_distribution<int> bit01(0, 1);
        for(size_t j=0;j<right_cols;++j) if(bit01(rng)) R.set_bit(j);
    }

    // Random input bits x
    std::vector<int> x = random_bits(nvars, rng);

    
    (void)compute_via_traversal(M0, M1, R, x);
    (void)compute_via_cpu_matvec(M0, M1, R, x);

    // Timed runs 
    auto t0 = Clock::now();
    auto out_trav = compute_via_traversal(M0, M1, R, x);
    auto t1 = Clock::now();

    auto out_cpu  = compute_via_cpu_matvec(M0, M1, R, x);
    auto t2 = Clock::now();

    bool match = (out_trav == out_cpu);
    std::chrono::duration<double> dt_trav = t1 - t0;
    std::chrono::duration<double> dt_cpu  = t2 - t1;

    // Report 
    std::cout << "Chain: nvars=" << nvars
              << " left_rows=" << left_rows
              << " right_cols=" << right_cols << "\n";
    std::cout << "Traversal time: " << dt_trav.count() << " s\n";
    std::cout << "CPU matvec time: " << dt_cpu.count() << " s\n";
    std::cout << "Validation: " << (match ? "PASS" : "FAIL") << "\n";

    // Print a few bits for sanity
    std::cout << "out[0..15]: ";
    for(size_t i=0;i<std::min<size_t>(16, out_trav.size()); ++i)
        std::cout << int(out_trav[i]);
    std::cout << "\n";

    return match ? 0 : 1;
}
