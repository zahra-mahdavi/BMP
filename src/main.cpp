#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <cassert>
#include "bitpack.hpp"
    #include "traversal.hpp"
#include "contractor.hpp"

using namespace std;

static vector<int> random_bits(size_t n){
    vector<int> b(n);
    for(size_t i=0;i<n;++i) b[i] = rand()&1;
    return b;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    size_t nvars = 16;        // number of matrices in chain
    size_t left_rows = 256;   // rows of leftmost matrix
    size_t right_cols = 128;  // size of R

    // Allow quick overrides
    for(int i=1;i<argc;i++){
        string s = argv[i];
        auto eat = [&](const string& key, size_t& dst){
            if(s.rfind(key,0)==0){
                dst = stoull(s.substr(key.size()));
                return true;
            }
            return false;
        };
        if(eat("--nvars=", nvars)) continue;
        if(eat("--left=", left_rows)) continue;
        if(eat("--right=", right_cols)) continue;
    }

    srand(42);

    // Build random row-switch chain
    vector<PackedMatrix> M0, M1;
    gen_random_rowswitch_chain(nvars, left_rows, right_cols, M0, M1);

    // Random R
    PackedVector R(right_cols);
    for(size_t j=0;j<right_cols;++j) if(rand()&1) R.set_bit(j);

    // Random input
    vector<int> x = random_bits(nvars);

    // Compute via traversal
    auto t0 = chrono::high_resolution_clock::now();
    auto out_trav = compute_via_traversal(M0, M1, R, x);
    auto t1 = chrono::high_resolution_clock::now();

    // Compute via CPU matvec (boolean semiring)
    auto out_cpu = compute_via_cpu_matvec(M0, M1, R, x);
    auto t2 = chrono::high_resolution_clock::now();

    // Compare
    bool match = (out_trav == out_cpu);
    chrono::duration<double> dt_trav = t1 - t0;
    chrono::duration<double> dt_cpu  = t2 - t1;

    cout << "Chain: nvars=" << nvars
         << " left_rows=" << left_rows
         << " right_cols=" << right_cols << "\n";
    cout << "Traversal: " << dt_trav.count() << " s\n";
    cout << "CPU matvec: " << dt_cpu.count() << " s\n";
    cout << "Validation: " << (match ? "PASS" : "FAIL") << "\n";

    // Print a few bits for sanity
    cout << "out[0..15]: ";
    for(size_t i=0;i<min<size_t>(16,out_trav.size());++i) cout << int(out_trav[i]);
    cout << "\n";

    return 0;
}
