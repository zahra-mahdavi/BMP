#pragma once
#include <vector>
#include <cstdint>
#include <cassert>
#include <random>
#include "bitpack.hpp"
#include "traversal.hpp"


inline LocVec packed_to_loc(const PackedMatrix& A){
    LocVec loc(A.rows);
    for(size_t r=0;r<A.rows;++r){
        loc[r] = (uint32_t)A.row_switch_col(r);
    }
    return loc;
}

inline std::vector<uint8_t> packedvec_to_bits(const PackedVector& v){
    std::vector<uint8_t> bits(v.nbits, 0);
    for(size_t i=0;i<v.nbits;++i){
        bits[i] = v.get_bit(i) ? 1u : 0u;
    }
    return bits;
}


inline TraversalBMP build_traversal(const std::vector<PackedMatrix>& M0,
                                    const std::vector<PackedMatrix>& M1,
                                    const PackedVector& R){
    assert(M0.size() == M1.size());
    TraversalBMP T;
    T.M0_loc.reserve(M0.size());
    T.M1_loc.reserve(M1.size());
    for(size_t i=0;i<M0.size();++i){
        T.M0_loc.push_back(packed_to_loc(M0[i]));
        T.M1_loc.push_back(packed_to_loc(M1[i]));
    }
    T.R_bits = packedvec_to_bits(R);
    T.n_left_rows = (uint32_t)M0.front().rows;
    return T;
}


inline std::vector<uint8_t> compute_via_traversal(const std::vector<PackedMatrix>& M0,
                                                  const std::vector<PackedMatrix>& M1,
                                                  const PackedVector& R,
                                                  const std::vector<int>& x_bits){
    TraversalBMP T = build_traversal(M0,M1,R);
    return traverse_all(T, x_bits);
}

// Compute output via chained matrix-vector on CPU.
// v0 = R; for i = nvars-1..0: v = A_i^{x_i} * v
inline std::vector<uint8_t> compute_via_cpu_matvec(const std::vector<PackedMatrix>& M0,
                                                   const std::vector<PackedMatrix>& M1,
                                                   const PackedVector& R,
                                                   const std::vector<int>& x_bits){
    PackedVector v = R;
    select_and_apply(M0, M1, x_bits, v);
    // Convert to 0/1 bits (length = rows of leftmost matrix)
    auto out = packedvec_to_bits(v);
    return out;
}

// generate a consistent random row-switch BMP chain with varying dims.
inline void gen_random_rowswitch_chain(size_t nvars,
                                       size_t left_rows,
                                       size_t right_cols,
                                       std::vector<PackedMatrix>& M0,
                                       std::vector<PackedMatrix>& M1){
    // Random intermediate widths, ensuring chain consistency:
    std::vector<size_t> widths(nvars+1);
    widths[0] = left_rows;
    widths[nvars] = right_cols;
    for(size_t i=1;i<nvars;++i){
        widths[i] = std::max<size_t>(1, (widths[i-1] + widths[i-1]/3)); // simple growth
    }
    
    M0.resize(nvars);
    M1.resize(nvars);
    for(size_t i=0;i<nvars;++i){
        M0[i] = random_row_switch_matrix(widths[i], widths[i+1]);
        M1[i] = random_row_switch_matrix(widths[i], widths[i+1]);
    }
}

