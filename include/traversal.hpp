#pragma once
#include <vector>
#include <cstdint>
#include <cassert>


using LocVec = std::vector<uint32_t>;

struct TraversalBMP {
    
    std::vector<LocVec> M0_loc;
    std::vector<LocVec> M1_loc;
    std::vector<uint8_t> R_bits; 
    uint32_t n_left_rows = 0;    

    inline size_t nvars() const { return M0_loc.size(); }
};

// Traverse a single output bit: start at left row r0, walk right-to-left following loc vectors.
inline uint8_t traverse_one(const TraversalBMP& T, const std::vector<int>& x_bits, uint32_t r0){
    assert(T.M0_loc.size() == T.M1_loc.size());
    assert(x_bits.size() == T.M0_loc.size());
    uint32_t idx = r0;
    for(int i=(int)T.nvars()-1;i>=0;--i){
        const LocVec& loc = (x_bits[i] ? T.M1_loc[i] : T.M0_loc[i]);
        assert(idx < loc.size());
        idx = loc[idx];
    }
    assert(idx < T.R_bits.size());
    return T.R_bits[idx] & 1u;
}

// Traverse all outputs.
inline std::vector<uint8_t> traverse_all(const TraversalBMP& T, const std::vector<int>& x_bits){
    std::vector<uint8_t> out(T.n_left_rows, 0);
    for(uint32_t r=0;r<T.n_left_rows;++r){
        out[r] = traverse_one(T, x_bits, r);
    }
    return out;
}
