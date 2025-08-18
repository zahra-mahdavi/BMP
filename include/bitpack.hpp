#pragma once
#include <vector>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <algorithm>

// Simple bit-packed boolean vector & matrix for boolean semiring (OR-AND).
// Vector is stored as 64-bit words; bit i represents entry i.
struct PackedVector {
    std::vector<uint64_t> words;
    size_t nbits = 0;
    PackedVector() = default;
    explicit PackedVector(size_t bits): nbits(bits) {
        words.assign((bits + 63)/64, 0ULL);
    }
    inline void clear() { std::fill(words.begin(), words.end(), 0ULL); }
    inline uint64_t word(size_t k) const { return words[k]; }
    inline size_t nwords() const { return words.size(); }
    inline void set_bit(size_t i) {
        assert(i < nbits);
        words[i>>6] |= (1ULL << (i&63));
    }
    inline bool get_bit(size_t i) const {
        assert(i < nbits);
        return (words[i>>6] >> (i&63)) & 1ULL;
    }
};

// PackedMatrix: rows x cols, row-major by 64-bit word-chunks across columns.
struct PackedMatrix {
    size_t rows = 0, cols = 0;
    size_t words_per_row = 0;
    std::vector<uint64_t> data; // size = rows * words_per_row

    PackedMatrix() = default;
    PackedMatrix(size_t r, size_t c): rows(r), cols(c) {
        words_per_row = (cols + 63)/64;
        data.assign(rows * words_per_row, 0ULL);
    }
    inline uint64_t* row_ptr(size_t r) { return &data[r * words_per_row]; }
    inline const uint64_t* row_ptr(size_t r) const { return &data[r * words_per_row]; }

    // set A[r,c] = 1
    inline void set_one(size_t r, size_t c) {
        assert(r < rows && c < cols);
        auto p = row_ptr(r);
        p[c>>6] |= (1ULL << (c&63));
    }

    // For row-switch matrices, each row has exactly one 1; return its column index.
    // If multiple bits set, returns the first.
    inline size_t row_switch_col(size_t r) const {
        const uint64_t* p = row_ptr(r);
        for(size_t k=0;k<words_per_row;++k){
            uint64_t w = p[k];
            if(w){
                // return index of least-significant set bit in w
                unsigned tz = (unsigned)__builtin_ctzll(w);
                return (k<<6) + tz;
            }
        }
        // no 1 found -> default 0
        return 0;
    }
};

// y = A ⊗ x  under boolean semiring (AND then OR across columns).
// That is: y_i = OR_j (A_ij & x_j).
// Requires: x.nbits == A.cols and y.nbits == A.rows.
inline void bool_matvec_cpu(const PackedMatrix& A, const PackedVector& x, PackedVector& y){
    assert(x.nbits == A.cols);
    if(y.nbits != A.rows){
        y = PackedVector(A.rows);
    }
    y.clear();
    for(size_t i=0;i<A.rows;++i){
        const uint64_t* row = A.row_ptr(i);
        bool bit = false;
        for(size_t k=0;k<A.words_per_row;++k){
            uint64_t w = row[k] & x.words[k];
            if(w){ bit = true; break; }
        }
        if(bit){
            y.set_bit(i);
        }
    }
}

// Convenience: apply selected matrix from array M[i] depending on input bit b (0/1).
inline void select_and_apply(const std::vector<PackedMatrix>& M0,
                             const std::vector<PackedMatrix>& M1,
                             const std::vector<int>& x_bits,
                             PackedVector& v) {
    assert(M0.size() == M1.size());
    assert(x_bits.size() == M0.size());
    for(int i=(int)M0.size()-1;i>=0;--i){
        const PackedMatrix& A = (x_bits[i] ? M1[i] : M0[i]);
        PackedVector y(A.rows);
        bool_matvec_cpu(A, v, y);
        v = std::move(y);
    }
}

// Build a random row-switch matrix of size r x c (each row has exactly one 1).
inline PackedMatrix random_row_switch_matrix(size_t r, size_t c) {
    PackedMatrix A(r,c);
    for(size_t i=0;i<r;++i){
        size_t col = (size_t) (rand() % c);
        A.set_one(i,col);
    }
    return A;
}
