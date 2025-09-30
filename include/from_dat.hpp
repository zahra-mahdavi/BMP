#pragma once
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

#include "loader_dat.hpp"   // BMPPositional
#include "traversal.hpp"    // TraversalBMP, LocVec
#include "bitpack.hpp"      // PackedMatrix, PackedVector

// Build TraversalBMP from parsed DAT (pad per-layer rows to max_rows and clamp cols to [0, right_cols-1]).
inline TraversalBMP build_from_dat(const BMPPositional& bmp) {
    if (bmp.T.empty()) throw std::runtime_error("build_from_dat: bmp.T is empty");
    TraversalBMP out;
    const size_t L = bmp.T.size();
    const size_t right_cols = bmp.R.size();

    // Find maximum row count among layers
    size_t max_rows = 0;
    for (const auto& loc : bmp.T) max_rows = std::max(max_rows, loc.size());
    if (max_rows == 0) throw std::runtime_error("build_from_dat: all layers have zero rows");

    out.M0_loc.resize(L);
    out.M1_loc.resize(L);

    for (size_t l = 0; l < L; ++l) {
        const auto& src = bmp.T[l];
        const size_t rows = src.size();

        out.M0_loc[l].resize(max_rows);
        out.M1_loc[l].resize(max_rows);

        for (size_t r = 0; r < max_rows; ++r) {
            // If this layer is shorter, repeat its last valid column (or 0 if empty)
            size_t col = (r < rows) ? src[r] : (rows ? src.back() : 0);
            if (right_cols > 0 && col >= right_cols) col = right_cols - 1; // clamp into range
            out.M0_loc[l][r] = col;
            out.M1_loc[l][r] = col; // TODO: if you have distinct loc1 per layer, assign here instead
        }
    }

    // R_bits
    out.R_bits.resize(right_cols);
    for (size_t i = 0; i < right_cols; ++i) {
        out.R_bits[i] = static_cast<uint8_t>(bmp.R[i] & 1);
    }

    out.n_left_rows = static_cast<uint32_t>(max_rows);
    return out;
}

// Convert DAT -> Packed for CPU matvec baseline (safe per layer)
inline void dat_to_packed(const BMPPositional& bmp,
                          std::vector<PackedMatrix>& M0,
                          std::vector<PackedMatrix>& M1,
                          PackedVector& R) {
    const size_t L = bmp.T.size();
    if (L == 0) throw std::runtime_error("dat_to_packed: bmp.T is empty");

    const size_t right_cols = bmp.R.size();
    M0.resize(L);
    M1.resize(L);

    for (size_t l = 0; l < L; ++l) {
        const auto& loc = bmp.T[l];
        const size_t rows = loc.size(); // rows may differ per layer
        M0[l] = PackedMatrix(rows, right_cols);
        M1[l] = PackedMatrix(rows, right_cols);

        for (size_t r = 0; r < rows; ++r) {
            size_t col = loc[r];
            if (right_cols == 0) continue;               // nothing to set
            if (col >= right_cols) col = right_cols - 1; // clamp defensively
            M0[l].set_one(r, col);
            M1[l].set_one(r, col); // TODO: use distinct col1 if available
        }
    }

    R = PackedVector(right_cols);
    for (size_t i = 0; i < right_cols; ++i) {
        if (bmp.R[i] & 1) R.set_bit(i);
    }
}

/* ==== SAFE versions of layouts: compact (no global padding) and arena (contiguous)
   These ensure each layer has at least n_left_rows entries to be compatible with
   traversal that expects access up to n_left_rows on every layer. Additional entries
   (if needed) repeat the last valid index of that layer.
*/

inline TraversalBMP build_from_dat_compact(const BMPPositional& bmp) {
    if (bmp.T.empty()) throw std::runtime_error("build_from_dat_compact: bmp.T is empty");
    TraversalBMP out;
    const size_t L = bmp.T.size();
    const size_t right_cols = bmp.R.size();
    const size_t n_left_rows = bmp.T.front().size(); // rows on the leftmost matrix

    out.M0_loc.resize(L);
    out.M1_loc.resize(L);

    for (size_t l = 0; l < L; ++l) {
        const auto& src1 = bmp.T[l];
        const size_t rows = src1.size();
        const size_t tgt  = rows < n_left_rows ? n_left_rows : rows; // ensure compatible length

        out.M0_loc[l].resize(tgt);
        out.M1_loc[l].resize(tgt);

        // fill actual rows
        for (size_t r = 0; r < rows; ++r) {
            uint32_t col0 = (src1[r] > 0) ? (static_cast<uint32_t>(src1[r]) - 1u) : 0u;
            if (right_cols && col0 >= right_cols) col0 = static_cast<uint32_t>(right_cols - 1);
            out.M0_loc[l][r] = col0;
            out.M1_loc[l][r] = col0;
        }
        // pad tail (if needed) by repeating the last valid index
        if (rows < tgt) {
            uint32_t last = rows ? out.M0_loc[l][rows-1] : 0u;
            for (size_t r = rows; r < tgt; ++r) {
                out.M0_loc[l][r] = last;
                out.M1_loc[l][r] = last;
            }
        }
    }

    out.R_bits.resize(right_cols);
    for (size_t i = 0; i < right_cols; ++i) out.R_bits[i] = static_cast<uint8_t>(bmp.R[i] & 1);
    out.n_left_rows = static_cast<uint32_t>(n_left_rows);
    return out;
}

inline TraversalBMP build_from_dat_arena(const BMPPositional& bmp,
                                         std::vector<uint32_t>& arena,
                                         std::vector<std::pair<size_t,size_t>>& layer_ranges) {
    if (bmp.T.empty()) throw std::runtime_error("build_from_dat_arena: bmp.T is empty");
    TraversalBMP out;
    const size_t L = bmp.T.size();
    const size_t right_cols = bmp.R.size();
    const size_t n_left_rows = bmp.T.front().size();

    // compute total rows across layers (actual rows only)
    size_t total = 0;
    for (const auto& v : bmp.T) total += v.size();
    arena.resize(total);
    layer_ranges.resize(L);

    // fill arena with zero-based, clamped cols
    size_t off = 0;
    for (size_t l = 0; l < L; ++l) {
        const auto& src1 = bmp.T[l];
        const size_t rows = src1.size();
        layer_ranges[l] = {off, rows};
        for (size_t r = 0; r < rows; ++r) {
            uint32_t col0 = (src1[r] > 0) ? (static_cast<uint32_t>(src1[r]) - 1u) : 0u;
            if (right_cols && col0 >= right_cols) col0 = static_cast<uint32_t>(right_cols - 1);
            arena[off + r] = col0;
        }
        off += rows;
    }

    // expose per-layer vectors (phase 1). Also pad to n_left_rows if needed by repeating last.
    out.M0_loc.resize(L);
    out.M1_loc.resize(L);
    for (size_t l = 0; l < L; ++l) {
        auto [o, rows] = layer_ranges[l];
        out.M0_loc[l].assign(&arena[o], &arena[o] + rows);
        out.M1_loc[l]  = out.M0_loc[l];
        if (rows < n_left_rows) {
            uint32_t last = rows ? out.M0_loc[l].back() : 0u;
            out.M0_loc[l].insert(out.M0_loc[l].end(), n_left_rows - rows, last);
            out.M1_loc[l] = out.M0_loc[l];
        }
    }

    out.R_bits.resize(right_cols);
    for (size_t i = 0; i < right_cols; ++i) out.R_bits[i] = static_cast<uint8_t>(bmp.R[i] & 1);
    out.n_left_rows = static_cast<uint32_t>(n_left_rows);
    return out;
}
