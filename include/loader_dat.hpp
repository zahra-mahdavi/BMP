#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <algorithm>

// Loader for .dat format seen in bmp_105.dat:
// Line 1: (ignored) enumeration of initial rows (size N0)
// Line 2: R vector with 2 ints: "0 1"
// Then repeated pairs: for each layer k:
//   - Line A: row indices (ignored except for length)
//   - Line B: column indices (1-based) in [1 .. rows_{k+1}], length = rows_k
// We only need Line B of each pair for traversal by positional index.
struct BMPPositional {
    std::vector<std::vector<uint32_t>> T; // T[k][r] = c (1-based), r in [0..rows_k-1]
    std::vector<int> R;                   // size 2, values {0,1}
};

inline std::vector<int> split_ints(const std::string& s){
    std::vector<int> out;
    std::istringstream iss(s);
    int x;
    while(iss >> x) out.push_back(x);
    return out;
}

inline BMPPositional load_bmp_dat(const std::string& path){
    std::ifstream fin(path);
    if(!fin) throw std::runtime_error("Cannot open file: " + path);
    std::string line;
    std::vector<std::string> lines;
    while(std::getline(fin, line)){
        if(!line.empty() && line.back()=='\r') line.pop_back(); // handle CRLF
        if(!line.empty()) lines.push_back(line);
    }
    if(lines.size() < 3) throw std::runtime_error("Unexpected .dat format: too few lines");

    // Line 1: initial rows enumeration (ignored)
    // Line 2: R vector (expected 2 ints)
    std::vector<int> R = split_ints(lines[1]);
    if(R.size() != 2) throw std::runtime_error("Expected R of size 2 on line 2");

    // Build pairs (equal-length neighbors)
    std::vector<std::vector<uint32_t>> T;
    for(size_t i=2;i+1<lines.size();){
        auto A = split_ints(lines[i]);
        auto B = split_ints(lines[i+1]);
        if(A.size() == B.size()){
            // store B as 1-based uint32_t
            std::vector<uint32_t> Tb(B.size());
            for(size_t j=0;j<B.size();++j){
                int c = B[j];
                if(c < 1) throw std::runtime_error("Column index must be 1-based positive");
                Tb[j] = static_cast<uint32_t>(c);
            }
            T.push_back(std::move(Tb));
            i += 2;
        }else{
            // skip singleton
            ++i;
        }
    }

    // Optional sanity: last layer should map to 1..2
    auto& last = T.back();
    uint32_t maxv = *std::max_element(last.begin(), last.end());
    if(maxv > 2){
        // warn but don't fail
        // throw std::runtime_error("Last T has indices > 2; unexpected for R of size 2");
    }

    BMPPositional bmp;
    bmp.T = std::move(T);
    bmp.R = std::move(R);
    return bmp;
}

inline std::vector<uint8_t> traverse_all_outputs(const BMPPositional& bmp){
    // Start from leftmost layer 0: rows_0 = len(T[0])
    size_t n0 = bmp.T.front().size();
    std::vector<uint8_t> out(n0, 0);
    for(size_t r=0;r<n0;++r){
        uint32_t idx = static_cast<uint32_t>(r+1); // 1-based
        for(const auto& T : bmp.T){
            idx = T[idx-1]; // step forward
        }
        // Map to R
        uint8_t bit = 0;
        if(idx==1 || idx==2){
            bit = static_cast<uint8_t>(bmp.R[idx-1]);
        }else if(idx==0 || idx==1){
            bit = static_cast<uint8_t>(bmp.R[idx]);
        }else{
            bit = 0; // default
        }
        out[r] = bit;
    }
    return out;
}
