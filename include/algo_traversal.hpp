#pragma once
#include "algo_iface.hpp"
#include "from_dat.hpp"
#include "traversal.hpp"

// Direct traversal along location-vectors
struct AlgoTraversal final : IBMPAlgorithm {
    TraversalBMP T;  // built once from DAT
    std::string name() const override { return "traversal"; }
    void preprocess(const BMPPositional& bmp) override {
        T = build_from_dat(bmp);
    }
    std::vector<uint8_t> eval(const std::vector<int>& x_bits) const override {
        return traverse_all(T, x_bits);
    }
};
