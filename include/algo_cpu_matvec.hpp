#pragma once
#include "algo_iface.hpp"
#include "from_dat.hpp"
#include "contractor.hpp"

// Baseline: boolean matrix-vector chain on CPU
struct AlgoCpuMatvec final : IBMPAlgorithm {
    std::vector<PackedMatrix> M0, M1;
    PackedVector R;
    std::string name() const override { return "cpu_matvec"; }
    void preprocess(const BMPPositional& bmp) override {
        dat_to_packed(bmp, M0, M1, R);
    }
    std::vector<uint8_t> eval(const std::vector<int>& x_bits) const override {
        return compute_via_cpu_matvec(M0, M1, R, x_bits);
    }
};
