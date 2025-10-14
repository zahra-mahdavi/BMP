#pragma once
#include <string>
#include <vector>
#include "loader_dat.hpp"

// Interface for BMP algorithms (Boolean Matrix Product traversal / eval variants)
struct IBMPAlgorithm {
    virtual ~IBMPAlgorithm() = default;
    virtual std::string name() const = 0;
    virtual void preprocess(const BMPPositional& bmp) = 0;
    // Evaluate y = f(x_bits) returning n_left_rows bits (0/1) as uint8_t
    virtual std::vector<uint8_t> eval(const std::vector<int>& x_bits) const = 0;
};
