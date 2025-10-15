#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <filesystem>
#include <cstring>

#include "loader_dat.hpp"
#include "from_dat.hpp"
#include "traversal.hpp"
#include "bitpack.hpp"
#include "contractor.hpp"
#include "bench_utils.hpp"
#include "algo_iface.hpp"
#include "algo_traversal.hpp"
#include "algo_traversal_unpadded.hpp"
#include "algo_traversal_padded.hpp"
#include "algo_cpu_matvec.hpp"

static inline uint64_t compute_bmp_volume(const BMPPositional& bmp){
    uint64_t sum = 0;
    const size_t L = bmp.T.size();
    const size_t n_right = bmp.R.size();
    for(size_t l=0; l<L; ++l){
        const uint64_t rows = (uint64_t)bmp.T[l].size();
        const uint64_t cols = (uint64_t)((l+1 < L) ? bmp.T[l+1].size() : n_right);
        sum += rows * cols;
    }
    return sum;
}

struct Cli {
    int reps = 100;
    std::string csv = "cells_per_sec.csv";
    std::vector<std::string> algos;
    std::vector<std::string> paths;
};

static Cli parse_cli(int argc, char** argv){
    Cli c;
    for(int i=1;i<argc;++i){
        std::string a = argv[i];
        auto getv = [&](const char* pfx)->const char*{
            size_t n = std::strlen(pfx);
            if(a.rfind(pfx,0)==0) return a.c_str() + n;
            return nullptr;
        };
        if(auto v = getv("--reps=")) { c.reps = std::max(1, std::atoi(v)); continue; }
        if(auto v = getv("--csv="))  { c.csv = v; continue; }
        if(auto v = getv("--algos=")){
            std::string s = v;
            size_t p = 0;
            while(true){
                size_t q = s.find(',', p);
                std::string t = s.substr(p, q==std::string::npos ? q : q-p);
                if(!t.empty()) c.algos.push_back(t);
                if(q==std::string::npos) break;
                p = q + 1;
            }
            continue;
        }
        if(a=="-h" || a=="--help"){
            std::cout << "Usage: bench_cells [--reps=100] [--csv=out.csv] [--algos=traversal_unpadded,traversal_padded_64,cpu_matvec] files.dat...\n";
            std::exit(0);
        }
        c.paths.push_back(a);
    }
    if(c.paths.empty()){ std::cerr << "No .dat files provided.\n"; std::exit(1); }
    return c;
}

static std::unique_ptr<IBMPAlgorithm> make_algo(const std::string& tag){
    if(tag == "traversal"){ return std::unique_ptr<IBMPAlgorithm>(new AlgoTraversal()); }
    else if(tag == "traversal_unpadded"){ return std::unique_ptr<IBMPAlgorithm>(new AlgoTraversalUnpadded()); }
    else if(tag.rfind("traversal_padded_",0)==0){
        auto P = std::make_unique<AlgoTraversalPadded>();
        size_t align = 64;
        try { align = std::stoul(tag.substr(std::string("traversal_padded_").size())); } catch(...){}
        P->align = align; return P;
    } else if(tag == "cpu_matvec"){ return std::unique_ptr<IBMPAlgorithm>(new AlgoCpuMatvec()); }
    std::cerr << "Unknown algo tag: " << tag << "\n"; std::exit(2);
}

int main(int argc, char** argv){
    Cli cli = parse_cli(argc, argv);
    if(cli.algos.empty()){
        cli.algos = {"traversal_unpadded", "traversal_padded_64", "cpu_matvec"};
    }
    std::cout << "bench_cells: reps=" << cli.reps << " csv=" << cli.csv << " algos=";
    for(size_t i=0;i<cli.algos.size();++i){ if(i) std::cout << ","; std::cout << cli.algos[i]; }
    std::cout << " files=" << cli.paths.size() << "\n";

    const std::filesystem::path csv_path{cli.csv};
    bool need_header = True;
    if (std::filesystem::exists(csv_path)) need_header = (std::filesystem::file_size(csv_path) == 0);
    std::ofstream csv(cli.csv, std::ios::app);
    if(!csv){ std::cerr << "Cannot open CSV: " << cli.csv << "\n"; return 3; }
    csv.setf(std::ios::scientific); csv << std::setprecision(6);
    if(need_header){
        csv << "file,algo,rep,n_left,n_right,nvars,bmp_volume,eval_s,cells_per_sec\n";
    }

    std::mt19937_64 rng(0xBEEFBEEF);

    for(const auto& path : cli.paths){
        try{
            BMPPositional bmp = load_bmp_dat(path);
            const uint64_t bmp_vol = compute_bmp_volume(bmp);
            const size_t nvars = bmp.T.size();
            const size_t n_left = bmp.T.front().size();
            const size_t n_right = bmp.R.size();

            for(const auto& tag : cli.algos){
                auto A = make_algo(tag);
                A->preprocess(bmp);

                for(int r=0; r<cli.reps; ++r){
                    std::vector<int> x(nvars);
                    std::uniform_int_distribution<int> d(0,1);
                    for(size_t i=0;i<nvars;++i) x[i] = d(rng);

                    auto t0 = std::chrono::steady_clock::now();
                    auto y = A->eval(x);
                    (void)y;
                    auto t1 = std::chrono::steady_clock::now();
                    std::chrono::duration<double> dt = t1 - t0;
                    double eval_s = dt.count();
                    double cps = (eval_s > 0.0) ? (double)bmp_vol / eval_s : 0.0;

                    csv << path << "," << tag << "," << r << ","
                        << n_left << "," << n_right << "," << nvars << ","
                        << bmp_vol << "," << std::setprecision(9) << eval_s << ","
                        << std::setprecision(9) << cps << "\n";
                }
                std::cout << "  wrote reps for algo=" << tag << " file=" << path << "\n";
            }
        }catch(const std::exception& e){
            std::cerr << "ERROR file='" << path << "': " << e.what() << "\n";
        }
    }

    std::cout << "Done. CSV: " << std::filesystem::absolute(csv_path).string() << "\n";
    return 0;
}
