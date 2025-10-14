#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "loader_dat.hpp"
#include "algo_iface.hpp"
#include "algo_traversal.hpp"
#include "algo_cpu_matvec.hpp"

struct Timer {
    using Clock = std::chrono::steady_clock;
    Clock::time_point t0;
    Timer() : t0(Clock::now()) {}
    double seconds() const {
        std::chrono::duration<double> dt = Clock::now() - t0;
        return dt.count();
    }
};

static std::vector<int> random_bits(size_t n, std::mt19937_64& rng){
    std::vector<int> x(n);
    std::uniform_int_distribution<int> d(0,1);
    for(size_t i=0;i<n;++i) x[i] = d(rng);
    return x;
}

static void write_csv_header(std::ofstream& csv){
    csv << "file,algo,n_left,n_right,reps,eval_mean_s,eval_std_s,consensus\n";
}

static void summarize_row(std::ofstream& csv,
                          const std::string& file, const std::string& algo,
                          size_t n_left, size_t n_right,
                          const std::vector<double>& times, bool consensus){
    // mean / std
    double mean = 0.0;
    for(double t: times) mean += t;
    mean /= (times.empty()?1:times.size());
    double var = 0.0;
    for(double t: times) var += (t-mean)*(t-mean);
    var /= (times.size()>=2 ? (times.size()-1) : 1);
    double std = std::sqrt(var);
    csv << file << "," << algo << "," << n_left << "," << n_right << ","
        << times.size() << "," << std::setprecision(9) << mean << "," << std << ","
        << (consensus?1:0) << "\n";
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " <bmp_*.dat> [more .dat ...] "
                  << "[--reps N=8] [--csv out.csv]\n";
        return 1;
    }

    int reps = 8;
    std::string csv_path = "consensus_results.csv";
    std::vector<std::string> inputs;
    for(int i=1;i<argc;++i){
        std::string a = argv[i];
        if(a=="--reps" && i+1<argc){ reps = std::max(1, std::atoi(argv[++i])); continue; }
        if(a=="--csv"  && i+1<argc){ csv_path = argv[++i]; continue; }
        if(a.rfind("--",0)==0){ std::cerr << "Unknown flag: " << a << "\n"; return 2; }
        inputs.push_back(a);
    }
    if(inputs.empty()){
        std::cerr << "No .dat inputs provided\n";
        return 3;
    }

    std::ofstream csv(csv_path);
    if(!csv){ std::cerr << "Cannot open CSV for writing: " << csv_path << "\n"; return 4; }
    write_csv_header(csv);

    std::mt19937_64 rng(0xB00L1234ULL);

    for(const auto& path : inputs){
        try {
            std::cout << "==> Loading: " << path << "\n";
            Timer t_load;
            BMPPositional bmp = load_bmp_dat(path);
            double t_load_s = t_load.seconds();

            const size_t L = bmp.T.size();        // number of layers/variables
            const size_t n_left = bmp.T.front().size();
            const size_t n_right = 2;             // by format (R of size 2)

            // Build algorithms
            std::vector<std::unique_ptr<IBMPAlgorithm>> algos;
            algos.emplace_back(new AlgoTraversal());
            algos.emplace_back(new AlgoCpuMatvec());
            for(auto& a : algos) a->preprocess(bmp);

            // Run reps with same x for all algorithms per rep
            std::vector<std::vector<uint8_t>> last_outputs(algos.size());
            std::vector<std::vector<double>> times(algos.size());
            bool consensus_all = true;

            for(int r=0; r<reps; ++r){
                auto x = random_bits(L, rng);
                for(size_t ai=0; ai<algos.size(); ++ai){
                    Timer te;
                    auto out = algos[ai]->eval(x);
                    double te_s = te.seconds();
                    times[ai].push_back(te_s);
                    last_outputs[ai] = std::move(out);
                }
                // consensus check vs algo 0
                for(size_t ai=1; ai<algos.size(); ++ai){
                    if(last_outputs[ai] != last_outputs[0]){
                        consensus_all = false;
                        // print first mismatch index to stderr
                        for(size_t k=0;k<last_outputs[0].size();++k){
                            if(last_outputs[ai][k] != last_outputs[0][k]){
                                std::cerr << "MISMATCH file=" << path
                                          << " algo=" << algos[ai]->name()
                                          << " at row=" << k
                                          << " ref=" << int(last_outputs[0][k])
                                          << " got=" << int(last_outputs[ai][k]) << "\n";
                                break;
                            }
                        }
                    }
                }
            }

            // Write rows
            for(size_t ai=0; ai<algos.size(); ++ai){
                summarize_row(csv, path, algos[ai]->name(), n_left, n_right, times[ai], consensus_all);
                std::cout << "  " << algos[ai]->name() << ": "
                          << "mean(eval) = " << std::setprecision(6)
                          << (times[ai].empty()?0.0:
                              std::accumulate(times[ai].begin(), times[ai].end(), 0.0)/times[ai].size())
                          << " s\n";
            }

            std::cout << "  Load time: " << t_load_s << " s | consensus=" << (consensus_all? "YES":"NO") << "\n";
        } catch(const std::exception& e){
            std::cerr << "ERROR processing '" << path << "': " << e.what() << "\n";
        }
    }

    std::cout << "Done. CSV: " << csv_path << "\n";
    return 0;
}
