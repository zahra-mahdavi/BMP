
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <filesystem>
#include <chrono>
#include <cmath>

#include "loader_dat.hpp"          // BMPPositional, load_bmp_dat(...)
#include "traversal.hpp"           // TraversalBMP, traverse_all(...)
#include "bitpack.hpp"             // PackedMatrix, PackedVector, bool_matvec_cpu, select_and_apply
#include "contractor.hpp"
#include "from_dat.hpp"            // build_from_dat (padded) + added compact/arena builders appended
#include "stats.hpp"

using Clock = std::chrono::high_resolution_clock;
using dsec  = std::chrono::duration<double>;

struct Cli {
    std::string algo = "both";           // traversal | cpu_matvec | both
    int reps = 300;                      // per meeting: increase sampling
    std::string csv = "results.csv";
    std::string layout = "padded";       // padded | vecvec | arena
    std::vector<std::string> paths;      // .dat files
};

static Cli parse_cli(int argc, char** argv){
    Cli c;
    for(int i=1;i<argc;++i){
        std::string a = argv[i];
        auto getv = [&](const char* pfx)->const char*{
            auto n = std::strlen(pfx);
            return (a.rfind(pfx,0)==0) ? a.c_str()+n : nullptr;
        };
        if(const char* v = getv("--algo="))     c.algo = v;
        else if(const char* v = getv("--reps="))   c.reps = std::max(1, std::atoi(v));
        else if(const char* v = getv("--csv="))    c.csv = v;
        else if(const char* v = getv("--layout=")) c.layout = v;
        else if(a=="-h" || a=="--help"){
            std::cout << "Usage: bench_bmp [--algo=traversal|cpu_matvec|both] "
                         "[--layout=padded|vecvec|arena] [--reps=K] [--csv=out.csv] files.dat...\n";
            std::exit(0);
        } else {
            c.paths.push_back(a);
        }
    }
    if(c.paths.empty()){
        std::cerr << "No .dat files provided.\n";
        std::exit(1);
    }
    return c;
}

// Build PackedMatrix chain + R from BMPPositional (1-based columns -> 0-based)
static void build_packed_from_dat_local(const BMPPositional& bmp,
                                        std::vector<PackedMatrix>& M0,
                                        std::vector<PackedMatrix>& M1,
                                        PackedVector& R)
{
    if(bmp.T.empty()) throw std::runtime_error("build_packed_from_dat_local: empty T");
    const size_t L = bmp.T.size();
    const size_t right_cols = bmp.R.size();
    M0.resize(L);
    M1.resize(L);
    for(size_t l=0;l<L;++l){
        const auto& loc = bmp.T[l];
        const size_t rows = loc.size();
        M0[l] = PackedMatrix(rows, right_cols);
        M1[l] = PackedMatrix(rows, right_cols);
        for(size_t r=0;r<rows;++r){
            size_t col1 = (size_t)(loc[r] ? (loc[r]-1) : 0);
            if(right_cols && col1 >= right_cols) col1 = right_cols - 1;
            M0[l].set_one(r, col1);
            M1[l].set_one(r, col1); // if M1 differs, set accordingly
        }
    }
    R = PackedVector(right_cols);
    for(size_t i=0;i<right_cols;++i){
        if(bmp.R[i] & 1) R.set_bit(i);
    }
}

static void gen_random_bits(std::mt19937& rng, size_t n, std::vector<int>& bits){
    bits.resize(n);
    std::bernoulli_distribution B(0.5);
    for(size_t i=0;i<n;++i) bits[i] = B(rng) ? 1 : 0;
}

static void summarize(const std::vector<double>& v, double& mean, double& std, double& se){
    if(v.empty()){ mean=std=se=0; return; }
    double s=0; for(double x: v) s+=x;
    mean = s / v.size();
    double var=0; for(double x: v){ double d=x-mean; var += d*d; }
    var /= (v.size()>1 ? (v.size()-1) : 1);
    std = std::sqrt(var);
    se  = std / std::sqrt((double)v.size());
}

int main(int argc, char** argv){
    try{
        Cli cli = parse_cli(argc, argv);
        std::ofstream csv(cli.csv);
        csv << "file,algo,layout,n_left,n_right,bmp_volume,load_s,pre_s,eval_mean_s,eval_std_s,eval_se_s,match\n";

        std::mt19937 rng(42);

        for(const std::string& path : cli.paths){
            try{
                auto t0 = Clock::now();
                BMPPositional bmp = load_bmp_dat(path.c_str());
                auto t1 = Clock::now();

                // Build TraversalBMP based on layout
                TraversalBMP T;
                double pre_s = 0.0;
                {
                    auto t_pre0 = Clock::now();
                    if(cli.layout == "vecvec"){
                        T = build_from_dat_compact(bmp);
                    } else if(cli.layout == "arena"){
                        static std::vector<uint32_t> arena;
                        static std::vector<std::pair<size_t,size_t>> ranges;
                        T = build_from_dat_arena(bmp, arena, ranges);
                    } else {
                        T = build_from_dat(bmp); // padded
                    }
                    auto t_pre1 = Clock::now();
                    pre_s = std::chrono::duration_cast<dsec>(t_pre1 - t_pre0).count();
                }

                // Build packed baseline
                std::vector<PackedMatrix> M0, M1;
                PackedVector Rvec;
                build_packed_from_dat_local(bmp, M0, M1, Rvec);

                const size_t L = T.M0_loc.size();
                const size_t n_left = T.n_left_rows;
                const size_t n_right = T.R_bits.size();

                // Simple proxy for volume: total rows across layers
                size_t bmp_volume = 0;
                for(const auto& v : bmp.T) bmp_volume += v.size();

                std::vector<double> t_trav, t_mul;
                t_trav.reserve(cli.reps);
                t_mul.reserve(cli.reps);
                int match_total = 0;

                for(int rep=0; rep<cli.reps; ++rep){
                    std::vector<int> x_bits;
                    gen_random_bits(rng, L, x_bits);

                    // Traversal
                    auto tt0 = Clock::now();
                    auto out_trav = traverse_all(T, x_bits);
                    auto tt1 = Clock::now();
                    t_trav.push_back(std::chrono::duration_cast<dsec>(tt1 - tt0).count());

                    bool need_mul = (cli.algo=="cpu_matvec" || cli.algo=="both");
                    if(need_mul){
                        auto tm0 = Clock::now();
                        PackedVector v = Rvec;
                        select_and_apply(M0, M1, x_bits, v); // right-to-left
                        auto tm1 = Clock::now();
                        t_mul.push_back(std::chrono::duration_cast<dsec>(tm1 - tm0).count());

                        int ok = 1;
                        if(v.nbits != n_left) ok = 0;
                        else{
                            for(size_t i=0;i<n_left;++i){
                                uint8_t b = v.get_bit(i) ? 1 : 0;
                                if(b != out_trav[i]){ ok = 0; break; }
                            }
                        }
                        match_total += ok;
                    }
                }

                double load_s = std::chrono::duration_cast<dsec>(t1 - t0).count();

                auto write_row = [&](const std::string& algo, const std::vector<double>& ts, int match_flag){
                    double mean=0, std=0, se=0; summarize(ts, mean, std, se);
                    csv << std::filesystem::path(path).filename().string() << ","
                        << algo << "," << cli.layout << ","
                        << n_left << "," << n_right << "," << bmp_volume << ","
                        << std::fixed << std::setprecision(6)
                        << load_s << "," << pre_s << ","
                        << mean << "," << std << "," << se << ","
                        << match_flag << "\n";
                };

                if(cli.algo=="traversal"){
                    write_row("traversal", t_trav, 0);
                } else if(cli.algo=="cpu_matvec"){
                    write_row("cpu_matvec", t_mul, 0);
                } else {
                    write_row("traversal", t_trav, match_total==cli.reps ? 1 : 0);
                    write_row("cpu_matvec", t_mul, match_total==cli.reps ? 1 : 0);
                }

            }catch(const std::exception& e){
                std::cerr << "ERROR while processing '" << path << "': " << e.what() << "\n";
            }catch(...){
                std::cerr << "ERROR while processing '" << path << "': unknown exception\n";
            }
        }

        std::cout << "All done. CSV: " << std::filesystem::absolute(cli.csv).string() << "\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr << "FATAL: " << e.what() << "\n";
        return 2;
    }
}
