#include <iostream>
#include <chrono>
#include <fstream>
#include "loader_dat.hpp"

int main(int argc, char** argv){
    using Clock = std::chrono::steady_clock;

    if(argc < 2){
        std::cerr << "Usage: demo_load <path-to-bmp_XXX.dat> [more.dat ...]\n";
        return 1;
    }

    double total_load = 0.0;
    double total_trav = 0.0;

    
    std::ofstream csv("timings.csv");
    csv << "File,LoadTime,TraversalTime,Layers,FirstLayerSize,Zeros,Ones\n";

    for(int f=1; f<argc; ++f){
        std::string path = argv[f];
        std::cout << "==== File: " << path << " ====\n";

        auto t0 = Clock::now();
        BMPPositional bmp = load_bmp_dat(path);
        auto t1 = Clock::now();

        auto out = traverse_all_outputs(bmp);
        auto t2 = Clock::now();

        size_t ones=0, zeros=0;
        for(auto b: out) (b?ones:zeros)++;

        std::cout << "Layers: " << bmp.T.size() << "\n";
        std::cout << "First layer size: " << bmp.T.front().size() << "\n";
        std::cout << "Outputs: zeros=" << zeros << " ones=" << ones << "\n";
        std::cout << "First 32 bits: ";
        for(size_t i=0;i<std::min<size_t>(32,out.size());++i) std::cout << int(out[i]);
        std::cout << "\n";

        std::chrono::duration<double> dt_load = t1 - t0;
        std::chrono::duration<double> dt_trav = t2 - t1;
        std::cout << "Load time: " << dt_load.count() << " s\n";
        std::cout << "Traversal time: " << dt_trav.count() << " s\n\n";

        total_load += dt_load.count();
        total_trav += dt_trav.count();

        
        csv << path << "," 
            << dt_load.count() << "," 
            << dt_trav.count() << "," 
            << bmp.T.size() << "," 
            << bmp.T.front().size() << "," 
            << zeros << "," 
            << ones << "\n";
    }

    std::cout << "==== Summary for " << (argc-1) << " files ====\n";
    std::cout << "Total load time: " << total_load << " s\n";
    std::cout << "Total traversal time: " << total_trav << " s\n";
    std::cout << "Average load time per file: " << (total_load / (argc-1)) << " s\n";
    std::cout << "Average traversal time per file: " << (total_trav / (argc-1)) << " s\n";

    csv.close();
    return 0;
}
