#include <iostream>
#include "loader_dat.hpp"

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: demo_load <path-to-bmp_XXX.dat>\n";
        return 1;
    }
    BMPPositional bmp = load_bmp_dat(argv[1]);
    auto out = traverse_all_outputs(bmp);
    size_t ones=0, zeros=0;
    for(auto b: out) (b?ones:zeros)++;
    std::cout << "Layers: " << bmp.T.size() << "\n";
    std::cout << "First layer size: " << bmp.T.front().size() << "\n";
    std::cout << "Outputs: zeros=" << zeros << " ones=" << ones << "\n";
    std::cout << "First 32 bits: ";
    for(size_t i=0;i<std::min<size_t>(32,out.size());++i) std::cout << int(out[i]);
    std::cout << "\n";
    return 0;
}
