#include "conv_txt_to_vec_fn.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main()
{
    
    std::vector<double> pressures_dbl = conv_txt_to_vec_fn("/autoDMP/pump_sinewave/sinewave_txt_file_generation/press[mbar].txt");
    
    for (auto file_line : pressures_dbl)
        std::cout << file_line << std::endl;

    return 0;
}