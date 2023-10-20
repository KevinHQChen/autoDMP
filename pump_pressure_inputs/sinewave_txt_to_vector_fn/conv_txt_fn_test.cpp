#include "conv_txt_to_vec_fn.hpp"


//File used to test conv_txt_to_vec_fn

int main()
{
    
    std::vector<double> pressures_dbl = conv_txt_to_vec_fn("/autoDMP/pump_pressure_inputs/sinewave_txt_file_generation/press[mbar].txt");
    
    for (auto file_line : pressures_dbl)
        std::cout << file_line << std::endl;

    return 0;
}