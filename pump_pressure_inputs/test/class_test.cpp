#include "include/pump_pressure_inputs.hpp"

int main(){
    Pressure_Pump_Input test_input("/autoDMP/pump_pressure_inputs/sinewave_txt_file_generation/p_vs_t.txt");
    
    for (auto file_line : test_input.input)
        std::cout << file_line.Time << ','<< file_line.Pressure << std::endl;

    return 0;
}

