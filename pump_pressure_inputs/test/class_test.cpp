/*
Program testing the implementation of the Pressure_Pump_Input class, 
intakes a pressure and time txt file generated from to_text_file.py 
using the sinewave_fn.py function. Creates test_input object of class 
type Pressure_Pump_Input which has a member variable 'input' which is a 
vector of struct type Pressure_Pump_Input_Single_Value, where each entry
has two member variables: Pressure and Time equal to the pressure and time 
from the text file respectivley. These values are printed to console
using the for loop.
*/
#include "include/pump_pressure_inputs.hpp"

int main(){
    Pressure_Pump_Input test_input("/autoDMP/pump_pressure_inputs/sinewave_txt_file_generation/p_vs_t.txt");
    
    for (auto file_line : test_input.input)
        std::cout << file_line.Time << ','<< file_line.Pressure << std::endl;

    return 0;
}

