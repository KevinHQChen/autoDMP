#include "include/pump_pressure_inputs.hpp"

Pressure_Pump_Input::Pressure_Pump_Input(std::string filepath)
{
    std::ifstream file;
    file.open(filepath);
    
    double t, p;

    while(file >> t >> p)
    {
        input.push_back(Pressure_Pump_Input_Single_Value{t, p});
    }
    file.close();

}
    