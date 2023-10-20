#pragma once
#include "util/util.hpp"

class Pressure_Pump_Input{

    public:
    //public member variables
        double Time;
        double Pressure;
    //public member functions
        std::vector<Pressure_Pump_Input> text_to_vector(std::string filepath);
    


};