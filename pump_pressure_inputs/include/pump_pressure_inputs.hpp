#pragma once
#include "util/util.hpp"

struct Pressure_Pump_Input_Single_Value{

    public:
    //public member variables
        double Time;
        double Pressure;
}

class Pressure_Pump_Input{

    public:
    //public member variables
        std::vector<Pressure_Pump_Input_Single_Value> input;
    //public constructor
        Pressure_Pump_Input (std::string filepath);
}