#pragma once
#include "util/util.hpp"
/*
struct type definition, holds a time value and a pressure value, of
type double, in respective public member variables
*/

struct Pressure_Pump_Input_Single_Value{

    public:
    //public member variables
        double Time;
        double Pressure;
};

/*
class type definition, holds a vector of struct type Pressure_Pump_Input_
Single_Value as a member variable named 'input', also has a constructor
which intakes a filepath to a txt file containing time and pressure values
in seperate columns (output for the to_text_file.py program)
in the form of a string. See source file for constructor definition.
*/

class Pressure_Pump_Input{

    public:
    //public member variables
        std::vector<Pressure_Pump_Input_Single_Value> input;
    //public constructor
        Pressure_Pump_Input (std::string filepath);
};