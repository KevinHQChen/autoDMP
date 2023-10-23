#include "include/pump_pressure_inputs.hpp"
/*
Constructor definition, intakes a filepath leading to txt file containing 
time and pressure values in separate columns 
(output for the to_text_file.py program) in the form of a string.
Constructor has a while loop that checks the rows of the text file
one at a time and if there is data in the columns (separated by whitespace) 
to be read into the double variables t (col1) and p (col2), 
then it adds an entry to the member vector 'input' of type 
Pressure_Pump_Input_Single_Value, where the entry has its member variables 
Time set to t and Pressure set to p respectivley.
*/

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
    