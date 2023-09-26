#include "conv_txt_to_vec_fn.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//Function that intakes the location of a text file with each line of the text file being a number and the final line being a whitespace and converts
//into a double vector of the values

std::vector<double> conv_txt_to_vec_fn(std::string filepath)
{
    std::ifstream file;
    file.open(filepath);
    std::vector<std::string> vec_values;
    std::string line;

    //While loop that runs as long as end of file not reached, adds each line to vector 
    //(note variables defined in while coniditon or body are created and destroyed on each iteration)
    while (!file.eof())
    {
        getline(file,line);
        vec_values.push_back(line);
    }
    //close file
    file.close();
    
    // initialize variable representing the size of the vector as an integer 
    unsigned int size_vec_values_vec = vec_values.size();

    // declare double vector that will hold the values in the string vector as doubles
    std::vector<double> vec_values_dbl;

    // for loop that takes each element in the  string vector and converts it to a double and adds each element sequentially to the double vector
    // note that the last element of the string vector, due to the nature of the sinewave_txt_file_gen procedure, is an empty string, hence the -1 in the range
    for (unsigned int i = 0; i < (size_vec_values_vec-1) ; i++)
        {
            vec_values_dbl.push_back(std::stod(vec_values[i]));
        }

    return vec_values_dbl;
}