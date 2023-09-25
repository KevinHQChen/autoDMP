#include "conv_txt_to_vec_fn.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

std::vector<double> conv_txt_to_vec_fn(std::string filepath)
{
    std::ifstream file;
    file.open(filepath);
    std::vector<std::string> vec_values;
    std::string line;

    while (!file.eof())
    {
        getline(file,line);
        vec_values.push_back(line);
    }
    
    file.close();
    
    unsigned int size_vec_values_vec = vec_values.size();

    std::vector<double> vec_values_dbl;

    for (unsigned int i = 0; i < (size_vec_values_vec-1) ; i++)
        {
            vec_values_dbl.push_back(std::stod(vec_values[i]));
        }

    return vec_values_dbl;
}