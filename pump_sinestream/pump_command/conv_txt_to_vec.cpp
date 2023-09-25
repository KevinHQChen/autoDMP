#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main()
{
    std::ifstream file;
    file.open("/autoDMP/pump_sinestream/sinewave_txt_file_generation/p.txt");
    std::vector<std::string> pressures;
    std::string line;

    while (!file.eof())
    {
        getline(file,line);
        pressures.push_back(line);
    }
    
    file.close();
    
    unsigned int size_press_vec = pressures.size();

    std::vector<double> pressures_dbl;
    for (unsigned int i = 0; i < (size_press_vec-1) ; i++)
        {
            pressures_dbl.push_back(std::stod(pressures[i]));
        }
        

    for (auto file_line : pressures_dbl)
        std::cout << file_line << std::endl;

    return 0;
}