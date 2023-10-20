#include <iostream>
#include <fstream>
#include <vector>
#include <string>

struct p_v_t{
    double Time;
    double Pressure;
};

//File used as precursor to conv_txt_to_vec_fn to get the logic correct, same general code as the fn
int main()
{
    std::ifstream file;
    file.open("/autoDMP/pump_pressure_inputs/sinewave_txt_file_generation/p_vs_t.txt");
    std::vector<p_v_t> press_vs_time;
    double time;
    double pressure;

    while (file >> time >> pressure)
    {
        press_vs_time.push_back(p_v_t{time, pressure});
    }
    
    file.close();
    
    // unsigned int size_press_vec = pressures.size();

    // std::vector<double> pressures_dbl;
    // for (unsigned int i = 0; i < (size_press_vec-1) ; i++)
    //     {
    //         pressures_dbl.push_back(std::stod(pressures[i]));
    //     }
        

    for (auto file_line : press_vs_time)
        std::cout << file_line.Time << ','<< file_line.Pressure << std::endl;

    return 0;
}