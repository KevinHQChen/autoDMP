"""
Short script used to generate sinewave using sinewave_fn and produces 3 text files from output, one of both pressure and time , one of 
only pressure, and one of only time. The text files of only pressure or only time can be used with conv_txt_to_vec_fn.cpp to convert the
txt files into double vectors for use in setting pump pressures.
"""


import numpy as np
from sinewave_fn import sinewave
array = sinewave(50, .5, 1, 0.001, 3 )
np.savetxt("p_vs_t.txt", array)
np.savetxt("press[mbar].txt", array[:,1])
np.savetxt("time[s].txt", array [:,0])
