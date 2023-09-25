import numpy as np
from sinewave_fn import sinewave
array = sinewave(50, .5, 1, 0.001, 3 )
np.savetxt("p_vs_t.txt", array)
np.savetxt("press[mbar].txt", array[:,1])
np.savetxt("time[s].txt", array [:,0])
