import math as m
from matplotlib import pyplot as plt
import numpy as np

"""
Function that intakes desired offset pressure, amplitude of sinewave (as % of offset), desired frequency, desired timestep, and desired number of periods
and outputs a n x 2 numpy array with time[s] and pressure[mbar] columns representing the desired sinewave, notice that it is currently setup to show 
at least the number of periods specified, it is not exact and depending on the arguments may overshoot the desired period. Also plots the sinewave upon function
call so the user can see what kind of wave they generated
""" 


def sinewave(p_offset,p_amp_percentage, f, del_t, num_periods):
    omega =2*m.pi*f
    p_amp = p_amp_percentage*p_offset
    period = 1/f
    p_vs_t = np.zeros((round(period/del_t+1)*num_periods,2))

    for t in range(0,round(period/del_t+1)*num_periods,1):
        p_vs_t[t][0] = t*del_t
        p_vs_t[t][1] = p_offset + p_amp*m.sin(omega*t*del_t)

    plt.plot(p_vs_t[:,0], p_vs_t[:,1])
    plt.show()
    
    return p_vs_t

