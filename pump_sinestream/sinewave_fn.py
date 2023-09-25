import math as m
from matplotlib import pyplot as plt
import numpy as np

def sinewave(p_offset,p_amp_percentage, f, del_t, num_periods):
    omega =2*m.pi*f
    p_amp = p_amp_percentage*p_offset
    period = 1/f
    p_vs_t = np.zeros((round(period/del_t+1)*num_periods,2))

    for t in range(0,round(period/del_t+1)*num_periods,1):
        p_vs_t[t][0] = t*del_t
        p_vs_t[t][1] = p_offset + p_amp*m.sin(omega*t*del_t)
    
    return p_vs_t

