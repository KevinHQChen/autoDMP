# import traceback
# try:
#     import numpy
# except:
#     traceback.print_exc()
import numpy as np
from scipy.signal import max_len_seq, periodogram
from numpy.fft import fft, ifft, fftshift, fftfreq

def prbs(selCh, minVal, maxVal, order):
    """function that generates a PRBS sequence.
    Parameters
    ----------
    selCh : bool array
        Channels to generate PRBS sequence for.
    minVal : array_like
        Minimum value of PRBS sequence for each channel.
    maxVal : array_like
        Maximum value of PRBS sequence for each channel.
    order : int
        Order of PRBS sequence.
    Returns
    -------
    prbs : array_like
        PRBS sequence.
    """
    print(selCh)
    prbs = np.zeros((3, 2**order-1))

    print('Generating PRBS sequence...')
    for i in range(len(selCh)):
        if selCh[i]:
            prbs_1_period = max_len_seq(order)[0] * (maxVal[i] - minVal[i]) - minVal[i]
            prbs[i, :] = prbs_1_period
    return prbs
