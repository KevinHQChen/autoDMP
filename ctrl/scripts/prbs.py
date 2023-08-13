import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import max_len_seq, periodogram
from numpy.fft import fft, ifft, fftshift, fftfreq

# import matplotlib
# matplotlib.use('qtagg')

# import sysid
# p = sysid.prbs.prbs([True, True, True], [-1, -1, -1], [1, 1, 1], 10)

def prbs(minVal, maxVal, order):
    """
    function that generates a PRBS sequence for a single channel.
    Parameters
    ----------
    minVal : float_array_like
        Minimum value of PRBS sequence for each channel.
    maxVal : float_array_like
        Maximum value of PRBS sequence for each channel.
    order : int
        Order of PRBS sequence.
    Returns
    -------
    prbs : array_like
        PRBS sequence.
    """

    # select taps for each pump to ensure each maximum length sequence is uncorrelated
    # (max_len_seq skips first tap; it is always equal to the order)
    taps = [
        (7),           # 2 taps
        (9, 8, 5),     # 4 taps
        (9, 8, 7, 5, 4) # 6 taps
    ]
    print('Generating PRBS sequence...')

    prbs = []
    for i, tap in enumerate(taps):
        prbs_1d = max_len_seq(order, taps=tap)[0] * (maxVal[i] - minVal[i]) + minVal[i]
        prbs.append(prbs_1d)

    # convert to numpy array
    prbs = np.array(prbs)

    return prbs
