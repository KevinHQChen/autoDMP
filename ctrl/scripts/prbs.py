import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import max_len_seq, periodogram
from numpy.fft import fft, ifft, fftshift, fftfreq

# import matplotlib
# matplotlib.use('qtagg')

# import sysid
# p = sysid.prbs.prbs([True, True, True], [-1, -1, -1], [1, 1, 1], 10)


def prbs(selCh, minVal, maxVal, order, flip=True):
    """
    function that generates a PRBS sequence.
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
    numCh = selCh.count(True)

    prbs = np.zeros((3, 2**order-1))

    print('Generating PRBS sequence...')
    for i in range(len(selCh)):
        if selCh[i]:
            prbs_1_period = max_len_seq(order)[0] * (maxVal[i] - minVal[i]) + minVal[i]
            if flip:
                prbs[i, :] = -prbs_1_period
            else:
                prbs[i, :] = prbs_1_period

    firstCh = True
    if numCh > 1:
        flipped_prbs = np.copy(prbs)
        for i in range(len(selCh)):
            if selCh[i]:
                if firstCh:
                    firstCh = False
                else:
                    flipped_prbs[i, :] = -prbs[i, :]
        prbs = np.concatenate((prbs, flipped_prbs), axis=1)

    return prbs
