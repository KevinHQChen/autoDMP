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

# Nope, std_dev is basically 1, so this does nothing. By the way, I've tested my code with another random sequence as follows:
# u = []
# for i in range(1023):
#    n = random.randint(-1,1)
#    u.append(n)
# var = np.var(u)
# mean = np.mean(u)
# nu = u - mean
# corr = np.correlate(u, u, mode='full')[len(u)-1:]
# corr = corr / var / len(nu)
# plt.plot(corr)
# plt.show()

# And this gives me the expected peak at 0 lag, and close to 0 everywhere else. But I expected a similar result with the max_len_seq, and I'm not seeing that at all
