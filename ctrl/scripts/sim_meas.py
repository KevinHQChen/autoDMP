import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import matrix_rank as rank
import control as ct

def sim_meas(selCh, u, x, dt):
    """
    simulate a n-input, n-output unstable second order underdamped system
    Parameters
    ----------
    """
    numCh = selCh.count(True)

    # Define SISO second order + integrator transfer function
    zeta = 0.5
    wn = 10
    num = [wn**2]
    num1 = [wn**2, 2*zeta*wn]
    den = [1, 2*zeta*wn, wn**2, 0]

    if numCh == 1:
        sysc_tf = ct.tf(num, den)
    elif numCh == 2:
        sysc_tf = ct.tf([[num, num1],
                         [num1, num]],
                        [[den, den],
                         [den, den]])

    sysc = ct.minreal(ct.ss(sysc_tf))
    sysd = ct.sample_system(sysc, dt, method='zoh')

    # Simulate the system
    x = np.dot(sysd.A, x) + np.dot(sysd.B, u)
    y = np.dot(sysd.C, x) + np.dot(sysd.D, u)

    y_ = np.zeros((3, 1))

    for i in range(selCh):
        if selCh[i]:
            y_[i] = y if numCh == 1 # else y_[i] = y[i]






plt.figure()
plt.subplot(1, 2, 1)
plt.plot(t, u[:, 0], label='input 1')
plt.plot(t, y[:, 0], label='output 1')
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(t, u[:, 1], label='input 2')
plt.plot(t, y[:, 1], label='output 2')
plt.legend()
plt.show()
