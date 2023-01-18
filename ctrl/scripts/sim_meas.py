import matplotlib.pyplot as plt
import numpy as np
import control
from control import tf, ss, mixsyn, step_response

# Define the plant
dt = 0.025
den = [0.2, 1.2, 1]
gtf = tf([[[1], [1]],
          [[2, 1], [2]]],
         [[den, den],
          [den, den]])
gss = ss(gtf, dt=dt)

# Analyze the plant
print('Plant transfer function:')
print(gtf)
print('Plant state space:')
print(gss)

# Analyze step response
t, y = step_response(gss, T=np.arange(0, 10, dt), squeeze=True)
# Plot the step response on 4 subplots with legends
plt.figure()
plt.subplot(2, 2, 1)
plt.plot(t, y[0, 0, :], label='input 1 to output 1')
plt.legend()
plt.subplot(2, 2, 2)
plt.plot(t, y[0, 1, :], label='input 1 to output 2')
plt.legend()
plt.subplot(2, 2, 3)
plt.plot(t, y[1, 0, :], label='input 2 to output 1')
plt.legend()
plt.subplot(2, 2, 4)
plt.plot(t, y[1, 1, :], label='input 2 to output 2')
plt.legend()
plt.show()






# Define the controller



class DiscreteMIMOStateSpaceSystem:
    def __init__(self, transfer_functions):
        self.sys = self.create_mimo_sys(transfer_functions)

    def create_mimo_sys(self, transfer_functions):
        """Convert SISO transfer functions to state-space and concatenate to create MIMO system"""
        A = []
        B = []
        C = []
        D = []
        for i in range(len(transfer_functions)):
            for j in range(len(transfer_functions[i])):
                sys = control.tf2ss(transfer_functions[i][j])
                A.append(sys.A)
                B.append(sys.B)
                C.append(sys.C)
                D.append(sys.D)
        A = np.concatenate(A, axis=1)
        B = np.concatenate(B, axis=1)
        C = np.concatenate(C, axis=0)
        D = np.concatenate(D, axis=0)
        return control.StateSpace(A, B, C, D)

    def simulate(self, inputs, initial_state=None):
        t, y, _ = control.forced_response(self.sys, T=inputs[0], U=inputs[1], X0=initial_state)
        return t, y

# class Model:
#     def __init__(self):
#         self.A = [[-1, -2], [-3, -4]]
#         self.B = [[1, 1], [1, 1]]
#         self.C = [[1, 0], [0, 1]]
#         self.D = [[0, 0], [0, 0]]
#         self.sys = control.ss(self.A, self.B, self.C, self.D)

# mimo_sys = Model()

# def sim_meas(u, y):
#     return val + 1
