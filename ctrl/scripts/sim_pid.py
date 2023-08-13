import matplotlib.pyplot as plt
import numpy as np
import control as ct

# Define an unstable second order underdamped system
dt = 0.025
t = np.arange(0, 10, dt)
zeta = 0.5
wn = 10
num = [wn**2]
den = [1, 2*zeta*wn, wn**2]
sysc = ct.tf(num, den)
sys = ct.sample_system(sysc, dt, method='zoh')

# Analyze the system (time and frequency domain)
t, y = ct.step_response(sys, t)
f, mag, phase = ct.bode(sys, dB=True, Hz=True, omega=np.logspace(-2, 2, 1000))
# Plot the results
# Time domain
plt.figure()
plt.plot(t, y)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Step Response')
# Frequency domain
# Bode plot
plt.figure()
ct.bode_plot(sys)
plt.show()
# Pole-zero plot
plt.figure()
ct.pzmap(sys)
plt.show()
# Nyquist plot
plt.figure()
ct.nyquist_plot(sys)
plt.show()

# Define a PID controller
kp = 200
ki = 10
kd = 5
pidc = ct.tf([kd, kp, ki], [1, 0])
# pid = ct.sample_system(pid, dt, method='zoh')

# Connect the PID controller to the system
sys_clc = ct.feedback(pidc*sysc)
sys_cl = ct.sample_system(sys_clc, dt, method='zoh')
sys_cl_ss = ct.tf2ss(sys_cl)
min_sys = ct.minreal(sys_cl_ss)

# Analyze the closed loop system (time and frequency domain)
# Time domain
# Step response
t, y = ct.step_response(sys_cl, t)
plt.figure()
plt.plot(t, y)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Step Response')
# Frequency domain
# Bode plot
plt.figure()
ct.bode_plot(sys_cl)
plt.show()
# Pole-zero plot
plt.figure()
ct.pzmap(sys_cl)
plt.show()
# Nyquist plot
plt.figure()
ct.nyquist_plot(sys_cl)
plt.show()

# Simulate the closed loop system
# Define the input
# Step input
u = np.ones_like(t)
# Ramp input
u = t
# Sinusoidal input
u = np.sin(2*np.pi*t)
# Square input
u = np.sign(np.sin(2*np.pi*t))
# Define the initial conditions
x0 = [0, 0, 0]
yout = [0]

# Iteratively simulate the system
for i in range(0, len(t)):
    t_i = t[i]
    u_i = u[i]
    if i == 0:
        t_i = [0, dt]
        u_i = [0, u_i]
    else:
        # t_i = [t[i-1], t[i]]
        t_i = [round(t[i-1], 3), round(t[i], 3)]
        u_i = [u[i-1], u[i]]
    response = ct.forced_response(min_sys, t_i, u_i, x0)
    # x0 = response[1][-1]
    x0 = response[2].T[-1]
    yout = np.append(yout, response[1][-1])

# Plot the results
plt.figure()
plt.plot(t, yout[1:])
plt.plot(t, u)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.show()

# Simulate the system
# Time domain
# Step response
response = ct.forced_response(min_sys, t, u, x0)
t = response[0]
y = response[1]
plt.figure()
plt.plot(t, y)
plt.plot(t, u)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Step Response')

# Frequency domain
# Bode plot
plt.figure()
ct.bode_plot(sys_cl)
plt.show()
# Pole-zero plot
plt.figure()
ct.pzmap(sys_cl)
plt.show()
# Nyquist plot
plt.figure()
ct.nyquist_plot(sys_cl)
plt.show()

