import control

class Model:
    def __init__(self):
        self.A = [[-1, -2], [-3, -4]]
        self.B = [[1, 1], [1, 1]]
        self.C = [[1, 0], [0, 1]]
        self.D = [[0, 0], [0, 0]]
        self.sys = control.ss(self.A, self.B, self.C, self.D)

mimo_sys = Model()

def sim_meas(u, y):
    return val + 1
