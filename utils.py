import numpy as np

def Phi(dt, w_h):
    pass


def skew(a): # returns skew-symmetric matrix of column vector x
    a = a.flatten()
    a_x = np.array([[0, -a[3], a[2]],
                    [a[3], 0, -a[1]],
                    [-a[2], a[1], a[0]] ])
    return a_x
