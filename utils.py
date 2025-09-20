import numpy as np
from numba import njit

I3 = np.eye(3)

@njit
def Phi(dt, w_h, simple=False): # returns state transition matrix (discrete-time), i.e. Phi = Exp(F dt)
    e = np.norm(w_h)
    p = e * dt
    a = skew(w_h)
    sp = np.sin(p)
    cp = np.cos(p)

    if simple:
        Phi11 = I3 - a*dt
        Phi12 = -I3 * dt
    else:
        Phi11 = I3 - a/e * sp + a**2 / e**2 * (1-cp)
        Phi12 = -I3 * dt - a**2 / e**3 * (p - sp) + a/e**2 * (1 - cp)

    Phi21 = np.zeros((3,3))
    Phi22 = np.eye(3)

    Phi = np.block([[Phi11, Phi12],
                    [Phi21, Phi22]])

    return Phi

@njit
def skew(a): # returns skew-symmetric matrix of column vector x
    a = a.flatten()
    a_x = np.array([[0, -a[3], a[2]],
                    [a[3], 0, -a[1]],
                    [-a[2], a[1], a[0]] ])
    return a_x


def Q(sigma_v, sigma_u, dt): # Discrete time noise covariance matrix (Eq. 6.93)
    Q11 = (sigma_v**2 * dt + sigma_u**2 * dt**3 / 3) * I3
    Q12 = - sigma_u**2 * dt**2 * I3 / 2
    Q22 = sigma_u**2 * dt * I3
    Q = np.block([[Q11, Q12], [Q12, Q22]])
    return Q

@njit
def K(P, H, R): # Kalman gain
    # B = P H'
    # S = H B + R
    # K = B inv(S)
    # MOre efficient way to solve: Ax = B --> np.linalg.solve(A,B)
    # K S = B; S' K' = B'; K = np.linalg.solve(S', B').T

    B = P @ H.T
    S = H @ B + R
    K = np.linalg.solve(S.T, B.T).T
    return K

@njit
def P_meas(K, H, P, R, Joseph = True): # returns the update covariance matrix, after measurement
    KH = K @ H
    IKH = I3 - KH

    if Joseph: 
        return IKH @ P @ IKH.T + K @ R @ K.T
    else:
        return IKH @ P

@njit 
def P_prop(P, Phi, Q): # returns the update of the covariance matrix, after propagation
    return Phi @ P @ Phi.T + Q


def measurement_indices(t_max, dt, measurement_freq): # returns the indices of t at which we trigger a measurement event
    n_steps = int(t_max / dt)
    expected_times = np.arange(0, t_max, 1/measurement_freq)
    indices = np.ceil(expected_times / dt).astype(int)  # Convert expected times to indices (round up to ensure >= condition)
    indices = indices[indices < n_steps]  # Filter out indices beyond array bounds
    return set(indices)  


