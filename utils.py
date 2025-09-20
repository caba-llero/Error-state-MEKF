import numpy as np
from numba import njit

I3 = np.eye(3)

def Phi(dt, w_h, simple=False, min_e=1e-7): # returns state transition matrix (discrete-time), i.e. Phi = Exp(F dt)
    e = np.linalg.norm(w_h)
    if e < min_e: # to avoid numerical instabilities as we divide by e
        simple=True

    p = e * dt
    a = skew(w_h)
    sp = np.sin(p)
    cp = np.cos(p)

    if simple:
        Phi11 = I3 - a*dt
        Phi12 = -I3 * dt
    else:
        Phi11 = I3 - a/e * sp + a@a / e**2 * (1-cp)
        Phi12 = -I3 * dt - a@a / e**3 * (p - sp) + a/e**2 * (1 - cp)

    Phi21 = np.zeros((3,3))
    Phi22 = np.eye(3)

    Phi = np.block([[Phi11, Phi12],
                    [Phi21, Phi22]])

    return Phi

def skew(a): # returns skew-symmetric matrix of column vector x
    a = a.flatten()
    a_x = np.array([[0, -a[2], a[1]],
                    [a[2], 0, -a[0]],
                    [-a[1], a[0], 0] ])
    return a_x


def Q(sigma_v, sigma_u, dt): # Discrete time noise covariance matrix (Eq. 6.93)
    Q11 = (sigma_v**2 * dt + sigma_u**2 * dt**3 / 3) * I3
    Q12 = - sigma_u**2 * dt**2 * I3 / 2
    Q22 = sigma_u**2 * dt * I3
    Q = np.block([[Q11, Q12], [Q12, Q22]])
    return Q

def K(P, H, R): # Kalman gain
    # B = P H'
    # S = H B + R
    # K = B inv(S)
    # MOre efficient way to solve: Ax = B --> np.linalg.solve(A,B)
    # K S = B; S' K' = B'; K = np.linalg.solve(S', B').T

    B = P @ H.T
    S = H @ B + R
    K = np.linalg.solve(S.T, B.T).T
    K_Z = K[:3, :] # Z part of the gain
    K_B = K[3:, :] # B part of the gain
    return K, K_Z, K_B

def P_meas(K, H, P, R, Joseph = True): # returns the update covariance matrix, after measurement
    KH = K @ H
    IKH = np.eye(6) - KH

    if Joseph: 
        return IKH @ P @ IKH.T + K @ R @ K.T
    else:
        return IKH @ P

def P_prop(P, Phi, Q): # returns the update of the covariance matrix, after propagation
    return Phi @ P @ Phi.T + Q


def measurement_indices(t_max, dt, measurement_freq): # returns the indices of t at which we trigger a measurement event
    n_steps = int(t_max / dt)
    expected_times = np.arange(0, t_max, 1/measurement_freq)
    indices = np.ceil(expected_times / dt).astype(int)  # Convert expected times to indices (round up to ensure >= condition)
    indices = indices[indices < n_steps]  # Filter out indices beyond array bounds
    return set(indices)  

def quat_mul(q1, q2):     # Quaternion product q1 ⊗ q2
    x1,y1,z1,w1 = q1
    x2,y2,z2,w2 = q2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    return np.array([x,y,z,w])

def quat_propagate(q, w, dt, p_min = 1e-7): 
    # propagates initial quaternion q thru angular vel w and delta time dt
    # done with q <-- dq ⊗ q
    dz = w*dt
    p = np.linalg.norm(dz)
    if p < p_min: # return the same quaternion to avoid numerical instabilities
        return q
    else:
        e = dz / p
        dq = np.hstack((e * np.sin(p/2), np.cos(p/2)))
        q = quat_mul(dq, q)
        return q


def Xi(q):
    return np.array([ [q[3], -q[2], q[1]],
                      [q[2], q[3], -q[0]],
                      [-q[1], q[0], q[3]],
                      [-q[0], -q[1], -q[2]]  
                         ])


def startracker_meas(q_t, q_h, sigma, rng, n):
    Z_n = rng.normal(0, sigma, n).reshape(-1,1) # noise on each axis
    q_m = q_t.reshape(-1,1) + 0.5 * Xi(q_t) @ Z_n # small angle approximation of q_m = q_n ⊗ q_t
    q_m = q_m.flatten()
    dq_m = quat_mul(q_m, quat_inv(q_h))
    dZ_m = 2 * dq_m[:3] / dq_m[3]
    return dZ_m

def quat_inv(q):
    return np.array([-q[0], -q[1], -q[2], q[3]])