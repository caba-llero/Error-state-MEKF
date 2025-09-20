import numpy as np
import matplotlib.pyplot as plt 
import utils as u
import importlib
importlib.reload(u)

'''
Notation

_h: "hat" accent, denotes estimated variable
_g: ground truth value of a variable (unknown to estimator)
_0: initial value of a variable
_m: measured variable (after adding noise)
_l: logged variable, saved for analysis/plotting

dX: error of variable X, verifying X_t = X_h + dX

Z: rotation vector - (3,1)
q: attitude quaternion - (4,1)
B: gyro bias - (3,1)
w: angular velocity - (3,1)
s: standard deviations of error state, i.e. diag(P)**0.5 - (6,1)

'''

######### 
# Constants
I3 = np.eye(3)
O3 = np.zeros((3,3))
pi = np.pi
arcsec_to_rad = pi / (180 * 3600)
degh_to_rads = pi / (180 * 3600)

#########

# Initialization
t_max = 60*90 # maximum integration time [s]
dt = 0.1 # integration timestep, for ground truth computing [s]
sigma_startracker = 6 # isotropic accuracy of startracker for each axis [arcsec]
sigma_v = 10**0.5 * 1e-7 # gyro angle random walk factor  [rad/s / sqrt(Hz) = rad/sqrt(s)]
sigma_u = 10**0.5 * 1e-10 # rate random walk coefficienct  [rad / s^(3/2)]
freq_startracker = 2 # frequency of startracker measurements [Hz]
freq_gyro = 100 # frequency of gyro measurements [Hz]

rng_seed = 1

Joseph = True # use Joseph formula to update the covariance matrix after startracker measurement. False = use simple form (perhaps more numerically unstable)

# Define initial values for estimates
q_h_0 = np.array([0,0,0,1])
B_h_0 =  np.array([0,0,0]) 
Pq = (6 * arcsec_to_rad)**2 * I3 # initial attitude error vector covariance [rad^2]
Pb = (0.2 * arcsec_to_rad)**2 * I3 # initial gyro bias error covariance [(rad/2)^2]


# Ground truth initial values
B_t_0 = np.array([0.1,0.1,0.1]) * degh_to_rads  # start with 0.1 deg/h bias
q_t_0 = np.array([1,0,0,1]) / 2**0.5 # initial attitude

# Define true angular velocity
def w_t_fun(t):
    w1 = 0.1*np.sin(0.01*t)
    w2 = 0.1*np.sin(0.0085*t)
    w3 = 0.001*np.cos(0.0085*t)
    return np.vstack((w1, w2, w3))

### Automatic initialization (not user input)
H = np.hstack((I3, O3)) # H = [I_3 0_3x3]
R = I3 * (sigma_startracker*arcsec_to_rad)**2 
Q = u.Q(sigma_v, sigma_u, dt)

times = np.arange(0, t_max, dt)

idx_gyro = u.measurement_indices(t_max, dt, freq_gyro)
idx_star = u.measurement_indices(t_max, dt, freq_startracker)
idx_all = idx_gyro | idx_star  # union of both sets
timesteps = len(idx_all)  # log at every measurement event

n = 3 # number of elements in state (3 for dZ, 3 for dB)
w_t_l = w_t_fun(times) # ground truth angular velocity 
rng = np.random.default_rng(seed=rng_seed)

# Empty arrays for logging variables
s_l = np.empty((6,timesteps))
q_h_l = np.empty((4,timesteps))
B_h_l = np.empty((3,timesteps))
q_t_l = np.empty((4,timesteps))
Z_d_l = np.empty((3,timesteps))
B_t_l = np.empty((3,timesteps))
t_l = np.empty(timesteps)

# Set initial values
q_t = q_t_0
B_t = B_t_0
B_h = B_h_0
q_h = q_h_0
q_d = u.quat_mul(q_t, u.quat_inv(q_h))
Z_d = 2 * q_d[:3] / q_d[3]
w_t = w_t_l[:,0]
P = np.block([[Pq, O3],[O3, Pb]])


for idx in range(1, timesteps):
    # propagate ground truth of quaternion and bias
    q_t = u.quat_propagate(q_t, w_t, dt)   
    B_t = B_t + rng.normal(0, sigma_v*dt**0.5, n) 

    if idx in idx_gyro: # propagate estimate of angular velocity and cov matrix
        w_m = w_t + B_t + rng.normal(0, sigma_u, n) # angular velocity true measurement
        w_h = w_m - B_h # angular velocity estimate
        Phi = u.Phi(dt, w_h) # state transition matrix
        P = Phi @ P @ P.T + Q

    if idx in idx_star: # update estimate, propagate covariance, inject to global state
        dZ_m = u.startracker_meas(q_t, sigma_startracker, rng, n) # obtain startracker measurement in form of a error vector
        K, K_Z, K_B = u.K(P, H, R) # calculate Kalman gain
        P = u.P_meas(K, H, P, R, Joseph) # propagate covariance
        dB_h = K_Z @ dZ_m # obtain estimate of bias error
        dZ_h = K_Z @ dZ_m # obtain estiamte of error vector
        B_h = B_h + dB_h # inject estimate of bias error into estimate of bias
        q_h = q_h.reshape(-1,1) + 0.5 * u.Xi(q_h) @ dB_h.reshape(-1,1) # inject estimate of error vector into estimate of atttitude quaternion
        q_h = q_h / np.linalg.norm(q_h) # second step of injection: normalization
        q_h = q_h.flatten()
        q_d = u.quat_mul(q_t, u.quat_inv(q_h))
        Z_d = 2 * q_d[:3] / q_d[3]

    if idx in idx_all: # log measurement
        s = np.sqrt( np.diag(P) ) # obtain std of each component 
        s_l[:,idx] = s 
        t_l[idx] = times[idx]
        q_h_l[:,idx] = q_h
        q_t_l[:,idx] = q_t
        Z_d_l[:,idx] = Z_d
        B_h_l[:,idx] = B_h.flatten()
        B_t_l[:,idx] = B_t


## Calculate errors
dB = B_t_l - B_h_l # bias error

plt.plot(t_l, Z_d_l[0,:])







