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
pi = np.pi
arcsec_to_rad = pi / (180 * 3600)

#########

# Initialization
t_max = 100 # maximum integration time [s]
dt = 0.01 # integration timestep, for ground truth computing [s]
sigma_startracker = 6 # isotropic accuracy of startracker for each axis [arcsec]
sigma_v = 10**0.5 * 1e-7 # gyro angle random walk factor  [rad/s / sqrt(Hz) = rad/sqrt(s)]
sigma_u = 10**0.5 * 1e-10 # rate random walk coefficienct  [rad / s^(3/2)]
freq_startracker = 2 # frequency of startracker measurements [Hz]
freq_gyro = 100 # frequency of gyro measurements [Hz]

# Estimated initial values
q_h_0 = np.array([[0,0,0,1]]).T
B_h_0 = np.array([[0,0,0]]).T
Pq = (6 * arcsec_to_rad)**2 * I3 # initial attitude error vector covariance [rad^2]
Pb = (0.2 * arcsec_to_rad)**2 * I3 # initial gyro bias error covariance [(rad/2)^2]

dZ_h = np.array([[0,0,0]]).T
dB_h = np.array([[0,0,0]]).T 


# True initial values
B_t_0 = np.array([[0,0,0]]).T # start with no gyro bias
q_t_0 = np.array([[0,0,0,1]]).T # initial attitude

# Define true angular velocity
def w_t_fun(t):
    w1 = 0.1*np.sin(0.01*t)
    w2 = 0.1*np.sin(0.0085*t)
    w3 = 0.001*np.cos(0.0085*t)
    return np.hstack((w1, w2, w3))

### Automatic initialization (not user input)
H = np.hstack((I3, np.zeros((3,3)))) # H = [I_3 0_3x3]
R = I3 * (sigma_startracker*arcsec_to_rad)**2 
Q = u.Q(sigma_v, sigma_u, dt)

times = np.arange(0, t_max, dt)

idx_gyro = u.measurement_indices(t_max, dt, freq_gyro)
idx_star = u.measurement_indices(t_max, dt, freq_startracker)
idx_all = idx_gyro | idx_star  # union of both sets
timesteps = len(idx_all)  # log at every measurement event

n = 6 # number of elements in state (3 for dZ, 3 for dB)
w_t_l = w_t_fun(t)

# Empty arrays for logging variables
s_l = np.empty((6,timesteps))
q_h_l = np.empty((4,timesteps))
B_h_l = np.empty((3,timesteps))

q_t_l = np.empty((4,timesteps))
B_t_l = np.empty((3,timesteps))

for idx, t in enumerate(times):
    # propagate ground truth

    if idx in idx_gyro: # propagate estimate
        pass

    if idx in idx_star: # update estimate
        pass
    
    if idx in idx_all: # log measurement
        pass








