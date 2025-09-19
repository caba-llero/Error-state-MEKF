import numpy as np
import matplotlib.pyplot as plt 
import utils as u


'''
Notation

_h: "hat" accent, denotes estimated variable
_t: true value of a variable (unknown to estimator)
_0: initial value of a variable
dX: error of variable X, verifying X_t = X_h + dX
Z: rotation vector (3,1)
q: attitude quaternion (4,1)
B: gyro bias (3,1)
'''

######### 
# Constants
I3 = np.eye(3)
pi = np.pi

#########

# Initialization
dt = 0.01 # integration timestep [s]
sigma_startracker = 6 # isotropic accuracy of startracker for each axis [arcsec]
sigma_v = 10**0.5 * 1e-7 # gyro angle random walk factor  [rad/s / sqrt(Hz) = rad/sqrt(s)]
sigma_u = 10**0.5 * 1e-10 # rate random walk coefficienct  [rad / s^(3/2)]


# Estimated initial values
q_h = np.array([0,0,0,1])
Pq = (6/3600 * 180 / pi) * I3 # initial attitude error vector covariance [rad^2]
Pb = (0.2/3600 * 180 / pi) * I3 # initial gyro bias error covariance [(rad/2)^2]

dZ_h = np.array([[0,0,0]]).T
dB_h = np.array([[0,0,0]]).T 


# True initial values
B_0 = np.array([[0,0,0]]).T # start with no gyro bias


### Automatic initialization (not user input)
H = np.hstack((I3, np.zeros((3,3)))) # H = [I_3 0_3x3]
R = I3 * sigma_startracker**2 

# Discrete time noise covariance matrix (Eq. 6.93)
Q11 = (sigma_v**2 * dt + sigma_u**2 * dt**3 / 3) * I3
Q12 = - sigma_u**2 * dt**2 * I3 / 2
Q22 = sigma_u**2 * dt * I3
Q = np.block([[Q11, Q12], [Q12, Q22]])


q_t = q_0


