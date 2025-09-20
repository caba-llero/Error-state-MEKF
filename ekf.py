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
_d: difference / error between true and estimated value
_n: noise

dX: error of variable X, verifying X_t = X_h + dX

Z: rotation vector - (3,)
q: attitude quaternion - (4,)
B: gyro bias - (3,)
w: angular velocity - (3,)
s: standard deviations of error state, i.e. diag(P)**0.5 - (6,)

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
t_max = 20*90 # maximum integration time [s]
dt = 0.01 # integration timestep, for ground truth computing [s]
sigma_startracker = 6 # isotropic accuracy of startracker for each axis [arcsec]
sigma_v = 10**0.5 * 1e-6 # gyro angle random walk factor  [rad/s / sqrt(Hz) = rad/sqrt(s)]
sigma_u = 10**0.5 * 1e-9 # rate random walk coefficienct  [rad / s^(3/2)]
freq_startracker = 1 # frequency of startracker measurements [Hz]
freq_gyro = 10 # frequency of gyro measurements [Hz]
init_inaccuracy = 10
rng_seed = 1

Joseph = True # use Joseph formula to update the covariance matrix after startracker measurement. False = use simple form (perhaps more numerically unstable)
simple_Phi = False # use the simplified state transition matrix (small angle)

# Define initial values for estimates
B_h_0 =  np.array([0,0,0]) 
Pq = (6 * arcsec_to_rad)**2 * I3 # initial attitude error vector covariance [rad^2]
Pb = (0.2 * degh_to_rads)**2 * I3 # initial gyro bias error covariance [(rad/s)^2]


# Ground truth initial values
B_t_0 = np.array([0.1,0.1,0.1]) * degh_to_rads  # start with 0.1 deg/h bias
q_t_0 = np.array([1,0,0,1]) / 2**0.5 # initial attitude

# Define true angular velocity
def w_t_fun(t):
    w1 = 0.1*np.sin(0.01*t) * pi/180
    w2 = 0.1*np.sin(0.0085*t) * pi/180
    w3 = 0.1*np.cos(0.0085*t) * pi/180
    return np.vstack((w1, w2, w3))

### Automatic initialization (not user input)
H = np.hstack((I3, O3)) # H = [I_3 0_3x3]
R = I3 * (sigma_startracker*arcsec_to_rad)**2 
Q_init = u.Q(sigma_v, sigma_u, dt, I3)

times = np.arange(0, t_max, dt)

idx_gyro = u.measurement_indices(t_max, dt, freq_gyro)
idx_star = u.measurement_indices(t_max, dt, freq_startracker)
idx_all = idx_gyro | idx_star  # union of both sets
timesteps = len(idx_all)  # log at every measurement event

n = 3 # number of elements in state (3 for dZ, 3 for dB)
w_t_l = w_t_fun(times) # ground truth angular velocity 
rng = np.random.default_rng(seed=rng_seed)

# Initial estimate for the attitude is a very noisy measurement
Z_n = rng.normal(0, sigma_startracker*arcsec_to_rad*init_inaccuracy, n).reshape(-1,1) # noise on each axis
q_m_0 = q_t_0.reshape(-1,1) + 0.5 * u.Xi(q_t_0) @ Z_n # small angle approximation of q_m = q_n ⊗ q_t
q_m_0 = q_m_0.flatten()
q_h_0 = q_m_0

# Empty arrays for logging variables
s_l = np.empty((6,timesteps))
q_h_l = np.empty((4,timesteps))
B_h_l = np.empty((3,timesteps))
q_t_l = np.empty((4,timesteps))
Z_d_l = np.empty((3,timesteps))
B_t_l = np.empty((3,timesteps))
t_l = np.empty(timesteps)
G_l = np.empty(timesteps) # pointing error

# Set initial values
q_t = q_t_0
B_t = B_t_0
B_h = B_h_0
q_h = q_h_0
q_d = u.quat_mul(q_t, u.quat_inv(q_h))
Z_d = u.quat_to_rotvec(q_d)
G = np.linalg.norm(Z_d)
P = np.block([[Pq, O3],[O3, Pb]])


"""
Iterate over the full truth time grid and log only at measurement events.
Use a separate log index k (0..timesteps-1) so log arrays match event count.
"""
k = 0

# Optionally log initial state if an event occurs at t=0
if 0 in idx_all and k < timesteps:
    s = np.sqrt(np.diag(P))
    s_l[:,k] = s
    t_l[k] = times[0]
    q_h_l[:,k] = q_h
    q_t_l[:,k] = q_t
    Z_d_l[:,k] = Z_d
    B_h_l[:,k] = B_h.flatten()
    B_t_l[:,k] = B_t
    k += 1

last_gyro_i = 0
for i in range(1, len(times)):
    # propagate ground truth of quaternion and bias to next time step
    w_t = w_t_l[:, i-1]
    q_t = u.quat_propagate(q_t, w_t, dt)
    B_t = B_t + rng.normal(0, sigma_u*dt**0.5, n)

    # propagate estimate on gyro event
    if i in idx_gyro:
        # elapsed time since last gyro event
        dt_g = times[i] - times[last_gyro_i]
        if dt_g <= 0:
            dt_g = dt
        # use instantaneous true rate at event time for simulated measurement
        w_t_meas = w_t_l[:, i] if i < w_t_l.shape[1] else w_t_l[:, -1]
        w_m = w_t_meas + B_t + np.random.standard_normal(n) * (sigma_v/np.sqrt(dt_g))
        w_h = w_m - B_h
        Phi = u.Phi(dt_g, w_h, I3, simple_Phi)
        Qk = u.Q(sigma_v, sigma_u, dt_g, I3)
        P = u.P_prop(P, Phi, Qk)
        q_h = u.quat_propagate(q_h, w_h, dt_g)   # propagate estimated attitude with gyro measurement
        last_gyro_i = i

    # update on star tracker event
    if i in idx_star:
        dZ_m = u.startracker_meas(q_t, q_h, sigma_startracker*arcsec_to_rad, n)
        K, K_Z, K_B = u.K(P, H, R)
        P = u.P_meas(K, H, P, R, Joseph)
        dB_h = K_B @ dZ_m
        dZ_h = K_Z @ dZ_m
        B_h = B_h + dB_h
        # Inject estimated attitude error (rotation vector) via exact quaternion mapping
        theta = np.linalg.norm(dZ_h)
        if theta > 0:
            axis = dZ_h / theta
            dq_err = np.hstack((axis * np.sin(0.5*theta), np.cos(0.5*theta)))
        else:
            dq_err = np.array([0.0, 0.0, 0.0, 1.0])
        q_h = u.quat_mul(dq_err, q_h)
        q_h = q_h / np.linalg.norm(q_h)
        q_d = u.quat_mul(q_t, u.quat_inv(q_h))
        Z_d = u.quat_to_rotvec(q_d)
        G = np.linalg.norm(Z_d) # pointing error

    # log at measurement events
    if i in idx_all and k < timesteps:
        s = np.sqrt(np.diag(P))
        s_l[:,k] = s
        t_l[k] = times[i]
        G_l[k] = G
        q_h_l[:,k] = q_h
        q_t_l[:,k] = q_t
        Z_d_l[:,k] = Z_d
        B_h_l[:,k] = B_h.flatten()
        B_t_l[:,k] = B_t
        k += 1


## Calculate errors
B_d = B_t_l - B_h_l # bias error

# Plot pointing error
plt.figure(figsize=(10, 6))
plt.plot(t_l, G_l)
plt.title("Total Pointing Error")
plt.ylabel("Error (rad)")
plt.xlabel("Time (s)")
plt.grid(True)

# Create plots for error components with 3-sigma bounds
fig, axs = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
fig.suptitle('Attitude and Gyro Bias Estimation Errors with 3-Sigma Bounds', fontsize=16)

# Attitude errors (Z_d_l)
for i in range(3):
    ax = axs[0, i]
    ax.plot(t_l, Z_d_l[i, :], 'b', label=f'Error Axis {i+1}')
    ax.plot(t_l, 3 * s_l[i, :], 'r--', label='3-sigma')
    ax.plot(t_l, -3 * s_l[i, :], 'r--')
    ax.set_title(f'Attitude Error (Z_d) - Component {i+1}')
    ax.set_ylabel('Error (rad)')
    ax.grid(True)
    ax.legend()

# Bias errors (B_d)
for i in range(3):
    ax = axs[1, i]
    ax.plot(t_l, B_d[i, :], 'b', label=f'Error Axis {i+1}')
    ax.plot(t_l, 3 * s_l[i+3, :], 'r--', label='3-sigma')
    ax.plot(t_l, -3 * s_l[i+3, :], 'r--')
    ax.set_title(f'Gyro Bias Error (B_d) - Component {i+1}')
    ax.set_ylabel('Error (rad/s)')
    ax.grid(True)
    ax.legend()

for ax in axs.flat:
    ax.set_xlabel('Time (s)')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()







