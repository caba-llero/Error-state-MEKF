import numpy as np
import matplotlib.pyplot as plt 
import utils as u
import importlib
import os
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

def _env_bool(name: str, default: bool) -> bool:
    v = os.getenv(name)
    if v is None:
        return default
    return v.strip().lower() in ("1", "true", "yes", "y", "on")

def _env_float(name: str, default: float) -> float:
    try:
        v = os.getenv(name)
        return float(v) if v is not None else default
    except Exception:
        return default

def _env_int(name: str, default: int) -> int:
    try:
        v = os.getenv(name)
        return int(float(v)) if v is not None else default
    except Exception:
        return default

# Initialization (overridable via environment variables)
dt = _env_float("MEKF_DT", 0.01)  # integration timestep [s]
t_max = _env_float("MEKF_TMAX", 20*60)  # maximum integration time [s]
sigma_startracker = _env_float("MEKF_SIGMA_ST", 6)  # [arcsec]
sigma_v = _env_float("MEKF_SIGMA_V", 10**0.5 * 1e-6)  # [rad/s / sqrt(Hz)]
sigma_u = _env_float("MEKF_SIGMA_U", 10**0.5 * 1e-9)  # [rad / s^(3/2)]
freq_startracker = _env_float("MEKF_FREQ_STAR", 1)  # [Hz]
freq_gyro = _env_float("MEKF_FREQ_GYRO", 10)  # [Hz]
init_inaccuracy = _env_float("MEKF_INIT_INACC", 10)
rng_seed = _env_int("MEKF_SEED", 1)

# Simplification flags for performance optimization (overridable)
Joseph = _env_bool("MEKF_JOSEPH", True)
simple_Phi = _env_bool("MEKF_SIMPLE_PHI", True)
simple_quaternions = _env_bool("MEKF_SIMPLE_QUAT", False)
simple_attitude_update = _env_bool("MEKF_SIMPLE_ATT_UPD", True)
simple_bias_update = _env_bool("MEKF_SIMPLE_BIAS_UPD", False)
enable_plots = _env_bool("MEKF_PLOTS", True)

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

# Measurement event indices and fast boolean flags
idx_gyro = u.measurement_indices(t_max, dt, freq_gyro)
idx_star = u.measurement_indices(t_max, dt, freq_startracker)
is_gyro = np.zeros(len(times), dtype=bool)
is_star = np.zeros(len(times), dtype=bool)
if len(idx_gyro) > 0:
    is_gyro[list(idx_gyro)] = True
if len(idx_star) > 0:
    is_star[list(idx_star)] = True
is_event = is_gyro | is_star
timesteps = int(is_event.sum())  # log at every measurement event

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
q_d = u.quat_mul(q_t, u.quat_inv(q_h), simple=False)
Z_d = u.quat_to_rotvec(q_d)
G = np.linalg.norm(Z_d)
P = np.block([[Pq, O3],[O3, Pb]])


"""
Iterate over the full truth time grid and log only at measurement events.
Use a separate log index k (0..timesteps-1) so log arrays match event count.
"""
k = 0

# Calculate event counts directly from index sets (more efficient)
ground_truth_steps = len(times)  # Total number of ground truth propagation steps
gyro_events = len(idx_gyro)      # Total number of gyro measurement events
startracker_events = len(idx_star)  # Total number of startracker measurement events
total_measurement_events = int(is_event.sum())  # Total number of measurement events (gyro + startracker)

# Start timing
import time
warmup = _env_bool("MEKF_WARMUP", True)
if warmup:
    try:
        u.warmup_jit()
    except Exception:
        pass
start_time = time.perf_counter()

# Optionally log initial state if an event occurs at t=0
if is_event[0] and k < timesteps:
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
gyro_dt_nominal = 1.0 / freq_gyro if freq_gyro > 0 else dt
Q_nominal = u.Q(sigma_v, sigma_u, gyro_dt_nominal, I3)

for i in range(1, len(times)):
    # propagate ground truth of quaternion and bias to next time step
    w_t = w_t_l[:, i-1]
    q_t = u.quat_propagate(q_t, w_t, dt, simple=simple_quaternions)
    B_t = B_t + rng.normal(0, sigma_u*dt**0.5, n)

    # propagate estimate on gyro event
    if is_gyro[i]:
        # elapsed time since last gyro event
        dt_g = times[i] - times[last_gyro_i]
        if dt_g <= 0:
            dt_g = dt
        # use instantaneous true rate at event time for simulated measurement
        w_t_meas = w_t_l[:, i] if i < w_t_l.shape[1] else w_t_l[:, -1]
        w_m = w_t_meas + B_t + np.random.standard_normal(n) * (sigma_v/np.sqrt(dt_g))
        w_h = w_m - B_h
        Phi = u.Phi(dt_g, w_h, I3, simple_Phi)
        # Reuse Q for nominal gyro interval if applicable
        if abs(dt_g - gyro_dt_nominal) < 1e-12:
            Qk = Q_nominal
        else:
            Qk = u.Q(sigma_v, sigma_u, dt_g, I3)
        P = u.P_prop(P, Phi, Qk)
        q_h = u.quat_propagate(q_h, w_h, dt_g, simple=simple_quaternions)   # propagate estimated attitude with gyro measurement
        last_gyro_i = i

    # update on star tracker event
    if is_star[i]:
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
            if simple_attitude_update:
                # Simplified attitude update for small angles
                dq_err = np.array([dZ_h[0]*0.5, dZ_h[1]*0.5, dZ_h[2]*0.5, 1.0])
            else:
                dq_err = np.hstack((axis * np.sin(0.5*theta), np.cos(0.5*theta)))
        else:
            dq_err = np.array([0.0, 0.0, 0.0, 1.0])
        q_h = u.quat_mul(dq_err, q_h, simple=simple_quaternions)
        q_h = q_h / np.linalg.norm(q_h)
        q_d = u.quat_mul(q_t, u.quat_inv(q_h), simple=False)
        Z_d = u.quat_to_rotvec(q_d)
        G = np.linalg.norm(Z_d) # pointing error

    # log at measurement events
    if is_event[i] and k < timesteps:
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

# End timing and print benchmarks
end_time = time.perf_counter()
total_time = end_time - start_time

print("=" * 60)
print("EKF SIMULATION BENCHMARKS")
print("=" * 60)
print(f"Ground truth propagations: {ground_truth_steps}")
print(f"Gyro measurement events: {gyro_events}")
print(f"Startracker measurement events: {startracker_events}")
print(f"Total measurement events: {total_measurement_events}")
print("-" * 40)
print(f"Total computation time: {total_time:.4f} seconds")
print(f"Average time per event: {total_time/total_measurement_events*1000:.2f} ms/event")
print(f"Average time per ground truth step: {total_time/ground_truth_steps*1000:.2f} ms/step")
print("-" * 40)
print("Simplification settings:")
print(f"  Simple quaternions: {simple_quaternions}")
print(f"  Simple attitude update: {simple_attitude_update}")
print(f"  Simple bias update: {simple_bias_update}")
print(f"  Simple Phi matrix: {simple_Phi}")
print(f"  Joseph form: {Joseph}")
print("-" * 40)
print("Runtime configuration:")
print(f"  t_max: {t_max} s, dt: {dt} s")
print(f"  freq_gyro: {freq_gyro} Hz, freq_startracker: {freq_startracker} Hz")
print(f"  plots enabled: {enable_plots}")
print("=" * 60)

## Calculate errors
B_d = B_t_l - B_h_l # bias error

# Plot pointing error
if enable_plots:
    try:
        plt.figure(figsize=(10, 6))
        plt.plot(t_l, G_l)
        plt.title("Total Pointing Error (Optimized Version)")
        plt.ylabel("Error (rad)")
        plt.xlabel("Time (s)")
        plt.grid(True)
        plt.show()
    except:
        print("Could not create pointing error plot")

# Create plots for error components with 3-sigma bounds
if enable_plots:
    try:
        fig, axs = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
        fig.suptitle('Attitude and Gyro Bias Estimation Errors with 3-Sigma Bounds (Optimized)', fontsize=16)

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
    except Exception as e:
        print(f"Could not create error components plot: {e}")







