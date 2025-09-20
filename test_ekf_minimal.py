import numpy as np
import time
import utils as u
import importlib
importlib.reload(u)

print("Starting minimal EKF test...")

# Test parameters
t_max = 1  # Very short test
dt = 0.01

# Constants
I3 = np.eye(3)
O3 = np.zeros((3,3))
pi = np.pi
arcsec_to_rad = pi / (180 * 3600)

# Initialize
times = np.arange(0, t_max, dt)
idx_gyro = u.measurement_indices(t_max, dt, 10)
idx_star = u.measurement_indices(t_max, dt, 1)
idx_all = idx_gyro | idx_star
timesteps = len(idx_all)

n = 3
rng = np.random.default_rng(seed=1)

# Event counters for benchmarking
event_counts = {
    'ground_truth': 0,
    'gyro_measurements': 0,
    'startracker_measurements': 0,
    'total_measurements': 0
}

# Start timing
start_time = time.perf_counter()

print("Running simulation loop...")

# Simple test loop
for i in range(1, len(times)):
    event_counts['ground_truth'] += 1
    if i in idx_gyro:
        event_counts['gyro_measurements'] += 1
    if i in idx_star:
        event_counts['startracker_measurements'] += 1
    if i in idx_all:
        event_counts['total_measurements'] += 1

# End timing
end_time = time.perf_counter()
total_time = end_time - start_time

print("=" * 60)
print("MINIMAL EKF SIMULATION BENCHMARKS")
print("=" * 60)
print(f"Ground truth propagations: {event_counts['ground_truth']","}")
print(f"Gyro measurement events: {event_counts['gyro_measurements']","}")
print(f"Startracker measurement events: {event_counts['startracker_measurements']","}")
print(f"Total measurement events: {event_counts['total_measurements']","}")
print("-" * 40)
print(f"Total computation time: {total_time:.4f} seconds")
print(f"Average time per event: {total_time/event_counts['total_measurements']*1000:.2f} ms/event")
print("=" * 60)
print("Test completed successfully!")
