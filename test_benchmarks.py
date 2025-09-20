import numpy as np
import time
import utils as u
import importlib
importlib.reload(u)

# Test basic functionality
print("Testing basic benchmark functionality...")

# Test parameters
t_max = 10  # Very short test
dt = 0.01
freq_gyro = 10
freq_startracker = 1

# Constants
I3 = np.eye(3)
O3 = np.zeros((3,3))
pi = np.pi
arcsec_to_rad = pi / (180 * 3600)

# Initialize
times = np.arange(0, t_max, dt)
idx_gyro = u.measurement_indices(t_max, dt, freq_gyro)
idx_star = u.measurement_indices(t_max, dt, freq_startracker)
idx_all = idx_gyro | idx_star
timesteps = len(idx_all)

print(f"Test parameters:")
print(f"  Time steps: {len(times)}")
print(f"  Gyro events: {len(idx_gyro)}")
print(f"  Startracker events: {len(idx_star)}")
print(f"  Total measurement events: {timesteps}")

# Start timing
start_time = time.perf_counter()

# Simple test loop
for i in range(1, len(times)):
    # Simple operations
    a = np.random.random(3)
    b = np.sin(a)
    c = np.linalg.norm(b)

end_time = time.perf_counter()
total_time = end_time - start_time

print(f"Test timing:")
print(f"  Total time: {total_time:.4f} seconds")
print("Test completed successfully!")
