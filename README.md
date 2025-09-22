# Error state MEKF 

An error state multiplicative extended Kalman filter, specifically for spacecraft attitude and gyro bias estimation. Loosely based on the EKF algorithms found in:

> [1] Fundamentals of Spacecraft Attitude Determination and Control, by Markley & Crassidis (2014)
> 
> [2] Quaternion kinematics for the error-state Kalman filter, by Joan Solà (2017)

We assume that the sensors available are a startracker and gyros on each axis. The estimator synthesizes startracker and gyro measurements based on the inputted noise parameters. The gyro measurement model assumes ARW and RRW; the startracker measurement model applies isotropic Gaussian noise to each axis. The user provides a predefined ground truth angular velocity in fucntion of time. Ground truth is propagated at an independent timestep to the kinematics equations, driven by gyro measurements (gyro measurement for dynamic model replacement is used). The user can plot pointing accuracy plots, estimation versus ground truth (error) plots for the gyro bias and attitude. 

The full documentation and theory behind the algorithm can be found in https://error-state-mekf.readthedocs.io/en/latest/theory.html

This algorithm is used in Pytellite.org to simulate the attitude estimation in a satellite (software-in-the-loop style). Check out https://github.com/caba-llero/pytellite for more info.

Installation
------------

```bash
pip install git+https://github.com/caba-llero/error-state-mekf.git
```

Basic usage
-----------

```python
from mekf import MEKF

mekf = MEKF(inputs={
    "t_max": 1000,
    "dt": 0.01,
    "sigma_startracker": 6,
    "sigma_v": 1e-6,
    "sigma_u": 1e-9,
    "freq_startracker": 1,
    "freq_gyro": 20,
})

results = mekf.calculate()
mekf.plot_all()
```
