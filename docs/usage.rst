Usage
=====

Basic Usage
-----------

The error state MEKF can be used to estimate spacecraft attitude and gyro bias from star tracker and gyro measurements. Here's a basic example:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from error_state_mekf import MEKF

    # Create MEKF instance
    mekf = MEKF(
        sigma_startracker=6,  # star tracker accuracy [arcsec]
        sigma_v=1e-6,         # gyro angle random walk [rad/s/sqrt(Hz)]
        sigma_u=1e-9,         # gyro rate random walk [rad/s^1.5]
        freq_startracker=1,   # star tracker frequency [Hz]
        freq_gyro=10          # gyro frequency [Hz]
    )

    # Run simulation
    results = mekf.simulate()

    # Plot results
    mekf.plot_results(results)

The main simulation is contained in the `ekf.py` file, which demonstrates the complete attitude estimation process.

Configuration Parameters
-------------------------

**Sensor Parameters:**

* ``sigma_startracker``: Isotropic accuracy of star tracker measurements [arcsec]
* ``sigma_v``: Gyro angle random walk coefficient [rad/s / sqrt(Hz)]
* ``sigma_u``: Gyro rate random walk coefficient [rad/s^1.5]
* ``freq_startracker``: Frequency of star tracker measurements [Hz]
* ``freq_gyro``: Frequency of gyro measurements [Hz]

**Algorithm Parameters:**

* ``init_inaccuracy``: Initial attitude error multiplier for noisy initialization
* ``Joseph``: Use Joseph formula for covariance update (default: True)
* ``simple_Phi``: Use simplified state transition matrix (default: False)

**Initial Conditions:**

* ``B_h_0``: Initial gyro bias estimate [rad/s]
* ``q_h_0``: Initial attitude quaternion estimate
* ``Pq``: Initial attitude error covariance [rad²]
* ``Pb``: Initial gyro bias error covariance [(rad/s)²]
