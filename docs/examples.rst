Examples
========

Basic Simulation
----------------

The main example demonstrates a complete attitude estimation scenario:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import utils as u

    # Set up simulation parameters
    t_max = 20*90  # maximum integration time [s]
    dt = 0.01      # integration timestep [s]
    sigma_startracker = 6  # star tracker accuracy [arcsec]

    # Define true angular velocity profile
    def w_t_fun(t):
        w1 = 0.1*np.sin(0.01*t) * np.pi/180
        w2 = 0.1*np.sin(0.0085*t) * np.pi/180
        w3 = 0.1*np.cos(0.0085*t) * np.pi/180
        return np.vstack((w1, w2, w3))

    # Run simulation (see ekf.py for complete implementation)
    # The simulation generates:
    # - Attitude estimation errors
    # - Gyro bias estimation errors
    # - 3-sigma uncertainty bounds
    # - Pointing error analysis

Custom Scenarios
----------------

You can modify the simulation parameters to test different scenarios:

.. code-block:: python

    # High accuracy sensors
    sigma_startracker = 1  # 1 arcsec accuracy
    sigma_v = 0.1e-6      # Low ARW
    sigma_u = 0.1e-9      # Low RRW

    # High update rates
    freq_startracker = 4  # 4 Hz star tracker
    freq_gyro = 100       # 100 Hz gyros

    # Challenging dynamics
    def w_t_fun(t):
        # More aggressive maneuvers
        return np.array([
            0.5*np.sin(0.05*t),
            0.3*np.cos(0.03*t),
            0.4*np.sin(0.07*t)
        ]).T

Performance Analysis
--------------------

The simulation provides several performance metrics:

.. code-block:: python

    # Calculate RMS errors
    rms_attitude = np.sqrt(np.mean(Z_d_l**2, axis=1))
    rms_bias = np.sqrt(np.mean(B_d**2, axis=1))

    # Calculate 3-sigma consistency
    consistency_ratio = np.mean(np.abs(Z_d_l) < 3*s_l, axis=1)

    # Pointing error analysis
    pointing_error = np.linalg.norm(Z_d_l, axis=0)

See the ``ekf.py`` file for the complete example with plotting and analysis.
