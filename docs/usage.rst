Usage
=====

Basic Usage
-----------

Instantiate the filter by passing a single ``inputs`` dictionary. Then call ``calculate()`` and use the available plotting helpers.

.. code-block:: python

    from mekf import MEKF

    # Create MEKF instance with minimal inputs
    inputs = {
        "t_max": 1000,
        "dt": 0.01,
        "sigma_startracker": 6,
        "sigma_v": 1e-6,
        "sigma_u": 1e-9,
        "freq_startracker": 1,
        "freq_gyro": 20,
        "init_inaccuracy": 30,
        "rng_seed": 1,
        # Optional: initial conditions and custom truth rate function
        # "B_h_0": np.array([0.0, 0.0, 0.0]),
        # "B_t_0": np.array([0.1, 0.1, 0.1]) * (np.pi / (180 * 3600)),
        # "q_t_0": np.array([1, 0, 0, 1]) / (2**0.5),
        # "w_t_fun": lambda t: ...
    }

    mekf = MEKF(inputs=inputs)
    results = mekf.calculate()

    # Plotting helpers (call any subset)
    mekf.plot_pointing_error()
    mekf.plot_bias()
    mekf.plot_attitude()
    mekf.plot_errors_with_bounds()
    # Or all at once:
    # mekf.plot_all()

Parameters
----------

All previous top-level values from ``ekf.py`` are now provided inside the ``inputs`` dict:

- ``t_max`` (float)
- ``dt`` (float)
- ``sigma_startracker`` (float, arcsec)
- ``sigma_v`` (float, rad/s/sqrt(Hz))
- ``sigma_u`` (float, rad/s^1.5)
- ``freq_startracker`` (Hz)
- ``freq_gyro`` (Hz)
- ``init_inaccuracy`` (float)
- ``rng_seed`` (int)
- ``Joseph`` (bool)
- ``simple_Phi`` (bool)
- ``B_h_0`` (3,)
- ``Pq`` (3x3)
- ``Pb`` (3x3)
- ``B_t_0`` (3,)
- ``q_t_0`` (4,)
- ``w_t_fun`` (callable)
