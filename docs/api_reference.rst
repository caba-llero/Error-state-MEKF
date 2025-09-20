API Reference
=============

Core Functions
--------------

.. automodule:: utils
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: utils.Phi

.. autofunction:: utils.Q

.. autofunction:: utils.K

.. autofunction:: utils.P_meas

.. autofunction:: utils.P_prop

.. autofunction:: utils.quat_mul

.. autofunction:: utils.quat_propagate

.. autofunction:: utils.quat_inv

.. autofunction:: utils.quat_to_rotvec

.. autofunction:: utils.startracker_meas

Utility Functions
-----------------

.. autofunction:: utils.skew

.. autofunction:: utils.Xi

.. autofunction:: utils.measurement_indices

Main Simulation
---------------

The main simulation code is contained in ``ekf.py`` and demonstrates the complete MEKF implementation.

Key Variables
-------------

**Notation Convention:**

* ``_h``: Estimated variable (hat notation)
* ``_g``: Ground truth value
* ``_0``: Initial value
* ``_m``: Measured variable
* ``_l``: Logged variable for analysis
* ``_d``: Difference/error between true and estimated value
* ``_n``: Noise

**State Variables:**

* ``Z``: Rotation vector (3,) - attitude error state
* ``q``: Attitude quaternion (4,) - spacecraft orientation
* ``B``: Gyro bias (3,) - systematic gyro error
* ``w``: Angular velocity (3,) - rotation rate
* ``s``: Standard deviations of error state (6,)
