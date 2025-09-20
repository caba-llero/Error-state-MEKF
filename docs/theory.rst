Theory
======

Error State Multiplicative Extended Kalman Filter
-------------------------------------------------

The error state multiplicative extended Kalman filter (MEKF) is specifically designed for attitude estimation problems. Unlike traditional EKF, the MEKF operates on the error state in the tangent space of the attitude manifold.

Key Advantages:

1. **Multiplicative Error Model**: Attitude errors are represented as small rotations, ensuring the attitude estimate remains normalized
2. **Consistent Covariance**: The covariance matrix accurately represents uncertainty on the manifold
3. **Numerical Stability**: Better conditioning for attitude estimation problems

State Representation
--------------------

**Nominal State:**

* ``q̂``: Attitude quaternion estimate
* ``b̂``: Gyro bias estimate

**Error State:**

* ``δθ``: Attitude error rotation vector (3×1)
* ``δb``: Bias error vector (3×1)

The true attitude and bias are related to the estimates by:

.. math::

    q_t = δq ⊗ q̂
    b_t = b̂ + δb

where ``δq`` is the quaternion corresponding to the error rotation vector ``δθ``.

Measurement Models
------------------

**Gyro Measurements:**

.. math::

    ω_m = ω_t + b_t + n_v + n_u

where:
- ``ω_m``: Measured angular velocity
- ``ω_t``: True angular velocity
- ``b_t``: True gyro bias
- ``n_v``: Angle random walk noise
- ``n_u``: Rate random walk noise

**Star Tracker Measurements:**

.. math::

    q_m = q_t ⊗ q_n

where:
- ``q_m``: Measured attitude quaternion
- ``q_t``: True attitude quaternion
- ``q_n``: Measurement noise (small rotation)

The measurement model in error state form becomes:

.. math::

    δz_m = H δx + n

where ``H = [I₃ 0₃×₃]`` and ``δx = [δθ δb]ᵀ``.

State Transition
----------------

The discrete-time state transition matrix for the error state is:

.. math::

    Φ = \begin{bmatrix}
        Φ_{11} & Φ_{12} \\
        0_{3×3} & I_{3×3}
    \end{bmatrix}

where:

.. math::

    Φ_{11} = I_3 - \frac{\sin(ω‖dt)}{\|ω\|} [ω×] + \frac{1 - \cos(ω‖dt)}{\|ω\|^2} [ω×]^2
    Φ_{12} = -I_3 dt

for small angle approximation, or the full nonlinear form for better accuracy.

Covariance Propagation
----------------------

The error state covariance propagation follows the standard Kalman filter equations:

.. math::

    P_k^- = Φ P_{k-1}^+ Φ^T + Q
    P_k^+ = (I - K H) P_k^-

where ``Q`` is the discrete-time process noise covariance matrix.

References
----------

This implementation is based on:

[1] Markley, F. L., & Crassidis, J. L. (2014). Fundamentals of spacecraft attitude determination and control. Springer.

[2] Solà, J. (2017). Quaternion kinematics for the error-state Kalman filter. arXiv preprint arXiv:1711.02508.
