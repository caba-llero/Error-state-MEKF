Theory
======

This section explains exactly the algorithm implemented in the code (``ekf.py`` and ``utils.py``), using the same symbols and operations. Update equations are written as in-place assignments using arrows.

State and Error Representation
------------------------------

Nominal state variables used in the code:

- ``q``: Attitude quaternion (4,)
- ``B``: Gyro bias (3,)

Error state variables (implicitly represented):

- ``dZ``: Attitude error rotation vector (3,)
- ``dB``: Bias error vector (3,)

The attitude error between truth and estimate is computed as

.. math::

   q_d = q_t \otimes q_h^{-1}, \quad Z_d = \mathrm{rotvec}(q_d)

where ``rotvec`` is the exact quaternion-to-rotation-vector map implemented in ``quat_to_rotvec``.

Event Timeline and Measurements
-------------------------------

The simulation advances on a fine time grid with step ``dt``. Two event types trigger filter actions:

- Gyro events at ``freq_gyro``
- Star tracker events at ``freq_startracker``

Events are scheduled via ``measurement_indices``. At a gyro event, the elapsed time since the last gyro event is ``dt_g`` and is used consistently in propagation.

Process and Measurement Models (as implemented)
-----------------------------------------------

Truth propagation per base time step (not a filter update):

.. math::

   q_t \leftarrow \operatorname{quat\_propagate}(q_t,\, w_t,\, dt),\qquad
   B_t \leftarrow B_t + n,\; n \sim \mathcal{N}\!\big(0,\, \sigma_u^2\,dt\, I_3\big).

Gyro measurement at a gyro event (simulation of sensor):

.. math::

   w_m = w_t + B_t + v,\qquad v \sim \mathcal{N}\!\left(0,\, \tfrac{\sigma_v^2}{dt_g}\, I_3\right),\qquad
   w_h = w_m - B_h.

Star tracker measurement at a star event (small-angle synthesis in body frame):

.. math::

   q_m = q_t + \tfrac{1}{2}\, \Xi(q_t)\, Z_n,\quad Z_n \sim \mathcal{N}(0,\, \sigma_{st}^2 I_3),\qquad
   d q_m = \operatorname{normalize}(q_m) \otimes q_h^{-1},\qquad
   dZ_m = \operatorname{rotvec}(d q_m).

Here ``Xi``, ``quat_inv``, and ``quat_to_rotvec`` are the exact functions from ``utils.py``. The measurement matrix and noise are set as

.. math::

   H = [\,I_3\; 0\,],\qquad R = \sigma_{st}^2 I_3,\quad \sigma_{st} = \text{sigma\_startracker}\cdot \text{arcsec\_to\_rad}.

Error-State Propagation
-----------------------

At each gyro event, the error-state covariance is propagated with the discrete transition ``Phi(dt_g, w_h)`` and process noise ``Q(sigma_v, sigma_u, dt_g)``.

Transition matrix (``Phi``): the code computes

.. math::

   \Phi = \begin{bmatrix} \Phi_{11} & \Phi_{12} \\ 0 & I_3 \end{bmatrix}

with

.. math::

   p = \|w_h\|\,dt, \quad a = [w_h\times], \quad s_p = \sin p, \; c_p = \cos p.

Full form used when the rate is not tiny and ``simple_Phi`` is False:

.. math::

   \Phi_{11} = I_3 - \frac{a}{\|w_h\|} s_p + \frac{a^2}{\|w_h\|^2} (1-c_p), \quad
   \Phi_{12} = -I_3\,dt - \frac{a^2}{\|w_h\|^3}(p - s_p) + \frac{a}{\|w_h\|^2}(1-c_p).

Simple form used when ``\|w_h\|`` is very small or ``simple_Phi`` is True:

.. math::

   \Phi_{11} = I_3 - a\,dt, \quad \Phi_{12} = -I_3\,dt.

Process noise (``Q``) is block-structured as in the code:

.. math::

   Q = \begin{bmatrix} Q_{11} & Q_{12} \\ Q_{12} & Q_{22} \end{bmatrix},\quad
   Q_{11} = (\sigma_v^2 dt + \sigma_u^2 dt^3/3) I_3,\; Q_{12} = -(\sigma_u^2 dt^2/2) I_3,\; Q_{22} = (\sigma_u^2 dt) I_3.

Covariance propagation update performed at gyro events:

.. math::

   P \leftarrow \Phi\, P\, \Phi^T + Q.

Nominal attitude is propagated with the measured rate minus estimated bias:

.. math::

   q_h \leftarrow \operatorname{quat\_propagate}(q_h,\, w_h,\, dt_g).

Measurement Update (Star Tracker Events)
----------------------------------------

Kalman gain and innovation are computed directly from the code:

.. math::

   S = H P H^T + R, \quad K = P H^T S^{-1}

The gain is split into attitude and bias parts (rows 0:3 and 3:6):

.. math::

   dZ_h \leftarrow K_Z\, dZ_m,\qquad dB_h \leftarrow K_B\, dZ_m,\qquad B_h \leftarrow B_h + dB_h.

Attitude correction uses the exact quaternion mapping of ``dZ_h`` to a unit quaternion ``dq_err`` and applies it multiplicatively:

.. math::

   d q_{err} = \operatorname{rotvec\_to\_quat}(dZ_h),\qquad q_h \leftarrow d q_{err} \otimes q_h,\qquad q_h \leftarrow \frac{q_h}{\lVert q_h \rVert}.

Covariance measurement update uses the Joseph form when ``Joseph = True`` (default), otherwise a simple form:

.. math::

   I_{KH} = I_6 - K H,\qquad P \leftarrow I_{KH} P I_{KH}^T + K R K^T.

If ``Joseph = False``, the code applies

.. math::

   P \leftarrow (I_6 - K H) P.

Function Definitions Used
-------------------------

Quaternion to rotation vector (``rotvec``), matching ``utils.quat_to_rotvec``:

.. math::

   q = \begin{bmatrix} v \\ w \end{bmatrix},\; \lVert q \rVert = 1,\; v \in \mathbb{R}^3,\; w \in \mathbb{R},\qquad
   \theta = 2\,\arctan2(\lVert v \rVert,\, w),\qquad
   \operatorname{rotvec}(q) = \begin{cases}
     0, & \lVert v \rVert < \varepsilon,\\
     v\, \dfrac{\theta}{\lVert v \rVert}, & \text{otherwise.}
   \end{cases}

Rotation vector to quaternion (used in ``ekf.py`` when injecting ``dZ_h``):

.. math::

   \theta = \lVert dZ_h \rVert,\quad \text{if } \theta > 0:\; e = dZ_h/\theta,\; d q_{err} = \begin{bmatrix} e\,\sin(\tfrac{\theta}{2}) \\ \cos(\tfrac{\theta}{2}) \end{bmatrix};\quad \text{else } d q_{err} = \begin{bmatrix} 0 \\ 1 \end{bmatrix}.

Logged Metrics
--------------

After each event, the implementation logs

- ``s = sqrt(diag(P))`` (state standard deviations)
- ``Z_d`` and ``B_d = B_t - B_h`` (errors)
- ``G = ||Z_d||`` (pointing error)

All computations above correspond exactly to functions ``Phi``, ``Q``, ``K``, ``P_prop``, ``P_meas``, ``quat_propagate``, ``quat_to_rotvec``, and the event-driven loop in ``ekf.py``.
