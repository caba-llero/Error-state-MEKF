Theory
======

This section presents the error-state multiplicative EKF implemented on the code.

Please note that the documentation currently is bare bones, and the notation used might not be very rigorous; I will be adding more details and improving the documentation as the project progresses. 


Notation
--------

.. list-table:: Main Variables and Parameters
   :header-rows: 1

   * - Code
     - Math
     - Description
   * - ``q``
     - :math:`\mathbf{q}`
     - Attitude quaternion
   * - ``B``
     - :math:`\boldsymbol{\beta}`
     - Gyro bias vector
   * - ``Z``
     - :math:`\boldsymbol{\delta\theta}`
     - Attitude error rotation vector
   * - ``w``
     - :math:`\boldsymbol{\omega}`
     - Angular velocity vector
   * - ``P``
     - :math:`P`
     - Error covariance matrix
   * - ``Phi``
     - :math:`\Phi`
     - State transition matrix
   * - ``Q``
     - :math:`Q`
     - Discrete-time process noise covariance matrix
   * - ``K``
     - :math:`K`
     - Kalman gain matrix
   * - ``H``
     - :math:`H`
     - Measurement matrix
   * - ``R``
     - :math:`R`
     - Measurement noise covariance matrix
   * - ``S``
     - :math:`S`
     - Innovation covariance matrix
   * - ``sigma_v``
     - :math:`\sigma_v`
     - Gyro angle random walk coefficient [rad/s^{0.5}]
   * - ``sigma_u``
     - :math:`\sigma_u`
     - Gyro rate random walk coefficient [rad/s^{1.5}]
   * - ``sigma_startracker``
     - :math:`\sigma_{st}`
     - Star tracker measurement accuracy [arcsec]
   * - ``dt``
     - :math:`\Delta t`
     - Base time step for ground truth propagation
   * - ``dt_g``
     - :math:`\Delta t_g`
     - Time step between gyro measurements
   * - ``I3``
     - :math:`I_3`
     - Identity matrix (3x3)


.. list-table:: Subscripts and Notation Conventions
   :header-rows: 1

   * - Code
     - Math
     - Description
   * - ``X_h``
     - :math:`\hat{X}`
     - Estimated value of variable
   * - ``X_t``
     - :math:`X_t`
     - Ground truth value of variable
   * - ``X_m``
     - :math:`X_m`
     - Measured value of variable
   * - ``X_0``
     - :math:`X_0`
     - Initial value of variable
   * - ``X_l``
     - :math:`X_l`
     - Logged variable for analysis
   * - ``X_d``
     - :math:`\delta X`
     - Error between truth and estimate (``X_t - X_h``)
   * - ``X_n``
     - :math:`X_n`
     - Noise term added to variable


State and Error Representation
------------------------------

The nominal state consists of the estimates for the attitude quaternion and gyro bias, where :math:`(\hat{q}, \hat{\beta}) \in \mathbb{S}^3 \times \mathbb{R}^3`. However, the Kalman filter updates the error states, then injects these updates to the nominal state. The error quaternion is the relation between the ground truth quaternion, represented multiplicatevely:

.. math::

   \delta q = q_t \otimes \hat{q}^{-1} = \begin{bmatrix} \boldsymbol{\delta\epsilon} \\ \delta\eta \end{bmatrix}

We use the error rotation vector :math:`\boldsymbol{\delta\theta} \in \mathbb{R}^3` in our state instead of the error quaternion:

.. math::

   \boldsymbol{\delta\theta} = 2 \arctan2 \left( \|\boldsymbol{\delta\epsilon}\|, \delta\eta \right) \frac{\boldsymbol{\delta\epsilon}}{\|\boldsymbol{\delta\epsilon}\|}

For small angles,

.. math::

   \delta q \approx \begin{bmatrix} \frac{1}{2} \boldsymbol{\delta\theta} \\ 1 \end{bmatrix}

A gyro measurement is modelled with a moving bias and white noise (random walk):

.. math::

   \boldsymbol{\omega}_m(t) = \boldsymbol{\omega}_t(t) + \boldsymbol{\beta}_t(t) + \mathbf{v}(t)

.. math::

   \mathbf{v}(t) \sim \mathcal{N}\!\left(0, \sigma_v^2 I_3 \right)

The bias grows as a Weiner process, at a level related to the rate random walk coefficient:

.. math::

   \boldsymbol{\dot{\beta_t}}(t) = \boldsymbol{u}(t)

.. math::
   \boldsymbol{u}(t) \sim \mathcal{N}\!\left(0, \sigma_u^2 I_3 \right)

The gyro bias error is defined as:

.. math::

    \boldsymbol{\delta \beta} = \boldsymbol{\beta}_t - \boldsymbol{\hat{\beta}}


Quaternion math
---------------

Throughout the documentation, we use the following quaternion multiplication convention:

.. math::

   q_1 \otimes q_2 = \begin{bmatrix} q_{1x} q_{2x} - q_{1y} q_{2y} - q_{1z} q_{2z} + q_{1w} q_{2w} \\ q_{1x} q_{2y} + q_{1y} q_{2x} + q_{1z} q_{2w} - q_{1w} q_{2z} \\ q_{1x} q_{2z} - q_{1y} q_{2w} + q_{1z} q_{2x} + q_{1w} q_{2y} \\ q_{1x} q_{2w} - q_{1y} q_{2z} - q_{1z} q_{2y} + q_{1w} q_{2x} \end{bmatrix}

In the code we use the function ``quat_mul`` in ``utils.py`` to perform this multiplication. The quaternion inverse is defined as:

.. math::

   q^{-1} = \begin{bmatrix} -q_x \\ -q_y \\ -q_z \\ q_w \end{bmatrix}


Initialization
--------------

The user provides an initial gyro bias :math:`\boldsymbol{\beta_0}` and initial attitude quaternion :math:`\boldsymbol{q_0}`. Additionally, the user provides initial values for the Kalman filter estimator: the attitude error covariance :math:`P_q` and the gyro bias error covariance :math:`P_b`; the initial estimate of the attitude error :math:`\boldsymbol{\delta\theta_0}` and the initial estimate of the gyro bias error :math:`\boldsymbol{\delta\beta_0}`.


Ground truth update
------------------------

The user inputs a desired ground truth angular velocity :math:`\boldsymbol{\omega_t}(t)` with the function ``w_t_fun`` in ``efk.py``. The user also inputs an initial gyro bias :math:`\boldsymbol{\beta_0}`. The ground truth quaternion and gyro bias are propagated at a period of :math:`\Delta t`.

Define :math:`\varphi = \|\boldsymbol{\omega}_t\| \Delta t`. The quaternion increment associated with the rotation :math:`\boldsymbol{\omega}_t \Delta t` is:

.. math::

   \Delta q = \begin{bmatrix} \mathbf{e} \sin(\frac{\varphi}{2}) \\ \cos(\frac{\varphi}{2}) \end{bmatrix}

.. math::

   \mathbf{e} = \frac{\boldsymbol{\omega}_t}{ \|\boldsymbol{\omega}_t\| }

This assumes that the angular velocity is constant throughout this timestep. The ground truth update is:

.. math::
   q \leftarrow \Delta q \otimes q
   


In a discrete step, the bias is updated as:

.. math::

    \boldsymbol{\beta_t} \leftarrow \boldsymbol{\beta_t} + \boldsymbol{u_\Delta}

   \mathbf{u_\Delta} \sim \mathcal{N}\!\left(0, \sigma_u^2 \Delta t I_3 \right)



Estimate propagation (gyro measurements)
---------------------------------------

Gyros are used for dynamic model replacement, i.e. I do not integrate the Euler rigid body equations. Instead, I use the gyro measurements to propagate the estimate. 

At gyro sampling instants separated by :math:`\Delta t_g`, the gyro provides a measurement which is synthesized from the ground truth angular velocity :math:`\boldsymbol{\omega}_t` and the gyro bias :math:`\boldsymbol{\beta}_t`.


.. math::

   \boldsymbol{\omega}_m = \boldsymbol{\omega}_t + \boldsymbol{\beta}_t + \mathbf{v}_\Delta

.. math::

   \mathbf{v}_\Delta \sim \mathcal{N}\!\left(0, \frac{\sigma_v^2}{\Delta t_g} I_3 \right)

The estimate of the angular velocity is:

.. math::

   \boldsymbol{\hat{\omega}} = \boldsymbol{\omega}_m - \hat{\boldsymbol{\beta}}


We propagate the estimated attitude quaternion the same way as the ground truth (see the section above).

The bias is constant in propagation:

.. math::

   \hat{\boldsymbol{\beta}} \leftarrow \hat{\boldsymbol{\beta}}

The covariance propagation in a gyro measurement step is discussed in the next two sections.

Linearized Error-State Propagation
----------------------------------

Let :math:`\boldsymbol{\delta x} = \begin{bmatrix} \boldsymbol{\delta\theta} \\ \delta\boldsymbol{\beta} \end{bmatrix} \in \mathbb{R}^6`. For a step :math:`\Delta t_g` with input :math:`\boldsymbol{\hat{\omega}}`, the first-order discrete transition is:

.. math::

   \boldsymbol{\delta x}  \leftarrow \Phi \boldsymbol{\delta x}

.. math::

   \Phi = \begin{bmatrix} \Phi_{11} & \Phi_{12} \\ 0 & I_3 \end{bmatrix}

With :math:`\varphi = \|\boldsymbol{\hat{\omega}}\| \Delta t_g`, :math:`s = \sin \varphi`, :math:`c = \cos \varphi`, and \boldsymbol{\hat{\omega}}_\times is the skew-symmetric matrix of :math:`\boldsymbol{\hat{\omega}}`,

.. math::

   \Phi_{11} = I_3 - \frac{\boldsymbol{\hat{\omega}}_\times}{\|\boldsymbol{\hat{\omega}}\|} s + \frac{\boldsymbol{\hat{\omega}}_\times^2}{\|\boldsymbol{\hat{\omega}}\|^2} (1 - c)

.. math::

   \Phi_{12} = -I_3 \Delta t_g - \frac{\boldsymbol{\hat{\omega}}_\times^2}{\|\boldsymbol{\hat{\omega}}\|^3} (\varphi - s) + \frac{\boldsymbol{\hat{\omega}}_\times}{\|\boldsymbol{\hat{\omega}}\|^2} (1 - c)

For :math:`\|\boldsymbol{\hat{\omega}}\| \Delta t_g \ll 1`, the approximation

.. math::

   \Phi_{11} \approx I_3 - \boldsymbol{\hat{\omega}}_\times \Delta t_g

.. math::

   \Phi_{12} \approx -I_3 \Delta t_g

can be used (user toggleable).


Process Noise Discretization
----------------------------

With gyro angle random walk density :math:`\sigma_v^2` and bias random walk density :math:`\sigma_u^2`, the discrete process covariance over :math:`\Delta t_g` is

.. math::

   Q = \begin{bmatrix} Q_{11} & Q_{12} \\ Q_{12} & Q_{22} \end{bmatrix}

.. math::

   Q_{11} = \left( \sigma_v^2 \Delta t_g + \frac{\sigma_u^2 \Delta t_g^3}{3} \right) I_3

.. math::

   Q_{12} = -\frac{\sigma_u^2 \Delta t_g^2}{2} I_3

.. math::

   Q_{22} = \sigma_u^2 \Delta t_g I_3

Finally, covaraince is propagated as:

.. math::

   P \leftarrow \Phi P \Phi^T + Q


Star tracker measurement synthesis
----------------------------------

At star tracker measurement events, a quaternion measurement :math:`q_m` of attitude is available. This is synthesized from the ground truth quaternion :math:`q_t` with white noise. This is done by synthesizing a R3 vector with white noise and then converting it to a quaternion.

.. math::
    \boldsymbol{\theta}_n \sim \mathcal{N}(0, \sigma_{st}^2 I_3)

.. math::
    q_m = q_t \otimes \mathbf{q}_n \approx q_t + \frac{1}{2} \Xi(q_t) \boldsymbol{\theta}_n


Here, we used small angle approximation (as the startracker measurement error is in the order of arcseconds), and we used the quaternion function :math:`\Xi` defined as:

.. math::
    \Xi(q) = \begin{bmatrix} q_w & -q_z & q_y \\ q_z & q_w & -q_x \\ -q_y & q_x & q_w \end{bmatrix}

Now, the estimated error quaternion for this measurement is:

.. math::
    \delta q_m = q_m \otimes \hat{q}^{-1} 

We then convert this to a rotation vector, explained in the first section. 

.. math::

    \boldsymbol{\delta q_m} \mapsto \boldsymbol{\delta\theta}_m

The observation model matrix is:

.. math::

    H = \begin{bmatrix} I_3 & 0 \end{bmatrix}

The measurement noise covariance matrix is the startracker measurement accuracy (assumed isotropic):

.. math::

    R = \sigma_{st}^2 I_3

We discuss the measurement update and injection in the next section, carried out in star tracker measurement events.


Measurement update and injection
--------------------------------

Compute the innovation covariance and gain, then the correction:

.. math::

   S \leftarrow H P H^T + R

.. math::

   K \leftarrow P H^T S^{-1}

.. math::

   \boldsymbol{\delta \hat x} \leftarrow K \boldsymbol{\delta\theta}_m

Split :math:`\boldsymbol{\delta \hat x} = \begin{bmatrix} \hat{\boldsymbol{\delta\theta}} \\ \hat{\boldsymbol{\delta\beta}} \end{bmatrix}` and inject into the global state:

.. math::

   \hat{\boldsymbol{\beta}} \leftarrow \hat{\boldsymbol{\beta}} + \widehat{\delta\boldsymbol{\beta}}

Map the estimated error vector to a quaternion:

.. math::

   \hat{\boldsymbol{\delta \theta}} \mapsto \hat{\boldsymbol{\delta q}}

and update the attitude multiplicatively:

.. math::

   \hat{q} \leftarrow \hat{\boldsymbol{\delta q}} \otimes \hat{q}


Update the covariance (the Joseph form is used by default; a simple form is also available):


.. math::

   P \leftarrow (I_6 - K H) P (I_6 - K H)^T + K R K^T
