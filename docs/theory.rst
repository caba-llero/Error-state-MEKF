Theory
======

This section presents the error-state multiplicative EKF implemented on the code.


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

   \delta q = q_t \otimes \hat{q}^{-1} = \begin{bmatrix} \boldsymbol{\epsilon} \\ \eta \end{bmatrix}

We use the error rotation vector :math:`\boldsymbol{\delta\theta} \in \mathbb{R}^3` in our state instead of the error quaternion:

.. math::

   \boldsymbol{\delta\theta} = 2 \arctan2 \left( \|\boldsymbol{\epsilon}\|, \eta \right) \frac{\boldsymbol{\epsilon}}{\|\boldsymbol{\epsilon}\|}

For small angles,

.. math::

   \delta q \approx \begin{bmatrix} \frac{1}{2} \boldsymbol{\delta\theta} \\ 1 \end{bmatrix}

A true gyro measurement is modelled as:

.. math::

   \boldsymbol{\omega}_m = \boldsymbol{\omega}_t + \boldsymbol{\beta}_t + \mathbf{v}_\Delta

.. math::

   \mathbf{v}_Delta \sim \mathcal{N}\!\left(0, \frac{\sigma_v^2}{\Delta t_g} I_3 \right)

Here, :math:`\Delta t_g` is the gyro sampling period. The gyro bias error is defined as:

.. math::

    \boldsymbol{\delta \beta} = \boldsymbol{\beta}_t - \boldsymbol{\hat{\beta}}



Ground truth propagation
------------------------

The user inputs a desired ground truth angular velocity :math:`\omega_t(t)` with the function ``w_t_fun`` in ``efk.py``. The ground truth quaternion is then propagated at a period of :math:`\Delta t`.

Define :math:`\varphi = \|\boldsymbol{\omega}_t\| \Delta t`. The quaternion increment associated with :math:`\boldsymbol{\omega}_h \Delta t_g` is:

.. math::

   \Delta q = \begin{bmatrix} \mathbf{e} \sin(\frac{\varphi}{2}) \\ \cos(\frac{\varphi}{2}) \end{bmatrix}

.. math::

   \mathbf{e} = \frac{\boldsymbol{\omega}_h}{ \|\boldsymbol{\omega}_h\| }

This assumes that the angular velocity is constant throughout this timestep. The ground truth update is then:

.. math::

   \hat{q} \leftarrow \Delta q \otimes \hat{q}

.. math::

   \hat{q} \leftarrow \frac{\hat{q}}{\|\hat{q}\|}


The bias grows as a Weiner process, at a level related to the rate random walk coefficient:

.. math::

   \boldsymbol{\dot{\beta_t}}(t) = \boldsymbol{u}

In a discrete step, this is done in the code with:

.. math::

    \boldsymbol{\beta_t} \leftarrow \boldsymbol{\beta_t} + \boldsymbol{u_Delta}

   \mathbf{u_\Delta} \sim \mathcal{N}\!\left(0, \sigma_u^2 \Delta t I_3 \right)




Discrete-Time Propagation (Gyro Events)
---------------------------------------

At gyro sampling instants separated by :math:`\Delta t_g`, the gyro provides

.. math::

   \boldsymbol{\omega}_m = \boldsymbol{\omega}_t + \boldsymbol{\beta}_t + \mathbf{v}_k

.. math::

   \mathbf{v}_k \sim \mathcal{N}\!\left(0, \frac{\sigma_v^2}{\Delta t_g} I_3 \right)

The propagation uses the bias-compensated rate

.. math::

   \boldsymbol{\omega}_h = \boldsymbol{\omega}_m - \hat{\boldsymbol{\beta}}

Define :math:`\varphi = \|\boldsymbol{\omega}_h\| \Delta t_g` and :math:`[\boldsymbol{\omega}_h]_\times` such that :math:`[\boldsymbol{\omega}_h]_\times \mathbf{x} = \boldsymbol{\omega}_h \times \mathbf{x}` for :math:`\mathbf{x} \in \mathbb{R}^3`. The quaternion increment associated with :math:`\boldsymbol{\omega}_h \Delta t_g` is

.. math::

   d q_g = \begin{bmatrix} \mathbf{e} \sin(\frac{\varphi}{2}) \\ \cos(\frac{\varphi}{2}) \end{bmatrix}

.. math::

   \mathbf{e} = \begin{cases} \boldsymbol{\omega}_h / \|\boldsymbol{\omega}_h\| & \|\boldsymbol{\omega}_h\| > 0 \\ 0 & \text{otherwise} \end{cases}

The nominal attitude update is

.. math::

   \hat{q} \leftarrow d q_g \otimes \hat{q}

.. math::

   \hat{q} \leftarrow \frac{\hat{q}}{\|\hat{q}\|}

The bias is constant in propagation:

.. math::

   \hat{\boldsymbol{\beta}} \leftarrow \hat{\boldsymbol{\beta}}


Linearized Error-State Propagation
----------------------------------

Let :math:`\boldsymbol{\delta x} = \begin{bmatrix} \boldsymbol{\delta\theta} \\ \delta\boldsymbol{\beta} \end{bmatrix} \in \mathbb{R}^6`. For a step :math:`\Delta t_g` with input :math:`\boldsymbol{\omega}_h`, the first-order discrete transition is

.. math::

   \boldsymbol{\delta x}_{k+1} = \Phi_k \boldsymbol{\delta x}_k + \mathbf{w}_k

.. math::

   \Phi_k = \begin{bmatrix} \Phi_{11} & \Phi_{12} \\ 0 & I_3 \end{bmatrix}

With :math:`\mathbf{a} = [\boldsymbol{\omega}_h]_\times`, :math:`\varphi = \|\boldsymbol{\omega}_h\| \Delta t_g`, :math:`s = \sin \varphi`, :math:`c = \cos \varphi`,

.. math::

   \Phi_{11} = I_3 - \frac{\mathbf{a}}{\|\boldsymbol{\omega}_h\|} s + \frac{\mathbf{a}^2}{\|\boldsymbol{\omega}_h\|^2} (1 - c)

.. math::

   \Phi_{12} = -I_3 \Delta t_g - \frac{\mathbf{a}^2}{\|\boldsymbol{\omega}_h\|^3} (\varphi - s) + \frac{\mathbf{a}}{\|\boldsymbol{\omega}_h\|^2} (1 - c)

For :math:`\|\boldsymbol{\omega}_h\| \Delta t_g \ll 1`, the approximation

.. math::

   \Phi_{11} \approx I_3 - \mathbf{a} \Delta t_g

.. math::

   \Phi_{12} \approx -I_3 \Delta t_g

is used.


Process Noise Discretization
----------------------------

With gyro angle random walk density :math:`\sigma_v^2` and bias random walk density :math:`\sigma_u^2`, the discrete process covariance over :math:`\Delta t_g` is

.. math::

   Q_k = \begin{bmatrix} Q_{11} & Q_{12} \\ Q_{12} & Q_{22} \end{bmatrix}

.. math::

   Q_{11} = \left( \sigma_v^2 \Delta t_g + \frac{\sigma_u^2 \Delta t_g^3}{3} \right) I_3

.. math::

   Q_{12} = -\frac{\sigma_u^2 \Delta t_g^2}{2} I_3

.. math::

   Q_{22} = \sigma_u^2 \Delta t_g I_3

Covariance propagation:

.. math::

   P \leftarrow \Phi_k P \Phi_k^T + Q_k


Measurement Model (Star Tracker Events)
---------------------------------------

At star tracker instants, a quaternion measurement :math:`q_m` of attitude is available. The innovation quaternion is

.. math::

   \delta q_m = q_m \otimes \hat{q}^{-1} = \begin{bmatrix} \boldsymbol{\epsilon}_m \\ \eta_m \end{bmatrix}

.. math::

   \|\delta q_m\| = 1

The corresponding rotation-vector innovation is

.. math::

   \mathbf{z}_k = \begin{cases} 0 & \|\boldsymbol{\epsilon}_m\| = 0 \\ 2 \arctan2 \left( \|\boldsymbol{\epsilon}_m\|, \eta_m \right) \frac{\boldsymbol{\epsilon}_m}{\|\boldsymbol{\epsilon}_m\|} & \text{otherwise} \end{cases}

Under the small-angle assumption,

.. math::

   \mathbf{z}_k = H \boldsymbol{\delta x}_k + \mathbf{n}_k

.. math::

   H = \begin{bmatrix} I_3 & 0 \end{bmatrix}

.. math::

   \mathbf{n}_k \sim \mathcal{N}(0, R)

.. math::

   R = \sigma_{st}^2 I_3


Measurement Update and Injection
--------------------------------

Compute the innovation covariance and gain, then the correction:

.. math::

   S_k = H P H^T + R

.. math::

   K_k = P H^T S_k^{-1}

.. math::

   \boldsymbol{\delta \hat x}_k = K_k \mathbf{z}_k

Split :math:`\boldsymbol{\delta \hat x}_k = \begin{bmatrix} \hat{\boldsymbol{\delta\theta}} \\ \widehat{\delta\boldsymbol{\beta}} \end{bmatrix}` and inject into the nominal state:

.. math::

   \hat{\boldsymbol{\beta}} \leftarrow \hat{\boldsymbol{\beta}} + \widehat{\delta\boldsymbol{\beta}}

Let :math:`\alpha = \|\hat{\boldsymbol{\delta\theta}}\|` and :math:`\mathbf{e} = \hat{\boldsymbol{\delta\theta}} / \alpha` if :math:`\alpha > 0`. Define

.. math::

   d q_{\text{err}} = \begin{bmatrix} \mathbf{e} \sin(\frac{\alpha}{2}) \\ \cos(\frac{\alpha}{2}) \end{bmatrix}

( or :math:`\begin{bmatrix} 0 \\ 1 \end{bmatrix}` if :math:`\alpha = 0` ),

and update the attitude multiplicatively:

.. math::

   \hat{q} \leftarrow d q_{\text{err}} \otimes \hat{q}

.. math::

   \hat{q} \leftarrow \frac{\hat{q}}{\|\hat{q}\|}

Update the covariance (Joseph form):

.. math::

   I_{KH} = I_6 - K_k H

.. math::

   P \leftarrow I_{KH} P I_{KH}^T + K_k R K_k^T


Timing of Updates
-----------------

Propagation occurs at gyro sampling instants :math:`\Delta t_g`. Measurement updates occur when a star tracker observation is available. Between observations, only propagation is performed.
