# Error state MEKF 

An error state multiplicative extended Kalman filter, specifically for attitude and gyro bias estimation. Loosely based on the EKF algorithms found in:

[1] Fundamentals of Spacecraft Attitude Determination and Control, by Markley & Crassidis (2014)
[2] Quaternion kinematics for the error-state Kalman filter, by Joan Solà (2017)

We assume that the sensors available are a startracker and gyros on each axis. The estimator synthesizes startracker and gyro measurements based on the inputted noise parameters. The gyro measurement model assumes ARW and RRW; the startracker measurement model applies isotropic Gaussian noise to each axis.

