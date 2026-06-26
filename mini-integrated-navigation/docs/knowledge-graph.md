# Knowledge Graph — mini-integrated-navigation

## L1: Definitions (Complete)
- Navigation coordinate frames (ECI, ECEF, NED, ENU, Body)
- WGS84 geodetic position (latitude, longitude, altitude)
- 3D vector and INS error state (15-state psi-angle)
- IMU measurement (gyro, accel, timestamp)
- IMU error parameters (bias, scale factor, noise PSD, misalignment)
- GNSS satellite measurement (pseudorange, Doppler, carrier phase)
- GNSS position solution (DOP, clock bias)
- Full navigation solution (position, velocity, attitude, biases)
- Quaternion, DCM, Euler angle, axis-angle representations
- Kalman filter state (x, P, F, H, Q, R, K)

## L2: Core Concepts (Complete)
- Strapdown inertial navigation (INS mechanization)
- GNSS positioning via pseudorange multilateration
- Sensor fusion via Kalman filtering
- Dead reckoning from IMU integration
- Coordinate frame transformations
- IMU error compensation (bias, scale factor, misalignment)
- GNSS error sources (ionosphere, troposphere, satellite clock)
- Complementary filtering for attitude estimation
- ZUPT (zero-velocity update) for pedestrian navigation
- Non-holonomic constraints for land vehicles

## L3: Mathematical Structures (Complete)
- SO(3) rotation group: quaternion algebra, DCM operations
- Linear algebra: matrix multiply, transpose, Cholesky decomposition
- Quaternion exponential/logarithm maps
- Skew-symmetric (cross-product) matrix
- Spherical linear interpolation (SLERP) on S^3
- WGS84 ellipsoid geometry (meridian/transverse radius)
- Normal gravity model (Somigliana formula)
- Attitude kinematics (quaternion differential equation)
- Gauss-Newton iteration for GNSS positioning
- ECEF/geodetic coordinate conversions

## L4: Fundamental Laws (Complete)
- Kalman filter optimality (linear-Gaussian, MMSE)
- Information filter (canonical form, Y = P^{-1})
- Extended Kalman filter (linearization about estimate)
- Error-state (indirect) Kalman filter for INS
- INS error dynamics (psi-angle formulation, Goshen-Meskin 1992)
- GPS satellite orbit model (Kepler + perturbations)
- Klobuchar ionospheric delay model
- Saastamoinen tropospheric delay model
- Covariance intersection (consistent fusion, unknown correlation)

## L5: Algorithms/Methods (Complete)
- Linear Kalman filter (predict + update + Joseph form)
- Scalar sequential measurement update
- EKF for nonlinear systems
- Information filter (predict/update in information space)
- Quaternion-based attitude representation
- TRIAD algorithm for attitude determination
- Madgwick AHRS filter (gradient descent)
- Mahony AHRS filter (explicit complementary, PI)
- Complementary filter (gyro + accel/mag blend)
- IMU strapdown mechanization
- 2-sample coning and sculling compensation
- Coarse INS alignment (leveling + gyrocompassing)
- ZUPT static detection
- Weighted Least Squares GNSS PVT
- RAIM residual-based fault detection
- RAIM Fault Detection and Exclusion (FDE)
- Solution separation integrity monitoring
- Normalized innovation squared (NIS) test
- Federated Kalman filter (Carlson, information sharing)
- Covariance intersection (optimal omega via trace minimization)
- Moving average filter for IMU data
- Outlier detection (MAD-based)
- Lever arm compensation (GNSS antenna to IMU)
- Navigation solution interpolation (SLERP for attitude)

## L6: Canonical Problems (Complete)
- Loosely coupled INS/GNSS (15-state EKF, position/velocity updates)
- Tightly coupled INS/GNSS (17-state EKF, pseudorange/Doppler updates)
- AHRS using IMU + magnetometer (Madgwick, Mahony, complementary)
- Standalone GNSS PVT (WLS with iono/tropo corrections)
- ZUPT-aided pedestrian navigation
- Land vehicle navigation with NHC

## L7: Applications (Partial+)
- UAV navigation (loose/tight INS/GNSS integration)
- Automotive dead reckoning (NHC + GNSS)
- Smartphone sensor fusion (AHRS + GNSS)
- Aviation integrity monitoring (RAIM)

## L8: Advanced Topics (Partial+)
- Federated (decentralized) Kalman filtering
- Covariance intersection for multi-sensor fusion
- Information filtering for efficient multi-rate fusion
- Solution separation for multi-constellation GNSS

## L9: Research Frontiers (Partial)
- AI-based sensor fusion (documented, not implemented)
- Cooperative navigation (documented)
- Quantum navigation sensors (documented)
