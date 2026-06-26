# Course Tree — mini-integrated-navigation

## Prerequisite Knowledge Dependencies

```
mini-integrated-navigation
│
├── mini-signal-system-theory (required)
│   ├── Fourier/Laplace transforms
│   ├── Linear systems theory
│   └── Random processes (power spectral density)
│
├── mini-digital-signal-process (required)
│   ├── Discrete-time signals
│   ├── Digital filter design
│   └── Spectral analysis
│
├── mini-communication-principle (required)
│   ├── Modulation/demodulation
│   ├── Spread spectrum (CDMA for GPS)
│   └── Signal detection theory
│
├── mini-control-automation (required)
│   ├── State-space representation
│   ├── Observability and controllability
│   └── Optimal estimation (Kalman filter)
│
├── mini-sensor-measurement (recommended)
│   ├── Inertial sensor principles
│   ├── Sensor error modeling
│   └── Allan variance analysis
│
├── mini-electromagnetic-wave (recommended)
│   ├── Radio wave propagation
│   ├── Ionospheric/tropospheric effects
│   └── Antenna theory (GNSS antenna)
│
└── mini-wireless-mobile-comm (optional)
    ├── Multipath fading
    ├── Diversity techniques
    └── RF front-end considerations
```

## Internal Module Dependency Tree

```
nav_common.h (no dependencies)
│
├── nav_rotation.h → nav_common.h
│   └── nav_rotation.c
│       ├── quaternion algebra
│       ├── DCM operations
│       ├── Euler conversions
│       └── TRIAD algorithm
│
├── nav_kalman.h → nav_common.h
│   └── nav_kalman.c
│       ├── matrix utilities (Cholesky, inverse)
│       ├── Linear Kalman filter
│       ├── Extended Kalman filter
│       ├── Information filter
│       └── Covariance intersection
│
├── nav_imu.h → nav_common.h, nav_rotation.h
│   └── nav_imu.c
│       ├── IMU error compensation
│       ├── INS mechanization (1-sample, 2-sample)
│       ├── Coarse alignment
│       ├── ZUPT
│       ├── NHC
│       └── INS error propagation
│
├── nav_gnss.h → nav_common.h
│   └── nav_gnss.c
│       ├── Satellite position (Kepler + perturbations)
│       ├── Azimuth/elevation
│       ├── Klobuchar iono model
│       ├── Saastamoinen tropo model
│       ├── WLS PVT solution
│       └── Coordinate conversions
│
├── nav_integration.h → all above
│   ├── nav_ins_gnss_loose.c
│   │   └── 15-state EKF, position/velocity updates
│   ├── nav_ins_gnss_tight.c
│   │   └── 17-state EKF, pseudorange/Doppler updates
│   └── lever arm compensation
│
├── nav_ahrs.c → nav_rotation.h, nav_imu.h
│   ├── Madgwick filter
│   ├── Mahony filter
│   └── Complementary filter
│
├── nav_integrity.c → nav_common.h, nav_gnss.h, nav_kalman.h
│   ├── RAIM residual detection
│   ├── RAIM FDE
│   ├── NIS test
│   └── Solution separation
│
├── nav_federated.c → nav_kalman.h
│   └── Federated Kalman filter with information sharing
│
└── nav_utils.c → nav_common.h, nav_rotation.h
    ├── Moving average filter
    ├── Outlier detection (MAD)
    ├── Solution interpolation
    ├── Haversine distance
    ├── Course over ground
    └── Speed over ground
```
