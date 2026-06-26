# Course Tree — mini-indoor-positioning

## Prerequisite Dependency Graph

```
L1: Definitions
├── Coordinate Systems (ENU, geodetic)
├── Distance Metrics (Euclidean 2D/3D)
├── RSSI, ToF, TDoA, AoA concepts
├── Position/Velocity/Attitude representation
└── Error Metrics (CEP, RMS, DRMS)

    ↓ requires

L2: Core Concepts
├── Trilateration (distance intersection)
├── Multilateration (hyperbolic intersection)
├── Fingerprint Positioning (pattern matching)
├── Dead Reckoning (inertial integration)
├── Sensor Fusion (multi-source estimation)
└── Map Matching (spatial constraints)

    ↓ requires

L3: Mathematical Structures
├── Quaternion Algebra (rotation representation)
├── Covariance Matrices (uncertainty)
├── Coordinate Transforms (geodetic ↔ ENU)
├── Signal Distance Metrics (Euclidean, Manhattan, Cosine)
└── DOP Geometry (linearized positioning)

    ↓ requires

L4: Fundamental Laws
├── Friis Path Loss (RSSI-distance relationship)
├── CRLB (estimation lower bound)
├── Error Propagation (uncertainty through nonlinear functions)
├── Kalman Optimality (MMSE for linear Gaussian)
└── Triangle Inequality (distance geometry)

    ↓ requires

L5: Algorithms & Methods
├── Trilateration: Linear LS, Gauss-Newton
├── TDOA: Chan algorithm, Taylor series
├── Fingerprint: NN, k-NN, WKNN, Probabilistic
├── INS: Strapdown mechanization, ZUPT
├── AHRS: Madgwick, Mahony filters
├── Kalman: Linear KF, EKF, UKF
├── Particle Filter: SIR resampling
├── AoA: Phase-based, triangulation, Stansfield
├── NLOS: RSSI-distance check, residual test
├── Robust: RANSAC, innovation outlier test
└── PDR: Step detection, stride estimation

    ↓ requires

L6: Canonical Problems
├── WiFi RSSI Trilateration (office, 4 APs)
├── IMU+UWB EKF Fusion (warehouse tracking)
├── Fingerprint Navigation (museum, 3 floors)
├── UWB Link Budget Analysis
├── Magnetic Field Positioning
└── Pedestrian Dead Reckoning

    ↓ requires

L7: Applications
├── Museum Visitor Guidance (WiFi fingerprint, 24 survey points)
├── Warehouse Logistics (UWB+IMU, Toyota-style)
└── Office WiFi Positioning (RSSI, 20m×15m)

    ↓ extends to

L8: Advanced Topics
├── Unscented Kalman Filter (sigma-point propagation)
├── Particle Filter (sequential Monte Carlo)
├── Chan Two-Step TDOA (hyperbolic WLS)
├── NLOS Mitigation (weighted+RANSAC)
└── Madgwick MARG (9-DOF AHRS)

    ↓ extends to

L9: Research Frontiers
├── 5G NR cm-level Positioning
├── Deep Neural Radio Maps
├── Cooperative Multi-Agent SLAM
├── Quantum-Enhanced Positioning
├── Terahertz Indoor Localization
└── 6G RIS-Aided Positioning
```

## Compilation Order

1. `indoor_positioning.h/c` — core types, RSSI model, trilateration, coordinate transforms
2. `fingerprint_positioning.h/c` — radio map, matching algorithms (depends on 1)
3. `inertial_navigation.h/c` — IMU, INS, AHRS, PDR (depends on 1)
4. `sensor_fusion.h/c` — KF, EKF, UKF, PF, complementary (depends on 1)
5. `tof_tdoa_positioning.h/c` — TWR, TDOA, AoA, UWB (depends on 1)
6. `positioning_error.h/c` — DOP, CRLB, error analysis (depends on 1)
