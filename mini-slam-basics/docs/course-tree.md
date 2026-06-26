# Course Tree — mini-slam-basics

## Prerequisite Dependency Tree

```
SLAM Basics
│
├── Probability & Statistics ───────────────┐
│   ├── Bayes rule                          │
│   ├── Gaussian distributions              │
│   ├── Multivariate normal                 │
│   ├── Conditional independence            │
│   └── Chi-squared distribution            │
│                                           │
├── Linear Algebra ─────────────────────────┤
│   ├── Matrix multiplication               │
│   ├── Matrix inversion                    │
│   ├── Cholesky decomposition              │
│   ├── Singular Value Decomposition (SVD) │
│   └── Positive-definite matrices          │
│                                           │
├── State Estimation ───────────────────────┤
│   ├── Kalman Filter (linear)              │
│   ├── Extended Kalman Filter              │
│   ├── Particle Filter                     │
│   └── Information Filter                  │
│                                           │
├── Lie Groups & Geometry ──────────────────┤
│   ├── SO(2) rotation group                │
│   ├── SE(2) rigid transformations         │
│   ├── SE(3) 3D transformations            │
│   ├── Tangent space / Lie algebra         │
│   └── Manifold optimization               │
│                                           │
├── Optimization ───────────────────────────┤
│   ├── Nonlinear least squares             │
│   ├── Gauss-Newton method                 │
│   ├── Levenberg-Marquardt                 │
│   ├── Robust kernels (Huber, Cauchy)      │
│   └── Sparse linear solvers               │
│                                           │
├── Sensor Models ──────────────────────────┤
│   ├── Range-bearing sensor                │
│   ├── LiDAR / laser range finder          │
│   ├── Camera (mono, stereo, RGB-D)        │
│   ├── Odometry (wheel encoders)           │
│   └── IMU (accelerometer + gyroscope)     │
│                                           │
└── Motion Models ──────────────────────────┤
    ├── Velocity motion model               │
    ├── Odometry motion model               │
    ├── Differential drive kinematics       │
    └── Ackermann steering kinematics       │
                                             │
                    ┌────────────────────────┘
                    ▼
            SLAM Algorithms
            ┌──────┼──────┐
            ▼      ▼      ▼
          EKF   FastSLAM  Graph
          SLAM  1.0/2.0   SLAM
```

## What This Module Provides

| Component | Prerequisites Satisfied |
|-----------|------------------------|
| `slam_core.c` | SE(2) algebra, linear algebra, Cholesky |
| `slam_ekf.c` | EKF-SLAM: predict, update, augment, NEES/NIS |
| `slam_fastslam.c` | Particle filter + per-landmark EKFs |
| `slam_graph.c` | Gauss-Newton, LM, Huber kernel |
| `slam_sensor.c` | Motion + observation models, LiDAR |
| `slam_data_assoc.c` | NN, JCBB, ICP, loop detection |
| `slam_loop_closure.c` | Scan descriptors, loop pipeline |
| `slam_formal.lean` | Formal properties in Lean 4 |

## Downstream Modules

This module serves as a prerequisite for:
- `mini-gnss-gps` — Sensor fusion with GNSS measurements
- `mini-inertial-navigation` — IMU preintegration in SLAM
- `mini-integrated-navigation` — Multi-sensor fusion (EKF/UKF)
- `mini-indoor-positioning` — WiFi/BLE SLAM
- `mini-uwb-localization` — UWB range-based SLAM
