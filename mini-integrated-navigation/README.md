# mini-integrated-navigation

Integrated Navigation: INS/GNSS Kalman Filter Integration Module.

## Module Status: COMPLETE ✅

- L1-L6: Complete
- L7: Complete (4 applications: UAV, automotive, smartphone, aviation)
- L8: Complete (3 advanced topics: federated KF, CI, info filter)
- L9: Partial (documented, not implemented)

## Overview

This module provides a comprehensive implementation of integrated navigation
systems, combining Inertial Navigation Systems (INS) with Global Navigation
Satellite Systems (GNSS) via Kalman filtering. Covers loose, tight, and
deep coupling architectures, AHRS algorithms, integrity monitoring, and
multi-sensor fusion.

## File Structure

```
mini-integrated-navigation/
├── Makefile              # make test builds and runs tests
├── README.md             # This file
├── include/              # 6 header files
│   ├── nav_common.h      # Core types, constants, coordinate frames
│   ├── nav_rotation.h    # Quaternion, DCM, Euler, TRIAD
│   ├── nav_kalman.h      # Linear KF, EKF, Info filter, CI
│   ├── nav_imu.h         # INS mechanization, alignment, ZUPT
│   ├── nav_gnss.h        # GNSS PVT, sat position, atmospheric models
│   └── nav_integration.h # Loose/tight coupling architectures
├── src/                  # 10 implementation files
│   ├── nav_rotation.c    # Quaternion algebra, SLERP, TRIAD
│   ├── nav_kalman.c      # KF/EKF/info filter/CI implementations
│   ├── nav_imu.c         # Strapdown mechanization, coning/sculling
│   ├── nav_gnss.c        # WLS positioning, sat orbit, iono/tropo
│   ├── nav_ins_gnss_loose.c  # 15-state EKF loose coupling
│   ├── nav_ins_gnss_tight.c  # 17-state EKF tight coupling
│   ├── nav_ahrs.c        # Madgwick, Mahony, complementary filters
│   ├── nav_integrity.c   # RAIM, FDE, NIS, solution separation
│   ├── nav_federated.c   # Federated Kalman filter
│   └── nav_utils.c       # MA filter, outlier detection, Haversine
├── tests/
│   └── test_nav.c        # 17 tests covering all modules
├── examples/
│   ├── example_ahrs.c            # Madgwick AHRS demonstration
│   ├── example_loose_coupling.c  # Loose INS/GNSS simulation
│   └── example_gnss_pvt.c        # GNSS PVT and coordinate transforms
└── docs/
    ├── knowledge-graph.md    # L1-L9 knowledge coverage table
    ├── coverage-report.md    # Detailed per-level assessment
    ├── gap-report.md         # Missing knowledge points
    ├── course-alignment.md   # Nine-school curriculum mapping
    └── course-tree.md        # Prerequisite dependency tree
```

## Line Count

- include/ + src/: **3130 lines**
- 6 header files, 10 source files
- All files compile cleanly with -Wall -Wextra

## Core Definitions (L1)

- Coordinate frames: ECI, ECEF, NED, ENU, Body, Wander
- WGS84 ellipsoid: geodetic position, radii of curvature, normal gravity
- 15-state INS error vector (psi-angle formulation)
- IMU measurement and error model (IEEE Std 952/1293)
- GNSS satellite measurement (pseudorange, Doppler, carrier phase)
- GNSS solution: DOP metrics (GDOP, PDOP, HDOP, VDOP, TDOP)
- Attitude representations: quaternion (Hamilton), DCM, Euler, axis-angle
- Kalman filter structures: state, covariance, transition, measurement matrices

## Core Theorems (L4)

1. **Kalman Filter Optimality** — Under linear-Gaussian assumptions, the KF
   produces the minimum mean-square error (MMSE) state estimate.
   `x̂ = E[x|z]`, `P = E[(x-x̂)(x-x̂)ᵀ]`

2. **Information Filter Duality** — The information form `Y = P⁻¹`,
   `y = Y·x` provides linear update: `Y⁺ = Y + HᵀR⁻¹H`

3. **INS Error Dynamics (psi-angle)** — The 15-state INS error model:
   `δṙ = -ω_en×δr + δv`
   `δv̇ = C_bⁿ·δfᵇ + fⁿ×ψ - (2ω_ieⁿ+ω_enⁿ)×δv + δg`
   `ψ̇ = -ω_inⁿ×ψ - C_bⁿ·δω_ibᵇ`

4. **Covariance Intersection** — Consistent fusion under unknown correlation:
   `P_CI⁻¹ = ωP_A⁻¹ + (1-ω)P_B⁻¹`

5. **GPS Satellite Orbit** — Kepler equation `E = M + e·sin(E)` with
   second-harmonic perturbation corrections (IS-GPS-200)

6. **Klobuchar Ionospheric Model** — `I_z = 5ns + A·cos(2π(t-φ)/P)`

## Core Algorithms (L5)

| Algorithm | File | Complexity |
|-----------|------|------------|
| Linear Kalman Filter | nav_kalman.c | O(n³ + m³) |
| Extended Kalman Filter | nav_kalman.c | O(n³ + m³) |
| Information Filter | nav_kalman.c | O(n³ + n²m) |
| Sequential Scalar Update | nav_kalman.c | O(n²) per element |
| Cholesky Decomposition | nav_kalman.c | O(n³/6) |
| TRIAD Attitude Determination | nav_rotation.c | O(1) |
| SLERP Interpolation | nav_rotation.c | O(1) |
| Madgwick AHRS | nav_ahrs.c | O(1) per sample |
| Mahony AHRS | nav_ahrs.c | O(1) per sample |
| INS Mechanization | nav_imu.c | O(1) per IMU sample |
| 2-Sample Coning/Sculling | nav_imu.c | O(1) |
| WLS GNSS PVT | nav_gnss.c | O(Kn_iter·m³) |
| GPS Satellite Position | nav_gnss.c | O(Kn_iter) |
| RAIM Residual Detection | nav_integrity.c | O(m³) |
| RAIM FDE | nav_integrity.c | O(K·m⁴) |
| Federated Kalman Filter | nav_federated.c | O(L·n³) |
| Covariance Intersection | nav_kalman.c | O(K·n³) |

## Nine-School Curriculum Mapping

| School | Key Courses | Module Coverage |
|--------|------------|-----------------|
| MIT | 6.003, 6.450, 16.06 | KF, GNSS, state estimation |
| Stanford | EE102A, EE359, AA272 | Estimation, GPS, navigation |
| Berkeley | EE123, EE117, EECS206A | KF, GNSS, sensor fusion |
| Illinois | ECE 310, ECE 459 | Filtering, spread spectrum |
| Michigan | EECS 351, AEROSP 551 | KF, INS/GNSS integration |
| Georgia Tech | ECE 4270, AE 6530 | Estimation, INS mechanization |
| TU Munich | SP, Comm, Navigation | KF, GNSS, integration |
| ETH | 227-0427, 151-0569 | Filtering, navigation systems |
| Tsinghua | 信号与系统, 导航原理 | Signal theory, navigation |

## Build and Test

```bash
make all       # Build static library libnav.a
make test      # Build and run test suite (17 tests)
make examples  # Build example programs
make clean     # Remove build artifacts
```

## References

1. Groves, P.D. (2013). *Principles of GNSS, Inertial, and Multisensor
   Integrated Navigation Systems*. Artech House.
2. Jekeli, C. (2001). *Inertial Navigation Systems with Geodetic
   Applications*. De Gruyter.
3. Titterton, D.H. & Weston, J.L. (2004). *Strapdown Inertial Navigation
   Technology*. IET.
4. Misra, P. & Enge, P. (2011). *Global Positioning System: Signals,
   Measurements, and Performance*. Ganga-Jamuna Press.
5. Kalman, R.E. (1960). "A New Approach to Linear Filtering and Prediction
   Problems." *ASME Journal of Basic Engineering*, 82, 35-45.
6. Maybeck, P.S. (1979). *Stochastic Models, Estimation, and Control*.
   Academic Press.
7. Julier, S.J. & Uhlmann, J.K. (1997). "A Non-divergent Estimation
   Algorithm in the Presence of Unknown Correlations." *ACC*.
8. Madgwick, S.O.H. (2010). "An Efficient Orientation Filter for Inertial
   and Inertial/Magnetic Sensor Arrays." Technical Report.
9. Carlson, N.A. (1990). "Federated Square Root Filter for Decentralized
   Parallel Processes." *IEEE Trans. AES*, 26(3).
10. Goshen-Meskin, D. & Bar-Itzhack, I.Y. (1992). "Unified Approach to
    Inertial Navigation System Error Modeling." *JGCD*, 15(3).
