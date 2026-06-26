# mini-indoor-positioning

Indoor positioning and navigation module — comprehensive implementation of
wireless, inertial, and sensor-fusion based positioning algorithms.

## Module Status: COMPLETE ✅

- **L1-L6**: Complete (all core knowledge covered)
- **L7**: Complete — 3 applications (museum, warehouse, office)
- **L8**: Complete — 5 advanced topics (UKF, PF, Chan TDOA, NLOS, MARG AHRS)
- **L9**: Partial — documented, not implemented (research-stage)
- **include/ + src/ lines**: 6,327 (threshold: 3,000 ✅)
- **Score**: 17/18

## Quick Start

```bash
make test     # Build and run comprehensive test suite
make examples # Build example programs
make clean    # Remove build artifacts
```

## Core Definitions (L1)

| Type | Description |
|------|-------------|
| `position2d_t`, `position3d_t` | 2D/3D position in ENU coordinates |
| `velocity2d_t`, `velocity3d_t` | Velocity vectors |
| `quaternion_t` | Unit quaternion for attitude representation |
| `navigation_state_t` | Full nav state (position+velocity+attitude) |
| `measurement_t` | Generic sensor measurement |
| `path_loss_model_t` | Indoor RSSI path loss parameters |
| `radio_map_t` | WiFi fingerprint database |
| `imu_sample_t`, `imu_calibration_t` | IMU data and calibration |
| `kalman_filter_t`, `ekf_t`, `ukf_t` | Kalman filter variants |
| `particle_filter_t` | Particle filter (SMC) |
| `error_ellipse_t` | 2D positioning uncertainty |
| `dop_metrics_t` | Dilution of precision |

## Core Theorems & Laws (L4)

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Indoor Path Loss | RSSI(d) = RSSI₀ − 10n·log₁₀(d/d₀) | `rssi_to_distance()`, `distance_to_rssi()` |
| Cramér-Rao Lower Bound (RSSI) | CRLB = trace(FIM⁻¹) | `crlb_rssi_positioning()` |
| Cramér-Rao Lower Bound (TOF) | σd ≥ c/(2√2·π·B·√SNR) | `crlb_tof_positioning()`, `uwb_ranging_crlb()` |
| Error Propagation | Pxy = (HᵀWH)⁻¹ | `propagate_positioning_error()` |
| Kalman Optimality | MMSE for linear Gaussian systems | `kf_predict()`, `kf_update()` |
| DOP Relation | Position σ = σ₀ · DOP | `compute_dop()` |
| Triangle Inequality | |a−b| + |b−c| ≥ |a−c| | `distance_2d()`, `distance_3d()` |

## Core Algorithms (L5)

| # | Algorithm | Function | Complexity |
|---|-----------|----------|------------|
| 1 | Linear LS Trilateration (2D) | `trilateration_2d()` | O(N) |
| 2 | Gauss-Newton Trilateration (3D) | `trilateration_3d()` | O(K·N) |
| 3 | Chan TDOA Multilateration | `tdoa_multilateration()` | O(N²) |
| 4 | Taylor Series TDOA | `tdoa_taylor_series()` | O(K·N²) |
| 5 | Nearest-Neighbor Fingerprint | `fingerprint_match_nn()` | O(M·K) |
| 6 | k-NN Fingerprint | `fingerprint_match_knn()` | O(M·K+MlogM) |
| 7 | Weighted k-NN Fingerprint | `fingerprint_match_wknn()` | O(M·K+MlogM) |
| 8 | Probabilistic Fingerprint (Horus) | `fingerprint_match_probabilistic()` | O(M·K) |
| 9 | Strapdown INS | `ins_mechanize()` | O(1) |
| 10 | ZUPT Correction | `ins_apply_zupt()` | O(1) |
| 11 | Madgwick AHRS (IMU+MARG) | `madgwick_update_imu()`, `_marg()` | O(1) |
| 12 | Mahony AHRS | `mahony_update_imu()` | O(1) |
| 13 | Kalman Filter | `kf_predict()`, `kf_update()` | O(n³) |
| 14 | Extended Kalman Filter | `ekf_predict()`, `ekf_update()` | O(n³) |
| 15 | Unscented Kalman Filter | `ukf_predict()`, `ukf_update()` | O(n³) |
| 16 | Particle Filter (SIR) | `pf_predict/update/resample` | O(N) |
| 17 | TWR / SDS-TWR | `twr_compute_distance()` | O(1) |
| 18 | AoA Triangulation | `aoa_triangulate()` | O(1) |
| 19 | RANSAC Robust Positioning | `ransac_positioning()` | O(I·N) |
| 20 | NLOS Detection | `detect_nlos_rssi_distance()` | O(1) |
| 21 | Step Detection + PDR | `pdr_process_accel()` | O(1) |
| 22 | Allan Variance | `compute_allan_variance()` | O(N) |
| 23 | First-Path Detection (UWB) | `uwb_first_path_tof()` | O(N) |

## Canonical Problems (L6)

| Problem | Example | Description |
|---------|---------|-------------|
| WiFi RSSI Positioning | `example_rssi_trilateration.c` | 4 APs, 20m×15m office floorplan |
| IMU+UWB Sensor Fusion | `example_imu_fusion.c` | EKF tracking in 50m×30m warehouse |
| Fingerprint Navigation | `example_fingerprint_navigation.c` | 3-floor museum, 24 survey points, 6 APs |

## Nine-School Curriculum Mapping

| School | Key Courses | Module Coverage |
|--------|-------------|-----------------|
| **MIT** | 6.003, 6.450, 6.630 | Kalman, path loss, UWB propagation |
| **Stanford** | EE102A, EE359, EE267 | AHRS, fingerprinting, IMU orientation |
| **Berkeley** | EE16A/B, EE123 | Sensor calibration, adaptive filtering |
| **Illinois** | ECE 310, ECE 459 | Estimation theory, wireless positioning |
| **Michigan** | EECS 351, 455, 411 | Kalman, UWB, microwave links |
| **Georgia Tech** | ECE 4270, 6601 | Optimal filtering, localization |
| **TU Munich** | Signal Processing, Comm | Bayesian estimation, channel models |
| **ETH Zürich** | 227-0427, 0436, 0455 | Estimation, positioning, EM |
| **清华** | 信号与系统, 通信原理, 导航 | Kalman, fingerprint, INS |

## Directory Structure

```
mini-indoor-positioning/
├── README.md                         ← This file
├── Makefile                          ← make test/examples/clean
├── include/                          ← 6 header files
│   ├── indoor_positioning.h          ← Core types, trilateration, RSSI
│   ├── fingerprint_positioning.h     ← Radio map, matching algorithms
│   ├── inertial_navigation.h         ← IMU, INS, AHRS, PDR
│   ├── sensor_fusion.h              ← KF, EKF, UKF, PF
│   ├── tof_tdoa_positioning.h       ← TWR, TDOA, AoA, UWB
│   └── positioning_error.h          ← DOP, CRLB, error analysis
├── src/                              ← 7 source files (6 C + 1 Lean)
│   ├── indoor_positioning.c
│   ├── fingerprint_positioning.c
│   ├── inertial_navigation.c
│   ├── sensor_fusion.c
│   ├── tof_tdoa_positioning.c
│   ├── positioning_error.c
│   └── indoor_positioning.lean      ← Lean 4 formalization
├── tests/
│   └── test_positioning.c           ← 80+ test assertions
├── examples/
│   ├── example_rssi_trilateration.c
│   ├── example_imu_fusion.c
│   └── example_fingerprint_navigation.c
└── docs/
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-alignment.md
    └── course-tree.md
```

## Key Metrics

| Metric | Value |
|--------|-------|
| include/ files | 6 headers |
| src/ files | 6 C + 1 Lean = 7 |
| include/ + src/ lines | 6,327 |
| struct definitions | 25+ |
| implemented algorithms | 36 |
| end-to-end examples | 3 |
| knowledge levels covered | 8 Complete, 1 Partial |
| compilation | `gcc -std=c11 -Wall -Wextra` passes |
| Lean theorems | 10+ formal definitions + theorem statements |

## Reference Textbooks

- Tsui, "Fundamentals of Global Positioning System Receivers" (2005)
- Groves, "Principles of GNSS, Inertial, and Multisensor Navigation" (2013)
- Bahl & Padmanabhan, "RADAR: In-Building RF Localization" (2000)
- Youssef & Agrawala, "The Horus WLAN Location System" (2005)
- Madgwick, "Efficient Orientation Filter for IMU/MARG" (2010)
- Chan & Ho, "Efficient Estimator for Hyperbolic Location" (1994)
- Arulampalam et al., "Tutorial on Particle Filters" (2002)
- Dardari et al., "Ranging with UWB in Multipath" (2009)
- Kaplan & Hegarty, "Understanding GPS/GNSS" (2017)
- Molisch, "Wireless Communications" (2011)
