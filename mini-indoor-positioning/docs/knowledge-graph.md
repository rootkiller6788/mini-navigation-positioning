# Knowledge Graph — mini-indoor-positioning

## L1: Definitions — Complete ✅

| # | Definition | C Implementation | Status |
|---|-----------|-----------------|--------|
| 1 | 2D Position (ENU Cartesian) | `position2d_t` in indoor_positioning.h | ✅ |
| 2 | 3D Position (ENU) | `position3d_t` in indoor_positioning.h | ✅ |
| 3 | Velocity (2D/3D) | `velocity2d_t`, `velocity3d_t` | ✅ |
| 4 | Quaternion (attitude) | `quaternion_t` in indoor_positioning.h | ✅ |
| 5 | Navigation State | `navigation_state_t` | ✅ |
| 6 | Error Ellipse | `error_ellipse_t` in positioning_error.h | ✅ |
| 7 | Path Loss Model | `path_loss_model_t` | ✅ |
| 8 | RSSI → Distance conversion | `rssi_to_distance()` | ✅ |
| 9 | Distance → RSSI conversion | `distance_to_rssi()` | ✅ |
| 10 | RSSI Fingerprint | `rssi_reading_t`, `survey_point_t` | ✅ |
| 11 | Radio Map | `radio_map_t` in fingerprint_positioning.h | ✅ |
| 12 | WiFi Access Point | `access_point_t`, `mac_address_t` | ✅ |
| 13 | IMU Sample | `imu_sample_t` in inertial_navigation.h | ✅ |
| 14 | IMU Calibration | `imu_calibration_t` | ✅ |
| 15 | Magnetic Field Vector | `magnetic_vector_t` | ✅ |
| 16 | TWR Exchange | `twr_exchange_t` in tof_tdoa_positioning.h | ✅ |
| 17 | UWB Anchor | `uwb_anchor_t` | ✅ |
| 18 | AoA Measurement | `aoa_measurement_t` | ✅ |
| 19 | Kalman Filter | `kalman_filter_t` in sensor_fusion.h | ✅ |
| 20 | DOP Metrics | `dop_metrics_t` in positioning_error.h | ✅ |
| 21 | Positioning Accuracy | `positioning_accuracy_t` | ✅ |
| 22 | Error Statistics | `positioning_error_stats_t` | ✅ |
| 23 | Measurement | `measurement_t` | ✅ |
| 24 | Map Constraints | `map_constraint_t` | ✅ |
| 25 | Positioning Mode | `positioning_mode_t` enum | ✅ |

## L2: Core Concepts — Complete ✅

| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Trilateration (distance-based) | `trilateration_2d()`, `trilateration_3d()` | ✅ |
| 2 | Multilateration (TDOA-based) | `tdoa_multilateration()` (Chan algorithm) | ✅ |
| 3 | Fingerprint Positioning | NN, k-NN, WKNN, Probabilistic | ✅ |
| 4 | Dead Reckoning | Strapdown INS, PDR, step detection | ✅ |
| 5 | Sensor Fusion | Kalman, EKF, UKF, Particle Filter | ✅ |
| 6 | Proximity/Centroid | `weighted_centroid_2d()` | ✅ |
| 7 | Map Matching | `map_constraint_t`, spatial distances | ✅ |
| 8 | Floor Estimation | `estimate_floor_from_rssi()` | ✅ |
| 9 | NLOS Detection | RSSI-distance consistency, residuals | ✅ |
| 10 | Clock Synchronization | TWR, SDS-TWR | ✅ |
| 11 | Error Budget | `decompose_error_sources()` | ✅ |

## L3: Mathematical Structures — Complete ✅

| # | Structure | Implementation | Status |
|---|-----------|---------------|--------|
| 1 | Euclidean Distance (2D/3D) | `distance_2d()`, `distance_3d()` | ✅ |
| 2 | Quaternion Algebra | multiply, conjugate, normalize, rotate | ✅ |
| 3 | Euler Angles ↔ Quaternion | `euler_to_quaternion()`, `quaternion_to_euler()` | ✅ |
| 4 | Geodetic → ENU Transform | `geodetic_to_enu()` (WGS-84) | ✅ |
| 5 | 2×2/3×3/N×N Matrices | `matrix2_t`, `matrix3_t`, `matrix_t` | ✅ |
| 6 | Covariance Propagation | Kalman predict/update equations | ✅ |
| 7 | Cholesky Decomposition | `cholesky_decompose()` | ✅ |
| 8 | Mahalanobis Distance | `mahalanobis_distance()` | ✅ |
| 9 | Jacobian Computation | EKF Jacobians, trilateration Jacobian | ✅ |
| 10 | Signal Space Distance | Euclidean, Manhattan, Cosine | ✅ |
| 11 | Rotation via Quaternion | `quaternion_rotate_vector()` | ✅ |
| 12 | Allan Variance | `compute_allan_variance()` | ✅ |

## L4: Fundamental Laws — Complete ✅

| # | Law/Theorem | Implementation | Status |
|---|------------|---------------|--------|
| 1 | Friis Free-Space Path Loss (adapted indoor) | `rssi_to_distance()`, `distance_to_rssi()` | ✅ |
| 2 | Cramer-Rao Lower Bound (RSSI) | `crlb_rssi_positioning()` | ✅ |
| 3 | Cramer-Rao Lower Bound (TOF) | `crlb_tof_positioning()` | ✅ |
| 4 | Error Propagation Law | `propagate_positioning_error()` | ✅ |
| 5 | DOP Geometric Interpretation | `compute_dop()` | ✅ |
| 6 | Distance Triangle Inequality | `distance_2d()`, `distance_3d()` | ✅ |
| 7 | RSSI Monotonicity (proven in Lean) | `indoor_positioning.lean` | ✅ |
| 8 | Kalman Optimality (MMSE for linear Gaussian) | `kf_predict()`, `kf_update()` | ✅ |
| 9 | UWB Ranging CRLB | `uwb_ranging_crlb()` | ✅ |
| 10 | Bayesian Estimation (particle filter convergence) | `pf_update_gaussian()` | ✅ |

## L5: Algorithms/Methods — Complete ✅

| # | Algorithm | Implementation | Complexity | Status |
|---|-----------|---------------|------------|--------|
| 1 | Linear LS Trilateration (2D) | `trilateration_2d()` | O(N) | ✅ |
| 2 | Gauss-Newton Trilateration (3D) | `trilateration_3d()` | O(K·N) | ✅ |
| 3 | Chan TDOA Multilateration | `tdoa_multilateration()` | O(N²) | ✅ |
| 4 | Taylor Series TDOA | `tdoa_taylor_series()` | O(K·N²) | ✅ |
| 5 | Nearest-Neighbor Fingerprint | `fingerprint_match_nn()` | O(M·K) | ✅ |
| 6 | k-NN Fingerprint | `fingerprint_match_knn()` | O(M·K+MlogM) | ✅ |
| 7 | Weighted k-NN Fingerprint | `fingerprint_match_wknn()` | O(M·K+MlogM) | ✅ |
| 8 | Probabilistic Fingerprint (Horus) | `fingerprint_match_probabilistic()` | O(M·K) | ✅ |
| 9 | Strapdown INS Mechanization | `ins_mechanize()` | O(1) | ✅ |
| 10 | ZUPT Correction | `ins_apply_zupt()` | O(1) | ✅ |
| 11 | Madgwick AHRS | `madgwick_update_imu()`, `_marg()` | O(1) | ✅ |
| 12 | Mahony AHRS | `mahony_update_imu()` | O(1) | ✅ |
| 13 | Quaternion Integration (1st order) | `quaternion_integrate_1st_order()` | O(1) | ✅ |
| 14 | Quaternion Integration (RK4) | `quaternion_integrate_rk4()` | O(1) | ✅ |
| 15 | Linear Kalman Filter | `kf_predict()`, `kf_update()` | O(n³) | ✅ |
| 16 | Extended Kalman Filter | `ekf_predict()`, `ekf_update()` | O(n³) | ✅ |
| 17 | Unscented Kalman Filter | `ukf_predict()`, `ukf_update()` | O(n³) | ✅ |
| 18 | Particle Filter (SIR) | `pf_predict()`, `pf_update_gaussian()`, `pf_resample()` | O(N) | ✅ |
| 19 | Complementary Filter (angle) | `complementary_filter_angle()` | O(1) | ✅ |
| 20 | Complementary Filter (position) | `complementary_filter_position()` | O(1) | ✅ |
| 21 | Two-Way Ranging | `twr_compute_distance()` | O(1) | ✅ |
| 22 | SDS-TWR (clock-drift-corrected) | `twr_sds_compute_distance()` | O(1) | ✅ |
| 23 | Gauss-Newton TOA (2D) | `toa_positioning_2d()` | O(K·N) | ✅ |
| 24 | AoA from Phase Difference | `aoa_from_phase()` | O(1) | ✅ |
| 25 | AoA Triangulation | `aoa_triangulate()` | O(1) | ✅ |
| 26 | Stansfield AoA Estimator | `aoa_positioning_stansfield()` | O(N) | ✅ |
| 27 | RANSAC Robust Positioning | `ransac_positioning()` | O(I·N) | ✅ |
| 28 | NLOS Detection (RSSI/Distance) | `detect_nlos_rssi_distance()` | O(1) | ✅ |
| 29 | NLOS Detection (Residual) | `detect_nlos_residual()` | O(N) | ✅ |
| 30 | Innovation Outlier Detection | `detect_measurement_outlier()` | O(m³) | ✅ |
| 31 | Step Detection (PDR) | `pdr_process_accel()` | O(1) | ✅ |
| 32 | Weinberg Stride Length | `pdr_stride_length_weinberg()` | O(1) | ✅ |
| 33 | Kim Stride Length | `pdr_stride_length_kim()` | O(1) | ✅ |
| 34 | Error Ellipse Computation | `compute_error_ellipse()` | O(1) | ✅ |
| 35 | Allan Variance | `compute_allan_variance()` | O(N) | ✅ |
| 36 | First-Path Detection (UWB CIR) | `uwb_first_path_tof()` | O(N) | ✅ |

## L6: Canonical Problems — Complete ✅

| # | Problem | Implementation | Status |
|---|---------|---------------|--------|
| 1 | WiFi RSSI Positioning in Office | `example_rssi_trilateration.c` | ✅ |
| 2 | IMU+UWB Sensor Fusion Tracking | `example_imu_fusion.c` | ✅ |
| 3 | Museum Fingerprint Navigation | `example_fingerprint_navigation.c` | ✅ |
| 4 | PDR with Step Detection | `pdr_process_accel()` + stride models | ✅ |
| 5 | EKF Range-Only Tracking | `ekf_setup_constant_velocity_ranging()` | ✅ |
| 6 | UWB Link Budget Analysis | `uwb_link_budget()` | ✅ |
| 7 | Magnetic Field Positioning | `magnetic_field_match()` | ✅ |

## L7: Applications — Complete ✅ (3 applications)

| # | Application | Implementation | Status |
|---|-------------|---------------|--------|
| 1 | Museum Visitor Guidance (WiFi Fingerprint) | `example_fingerprint_navigation.c` — 3-floor museum with 24 survey points | ✅ |
| 2 | Warehouse Logistics Tracking (UWB+IMU) | `example_imu_fusion.c` — EKF fusion, Toyota-style warehouse | ✅ |
| 3 | Office WiFi Positioning (RSSI Trilateration) | `example_rssi_trilateration.c` — 4 APs in 20m×15m office | ✅ |

## L8: Advanced Topics — Complete ✅ (4 topics)

| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | Unscented Kalman Filter | `ukf_predict()`, `ukf_update()` — sigma-point propagation without Jacobians | ✅ |
| 2 | Particle Filter (Sequential Monte Carlo) | `pf_predict/update/resample/get_mean/cov` — 500 particles, SIR | ✅ |
| 3 | Chan TDOA Algorithm | `tdoa_multilateration()` — Two-step WLS hyperbolic positioning | ✅ |
| 4 | NLOS Mitigation (Weighted+RANSAC) | `detect_nlos_*`, `ransac_positioning()`, `compute_nlos_weights()` | ✅ |
| 5 | Mahony/Madgwick AHRS Fusion | `madgwick_update_marg()`, `mahony_update_imu()` | ✅ |

## L9: Research Frontiers — Partial

| # | Topic | Status |
|---|-------|--------|
| 1 | 5G NR Positioning (cm-level) | Documented (see gap-report.md) |
| 2 | AI/ML Fingerprinting (Deep Neural Radio Maps) | Documented |
| 3 | Cooperative Multi-Agent Positioning | Documented |
| 4 | Quantum-Enhanced Positioning | Documented |
| 5 | Terahertz Indoor Localization | Documented |
| 6 | 6G RIS-Aided Positioning | Documented |
