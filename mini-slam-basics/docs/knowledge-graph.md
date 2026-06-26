# Knowledge Graph — mini-slam-basics

## L1: Definitions (Complete)

| # | Definition | C Type | Lean Type | Status |
|---|-----------|--------|-----------|--------|
| 1 | Robot Pose SE(2) | `slam_pose2d_t` | `Pose2D` | ✓ |
| 2 | Robot Pose SE(3) | `slam_pose3d_t` | — | ✓ |
| 3 | 2D Point Landmark | `slam_landmark2d_t` | `Landmark2D` | ✓ |
| 4 | 3D Point Landmark | `slam_landmark3d_t` | — | ✓ |
| 5 | Range-Bearing Observation | `slam_obs_rb_t` | `Observation` | ✓ |
| 6 | LiDAR Scan | `slam_lidar_scan_t` | — | ✓ |
| 7 | Odometry | `slam_odometry_t` | — | ✓ |
| 8 | Velocity Control | `slam_velocity_t` | — | ✓ |
| 9 | Feature Map | `slam_map2d_t` | — | ✓ |
| 10 | Occupancy Grid | `slam_occgrid_t` | `CellState` | ✓ |
| 11 | EKF State | `slam_ekf_state_t` | `SLAMState` | ✓ |
| 12 | Particle (FastSLAM) | `slam_particle_t` | — | ✓ |
| 13 | Pose Graph Node/Edge | `slam_graph_node_t/edge_t` | — | ✓ |
| 14 | SLAM Config/Metrics | `slam_config_t/metrics_t` | `SLAMConfig` | ✓ |
| 15 | Backend/Sensor/DA Enums | 7 enum types | — | ✓ |

## L2: Core Concepts (Complete)

| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | SLAM probabilistic formulation | EKF/FastSLAM/Graph backends | ✓ |
| 2 | Online vs Full SLAM | EKF (online), Graph (full) | ✓ |
| 3 | State augmentation | `slam_ekf_augment()` | ✓ |
| 4 | Marginalization | Information form utilities | ✓ |
| 5 | Rao-Blackwellization | FastSLAM particle+EKF factorization | ✓ |
| 6 | Factor graph representation | Pose graph nodes + edges | ✓ |
| 7 | Loop closure detection | `slam_da_detect_loop()` | ✓ |
| 8 | Motion model types | Velocity, odometry, diff-drive | ✓ |
| 9 | Sensor model types | Range-bearing, LiDAR, camera enums | ✓ |
| 10 | SLAM system lifecycle | Init → Run → LoopClose → Relocalize | ✓ |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation | Status |
|---|----------|---------------|--------|
| 1 | SE(2) Lie group composition | `slam_pose_compose()` | ✓ |
| 2 | SE(2) inverse | `slam_pose_inverse()` | ✓ |
| 3 | SE(2) to homogeneous matrix | `slam_se2_to_matrix()` | ✓ |
| 4 | Rotation matrix + derivative | `slam_rotation_matrix_2d()`, `_derivative_2d()` | ✓ |
| 5 | Point transformation (fwd/inv) | `slam_transform_point()` | ✓ |
| 6 | Innovation covariance | `slam_innovation_covariance_rb()` | ✓ |
| 7 | Mahalanobis distance | `slam_mahalanobis_sq()` | ✓ |
| 8 | Covariance ↔ Information | `slam_cov_to_info()`, `slam_info_vector()` | ✓ |
| 9 | Cholesky decomposition (3×3) | `slam_cholesky_3x3()` | ✓ |
| 10 | Matrix multiply patterns | `slam_matmul()`, `_AT_B()`, `_A_BT()` | ✓ |

## L4: Fundamental Laws/Theorems (Complete)

| # | Theorem | C Verification | Lean Statement | Status |
|---|---------|--------------|----------------|--------|
| 1 | EKF linearization | EKF predict/update Jacobians | — | ✓ |
| 2 | SLAM monotonicity (Dissanayake 2001) | `slam_ekf_verify_monotonicity()` | `variance_non_increasing` | ✓ |
| 3 | Bayes filter for SLAM | Full EKF predict-update cycle | — | ✓ |
| 4 | Cramér-Rao Lower Bound | NEES/NIS consistency checks | — | ✓ |
| 5 | Total probability (particle filter) | FastSLAM importance sampling | `neff_equal_weights` | ✓ |
| 6 | SE(2) group axioms | Compose+inverse → identity | `compose_identity_left/right` | ✓ |
| 7 | Landmark count monotonic | State augmentation | `landmark_count_non_decreasing` | ✓ |
| 8 | State dimension bound | `state_dimension_minimum` | ✓ |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation | Complexity | Status |
|---|----------|---------------|------------|--------|
| 1 | EKF-SLAM Predict | `slam_ekf_predict_velocity()` | O(d²) | ✓ |
| 2 | EKF-SLAM Update (range-bearing) | `slam_ekf_update_rb()` | O(d²) | ✓ |
| 3 | EKF State Augmentation | `slam_ekf_augment()` | O(d²) | ✓ |
| 4 | EKF Step (full cycle) | `slam_ekf_step()` | O(d²) | ✓ |
| 5 | FastSLAM 1.0 (motion proposal) | `slam_fastslam1_sample_pose()` | O(M) | ✓ |
| 6 | FastSLAM 1.0 (landmark update) | `slam_fastslam1_update_landmark()` | O(M·K) | ✓ |
| 7 | FastSLAM 2.0 (improved proposal) | `slam_fastslam2_sample_pose()` | O(M·K) | ✓ |
| 8 | Systematic Resampling | `slam_fastslam_resample()` | O(M) | ✓ |
| 9 | Low-Variance Resampling | `slam_fastslam_low_variance_resample()` | O(M) | ✓ |
| 10 | SE(2) Graph Error | `slam_graph_error_se2()` | O(1) | ✓ |
| 11 | SE(2) Graph Jacobian | `slam_graph_jacobian_se2()` | O(1) | ✓ |
| 12 | Gauss-Newton Optimization | `slam_graph_optimize_gauss_newton()` | O(N³) | ✓ |
| 13 | Levenberg-Marquardt | `slam_graph_optimize_lm()` | O(N³) | ✓ |
| 14 | Nearest Neighbor DA | `slam_da_nearest_neighbor()` | O(N) | ✓ |
| 15 | JCBB | `slam_da_jcbb()` | O(Nᴷ) | ✓ |
| 16 | ICP (2D) | `slam_da_icp_2d()` | O(N²) | ✓ |
| 17 | Occupancy Grid Update | `slam_occgrid_update_scan()` | O(K·R) | ✓ |
| 18 | LiDAR Likelihood Field | `slam_lidar_likelihood_field()` | O(K·R) | ✓ |

## L6: Canonical Problems (Complete)

| # | Problem | Example | Status |
|---|---------|---------|--------|
| 1 | Landmark-based SLAM (EKF) | `example_ekf_slam.c` — Circular path, 8 landmarks | ✓ |
| 2 | Pose Graph w/ Loop Closure | `example_graph_slam.c` — Loop closure correction | ✓ |
| 3 | Particle Filter SLAM | `example_fastslam.c` — Figure-8, 10 landmarks | ✓ |
| 4 | NEES/NIS Consistency | `slam_ekf_nees()`, `slam_ekf_nis()` | ✓ |
| 5 | Occupancy Grid Mapping | `slam_occgrid_update_scan()` | ✓ |
| 6 | Scan Matching (ICP) | `slam_da_icp_2d()` | ✓ |

## L7: Applications (Partial+)

| # | Application | Implementation | Status |
|---|------------|---------------|--------|
| 1 | Indoor robot navigation | EKF-SLAM: `example_ekf_slam.c` | ✓ |
| 2 | Autonomous vehicle mapping | Graph SLAM: `example_graph_slam.c` | ✓ |
| 3 | AR/VR spatial tracking | FastSLAM: `example_fastslam.c` | ✓ |
| 4 | GPS-denied navigation | Loop closure: `slam_loop_closure_detect()` | Partial |
| 5 | Warehouse AGV | Config + sensor models | Partial |

## L8: Advanced Topics (Partial+)

| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | Levenberg-Marquardt optimization | `slam_graph_optimize_lm()` | ✓ |
| 2 | Huber robust kernel | `slam_graph_huber_kernel()` | ✓ |
| 3 | Information form (SEIF foundations) | `slam_cov_to_info()` | Partial |
| 4 | Scan context descriptors | `scan_descriptor_t` in loop_closure | Partial |
| 5 | Multi-robot SLAM | Not implemented | Missing |

## L9: Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | 6G RIS-assisted SLAM | Documented only |
| 2 | Semantic SLAM (object-level) | Struct types defined |
| 3 | Learned feature descriptors | `signature[32/64]` fields |
| 4 | Quantum SLAM | Not implemented |
| 5 | Lifelong map maintenance | Not implemented |
