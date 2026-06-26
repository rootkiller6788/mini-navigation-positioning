# mini-slam-basics

Simultaneous Localization and Mapping (SLAM) — the problem of building a map of an unknown environment while simultaneously tracking the robot's pose within that map.

**Reference Textbooks:**
- Thrun, Burgard & Fox (2005) *Probabilistic Robotics*, MIT Press
- Dissanayake, Newman, Clark et al. (2001) *A Solution to the SLAM Problem*, IEEE T-RO
- Montemerlo, Thrun et al. (2002) *FastSLAM*, AAAI
- Dellaert & Kaess (2006) *Square Root SAM*, IJRR

---

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Partial+ (3 applications: indoor nav, autonomous vehicle, AR/VR)
- **L8**: Partial+ (LM optimization, Huber kernel, scan descriptors)
- **L9**: Partial (documented, not implemented)

| Level | Status | Score |
|-------|--------|-------|
| L1 Definitions | Complete | 2 |
| L2 Core Concepts | Complete | 2 |
| L3 Math Structures | Complete | 2 |
| L4 Fundamental Laws | Complete | 2 |
| L5 Algorithms/Methods | Complete | 2 |
| L6 Canonical Problems | Complete | 2 |
| L7 Applications | Partial | 1 |
| L8 Advanced Topics | Partial | 1 |
| L9 Research Frontiers | Partial | 1 |
| **Total** | | **15/18** |

**Line Count**: include/ + src/ = **6171 lines** (threshold: 3000) ✅

---

## Core Definitions

### Robot Pose
| Type | Symbol | Definition |
|------|--------|-----------|
| SE(2) pose | (x, y, θ) | 2D rigid transformation: position + heading |
| SE(3) pose | (x, y, z, qw, qx, qy, qz) | 3D rigid transformation: position + quaternion |

### Landmarks & Observations
| Type | Symbol | Definition |
|------|--------|-----------|
| 2D Landmark | m_j = (m_jx, m_jy) | Point feature in world frame |
| Range-bearing | z = (r, φ) | Range r ≥ 0, bearing φ ∈ [−π, π) relative to robot heading |
| Odometry | (δ_rot1, δ_trans, δ_rot2) | Relative motion between poses |

### State Representation
| Backend | State Vector | Dimension |
|---------|-------------|-----------|
| EKF-SLAM | μ = [x, y, θ, m1x, m1y, ..., mNx, mNy]^T | 3 + 2N |
| FastSLAM | M particles × (pose + per-landmark EKFs) | M × (3 + 2K) |
| Graph SLAM | N pose nodes + E edges | 3N variables |

---

## Core Theorems

### EKF-SLAM Convergence (Dissanayake et al. 2001)
```
The determinant of any submatrix of the map covariance matrix
decreases monotonically as successive observations are made.
In the limit, landmark estimates become fully correlated.
```

### SLAM Probabilistic Formulation
```
p(x_t, m | z_{1:t}, u_{1:t})
  = η · p(z_t | x_t, m) · ∫ p(x_t | x_{t-1}, u_t) · p(x_{t-1}, m | z_{1:t-1}, u_{1:t-1}) dx_{t-1}
```

### FastSLAM Factorization (Rao-Blackwellization)
```
p(x_{1:t}, Θ | z_{1:t}, u_{1:t})
  = p(x_{1:t} | z_{1:t}, u_{1:t}) · Π_j p(θ_j | x_{1:t}, z_{1:t})
```

### Pose Graph Optimization
```
X* = argmin_X Σ_{(i,j)∈E} e_{ij}(X)^T · Ω_{ij} · e_{ij}(X)
where e_{ij} = z_{ij} ⊖ (x_i^{-1} ∘ x_j)
```

### Motion Model (Velocity)
```
For |ω| > ε:  x' = x + (v/ω)(sin(θ+ωΔt) − sin(θ))
              y' = y + (v/ω)(cos(θ) − cos(θ+ωΔt))
              θ' = θ + ωΔt
For |ω| ≤ ε:  x' = x + v·Δt·cos(θ)
              y' = y + v·Δt·sin(θ)
```

### Observation Model (Range-Bearing)
```
r̂ = √((m_jx − x)² + (m_jy − y)²)
φ̂ = atan2(m_jy − y, m_jx − x) − θ
```

---

## Core Algorithms

| Algorithm | Function | Complexity | Reference |
|-----------|----------|------------|-----------|
| EKF-SLAM Predict | `slam_ekf_predict_velocity()` | O(d²) | Smith et al. (1990) |
| EKF-SLAM Update | `slam_ekf_update_rb()` | O(d²) | Dissanayake et al. (2001) |
| State Augmentation | `slam_ekf_augment()` | O(d²) | Thrun Ch.10 |
| EKF Full Step | `slam_ekf_step()` | O(d²) | Thrun Ch.10 |
| FastSLAM 1.0 Sampling | `slam_fastslam1_sample_pose()` | O(M) | Montemerlo et al. (2002) |
| FastSLAM Landmark Update | `slam_fastslam1_update_landmark()` | O(M·K) | Montemerlo et al. (2002) |
| FastSLAM 2.0 Proposal | `slam_fastslam2_sample_pose()` | O(M·K) | Montemerlo & Thrun (2003) |
| Systematic Resampling | `slam_fastslam_resample()` | O(M) | Kitagawa (1996) |
| Gauss-Newton Graph Opt | `slam_graph_optimize_gauss_newton()` | O(N³) | Lu & Milios (1997) |
| Levenberg-Marquardt | `slam_graph_optimize_lm()` | O(N³) | Kümmerle et al. (2011) |
| Nearest Neighbor DA | `slam_da_nearest_neighbor()` | O(N) | Bar-Shalom (1988) |
| JCBB | `slam_da_jcbb()` | O(Nᴷ) | Neira & Tardos (2001) |
| ICP Scan Matching | `slam_da_icp_2d()` | O(N²) | Besl & McKay (1992) |
| Occupancy Grid Update | `slam_occgrid_update_scan()` | O(K·R) | Thrun Ch.9 |

---

## Classic Problems

1. **Landmark-based EKF-SLAM** — `example_ekf_slam.c`
   Robot drives circular path, maps 8 landmarks via range-bearing EKF-SLAM

2. **Pose Graph with Loop Closure** — `example_graph_slam.c`
   Odometry chain with loop closure constraint, Gauss-Newton optimization

3. **Particle Filter FastSLAM** — `example_fastslam.c`
   Figure-8 trajectory, 30 particles, 10 landmarks, systematic resampling

4. **Filter Consistency (NEES/NIS)** — tests/test_slam.c
   Normalized error checks against χ² distribution

5. **Occupancy Grid Mapping** — `slam_occgrid_update_scan()`
   LiDAR scan integration with Bresenham ray casting

6. **ICP Scan Matching** — `slam_da_icp_2d()`
   2D point cloud alignment for loop closure constraints

---

## Nine-School Course Mapping

| School | Courses | SLAM Topics |
|--------|---------|-------------|
| MIT | 6.141 Robotics, 6.438 Inference | EKF-SLAM, factor graphs, particle filters |
| Stanford | CS 225A, CS 226, AA 274 | EKF, FastSLAM, pose graph optimization |
| Berkeley | EE C106A/B, CS 287 | Motion/sensor models, graph SLAM, ICP |
| Illinois | ECE 470, ECE 550 | SLAM formulation, pose graphs |
| Michigan | ROB 530, EECS 568 | EKF-SLAM, occupancy grids, feature SLAM |
| Georgia Tech | CS 7630, ECE 6557 | Full pipeline, industrial SLAM |
| TU Munich | IN 2356, IN 2138 | EKF/Graph SLAM, LiDAR, ICP |
| ETH Zurich | 151-0854, 263-5902 | EKF/FastSLAM, visual SLAM |
| Tsinghua | 机器人学, 计算机视觉 | Motion models, Kalman, visual SLAM |

---

## Build & Test

```bash
make          # Build library and test binary
make test     # Build and run all tests (25/25 pass)
make examples # Build all 3 examples
make clean    # Remove build artifacts
```

---

## File Structure

```
mini-slam-basics/
├── Makefile                          # Build system
├── README.md                         # This file (COMPLETE ✅)
├── include/
│   ├── slam_types.h                  (471 lines) - Core types, enums, structs
│   ├── slam_ekf.h                    (298 lines) - EKF-SLAM API
│   ├── slam_fastslam.h              (289 lines) - FastSLAM particle filter API
│   ├── slam_graph.h                  (280 lines) - Pose graph optimization API
│   ├── slam_sensor.h                 (345 lines) - Sensor/motion models API
│   └── slam_data_assoc.h            (267 lines) - Data association API
├── src/
│   ├── slam_core.c                   (704 lines) - SE(2) algebra, linear algebra
│   ├── slam_ekf.c                    (666 lines) - EKF-SLAM implementation
│   ├── slam_fastslam.c              (717 lines) - FastSLAM 1.0/2.0
│   ├── slam_graph.c                  (558 lines) - Gauss-Newton, LM, Huber
│   ├── slam_sensor.c                 (444 lines) - Motion, observation, LiDAR
│   ├── slam_data_assoc.c            (443 lines) - NN, JCBB, ICP
│   ├── slam_loop_closure.c          (257 lines) - Loop detection pipeline
│   └── slam_formal.lean             (432 lines) - Lean 4 formal verification
├── tests/
│   └── test_slam.c                   - 25 test groups, all passing
├── examples/
│   ├── example_ekf_slam.c            - EKF-SLAM circular path demo
│   ├── example_graph_slam.c          - Pose graph loop closure demo
│   └── example_fastslam.c            - FastSLAM particle filter demo
├── docs/
│   ├── knowledge-graph.md            - L1-L9 knowledge coverage table
│   ├── coverage-report.md            - Level-by-level assessment
│   ├── gap-report.md                 - Missing items and priorities
│   ├── course-alignment.md           - Nine-school curriculum mapping
│   └── course-tree.md                - Prerequisite dependency tree
├── demos/
└── benches/
```

---

## References

1. Thrun, S., Burgard, W., & Fox, D. (2005). *Probabilistic Robotics*. MIT Press.
2. Dissanayake, M.W.M.G., Newman, P., Clark, S., Durrant-Whyte, H.F., & Csorba, M. (2001). A Solution to the Simultaneous Localization and Map Building (SLAM) Problem. *IEEE Trans. Robotics & Automation*, 17(3), 229-241.
3. Montemerlo, M., Thrun, S., Koller, D., & Wegbreit, B. (2002). FastSLAM: A Factored Solution to the Simultaneous Localization and Mapping Problem. *AAAI*.
4. Montemerlo, M. & Thrun, S. (2003). Simultaneous Localization and Mapping with Unknown Data Association Using FastSLAM. *ICRA*.
5. Dellaert, F. & Kaess, M. (2006). Square Root SAM: Simultaneous Localization and Mapping via Square Root Information Smoothing. *IJRR*, 25(12), 1181-1203.
6. Kümmerle, R., Grisetti, G., Strasdat, H., Konolige, K., & Burgard, W. (2011). g2o: A General Framework for Graph Optimization. *ICRA*.
7. Lu, F. & Milios, E. (1997). Globally Consistent Range Scan Alignment for Environment Mapping. *Autonomous Robots*, 4, 333-349.
8. Smith, R., Self, M., & Cheeseman, P. (1990). Estimating Uncertain Spatial Relationships in Robotics. *Autonomous Robot Vehicles*, Springer.
9. Bailey, T. & Durrant-Whyte, H. (2006). Simultaneous Localization and Mapping (SLAM): Part II. *IEEE Robotics & Automation Magazine*, 13(3), 108-117.
10. Neira, J. & Tardos, J.D. (2001). Data Association in Stochastic Mapping Using the Joint Compatibility Test. *IEEE Trans. Robotics & Automation*, 17(6), 890-897.
