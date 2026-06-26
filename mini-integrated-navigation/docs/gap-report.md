# Gap Report — mini-integrated-navigation

## Missing Knowledge Points

### L3: Math Structures
- None (complete coverage)

### L4: Fundamental Laws
- Schmidt-Kalman filter (consider covariance) — not implemented
- Sequential probability ratio test (SPRT) for integrity — not implemented

### L5: Algorithms
- Unscented Kalman Filter (UKF) — not implemented (sigma-point approach)
- Particle filter for non-Gaussian navigation — not implemented
- Multi-rate Kalman filtering with time-delayed measurements — not implemented

### L6: Canonical Problems
- Deeply coupled (ultra-tight) INS/GNSS integration — documented only
- Visual-inertial odometry (VIO) — not in scope
- RTK/PPP carrier-phase positioning — not implemented

### L7: Applications
- Marine (ship) navigation — not covered
- Spacecraft attitude determination — not covered
- Indoor positioning (Wifi/BLE fusion) — not in scope

### L8: Advanced Topics
- Gaussian sum filters — not implemented
- Interacting Multiple Model (IMM) estimation — not implemented
- Factor graph optimization (iSAM/GTSAM) — not implemented

### L9: Research Frontiers
- All items are documentation-only (no implementation required for COMPLETE)

## Priority for Future Work

| Priority | Item | Effort |
|----------|------|--------|
| High | UKF implementation | Medium |
| High | Schmidt-Kalman filter | Medium |
| Medium | Particle filter | Large |
| Medium | Deep coupling architecture | Large |
| Low | IMM estimation | Medium |
| Low | Factor graph SLAM interface | Large |
