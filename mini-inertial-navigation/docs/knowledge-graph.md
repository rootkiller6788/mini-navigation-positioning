# Knowledge Graph - mini-inertial-navigation

## L1: Definitions (Complete)
- 20 typedefs and enums covering all core INS types
- Coordinate frames (ECI, ECEF, NED, ENU, BODY)
- Vector/matrix types, position types, attitude types
- IMU sensor types, error model types, calibration types
- Kalman filter state type, GPS measurement types

## L2: Core Concepts (Complete)
- Inertial reference frames and coordinate transforms
- Earth rotation and Coriolis acceleration
- Quaternion algebra for attitude representation
- Strapdown velocity and position update equations
- Deterministic sensor errors (bias, scale factor, misalignment)
- IMU performance grades and applications

## L3: Math Structures (Complete)
- Quaternion group (SO(3) double cover)
- Quaternion kinematics ODE
- Skew-symmetric matrices (so(3) Lie algebra)
- Direction cosine matrix orthogonality
- Euler angle gimbal lock analysis
- Coriolis pseudo-force in rotating frame

## L4: Fundamental Laws (Complete)
- Newton laws in rotating NED frame
- Quaternion norm preservation theorem
- Strapdown velocity equation derivation
- Somigliana closed-form gravity (WGS84)
- Schuler period: T_s = 5067 s
- INS error propagation (cubic gyro dominance)
- Bortz coning error: pi*A^2 per cycle

## L5: Algorithms/Methods (Complete)
- Quaternion integration (Euler, exact, Bortz coning)
- Sculling compensation (2-sample and 4-sample)
- Strapdown mechanization (attitude+velocity+position)
- Level alignment from accelerometer readings
- Allan variance computation and decomposition
- 6-position static accelerometer calibration
- Multi-position gyroscope calibration (LS)
- Cross-axis misalignment calibration
- Thermal calibration (polynomial LS)
- Error-state Kalman filter (15-state)
- GPS measurement update (loosely coupled)
- ZUPT detection (SHOE) and Kalman update

## L6: Canonical Problems (Complete)
- AHRS: Complementary filter attitude estimation
- Strapdown INS: Full trajectory from IMU data
- Pedestrian ZUPT: Zero-velocity aided dead reckoning
- GPS/INS loosely coupled integration loop

## L7: Applications (Complete)
- AHRS for UAV/drone stabilization
- Boeing 747 navigation-grade INS accuracy
- Pedestrian indoor tracking (first responder)
- Smartphone orientation (consumer IMU)
- Automotive dead reckoning (GPS-denied)
- GNSS outage budget estimation

## L8: Advanced Topics (Partial)
- [Done] Particle filter for non-Gaussian INS errors
- [Done] Transfer alignment via velocity matching
- [TODO] Tightly coupled GPS/INS
- [TODO] In-motion alignment
- [TODO] Lever arm compensation

## L9: Research Frontiers (Partial)
- Chip-scale atomic IMUs (documented)
- Quantum inertial sensors (documented)
- AI-enhanced drift compensation (documented)
