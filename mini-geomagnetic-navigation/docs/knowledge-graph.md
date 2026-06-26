# Knowledge Graph -- mini-geomagnetic-navigation

## L1: Definitions

| Term | Definition | C Type | Location |
|------|-----------|--------|----------|
| Geodetic Coordinate | WGS84 (lat, lon, alt) | `GeodeticCoord` | `geomag_core.h` |
| ECEF Coordinate | Earth-Centered Earth-Fixed Cartesian | `ECEFCoord` | `geomag_core.h` |
| NED Coordinate | North-East-Down local tangent plane | `NEDCoord` | `geomag_core.h` |
| Magnetic Field Vector | B = (Bx, By, Bz) [nT] | `MagVector` | `geomag_core.h` |
| Magnetic Elements | D, I, F, H, X, Y, Z | `MagneticElements` | `geomag_core.h` |
| Gauss Coefficients | g_n^m, h_n^m for spherical harmonics | `GaussCoeff` | `geomag_core.h` |
| IGRF Model | International Geomagnetic Reference Field | `IGRFModel` | `geomag_core.h` |
| Magnetometer Types | Scalar, Triaxial, Gradiometer | `MagnetometerType` | `geomag_core.h` |
| Magnetic Anomaly | Observed - IGRF residual | `MagneticAnomaly` | `geomag_core.h` |
| Navigation Solution | Position, velocity, attitude | `NavSolution` | `geomag_core.h` |
| Magnetic Map | Gridded field data for MAGCOM | `MagneticMap` | `geomag_core.h` |
| Geomagnetic Activity | Kp, Dst, AE indices | `GeomagneticActivity` | `geomag_core.h` |
| Legendre State | Schmidt P_n^m computation state | `LegendreState` | `geomag_core.h` |
| Quaternion | Hamilton rotation quaternion | `Quaternion` | `geomag_math.h` |
| Kalman Filter | Linear KF state | `KalmanFilter` | `geomag_kalman.h` |
| Extended KF | EKF for nonlinear systems | `ExtendedKalmanFilter` | `geomag_kalman.h` |

## L2: Core Concepts

| Concept | Description | Implementation |
|---------|-------------|----------------|
| Spherical Harmonic Expansion | Scalar potential V(r,theta,phi) | `igrf_compute_field()` |
| Magnetic Compass | Heading from horizontal B components | `mag_heading_from_triaxial()` |
| Magnetic Map Matching | MAGCOM correlation positioning | `magcom_correlation_match()` |
| Magnetic Anomaly Detection | Threshold/OBF detectors | `mad_threshold_detect()` |
| Sensor Calibration | Hard-iron, soft-iron correction | `magnetometer_calibrate()` |
| Magnetic Gradiometry | Spatial gradient measurement | `compute_magnetic_gradient()` |
| INS/MAG Integration | EKF fusion of inertial + magnetic | `ins_mag_ekf_update()` |
| Geomagnetic Secular Variation | Time rate of change of main field | `igrf_compute_secular_variation()` |

## L3: Mathematical Structures

| Structure | Description | APIs |
|-----------|-------------|------|
| 3x3 Matrix Algebra | multiply, inverse, transpose, det | `mat3x3_*()` |
| Vector Algebra | dot, cross, norm, normalize | `vec3_*()` |
| Quaternion Algebra | multiply, conjugate, normalize, SLERP | `quat_*()` |
| Spherical Trigonometry | Great-circle distance, bearing | `great_circle_*()` |
| Rotation Matrices | ECEF-to-NED, quaternion-to-DCM | `ecef_to_ned_rotation()`, `quat_to_dcm()` |
| Bilinear Interpolation | 2D grid interpolation | `bilinear_interpolate_2d()` |

## L4: Fundamental Laws

| Law/Theorem | Statement | Verification |
|-------------|-----------|-------------|
| Laplace Equation (magnetostatics) | del^2 V = 0 above Earth surface | Basis for spherical harmonic expansion in IGRF |
| Maxwell: div B = 0 | Magnetic field is solenoidal | Gradient tensor trace = 0 check |
| IGRF Spherical Harmonic Synthesis | B = -grad(V) with V as SH expansion | `igrf_compute_field()` implementation |
| Magnetic Dipole Approximation | n=1 terms ~90% of field | `compute_dipole_field()` |
| Lowes-Mauersberger Spectrum | Spatial power spectrum of geomagnetic field | `compute_field_energy()` |
| Dipole Moment Equation | m = (4*pi*a^3/mu_0) * sqrt(g10^2+g11^2+h11^2) | `compute_dipole_moment()` |

## L5: Algorithms/Methods

| Algorithm | Complexity | Implementation |
|-----------|------------|----------------|
| Schmidt Legendre Recurrence | O(nmax^2) | `compute_schmidt_legendre()` |
| IGRF Field Synthesis | O(nmax^2) | `igrf_compute_field()` |
| Secular Variation Prediction | O(nmax^2) | `igrf_predict_field()` |
| Ellipsoid Fitting Calibration | O(N + 9^3) | `magnetometer_calibrate_ls()` |
| OBF Magnetic Anomaly Detection | O(N) | `mad_obf_detect()` |
| NCC Map Matching | O(N_search * N_points) | `magcom_correlation_match()` |
| MAD Map Matching | O(N_search * N_points) | `magcom_mad_match()` |
| Linear Kalman Filter | O(n^3 + n^2*m) | `kalman_predict/update()` |
| Extended Kalman Filter | O(n^3) | `ekf_predict/update()` |
| Golden-Section Search | O(log((b-a)/tol)) | `golden_section_search()` |
| Gradient Descent 2D | O(maxiter) | `gradient_descent_2d()` |
| Dipole Moment Estimation | O(N) | `estimate_dipole_moment()` |

## L6: Canonical Problems

| Problem | Solution | Example |
|---------|----------|---------|
| Compute B at any point on Earth | `igrf_compute_field()` | IGRF positioning demo |
| Magnetic Compass Correction | Declination from IGRF | Compass nav demo |
| Map-based Localization | MAGCOM NCC/MAD matching | MAGCOM demo |
| INS Drift Correction | Magnetic gradient EKF update | EKF integration |
| Magnetic Anomaly Detection | Threshold + OBF detectors | Sensor tests |
| Magnetic Pole Location | Iterative gradient search | `compute_magnetic_poles()` |
| Dipole Moment from Anomalies | Linear least squares | `estimate_dipole_moment()` |

## L7: Applications

| Application | Description | Status |
|-------------|-------------|--------|
| Underwater Vehicle Navigation | MAGCOM + INS for AUV | Partial (map matching + IGRF) |
| GPS-Denied Positioning | Magnetic field as navigation signal | Partial (IGRF single-point) |
| Aircraft Magnetic Backup | Compass + IGRF declination | Complete (compass nav demo) |
| Magnetic Survey | Anomaly detection + mapping | Complete (MAD + gradient) |

## L8: Advanced Topics

| Topic | Description | Status |
|-------|-------------|--------|
| Bayesian Magnetic SLAM | Simultaneous localization and magnetic mapping | Documented |
| Gradient Tensor Navigation | Full tensor invariants for positioning | Partial (NSS computation) |
| AI-Enhanced Map Matching | Machine learning for magnetic features | Documented |
| Multi-Sensor Fusion | EKF with magnetometer + IMU + GNSS | Partial (INS/MAG EKF) |

## L9: Research Frontiers

| Topic | Description | Status |
|-------|-------------|--------|
| Quantum Magnetometer Navigation | NV-diamond, SERF magnetometers for cm-level | Documented |
| 6G Magnetic Positioning | Magnetic MIMO for indoor localization | Documented |
| Swarm Magnetic Survey | Multi-UAV cooperative magnetic mapping | Documented |
