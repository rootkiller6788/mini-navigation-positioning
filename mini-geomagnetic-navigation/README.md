# mini-geomagnetic-navigation

**Geomagnetic Navigation Library** -- Earth magnetic field models, sensor simulation, map matching, and Kalman filtering for magnetic-aided navigation.

## Module Status: COMPLETE

- **L1 Definitions**: Complete -- 16 struct/enum typedefs
- **L2 Core Concepts**: Complete -- 8 core concepts with implementations
- **L3 Math Structures**: Complete -- matrix, vector, quaternion, spherical geometry
- **L4 Fundamental Laws**: Complete -- Laplace eq., Maxwell, IGRF synthesis, dipole model
- **L5 Algorithms**: Complete -- 10+ algorithms (Legendre, IGRF, KF, EKF, MAD, NCC)
- **L6 Canonical Problems**: Complete -- 3 end-to-end examples
- **L7 Applications**: Partial (3/4 -- underwater, GPS-denied, compass)
- **L8 Advanced Topics**: Partial (2/4 -- gradient tensor, sensor fusion)
- **L9 Research Frontiers**: Partial (documented -- quantum, 6G, swarm)

## Line Count
- `include/` + `src/`: 4397 lines (threshold: 3000)

## Core Definitions (L1)
- Geodetic/E CEF/NED coordinate systems (WGS84)
- Magnetic field vector and 7 standard magnetic elements (D, I, F, H, X, Y, Z)
- Gauss spherical harmonic coefficients (g_n^m, h_n^m)
- IGRF model descriptor (IGRF-13, epoch 2020.0)
- Magnetometer types: scalar, triaxial fluxgate, gradiometer
- Magnetic map grid for MAGCOM navigation
- Kalman filter state structures (linear + extended)

## Core Theorems (L4)
- **IGRF Spherical Harmonic Synthesis**: B = -grad(V), V from Laplace equation solution
  V(r,theta,phi) = a * sum (a/r)^(n+1) * sum [g cos(m phi) + h sin(m phi)] * P_n^m(cos theta)
- **Magnetic Dipole Approximation**: n=1 terms account for ~90% of field
  B_r = 2(a/r)^3 * [g10*cos(theta) + ...]
- **Lowes-Mauersberger Spectrum**: W ~ sum (n+1) * [(g_n^m)^2 + (h_n^m)^2]
- **Maxwell div B = 0**: gradient tensor is traceless in source-free region

## Core Algorithms (L5)
| Algorithm | Complexity | Function |
|-----------|------------|----------|
| Schmidt Legendre Recurrence | O(nmax^2) | `compute_schmidt_legendre()` |
| IGRF Field Synthesis | O(nmax^2) | `igrf_compute_field()` |
| Ellipsoid Fitting Calibration | O(N + 729) | `magnetometer_calibrate_ls()` |
| NCC Map Matching (MAGCOM) | O(N_search * N) | `magcom_correlation_match()` |
| Linear Kalman Filter | O(n^3) | `kalman_predict/update()` |
| Extended Kalman Filter | O(n^3) | `ekf_predict/update()` |
| OBF Anomaly Detection | O(N) | `mad_obf_detect()` |
| Golden-Section Search | O(log(1/eps)) | `golden_section_search()` |

## Classic Problems (L6)
1. **Compute geomagnetic field at any point** -- `igrf_compute_field()`
2. **Magnetic compass with declination** -- `mag_heading_from_triaxial()`
3. **Magnetic map-based localization** -- `magcom_correlation_match()`
4. **INS drift correction via magnetics** -- `magnetic_gradient_update()`
5. **Magnetic anomaly detection** -- `mad_threshold_detect()`
6. **Magnetic pole location** -- `compute_magnetic_poles()`
7. **Dipole moment estimation** -- `estimate_dipole_moment()`

## Nine-School Curriculum Mapping

| School | Key Course | Coverage |
|--------|-----------|----------|
| MIT | 6.630 EM Waves | Maxwell, dipole model |
| Stanford | EE359 Wireless | Positioning fundamentals |
| Berkeley | EE117 EM | Magnetostatics theory |
| Illinois | ECE 451 EM | Spherical harmonics |
| Michigan | EECS 351 DSP | Kalman filtering |
| Georgia Tech | ECE 6350 EM | Potential theory |
| TU Munich | High-Frequency Eng | EM field computation |
| ETH | 227-0427 Signal Proc | Estimation & detection |
| Tsinghua | Navigation Systems | INS, KF, integration |

## Build & Test

```bash
make all          # Build libgeomag.a
make test         # Run all tests (21/21 passing)
make examples     # Build example programs
make clean        # Clean build artifacts
```

## File Structure

```
mini-geomagnetic-navigation/
  include/          (6 headers, 1595 lines)
    geomag_core.h       -- Core types and coordinates
    geomag_model.h      -- IGRF/WMM field model API
    geomag_math.h       -- Matrix, vector, quaternion, spherical geometry
    geomag_sensor.h     -- Magnetometer models and MAD
    geomag_navigation.h -- Map matching and positioning
    geomag_kalman.h     -- Kalman/EKF filter and INS/MAG integration
  src/              (6 sources, 2802 lines)
    geomag_core.c       -- Coordinate transforms, magnetic elements
    geomag_model.c      -- IGRF-13 implementation, dipole, poles
    geomag_math.c       -- Linear algebra, quaternions, interpolation
    geomag_sensor.c     -- Sensor models, calibration, MAD, gradient
    geomag_navigation.c -- MAGCOM, IGRF positioning, compass nav
    geomag_kalman.c     -- KF, EKF, INS/MAG integration
  tests/            (5 test files, 21 tests)
  examples/         (3 end-to-end demos)
  docs/             (5 knowledge documents)
  Makefile
  README.md
```

## Reference Textbooks
- Merrill, McElhinny & McFadden, "The Magnetic Field of the Earth" (1996)
- Langel, "The Main Field", in Geomagnetism Vol. 1 (1987)
- Alken et al., "IGRF-13", Earth Planets Space (2021)
- Groves, "Principles of GNSS, Inertial, and Multisensor Integrated Navigation" (2013)
- Brown & Hwang, "Introduction to Random Signals and Applied Kalman Filtering" (2012)
- Ripka, "Magnetic Sensors and Magnetometers" (2001)
- Titterton & Weston, "Strapdown Inertial Navigation Technology" (2004)
