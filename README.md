# Mini Navigation & Positioning

A collection of **from-scratch, zero-dependency C implementations** of navigation, positioning, and timing algorithms for autonomous systems. Each module maps to MIT, Stanford, and other top-tier university courses, bridging theory and practice by translating textbook equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-geomagnetic-navigation](mini-geomagnetic-navigation/) | IGRF/WMM spherical harmonic field models, magnetic map matching (MAGCOM), gradient navigation, magnetometer hard-iron/soft-iron calibration, Kalman filter mag-aided INS | MIT 16.687, Stanford AA272C |
| [mini-gnss-gps](mini-gnss-gps/) | ECEF/geodetic coordinate transforms, pseudorange correction, Bancroft direct solver, iterative least-squares PVT, Hatch carrier smoothing, Klobuchar ionosphere, Saastamoinen troposphere, DOP, C/A code generation, Doppler, acquisition search | Stanford AA272C, MIT 16.687, Misra & Enge |
| [mini-indoor-positioning](mini-indoor-positioning/) | WiFi/BLE/Magnetic RSSI fingerprinting, pedestrial dead reckoning (PDR), IMU-based inertial navigation, ToF/TDoA/AOA time-based positioning, linear/extended Kalman filter sensor fusion, error metrics (DOP, CEP, RMSE) | CMU 16-833, Stanford CS225A, MIT 6.882 |
| [mini-inertial-navigation](mini-inertial-navigation/) | Quaternion/DCM/Euler attitude representations, IMU calibration (bias, scale factor, misalignment), strapdown mechanization equations, Schuler oscillation, Allan variance error characterization, GNSS/INS loose integration | MIT 16.485, Stanford AA272C, Titterton & Weston |
| [mini-integrated-navigation](mini-integrated-navigation/) | GNSS WLS position solution, satellite ephemeris computation, IMU mechanization, loosely/tightly/deeply coupled GNSS+INS architectures, Kalman filter framework, rotation representations | Stanford AA272C, MIT 16.687, Groves (2013) |
| [mini-slam-basics](mini-slam-basics/) | EKF-SLAM with landmark initialization, FastSLAM (Rao-Blackwellized particle filter), graph-based pose graph optimization, data association (Mahalanobis, JCBB), range-bearing sensor models, SE(2)/SE(3) Lie group poses | MIT 16.485, Stanford CS231A, Thrun et al. |
| [mini-timing-sync](mini-timing-sync/) | Allan variance and frequency stability, clock models (offset/skew/drift/aging), NTP client (RFC 5905), PTP engine (IEEE 1588-2019), phase/frequency detectors, digital PLL, time transfer (two-way, common-view, all-in-view) | NIST, IEEE 1588, MIT 6.241J |
| [mini-uwb-localization](mini-uwb-localization/) | IEEE 802.15.4a channel models (LOS/NLOS residential/office/industrial), Saleh-Valenzuela multipath, two-way ranging (TWR), TOA positioning, NLOS detection and mitigation, particle-filter tracking, self-contained linear algebra | IEEE 802.15.4a, Stanford CS444M, MIT 6.882 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and reference annotations
- **Practical demos** — GNSS receivers, INS mechanization, SLAM robot, NTP time server, UWB localizer

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-geomagnetic-navigation
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-navigation-positioning/
├── mini-geomagnetic-navigation/   # IGRF/WMM spherical harmonics, MAGCOM, mag-aided Kalman
├── mini-gnss-gps/                 # Pseudorange, Bancroft PVT, iono/tropo, C/A code, DOP
├── mini-indoor-positioning/       # WiFi/BLE fingerprinting, PDR, ToF/TDoA sensor fusion
├── mini-inertial-navigation/      # Quaternion/DCM, strapdown mechanization, IMU calibration
├── mini-integrated-navigation/    # GNSS+INS coupling architectures, Kalman filter framework
├── mini-slam-basics/              # EKF-SLAM, FastSLAM, pose graph optimization, data association
├── mini-timing-sync/              # NTP, PTP, Allan variance, clock models, digital PLL
└── mini-uwb-localization/         # IEEE 802.15.4a channel, TWR ranging, NLOS mitigation
```

## License

MIT
