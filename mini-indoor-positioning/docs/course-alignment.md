# Course Alignment — mini-indoor-positioning

## Nine-School Curriculum Mapping

### MIT
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| 6.003 Signal Processing | Estimation theory, Kalman filtering | `sensor_fusion.c` — KF, EKF, UKF |
| 6.450 Digital Communications | RSSI models, channel characterization | `indoor_positioning.c` — path loss model |
| 6.630 EM Waves | Wave propagation, UWB | `tof_tdoa_positioning.c` — TOF/AoA |

### Stanford
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| EE102A Signal Processing | Fourier analysis, filtering | AHRS complementary filters |
| EE359 Wireless Communications | Indoor propagation, fingerprinting | `fingerprint_positioning.c` |
| EE267 Virtual Reality | IMU orientation, sensor fusion | `inertial_navigation.c` — Madgwick/Mahony |

### UC Berkeley
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| EE16A/B Circuits | Sensor interfacing | IMU calibration model |
| EE123 Digital Signal Processing | Adaptive filtering, LMS | Kalman filter implementation |
| EE117 EM Waves | Antenna arrays, AoA | `tof_tdoa_positioning.c` — AoA |

### Illinois (UIUC)
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| ECE 310 DSP | Estimation, spectral analysis | Allan variance |
| ECE 459 Communications | Wireless positioning fundamentals | RSSI/TDOA positioning |

### Michigan
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| EECS 351 DSP | Kalman filtering | `sensor_fusion.c` |
| EECS 455 Communications | UWB, ranging | `tof_tdoa_positioning.c` |
| EECS 411 Microwave | Link budget | `uwb_link_budget()` |

### Georgia Tech
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| ECE 4270 DSP | Optimal filtering | EKF, UKF |
| ECE 6601 Communications | Wireless localization | Trilateration, multilateration |

### TU Munich
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| Signal Processing | Bayesian estimation | Particle filter |
| Communications | Indoor channel models | Path loss model |
| High-Frequency Engineering | UWB technology | UWB ranging CRLB |

### ETH Zürich
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| 227-0427 Signal Processing | Estimation, detection | Error analysis, outlier detection |
| 227-0436 Communications | Positioning systems | Fingerprint, trilateration |
| 227-0455 EM Waves | Antenna, propagation | AoA, link budget |

### 清华 (Tsinghua)
| Course | Topic | Module Implementation |
|--------|-------|----------------------|
| 信号与系统 | 估计理论 | Kalman filter, error propagation |
| 通信原理 | 室内定位 | RSSI fingerprint, TDOA |
| 导航与定位 | 惯性导航 | Strapdown INS, PDR |
| 电磁场 | 天线阵列 | AoA estimation |

## Reference Textbooks

| Textbook | Author(s) | Module Coverage |
|----------|-----------|-----------------|
| Fundamentals of GPS Receivers | Tsui (2005) | Coordinate transforms, DOP |
| Principles of GNSS, Inertial, and Multisensor Navigation | Groves (2013) | INS, Kalman, sensor fusion |
| Wireless Communications | Molisch (2011) | Path loss, channel models |
| Digital Communications | Proakis & Salehi (2008) | Signal processing foundations |
| Antenna Theory | Balanis (2016) | AoA, array processing |
| RADAR: In-Building RF Localization | Bahl & Padmanabhan (2000) | WiFi fingerprinting |
| The Horus WLAN Location System | Youssef & Agrawala (2005) | Probabilistic fingerprinting |
