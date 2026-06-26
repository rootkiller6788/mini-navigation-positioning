# Knowledge Graph — mini-timing-sync

## L1: Definitions (Complete)
| # | Definition | C Type | Lean Type |
|---|-----------|--------|-----------|
| 1 | Timestamp (seconds + nanoseconds) | `Timestamp` | `ClockTime` |
| 2 | Time offset (slave - master) | `TimeOffset` | `TimeOffset` |
| 3 | Clock skew (frequency offset) | `ClockSkew` | - |
| 4 | Clock drift (aging rate) | `ClockDrift` | - |
| 5 | Clock state vector | `ClockState` | `ClockState` |
| 6 | PTP port state | `PtpPortState` | `PtpPortState` |
| 7 | Sync status | `SyncStatus` | `SyncStatus` |
| 8 | Holdover configuration | `HoldoverConfig` | - |
| 9 | NTP stratum | `NtpStratum` | - |
| 10 | Clock quality | `ClockQuality` | - |
| 11 | PTP timestamps {t1,t2,t3,t4} | `PtpTimestamps` | - |
| 12 | NTP packet header | `NtpPacket` | - |
| 13 | Oscillator specification | `OscillatorSpec` | - |
| 14 | Clock parameters (x0,y0,D,A) | `ClockParameters` | - |
| 15 | Phase error | `PhaseDetectorOutput` | - |
| 16 | Allan variance result | `AllanResult` | - |
| 17 | Noise type classification | `NoiseType` | - |
| 18 | Time transfer link | `TimeTransferLink` | - |
| 19 | Asymmetry model | `AsymmetryModel` | - |

## L2: Core Concepts (Complete)
| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Two-way time transfer | `timing_compute_offset_delay()` |
| 2 | PI clock servo/discipline | `pi_servo_update()` |
| 3 | Holdover management | `holdover_should_enter()`, `holdover_can_exit()` |
| 4 | BMCA (Best Master Clock) | `bmca_compare_datasets()` |
| 5 | NTP clock filter | `ntp_clock_filter()` |
| 6 | DPLL phase locking | `dpll_update()`, `dpll_lock_detect()` |
| 7 | PLL loop dynamics | `pll_compute_dynamics()` |
| 8 | Timestamp validation | `ptp_validate_timestamps()` |

## L3: Mathematical Structures (Complete)
| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Clock state-space model (3-state) | `clock_state_transition()` |
| 2 | Kalman filter on clock state | `clock_kalman_predict()`, `clock_kalman_update()` |
| 3 | Least-squares regression | `clock_linear_regression()` |
| 4 | Polynomial clock model | `clock_model_time_error()` |
| 5 | Power-law spectral density | `NoiseCoefficients` |
| 6 | Marzullo interval endpoints | `MarzulloEndpoint` (in ntp_client.c) |

## L4: Fundamental Laws (Complete)
| # | Law/Theorem | Formula | Code | Lean |
|---|------------|---------|------|------|
| 1 | Two-way offset equation | offset=((t2-t1)-(t4-t3))/2 | `timing_compute_offset_delay()` | `two_way_offset_symmetry` |
| 2 | Two-way delay equation | delay=((t2-t1)+(t4-t3))/2 | same function | `two_way_delay_symmetry` |
| 3 | NTP offset equation | offset=((T2-T1)+(T3-T4))/2 | `timing_ntp_offset_delay()` | - |
| 4 | Allan variance formula | sigma_y^2=(1/2M)*sum(diffs^2) | `allan_variance_compute()` | `allan_variance_nonneg` |
| 5 | Clock offset evolution | x(t)=x0+y0*t+D*t^2/2 | `clock_model_time_error()` | `clock_offset_linear_evolution` |
| 6 | PTP causality | t1<t2, t3<t4 | `ptp_validate_timestamps()` | `ptp_causality` |
| 7 | Sagnac effect | dt=(2*w_E*A)/c^2 | `sagnac_correction_ns()` | - |
| 8 | Hadamard drift insensitivity | delta^2 annihilates linear | `hadamard_variance()` | `hadamard_drift_insensitivity` |

## L5: Algorithms/Methods (Complete)
| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Overlapping Allan variance | `allan_variance_compute()` |
| 2 | Modified Allan variance | `modified_allan_variance()` |
| 3 | Hadamard variance | `hadamard_variance()` |
| 4 | Time deviation (TDEV) | `time_deviation_compute()` |
| 5 | Noise coefficient fitting | `allan_fit_noise_coefficients()` |
| 6 | Noise type identification | `allan_identify_noise_type()` |
| 7 | BMCA comparison | `bmca_compare_datasets()` |
| 8 | Marzullo intersection | `ntp_intersection_algorithm()` |
| 9 | NTP clock combining | `ntp_combine_clocks()` |
| 10 | NTP clock filter | `ntp_clock_filter()` |
| 11 | Common-view GPS time transfer | `common_view_transfer()` |
| 12 | All-in-view GPS time transfer | `all_in_view_transfer()` |
| 13 | Kalman predict/update | `clock_kalman_predict()`, `clock_kalman_update()` |
| 14 | Least-squares clock regression | `clock_linear_regression()` |
| 15 | Phase detectors (linear, XOR, PFD, Hogge) | various in phase_detector.c |

## L6: Canonical Problems (Complete)
| # | Problem | Example |
|---|---------|---------|
| 1 | PTP master-slave synchronization | `example_ptp_sync.c` |
| 2 | Oscillator stability characterization | `example_allan_variance.c` |
| 3 | NTP client-server sync | `example_ntp_time_transfer.c` |
| 4 | Holdover management during GNSS outage | `holdover_estimate_uncertainty()` |
| 5 | Cycle slip detection | `detect_cycle_slip()` |
| 6 | PLL lock detection | `dpll_lock_detect()` |

## L7: Applications (Complete - 4 applications)
| # | Application | Implementation |
|---|------------|---------------|
| 1 | 5G fronthaul timing (ITU-T G.8271.1) | `ptp_5g_fronthaul_check()` |
| 2 | IEC 61850 power grid synchrophasor | `ptp_power_grid_check()` |
| 3 | GPSDO stability specification | `gpsdo_stability_check()` |
| 4 | MiFID II financial timestamping | `ntp_mifid_compliance()` |

## L8: Advanced Topics (Partial - 3/5)
| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Three-cornered hat method | `three_cornered_hat()` |
| 2 | White Rabbit sub-ns timing | `white_rabbit_link_model()`, `white_rabbit_phase_track()` |
| 3 | TWSTFT delay modeling | `twstft_delay_model()` |

## L9: Research Frontiers (Partial - documented)
| # | Topic | Status |
|---|-------|--------|
| 1 | Quantum clock synchronization | Documented in gap-report |
| 2 | Optical clock networks | Documented |
| 3 | Chip-scale atomic clocks for timing | Documented |
