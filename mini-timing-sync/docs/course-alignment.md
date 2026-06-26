# Course Alignment — mini-timing-sync

## MIT
| Course | Topic | Mapping |
|--------|-------|---------|
| 6.450 Digital Comm | Symbol timing recovery, PLL | phase_detector.c |
| 6.003 Signal Processing | Kalman filtering, estimation | timing_sync.c Kalman functions |
| 6.082 Computer Networks | NTP, time synchronization | ntp_client.c |

## Stanford
| Course | Topic | Mapping |
|--------|-------|---------|
| EE359 Wireless Comm | Frequency/timing offset estimation | clock_model.c |
| EE102A Signal Processing | Allan variance, clock modeling | allan_variance.c |
| EE384S Digital Comm | PTP, IEEE 1588 | ptp_engine.c |

## Berkeley
| Course | Topic | Mapping |
|--------|-------|---------|
| EE123 DSP | PLL, phase detectors, DPLL | phase_detector.c |
| EE16A/B Circuits | Oscillator phase noise | allan_variance.c |
| EE290C Advanced Comm | Network synchronization | time_transfer.c |

## ETH Zurich
| Course | Topic | Mapping |
|--------|-------|---------|
| 227-0436 Comm | Clock recovery, synchronization | timing_sync.c |
| 227-0455 EM | Sagnac effect, TWSTFT | time_transfer.c |
| 227-0427 Signal Processing | Kalman filters, estimation | timing_sync.c |

## TU Munich
| Course | Topic | Mapping |
|--------|-------|---------|
| High-Frequency Eng | Oscillator characterization | allan_variance.c |
| Communications | Synchronization in wireless | clock_model.c |
| Signal Processing | PLL/DLL design | phase_detector.c |

## CMU
| Course | Topic | Mapping |
|--------|-------|---------|
| 18-345 Telecom Networks | NTP, PTP, synchronization | ntp_client.c, ptp_engine.c |
| 15-424 CPS | Clock invariants, formal proofs | timing_sync.lean |

## Tsinghua
| Course | Topic | Mapping |
|--------|-------|---------|
| Communication Principles | Clock sync, timing recovery | timing_sync.c, phase_detector.c |
| Signal and Systems | Frequency stability, Allan variance | allan_variance.c |

## Cambridge
| Course | Topic | Mapping |
|--------|-------|---------|
| Digital Comm | Symbol timing, PLL | phase_detector.c |
| Time and Frequency | Atomic clocks, stability analysis | allan_variance.c |
