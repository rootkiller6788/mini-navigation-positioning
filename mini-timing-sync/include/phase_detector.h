/**
 * @file phase_detector.h
 * @brief Phase and frequency detectors, digital PLL components
 *
 * L1: PhaseError, FrequencyError, PhaseDetectorType
 * L2: PLL lock detection, phase-frequency detector operation
 * L3: Transfer functions in s-domain and z-domain
 * L4: PLL loop dynamics, lock range, pull-in range
 * L5: Digital phase detector algorithms, loop filter design
 * L6: PLL locking transient analysis
 *
 * Reference: Gardner, F.M. "Phaselock Techniques" (2005)
 *            Best, R.E. "Phase-Locked Loops" (2007)
 * Course: MIT 6.450, Stanford EE359, Berkeley EE123
 */

#ifndef PHASE_DETECTOR_H
#define PHASE_DETECTOR_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* L1: Phase detector type enumeration */
typedef enum {
    PD_TYPE_LINEAR      = 0,  /* Multiplier/linear detector */
    PD_TYPE_XOR         = 1,  /* XOR gate phase detector */
    PD_TYPE_PFD         = 2,  /* Phase-Frequency Detector (tri-state) */
    PD_TYPE_HOGGE       = 3,  /* Hogge phase detector (clock/data) */
    PD_TYPE_ALEXANDER   = 4,  /* Alexander (bang-bang) phase detector */
    PD_TYPE_TANLOCK     = 5   /* Tanlock phase detector */
} PhaseDetectorType;

/* L1: Phase detector measurement result */
typedef struct {
    double phase_error_rad;    /* Phase error in radians */
    double frequency_error_hz; /* Frequency error in Hz */
    double lock_indicator;     /* Lock confidence [0, 1] */
    int    cycle_slip_detected; /* 1 if cycle slip detected */
} PhaseDetectorOutput;

/* L2: Digital PLL loop filter coefficients */
typedef struct {
    double alpha;        /* Proportional path coefficient */
    double beta;         /* Integral path coefficient */
    double gamma;        /* Double-integral path (type-3 PLL) */
    double T_update_s;   /* Update interval [seconds] */
} LoopFilterCoeffs;

/* L2: Digital PLL state */
typedef struct {
    double phase_accumulator;    /* NCO phase accumulator */
    double frequency_word;       /* NCO frequency control word */
    double integral_path;        /* Loop filter integrator state */
    double double_integral_path; /* Type-3 loop filter second integrator */
    double lock_metric;          /* Lock detection metric */
    int    lock_count;           /* Consecutive lock samples */
} DpllState;

/* L3: PLL transfer function parameters */
typedef struct {
    double natural_freq_hz;   /* Natural frequency omega_n [Hz] */
    double damping_factor;    /* Damping factor zeta */
    double loop_bandwidth_hz; /* Loop noise bandwidth BL [Hz] */
    double phase_margin_deg;  /* Phase margin [degrees] */
} PllDynamics;

/**
 * L4: Compute PLL natural frequency and damping factor from
 * loop filter coefficients.
 *
 * For a type-2 PLL with PI loop filter:
 *   omega_n = sqrt(K * beta / T)
 *   zeta    = (alpha / 2) * sqrt(T / (K * beta))
 *
 * Where K is the combined VCO + phase detector gain.
 *
 * @param coeffs  Loop filter coefficients
 * @param K_gain  Combined detector + VCO gain [rad/s per rad]
 * @param dyn     [out] PLL dynamics parameters
 */
void pll_compute_dynamics(const LoopFilterCoeffs *coeffs, double K_gain,
                          PllDynamics *dyn);

/**
 * L4: Compute PLL lock range (hold-in range).
 *
 * For type-2 PLL: hold-in range is theoretically infinite
 * (limited only by VCO tuning range).
 * For type-1 PLL: lock_range = K_vco * K_pd * pi/2
 *
 * @param pd_gain   Phase detector gain [V/rad or LSB/rad]
 * @param vco_gain  VCO gain [Hz/V or Hz/LSB]
 * @param pd_type   Type of phase detector
 * @return Lock range in Hz (0 = infinite for type-2)
 */
double pll_lock_range(double pd_gain, double vco_gain,
                      PhaseDetectorType pd_type);

/**
 * L5: Compute phase error from two sinusoidal signals.
 *
 * For two signals: s1(t) = A1*sin(wt + phi1), s2(t) = A2*sin(wt + phi2)
 * The product after lowpass filtering gives:
 *   V_out = (A1*A2/2) * sin(phi1 - phi2) ~= K_pd * (phi1 - phi2) for small errors
 *
 * @param sample1  First signal sample
 * @param sample2  Second signal sample
 * @param K_pd     Phase detector gain
 * @return Phase error in radians
 */
double phase_error_from_samples(double sample1, double sample2, double K_pd);

/**
 * L5: XOR phase detector characteristic.
 *
 * XOR output duty cycle encodes phase difference.
 * For signals with 50% duty cycle: V_avg = V_high * |delta_phi| / pi
 * Range: [-pi/2, pi/2] for square waves.
 *
 * @param sig1_level  0 or 1 for first signal
 * @param sig2_level  0 or 1 for second signal
 * @param v_high      Logic high voltage
 * @return Analog phase error voltage
 */
double xor_phase_detector(int sig1_level, int sig2_level, double v_high);

/**
 * L5: Phase-Frequency Detector (PFD) output.
 *
 * PFD produces UP/DN pulses. Net output:
 *   Q = (UP pulse width - DN pulse width) / T_ref
 *
 * Range: [-2*pi, 2*pi], enabling frequency discrimination.
 *
 * @param ref_edge_time   Time of reference edge [s]
 * @param div_edge_time   Time of divided VCO edge [s]
 * @param T_ref           Reference period [s]
 * @return Normalized phase error output
 */
double pfd_output(double ref_edge_time, double div_edge_time, double T_ref);

/**
 * L5: Hogge phase detector for clock-and-data recovery (CDR).
 *
 * Uses two DFFs clocked on rising and falling data edges.
 * Output = (A - B) where A samples at data transition center
 * and B samples at data edge for timing error.
 *
 * @param data_sample     Current data sample value
 * @param prev_data       Previous data sample
 * @param edge_sample     Sample at nominal edge position
 * @return Phase error direction (+1 = early, -1 = late, 0 = no transition)
 */
int hogge_phase_detect(double data_sample, double prev_data,
                       double edge_sample);

/**
 * L6: Lock detection algorithm for digital PLL.
 *
 * Monitors phase error variance over a moving window.
 * Lock declared when variance < threshold for N_consecutive samples.
 *
 * @param phase_error_rad     Current phase error [rad]
 * @param lock_threshold_rad2 Variance threshold for lock [rad^2]
 * @param lock_count_needed   Required consecutive samples below threshold
 * @param dpll                DPLL state (lock_metric and lock_count updated)
 * @return 1 if PLL is locked, 0 if not
 */
int dpll_lock_detect(double phase_error_rad, double lock_threshold_rad2,
                     int lock_count_needed, DpllState *dpll);

/**
 * L6: Digital PLL update iteration.
 *
 * Full update cycle:
 * 1. Phase error detection
 * 2. Loop filtering (PI/lead-lag)
 * 3. NCO frequency/phase update
 *
 * @param phase_error_rad  Measured phase error
 * @param coeffs           Loop filter coefficients
 * @param dpll             DPLL state (updated in place)
 * @param freq_out_hz      [out] Current output frequency estimate
 */
void dpll_update(double phase_error_rad, const LoopFilterCoeffs *coeffs,
                 DpllState *dpll, double *freq_out_hz);

/**
 * Initialize DPLL state for a given nominal frequency.
 *
 * @param dpll             State to initialize
 * @param nominal_freq_hz  Nominal (free-running) frequency [Hz]
 * @param sampling_freq_hz NCO clock frequency [Hz]
 * @param nco_bit_width    NCO accumulator bit width
 */
void dpll_init(DpllState *dpll, double nominal_freq_hz,
               double sampling_freq_hz, int nco_bit_width);

/**
 * L6: Cycle slip detection.
 *
 * Detects when phase error jumps by ~2*pi, indicating a cycle slip.
 * Uses a simple threshold on the phase error difference.
 *
 * @param current_phase_error  Current phase error [rad]
 * @param prev_phase_error     Previous phase error [rad]
 * @param threshold_rad        Detection threshold (typically ~pi) [rad]
 * @return 1 if cycle slip detected, 0 otherwise
 */
int detect_cycle_slip(double current_phase_error, double prev_phase_error,
                      double threshold_rad);

#ifdef __cplusplus
}
#endif
#endif
