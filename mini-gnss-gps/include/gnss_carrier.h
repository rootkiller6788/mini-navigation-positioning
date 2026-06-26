#ifndef GNSS_CARRIER_H
#define GNSS_CARRIER_H
#include "gnss_common.h"
#include "gnss_pseudorange.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Carrier Phase Measurement
 *
 * Φ = ρ + c·(dt_r - dt_s) + λ·N - I + T + M_φ + ε_φ
 *
 * Φ   : carrier phase measurement [m]
 * λ·N : integer ambiguity term (N cycles × λ meters/cycle)
 * I   : ionospheric advance (negative sign — carrier is advanced)
 *       (Note: code is delayed by +I, carrier is advanced by -I)
 *
 * The ambiguity N is initially unknown and must be resolved for cm-level
 * positioning.
 * ========================================================================= */

/** Single-epoch carrier phase measurement for one satellite */
typedef struct {
    gnss_satid_t   sat_id;
    gnss_gpstime_t time_rx;
    double         carrier_phase;    /* Accumulated phase [cycles] */
    double         phase_range;      /* λ × cycles [m] */
    double         signal_strength;  /* C/N₀ [dB-Hz] */
    int            valid;            /* 1 = no cycle slip detected */
    int            cycle_slip_flag;  /* 1 = slip detected this epoch */
} gnss_carrier_phase_t;

/* -------------------------------------------------------------------------
 * L2: Integer Ambiguity
 * ------------------------------------------------------------------------- */

/** Integer ambiguity state for one satellite */
typedef struct {
    int    prn;
    double float_ambiguity;    /* Float estimate [cycles] */
    int    fixed_ambiguity;    /* Fixed integer [cycles] (-1 = not fixed) */
    double fix_ratio;          /* LAMBDA ratio test value */
    int    is_fixed;           /* 1 = successfully fixed to integer */
} gnss_ambiguity_t;

/** Double-differenced ambiguity set (for relative positioning) */
typedef struct {
    int     n_dd;              /* Number of double-difference pairs */
    int     ref_prn;           /* Reference satellite PRN */
    double *float_ambiguities; /* [n_dd] float estimates [cycles] */
    int    *fixed_ambiguities; /* [n_dd] fixed integers [cycles] */
    double  ratio;             /* LAMBDA ratio test */
    int     all_fixed;         /* 1 = all ambig fixed successfully */
} gnss_ambiguity_dd_t;

/* -------------------------------------------------------------------------
 * L2: Cycle Slip Detection
 * ------------------------------------------------------------------------- */

/** Cycle slip detector configuration */
typedef struct {
    double jump_threshold;     /* Minimum jump to flag as slip [cycles] (typ. 3) */
    double iqr_factor;         /* IQR multiplier for outlier detection */
    int    window_size;        /* Moving-window size for statistics */
} gnss_slip_detector_config_t;

/** Cycle slip detection status per satellite */
typedef struct {
    int     n_slips_detected;
    int    *slip_epochs;       /* Epoch indices where slips occurred */
    int     current_slip;      /* 1 = slip at current epoch */
    double  metric;            /* Detection metric value (e.g., GF diff) */
    double  threshold;         /* Current threshold */
} gnss_slip_detector_t;

/* -------------------------------------------------------------------------
 * L5: Carrier Smoothing — Hatch Filter (1982)
 * ------------------------------------------------------------------------- */

/**
 * @brief Hatch filter state for carrier-smoothed pseudorange
 *
 * Update equation:
 *   P_sm[n] = (1/N)·P[n] + (1 - 1/N)·(P_sm[n-1] + Φ[n] - Φ[n-1])
 *
 * Equivalently:
 *   P_sm[n] = Φ[n] + (P_sm[n-1] - Φ[n-1]) + (1/N)·(P[n] - P[n-1] - (Φ[n] - Φ[n-1]))
 *
 * where:
 *   P[n]    = current pseudorange measurement
 *   Φ[n]    = current carrier phase (in meters, λ × cycles)
 *   N       = smoothing window (typ. 100 for GPS L1, 200 epochs)
 *   P_sm[n] = smoothed pseudorange output
 *
 * Key property: reduces pseudorange noise by √N, but ionospheric
 * divergence limits the practical smoothing window.
 *
 * Reference: Hatch, R. (1982). "The Synergism of GPS Code and Carrier
 * Measurements." Proc. 3rd Int'l Geodetic Symp. on Satellite Doppler
 * Positioning.
 */
typedef struct {
    int     window_width;     /* Smoothing window N (epochs) */
    int     n_samples;        /* Number of samples accumulated */
    double  prev_p_sm;        /* P_sm[n-1] [m] */
    double  prev_phi;         /* Φ[n-1] [m] */
    int     initialized;      /* 1 = has previous values */
    int     reset_on_slip;    /* 1 = auto-reset on cycle slip */
} gnss_hatch_filter_t;

/* -------------------------------------------------------------------------
 * L5: Ionosphere-Free Linear Combination (dual-frequency)
 * ------------------------------------------------------------------------- */

/**
 * @brief Ionosphere-free combination coefficients
 *
 * For dual-frequency receivers (L1/L2):
 *   P_IF = (f₁²·P₁ - f₂²·P₂) / (f₁² - f₂²)
 *   Φ_IF = (f₁²·Φ₁ - f₂²·Φ₂) / (f₁² - f₂²)
 *
 * This eliminates first-order ionospheric delay (99.9%).
 * Trade-off: noise amplified by factor ≈ 3 for GPS L1/L2.
 */
typedef struct {
    double f1, f2;           /* Carrier frequencies [Hz] */
    double coeff_p1;         /* f₁² / (f₁² - f₂²) */
    double coeff_p2;         /* -f₂² / (f₁² - f₂²) */
    double noise_factor;     /* √(coeff_p1² + coeff_p2²) */
} gnss_ionofree_combo_t;

/** Geometry-free combination (detects ionospheric TEC) */
typedef struct {
    double f1, f2;           /* Carrier frequencies [Hz] */
    double coeff_phi1;       /* coefficient for Φ₁ */
    double coeff_phi2;       /* coefficient for Φ₂ */
} gnss_geomfree_combo_t;

/** Melbourne-Wübbena combination (wide-lane ambiguity resolution) */
typedef struct {
    double f1, f2;           /* Carrier frequencies [Hz] */
    double coeff_p;          /* Narrow-lane code coefficient */
    double coeff_phi;        /* Wide-lane phase coefficient */
    double wide_lane_lambda; /* Wide-lane wavelength ≈ 0.86 m (GPS L1/L2) */
} gnss_mw_combo_t;

/* -------------------------------------------------------------------------
 * L6: Differential GPS (DGPS) corrections
 * ------------------------------------------------------------------------- */

/** Pseudorange correction from a reference station */
typedef struct {
    int32_t base_prn;
    double  prc;             /* Pseudorange correction [m] */
    double  rrc;             /* Range-rate correction [m/s] */
    double  iod;             /* Issue of data */
    double  udre;            /* User Differential Range Error [m] */
} gnss_dgps_correction_t;

/* -------------------------------------------------------------------------
 * API: Carrier phase processing
 * ------------------------------------------------------------------------- */

/**
 * @brief Initialize Hatch carrier smoothing filter
 *
 * @param f         Filter state
 * @param window_n  Smoothing window width (epochs)
 * @param reset_on_slip  Auto-reset on cycle slip detection
 */
void gnss_hatch_init(gnss_hatch_filter_t *f, int window_n, int reset_on_slip);

/**
 * @brief Apply one epoch of Hatch carrier smoothing
 *
 * @param f        Filter state
 * @param pseudorange  Current pseudorange [m]
 * @param phase_range  Current carrier phase in meters (λ·cycles) [m]
 * @param cycle_slip   1 = cycle slip detected → reset filter
 * @return Smoothed pseudorange [m]
 */
double gnss_hatch_smooth(gnss_hatch_filter_t *f,
                          double pseudorange, double phase_range,
                          int cycle_slip);

/** @brief Reset Hatch filter (e.g., after cycle slip or satellite change) */
void gnss_hatch_reset(gnss_hatch_filter_t *f);

/* -------------------------------------------------------------------------
 * API: Cycle slip detection
 * ------------------------------------------------------------------------- */

/**
 * @brief Geometry-free cycle slip detection
 *
 * GF = λ₁·Φ₁ - λ₂·Φ₂ = (1 - f₁²/f₂²)·I + (λ₁·N₁ - λ₂·N₂)
 *
 * GF changes only when ionospheric TEC changes or cycle slips occur.
 * TEC changes slowly (~0.01 cycle/s), so large jumps indicate cycle slips.
 *
 * Returns absolute GF change from previous epoch; > threshold → cycle slip.
 */
double gnss_slip_detect_gf(double phi1_cycles, double phi2_cycles,
                            double prev_phi1, double prev_phi2,
                            double lambda1, double lambda2);

/**
 * @brief Cycle slip detection from phase-rate comparison
 *
 * Compares current Doppler-derived range-rate with phase change rate.
 * Large discrepancy indicates cycle slip.
 *
 * @return detection metric (cycles); > threshold → slip
 */
double gnss_slip_detect_doppler(double delta_phi_cycles,
                                 double doppler_hz, double dt);

/**
 * @brief Initialize a cycle-slip detector
 */
void gnss_slip_detector_init(gnss_slip_detector_t *detector);

/**
 * @brief Run full cycle-slip detection suite
 *
 * @return 1 if slip detected, 0 if clean
 */
int gnss_slip_detect(gnss_slip_detector_t *detector,
                      double phi1, double phi2,
                      double prev_phi1, double prev_phi2,
                      double doppler_hz, double dt,
                      const gnss_slip_detector_config_t *cfg);

/* -------------------------------------------------------------------------
 * API: Linear combinations
 * ------------------------------------------------------------------------- */

gnss_ionofree_combo_t gnss_ionofree_combination(double f1, double f2);
double gnss_ionofree_code(double P1, double P2, const gnss_ionofree_combo_t *ic);
double gnss_ionofree_phase(double Phi1, double Phi2,
                            const gnss_ionofree_combo_t *ic);
double gnss_ionofree_noise_factor(const gnss_ionofree_combo_t *ic);

gnss_geomfree_combo_t gnss_geomfree_combination(double f1, double f2);
double gnss_geomfree_phase(double Phi1_m, double Phi2_m,
                            const gnss_geomfree_combo_t *gc);

gnss_mw_combo_t gnss_mw_combination(double f1, double f2);
double gnss_mw_wide_lane(double P_m, double Phi1_m, double Phi2_m,
                          const gnss_mw_combo_t *mc);

/* -------------------------------------------------------------------------
 * API: DGPS
 * ------------------------------------------------------------------------- */

/**
 * @brief Apply DGPS correction to a pseudorange
 *
 * P_corrected = P_raw + PRC + RRC · (t - t0)
 */
double gnss_dgps_apply(double pseudorange, const gnss_dgps_correction_t *corr,
                        double time_of_applicability);

/**
 * @brief Compute position-domain DGPS correction from reference station
 */
int gnss_dgps_position_correction(gnss_ecef_t ref_known,
                                   gnss_ecef_t ref_computed,
                                   double *dx, double *dy, double *dz);

#ifdef __cplusplus
}
#endif
#endif /* GNSS_CARRIER_H */
