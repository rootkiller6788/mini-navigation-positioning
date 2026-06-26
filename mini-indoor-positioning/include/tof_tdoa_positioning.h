/**
 * @file tof_tdoa_positioning.h
 * @brief Time-based indoor positioning: ToF, TDoA, UWB, AOA
 *
 * Knowledge Coverage:
 *   L1 - Definitions: Time of Flight (ToF), Time Difference of Arrival
 *        (TDoA), Angle of Arrival (AoA), Two-Way Ranging (TWR),
 *        Symmetric Double-Sided TWR (SDS-TWR)
 *   L2 - Core Concepts: hyperbolic positioning, circular positioning,
 *        clock synchronization, timestamp exchange protocol
 *   L3 - Mathematical Structures: hyperbolic system of equations,
 *        direction finding, interferometry basics
 *   L4 - Fundamental Laws: constant speed of light (signal propagation),
 *        clock drift models, Cramer-Rao bound for TOA/TDOA
 *   L5 - Algorithms: linear least squares for TOA, Chan's algorithm
 *        for TDOA, Taylor series iterative TDOA, MUSIC for AoA
 *   L6 - Canonical Problems: UWB indoor positioning (Decawave DW1000),
 *        acoustic chirp ToF, 5G NR positioning
 *
 * Reference: Dardari et al., "Ranging With Ultrawide Bandwidth Signals
 *            in Multipath Environments," Proc. IEEE, 2009.
 *            IEEE Std 802.15.4z-2020 — UWB enhanced ranging
 *            Chan & Ho (1994) — TDOA algorithm
 *
 * Course Alignment:
 *   - MIT 6.630 (EM Waves) — wave propagation
 *   - Stanford EE359 (Wireless) — UWB tech
 *   - ETH 227-0455 (EM)
 *   - 清华 通信原理 — 扩频/UWB
 */

#ifndef TOF_TDOA_POSITIONING_H
#define TOF_TDOA_POSITIONING_H

#include "indoor_positioning.h"

/* ============================================================================
 * L1 - Definitions: Time-of-Flight and UWB
 * ============================================================================ */

/** Speed of light in vacuum (m/s) */
#define SPEED_OF_LIGHT_MPS 299792458.0

/** Maximum UWB ranging distance (meters) */
#define UWB_MAX_RANGE 200.0

/** UWB channel definitions (IEEE 802.15.4) */
#define UWB_CHANNEL_1 3494.4   /**< Center frequency MHz, BW 499.2 MHz */
#define UWB_CHANNEL_2 3993.6
#define UWB_CHANNEL_3 4492.8
#define UWB_CHANNEL_5 6489.6

/**
 * @brief Time of Flight measurement
 *
 * L1 Definition: ToF is the signal propagation time from
 * transmitter to receiver. Distance = ToF * c.
 *
 * Two-way ranging (TWR) avoids clock synchronization by
 * round-trip timing.
 */
typedef struct {
    uint64_t t_sp;         /**< Transmit timestamp of poll message (ticks) */
    uint64_t t_rp;         /**< Receive timestamp of poll message at responder */
    uint64_t t_sr;         /**< Transmit timestamp of response message */
    uint64_t t_rr;         /**< Receive timestamp of response message at initiator */
    uint64_t t_sf;         /**< Transmit timestamp of final message (SDS-TWR) */
    uint64_t t_rf;         /**< Receive timestamp of final message */
    double   tick_period_s; /**< Clock tick period in seconds */
    int      is_ss_twr;    /**< 1 if Symmetric Double-Sided TWR */
} twr_exchange_t;

/**
 * @brief Compute distance from two-way ranging exchange
 *
 * Standard TWR: d = (t_rr - t_sp - (t_sr - t_rp)) * c / 2
 * SDS-TWR: d = (T_round1 * T_round2 - T_reply1 * T_reply2) /
 *               (T_round1 + T_round2 + T_reply1 + T_reply2) * c
 *
 * @param twr TWR timestamp exchange data
 * @return Estimated distance in meters, or -1 on error
 *
 * L5: Two-Way Ranging distance estimation.
 * Complexity: O(1)
 */
double twr_compute_distance(const twr_exchange_t *twr);

/**
 * @brief Correct TWR distance for clock drift using SDS-TWR
 *
 * SDS-TWR eliminates first-order clock drift errors.
 *
 * @param twr SDS-TWR exchange data
 * @return Clock-drift-corrected distance in meters
 *
 * L5: Symmetric Double-Sided TWR.
 * Reference: IEEE 802.15.4z-2020
 */
double twr_sds_compute_distance(const twr_exchange_t *twr);

/* ============================================================================
 * L2 - UWB Positioning System Models
 * ============================================================================ */

/**
 * @brief UWB anchor (fixed reference node)
 *
 * L1 Definition: A UWB anchor is a fixed-position device that
 * participates in ranging exchanges with mobile tags.
 */
typedef struct {
    position3d_t position;       /**< Surveyed anchor position */
    uint16_t     anchor_id;      /**< Unique anchor identifier */
    double       antenna_delay;  /**< Antenna group delay in seconds */
    int          active;         /**< 1 if anchor is operational */
} uwb_anchor_t;

/**
 * @brief UWB ranging measurement
 */
typedef struct {
    uint16_t anchor_id;       /**< Anchor ID */
    double   distance;         /**< Measured distance in meters */
    double   distance_std;     /**< Distance uncertainty std dev */
    double   rssi;             /**< RSSI of received UWB signal */
    uint64_t timestamp_us;     /**< Measurement timestamp */
    uint8_t  los_flag;         /**< 1 = line-of-sight, 0 = NLOS */
    double   los_confidence;   /**< LOS probability [0, 1] */
} uwb_ranging_t;

/* ============================================================================
 * L5 - ToF/TDoA Positioning Algorithms
 * ============================================================================ */

/**
 * @brief Solve 2D position from TOA measurements using iterative
 * Gauss-Newton non-linear least squares.
 *
 * Minimizes: sum_i (||p_user - p_anchor_i|| - d_i)^2
 *
 * @param anchor_positions Known anchor positions
 * @param toa_distances Distance measurements (TOA * c)
 * @param n_anchors Number of anchors (>= 3 for 2D)
 * @param initial_guess Starting position estimate
 * @param[out] result Final position estimate
 * @param max_iter Maximum Gauss-Newton iterations
 * @param tol Convergence tolerance (meters)
 * @return 0 on success, -1 on no convergence
 *
 * L5: Non-linear least squares (Gauss-Newton) for TOA positioning.
 * Complexity: O(K*N^2) for K iterations, N anchors
 */
int toa_positioning_2d(const position2d_t *anchor_positions,
                       const double *toa_distances,
                       int n_anchors,
                       const position2d_t *initial_guess,
                       position2d_t *result,
                       int max_iter, double tol);

/**
 * @brief Taylor series iterative TDOA positioning
 *
 * Linearizes the hyperbolic TDOA equations around an initial
 * guess and iteratively solves for the position correction.
 *
 * @param anchors Anchor positions (index 0 is reference)
 * @param tdoa_s TDOA measurements in seconds relative to anchor[0]
 * @param n_anchors Total anchors (>= 4)
 * @param speed Signal propagation speed in m/s
 * @param initial_guess Starting position
 * @param[out] result Final estimate
 * @param max_iter Maximum iterations
 * @param tol Convergence tolerance
 * @return 0 on success
 *
 * L5: Taylor series TDOA positioning.
 * Reference: Foy, "Position-location solutions by Taylor-series
 * estimation," IEEE Trans. AES, 1976.
 *
 * Complexity: O(K*N^2) for K iterations
 */
int tdoa_taylor_series(const position3d_t *anchors,
                       const double *tdoa_s,
                       int n_anchors,
                       double speed,
                       const position3d_t *initial_guess,
                       position3d_t *result,
                       int max_iter, double tol);

/* ============================================================================
 * L5 - Non-Line-of-Sight (NLOS) Detection & Mitigation
 * ============================================================================ */

/**
 * @brief NLOS detection using RSSI vs distance consistency check
 *
 * If measured RSSI is significantly lower than expected for the
 * measured distance, the signal likely traversed an NLOS path.
 *
 * @param measured_distance Measured ToF distance in meters
 * @param measured_rssi Measured RSSI in dBm
 * @param model Path loss model parameters
 * @param n_sigma_thresh Number of sigma threshold for NLOS flag
 * @return 1 if NLOS detected, 0 if LOS
 *
 * Complexity: O(1)
 */
int detect_nlos_rssi_distance(double measured_distance, double measured_rssi,
                              const path_loss_model_t *model,
                              double n_sigma_thresh);

/**
 * @brief NLOS mitigation by de-weighting suspected NLOS measurements
 *
 * Applies a weight to each ranging measurement based on LOS confidence.
 * w_i = los_confidence_i (0-1)
 *
 * Used in weighted least-squares positioning.
 *
 * @param range_measurements Array of ranging measurements
 * @param n_measurements Number of measurements
 * @param[out] weights Computed weights for each measurement
 *
 * Complexity: O(N)
 */
void compute_nlos_weights(const uwb_ranging_t *range_measurements,
                          int n_measurements,
                          double *weights);

/**
 * @brief NLOS identification via range residual test
 *
 * After initial positioning, compute residuals r_i = |d_meas_i - d_pred_i|.
 * Measurements with large residuals are likely NLOS.
 *
 * @param anchor_positions Anchor positions
 * @param distances Measured distances
 * @param n_anchors Number of anchors
 * @param estimated_pos Current position estimate
 * @param[out] residuals Residual for each measurement
 * @param threshold Residual threshold for NLOS flag (meters)
 * @param[out] nlos_flags 1 for each NLOS measurement
 * @return Number of measurements flagged as NLOS
 *
 * L5: Residual-based NLOS detection.
 * Complexity: O(N)
 */
int detect_nlos_residual(const position2d_t *anchor_positions,
                         const double *distances,
                         int n_anchors,
                         const position2d_t *estimated_pos,
                         double *residuals,
                         double threshold,
                         int *nlos_flags);

/* ============================================================================
 * L5 - Angle of Arrival (AoA) Estimation
 * ============================================================================ */

/**
 * @brief AoA measurement from a 2-element antenna array
 *
 * L1 Definition: Angle of Arrival is the direction from which
 * a signal arrives at the receiver. For a 2-element array with
 * spacing d, the phase difference is:
 *   delta_phi = 2*pi*d*sin(theta)/lambda
 *   theta = arcsin(delta_phi * lambda / (2*pi*d))
 */
typedef struct {
    double antenna_spacing_m;     /**< Antenna element spacing in meters */
    double wavelength_m;          /**< Signal wavelength in meters */
    double phase_difference_rad;  /**< Measured phase difference */
    double snr;                   /**< Signal-to-noise ratio */
    double angle_std;             /**< Angle uncertainty standard deviation */
} aoa_measurement_t;

/**
 * @brief Estimate AoA from phase difference measurement
 *
 * @param meas AoA measurement data
 * @return Estimated angle in radians (from array broadside), or NaN on error
 *
 * Complexity: O(1)
 */
double aoa_from_phase(const aoa_measurement_t *meas);

/**
 * @brief Compute position from two AoA measurements (triangulation)
 *
 * Given AoA from two known anchor positions, compute intersection
 * of the two lines of bearing.
 *
 * @param anchor1 First anchor position
 * @param anchor2 Second anchor position
 * @param aoa1 AoA at anchor 1 in radians (from east, CCW)
 * @param aoa2 AoA at anchor 2 in radians (from east, CCW)
 * @param[out] result Intersection position
 * @return 0 on success, -1 if lines are parallel
 *
 * L5: Triangulation from two angle measurements.
 * Complexity: O(1)
 */
int aoa_triangulate(const position2d_t *anchor1, const position2d_t *anchor2,
                    double aoa1, double aoa2,
                    position2d_t *result);

/**
 * @brief Position estimation from multiple AoA measurements using
 * maximum likelihood (Stansfield estimator)
 *
 * Minimizes: sum_i (theta_i_measured - theta_i_model)^2
 *
 * @param anchor_positions Anchor positions
 * @param aoa_measurements AoA in radians from each anchor
 * @param n_anchors Number of anchors
 * @param[out] result Estimated position
 * @return 0 on success
 *
 * L5: Stansfield (pseudolinear) AoA positioning estimator.
 * Reference: Pages-Zamora et al., "Closed-form solution for positioning
 * based on angle of arrival measurements," IEEE PIMRC, 2002.
 *
 * Complexity: O(N)
 */
int aoa_positioning_stansfield(const position2d_t *anchor_positions,
                               const double *aoa_measurements,
                               int n_anchors,
                               position2d_t *result);

/* ============================================================================
 * L6 - UWB Ranging Real-World Considerations
 * ============================================================================ */

/**
 * @brief Compute UWB link budget
 *
 * Checks if the received signal power exceeds the receiver sensitivity
 * for reliable ranging.
 *
 * @param tx_power_dbm Transmit power in dBm
 * @param tx_antenna_gain_dbi Transmit antenna gain
 * @param rx_antenna_gain_dbi Receive antenna gain
 * @param distance Distance in meters
 * @param frequency_mhz Center frequency in MHz
 * @param rx_sensitivity_dbm Receiver sensitivity in dBm
 * @return Link margin in dB (positive = link works)
 *
 * L6: UWB link budget analysis.
 * Complexity: O(1)
 */
double uwb_link_budget(double tx_power_dbm, double tx_antenna_gain_dbi,
                       double rx_antenna_gain_dbi, double distance,
                       double frequency_mhz, double rx_sensitivity_dbm);

/**
 * @brief Estimate UWB ranging precision (Cramer-Rao lower bound)
 *
 * CRLB for TOA estimation with UWB signal of bandwidth B:
 *   sigma_d >= c / (2*sqrt(2)*pi*B*SNR)
 *
 * @param bandwidth_hz Signal bandwidth in Hz
 * @param snr_linear Signal-to-noise ratio (linear, not dB)
 * @return Minimum achievable range error std dev in meters
 *
 * L4: CRLB for UWB TOA ranging.
 * Complexity: O(1)
 */
double uwb_ranging_crlb(double bandwidth_hz, double snr_linear);

/**
 * @brief First-path detection for multipath mitigation in UWB
 *
 * UWB's fine time resolution enables separating the direct path
 * from multipath reflections. This function estimates the ToF of
 * the first (direct) arrival path given a channel impulse response.
 *
 * @param cir_magnitude Channel impulse response magnitude samples
 * @param cir_sample_period_s Time between CIR samples in seconds
 * @param n_samples CIR length
 * @param noise_threshold Detection threshold relative to noise floor
 * @return First path ToF in seconds, or -1 if undetected
 *
 * L6: First path detection in UWB CIR.
 * Complexity: O(N)
 */
double uwb_first_path_tof(const double *cir_magnitude,
                          double cir_sample_period_s,
                          int n_samples,
                          double noise_threshold);

#endif /* TOF_TDOA_POSITIONING_H */
