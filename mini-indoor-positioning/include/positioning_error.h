/**
 * @file positioning_error.h
 * @brief Error analysis and accuracy metrics for indoor positioning
 *
 * Knowledge Coverage:
 *   L1 - Definitions: CEP, DRMS, 2DRMS, R95, RMSE, bias, precision,
 *        availability, integrity, continuity
 *   L2 - Core Concepts: error budget, dilution of precision (DOP),
 *        sensitivity analysis, error propagation
 *   L3 - Mathematical Structures: covariance propagation, Jacobian
 *        of positioning equations, Fisher information matrix
 *   L4 - Fundamental Laws: Cramer-Rao lower bound (CRLB) for
 *        positioning, error propagation law
 *   L5 - Algorithms: error ellipse computation, CEP estimation,
 *        confidence interval computation, outlier detection
 *
 * Reference: Groves (2013), Ch.14 - "Performance Assessment"
 *            Kaplan & Hegarty, "Understanding GPS/GNSS: Principles
 *            and Applications" (2017), Ch.7 - "Performance"
 *
 * Course Alignment:
 *   - MIT 6.003 (Signal Processing) — estimation error
 *   - Stanford EE363 (Linear Dynamical Systems) — covariance analysis
 *   - ETH 227-0427 (Signal Processing)
 *   - 清华 信号与系统 — 估计误差分析
 */

#ifndef POSITIONING_ERROR_H
#define POSITIONING_ERROR_H

#include "indoor_positioning.h"
#include <stdint.h>

/* ============================================================================
 * L1 - Definitions: Error Metrics
 * ============================================================================ */

/**
 * @brief Dilution of Precision (DOP) metrics
 *
 * L1 Definition: DOP quantifies the geometric contribution of
 * anchor/satellite geometry to positioning error. Lower = better.
 *
 * GDOP = Geometric DOP (3D position + time)
 * PDOP = Position DOP (3D position)
 * HDOP = Horizontal DOP (2D position)
 * VDOP = Vertical DOP (height)
 * TDOP = Time DOP (clock bias)
 *
 * Reference: Kaplan & Hegarty (2017), Section 7.3
 */
typedef struct {
    double gdop;  /**< Geometric dilution of precision */
    double pdop;  /**< Position dilution of precision (3D) */
    double hdop;  /**< Horizontal dilution of precision (2D) */
    double vdop;  /**< Vertical dilution of precision */
    double tdop;  /**< Time dilution of precision */
} dop_metrics_t;

/**
 * @brief Comprehensive positioning error statistics
 *
 * L1 Definition: Full error characterization for a positioning
 * system, including accuracy (mean error), precision (std dev),
 * and reliability metrics.
 */
typedef struct {
    double mean_error_x;       /**< Mean error in X (east) direction */
    double mean_error_y;       /**< Mean error in Y (north) direction */
    double mean_error_z;       /**< Mean error in Z (up) direction */
    double std_dev_x;          /**< Standard deviation in X */
    double std_dev_y;          /**< Standard deviation in Y */
    double std_dev_z;          /**< Standard deviation in Z */
    double rms_2d;             /**< Root-mean-square 2D error */
    double rms_3d;             /**< Root-mean-square 3D error */
    double cep50;              /**< Circular error probable 50% */
    double cep95;              /**< Circular error probable 95% */
    double sep50;              /**< Spherical error probable 50% */
    double sep95;              /**< Spherical error probable 95% */
    double drms;               /**< Distance RMS */
    double r95;                /**< 95% confidence radius */
    double cep99;              /**< Circular error probable 99% */
    double max_error;          /**< Maximum observed error */
    double availability;       /**< Percentage of time position available */
    int    n_samples;          /**< Number of samples in computation */
    int    n_outliers;         /**< Number of outliers detected */
} positioning_error_stats_t;

/* ============================================================================
 * L5 - Error Metric Computation
 * ============================================================================ */

/**
 * @brief Compute DOP from geometry matrix G
 *
 * Given N anchors in 3D, form the geometry matrix G (N x 4):
 *   G_i = [cos(az_i)*cos(el_i), sin(az_i)*cos(el_i), sin(el_i), 1]
 *
 * Then Q = (G^T * G)^{-1} and DOPs = sqrt(diag(Q))
 *
 * @param anchor_positions Array of anchor positions (at least 4)
 * @param user_position User/receiver position
 * @param n_anchors Number of anchors (>=4)
 * @param[out] dop Computed DOP metrics
 * @return 0 on success, -1 if G^T*G is singular
 *
 * L5: DOP computation from geometry matrix.
 * Complexity: O(N) for N anchors
 */
int compute_dop(const position3d_t *anchor_positions,
                const position3d_t *user_position,
                int n_anchors,
                dop_metrics_t *dop);

/**
 * @brief Compute positioning error statistics from truth vs estimate data
 *
 * @param truth Array of ground-truth positions
 * @param estimate Array of estimated positions
 * @param n Number of samples
 * @param[out] stats Computed error statistics
 *
 * L5: Aggregate positioning error statistics.
 * Complexity: O(N)
 */
void compute_error_stats(const position3d_t *truth,
                         const position3d_t *estimate,
                         int n,
                         positioning_error_stats_t *stats);

/**
 * @brief Compute 2D error ellipse from covariance matrix
 *
 * Given the 2x2 covariance submatrix for (x,y), compute the
 * error ellipse parameters at a given confidence level.
 *
 * The error ellipse semi-axes are:
 *   a = sqrt(chi2_inv * lambda_max)
 *   b = sqrt(chi2_inv * lambda_min)
 * where lambda are eigenvalues of the 2x2 covariance.
 *
 * @param cov_xx Variance in x
 * @param cov_yy Variance in y
 * @param cov_xy Covariance between x,y
 * @param confidence Confidence level (0 to 1, e.g., 0.95)
 * @param[out] ellipse Error ellipse parameters
 *
 * L5: Error ellipse computation from Gaussian covariance.
 * Reference: Johnson & Wichern, "Applied Multivariate Statistical
 * Analysis" (2007), Ch.4.
 *
 * Complexity: O(1)
 */
void compute_error_ellipse(double cov_xx, double cov_yy, double cov_xy,
                           double confidence, error_ellipse_t *ellipse);

/**
 * @brief Chi-squared inverse (percent point function) approximation
 *
 * Computes the approximate value x such that P(chi2_k <= x) = p
 * for k degrees of freedom.
 *
 * Uses Wilson-Hilferty transformation for k >= 1.
 *
 * @param p Probability (0 to 1)
 * @param k Degrees of freedom
 * @return Approximate chi2 inverse value
 *
 * Complexity: O(1)
 */
double chi2_inv(double p, int k);

/* ============================================================================
 * L3/L4 - Cramer-Rao Lower Bound for Positioning
 * ============================================================================ */

/**
 * @brief Compute the Cramer-Rao Lower Bound (CRLB) for RSSI-based
 * 2D positioning with N anchors.
 *
 * L4 Law: The CRLB provides the theoretical minimum variance
 * achievable by any unbiased estimator of position.
 *
 * For RSSI positioning, the Fisher Information Matrix (FIM) is:
 *   FIM = sum_i (10*n/(ln(10)*sigma_i*d_i))^2 * G_i * G_i^T
 * where G_i is the unit vector from anchor i to user.
 *
 * CRLB = trace(FIM^{-1})
 *
 * @param anchor_positions Anchor positions
 * @param distances True distances to anchors
 * @param rssi_std RSSI measurement std dev per anchor
 * @param path_loss_exp Path loss exponent
 * @param n_anchors Number of anchors
 * @return CRLB in meters^2 (sum of position variances)
 *
 * L4: CRLB for RSSI positioning.
 * Reference: Patwari et al., "Relative location estimation in wireless
 * sensor networks," IEEE Trans. Signal Processing, 2003.
 *
 * Complexity: O(N)
 */
double crlb_rssi_positioning(const position2d_t *anchor_positions,
                             const double *distances,
                             const double *rssi_std,
                             double path_loss_exp,
                             int n_anchors);

/**
 * @brief Compute CRLB for TOF-based positioning
 *
 * For TOF, the FIM is:
 *   FIM = sum_i (1/sigma_i^2) * G_i * G_i^T
 *
 * @param anchor_positions Anchor positions
 * @param noise_std TOF measurement noise std dev per anchor (meters)
 * @param n_anchors Number of anchors
 * @return CRLB in meters^2
 *
 * L4: CRLB for TOF positioning.
 * Complexity: O(N)
 */
double crlb_tof_positioning(const position2d_t *anchor_positions,
                            const double *noise_std,
                            int n_anchors);

/**
 * @brief Propagate measurement uncertainty through trilateration
 *
 * Given anchor geometry and distance measurement uncertainties,
 * compute the resulting position uncertainty covariance using
 * first-order error propagation.
 *
 * P_xy = (H^T * W * H)^{-1}
 * where H is the Jacobian of distance equations and W = diag(1/sigma_i^2)
 *
 * @param anchor_positions Anchor positions
 * @param distance_std Measurement standard deviations
 * @param user_pos User position (linearization point)
 * @param n_anchors Number of anchors
 * @param[out] cov_xx Output: X variance
 * @param[out] cov_yy Output: Y variance
 * @param[out] cov_xy Output: XY covariance
 *
 * L4: Error propagation through positioning equations.
 * Complexity: O(N)
 */
void propagate_positioning_error(const position2d_t *anchor_positions,
                                 const double *distance_std,
                                 const position2d_t *user_pos,
                                 int n_anchors,
                                 double *cov_xx, double *cov_yy,
                                 double *cov_xy);

/* ============================================================================
 * L5 - Outlier Detection and Integrity
 * ============================================================================ */

/**
 * @brief Detect outlier measurements using chi-squared test on
 * innovation sequence of a Kalman filter.
 *
 * If normalized innovation squared (NIS) exceeds chi2 threshold,
 * the measurement is flagged as an outlier.
 *
 * @param innovation Measurement innovation vector
 * @param innovation_cov Innovation covariance matrix (n_measure x n_measure)
 * @param n_measure Measurement dimension
 * @param confidence Confidence level (e.g., 0.99)
 * @return 1 if outlier detected, 0 if nominal
 *
 * L5: Innovation-based outlier detection.
 * Complexity: O(m^3) for m = n_measure
 */
int detect_measurement_outlier(const double *innovation,
                               const double *innovation_cov,
                               int n_measure,
                               double confidence);

/**
 * @brief RANSAC-based robust position estimation
 *
 * Random Sample Consensus: randomly selects subsets of measurements,
 * estimates position, and selects the estimate consistent with
 * the most measurements.
 *
 * @param distances Distance measurements
 * @param anchor_positions Anchor positions
 * @param n_anchors Number of anchors
 * @param n_subsets Minimum subset size (3 for 2D)
 * @param n_iterations Number of RANSAC iterations
 * @param inlier_threshold Distance threshold for inliers (meters)
 * @param[out] result Robust position estimate
 * @return Number of inliers for the best model
 *
 * L5: RANSAC for robust positioning against outliers.
 * Complexity: O(I * N) for I iterations, N anchors
 */
int ransac_positioning(const double *distances,
                       const position2d_t *anchor_positions,
                       int n_anchors,
                       int n_subsets,
                       int n_iterations,
                       double inlier_threshold,
                       position2d_t *result);

/* ============================================================================
 * L2 - Error Budget Analysis
 * ============================================================================ */

/**
 * @brief Analyze error budget by decomposing total error into sources
 *
 * @param truth Ground-truth positions
 * @param estimate Estimated positions
 * @param n Number of samples
 * @param[out] bias_component Computed bias (systematic error)
 * @param[out] noise_component Noise (random error) std dev
 * @param[out] drift_component Time-correlated error drift rate
 * @param[out] total_rms Total RMS error
 *
 * Total error = bias + noise + drift + multipath...
 *
 * L2: Error budget decomposition.
 * Complexity: O(N)
 */
void decompose_error_sources(const position3d_t *truth,
                             const position3d_t *estimate,
                             int n,
                             double *bias_component,
                             double *noise_component,
                             double *drift_component,
                             double *total_rms);

/**
 * @brief Compute the Allan variance for IMU noise characterization
 *
 * Allan variance separates noise sources (angle random walk,
 * bias instability, rate random walk) by analyzing the IMU
 * signal at different averaging times.
 *
 * @param signal Time-series signal (gyro or accel readings for one axis)
 * @param n_samples Number of samples
 * @param sample_rate_hz Sample rate
 * @param[out] tau Array of averaging times (output, pre-allocated, size m)
 * @param[out] allan_var Array of Allan variance values
 * @param m Number of clusters to compute
 *
 * L5: Allan variance for IMU noise characterization.
 * Reference: IEEE Std 952-1997, "Standard Specification Format Guide
 * and Test Procedure for Single-Axis Interferometric Fiber Optic Gyros"
 */
void compute_allan_variance(const double *signal, int n_samples,
                            double sample_rate_hz,
                            double *tau, double *allan_var, int m);

/* ============================================================================
 * L1 - Reliability Metrics
 * ============================================================================ */

/**
 * @brief Compute system availability (percentage of time valid
 * positioning solution is available)
 *
 * @param n_valid Number of epochs with valid solution
 * @param n_total Total number of epochs
 * @return Availability ratio [0, 1]
 */
static inline double compute_availability(int n_valid, int n_total) {
    if (n_total <= 0) return 0.0;
    return (double)n_valid / (double)n_total;
}

/**
 * @brief Compute integrity risk (probability of undetected error
 * exceeding alert limit)
 *
 * @param error_samples Array of positioning errors
 * @param n_samples Number of samples
 * @param alert_limit Maximum allowable error before alert
 * @return Integrity risk [0, 1]
 */
double compute_integrity_risk(const double *error_samples, int n_samples,
                              double alert_limit);

/**
 * @brief Compute continuity risk (probability of service interruption
 * after initiation)
 *
 * @param epoch_durations Duration of each continuous positioning segment
 * @param n_segments Number of segments
 * @param interruption_threshold Minimum interruption considered failure
 * @return Continuity risk [0, 1]
 */
double compute_continuity_risk(const double *epoch_durations,
                               int n_segments,
                               double interruption_threshold);

#endif /* POSITIONING_ERROR_H */
