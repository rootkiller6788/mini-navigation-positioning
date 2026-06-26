/**
 * geomag_sensor.h -- Magnetometer Sensor Models & Calibration
 *
 * L2: Core concepts of magnetic field measurement
 * L5: Sensor calibration algorithms (hard-iron, soft-iron, alignment)
 * L6: Magnetic anomaly detection (MAD)
 *
 * Magnetometer types modeled:
 *   1. Scalar (total field): optically pumped (Cs, K, He-4), Overhauser, proton precession
 *      - Measures |B| only, insensitive to orientation
 *      - Precision: 0.001-1 nT, noise floor ~0.001 nT/sqrt(Hz)
 *   2. Triaxial fluxgate: measures Bx, By, Bz in sensor body frame
 *      - Precision: 0.1-10 nT, bandwidth DC-100 Hz
 *      - Susceptible to hard-iron (permanent magnet bias) and soft-iron (induced) errors
 *   3. Gradiometer: differential measurement suppressing common-mode
 *      - Two sensors at known baseline (typically 0.5-5 m)
 *      - Suppresses diurnal variation, distant magnetic disturbances
 *
 * Reference:
 *   Ripka, "Magnetic Sensors and Magnetometers" (2001)
 *   Gebre-Egziabher et al., "Calibration of Strapdown Magnetometers",
 *     IEEE Trans. Aerospace and Electronic Systems (2006)
 *   Munschy et al., "Scalar and vector magnetometers for UAV", Sensors (2020)
 */

#ifndef GEOMAG_SENSOR_H
#define GEOMAG_SENSOR_H

#include "geomag_core.h"

/* ========================================================================
 * L2: Scalar magnetometer model
 * ======================================================================== */

/**
 * L2: Model a scalar magnetometer reading.
 *
 * Reading = |B_true| + bias + noise
 *
 * Bias sources: heading error (few nT), temperature drift
 * Noise model: Gaussian white noise with given standard deviation
 *
 * @param B_true    True magnetic field vector [nT]
 * @param bias      Sensor bias [nT]
 * @param noise_std Measurement noise standard deviation [nT]
 * @param reading   Output scalar reading [nT]
 */
void scalar_magnetometer_read(const MagVector *B_true, double bias,
                               double noise_std, double *reading);

/**
 * L2: Model a triaxial fluxgate magnetometer reading.
 *
 * B_raw = S * M * (B_body + B_hard) + w
 *
 * where:
 *   B_body = true field in body frame
 *   B_hard = 3x1 hard-iron bias vector
 *   S = 3x3 diagonal scale factor matrix
 *   M = 3x3 soft-iron + misalignment matrix
 *   w = Gaussian white noise [nT]
 *
 * @param B_body    True field in body frame [nT]
 * @param bias      Hard-iron bias [nT]
 * @param scale     3x3 scale factor (diagonal)
 * @param soft_iron 3x3 soft-iron/misalignment matrix (row-major)
 * @param noise_std Noise std per axis [nT]
 * @param reading   Output: raw triaxial reading [nT]
 */
void triaxial_magnetometer_read(const MagVector *B_body,
                                 const double bias[3],
                                 const double scale[9],
                                 const double soft_iron[9],
                                 double noise_std,
                                 double reading[3]);

/**
 * L5: Two-step magnetometer calibration.
 *
 * Step 1 (hard-iron): estimate bias as center of measurement ellipsoid.
 *   Given N measurements on a sphere, the hard-iron offset is the
 *   centroid of the measurement cloud.
 *
 * Step 2 (soft-iron): estimate scale factors and non-orthogonality.
 *   Fits an ellipsoid to measurements, computes correction matrix
 *   to map ellipsoid back to sphere.
 *
 * Algorithm: Iterative least-squares ellipsoid fitting (Li et al., 2005)
 *
 * @param measurements  Array of N triaxial measurements [nT]
 * @param N             Number of measurements
 * @param bias_out      Output hard-iron bias [nT]
 * @param scale_out     Output 3x3 combined correction matrix (row-major)
 * @return 0 on success, -1 if insufficient data (N<9)
 */
int magnetometer_calibrate_ls(const double measurements[][3], int N,
                               double bias_out[3], double scale_out[9]);

/**
 * L5: Compute calibration quality metric (sphericity).
 *
 * After calibration, ideal measurements should lie on a sphere of radius |B|.
 * This function computes the RMS deviation from that sphere.
 *
 * sphericity_error = sqrt( (1/N) * sum (|B_cal_i| - |B_ref|)^2 )
 *
 * @param calibrated  Calibrated triaxial measurements [nT]
 * @param N           Number of measurements
 * @param B_ref       Reference total field magnitude [nT]
 * @return RMS sphericity error [nT]
 */
double calibration_sphericity_error(const double calibrated[][3], int N,
                                     double B_ref);

/* ========================================================================
 * L5: Magnetic gradient measurement
 * ======================================================================== */

/**
 * L5: Compute magnetic field gradient from gradiometer measurements.
 *
 * A gradiometer consists of two triaxial sensors separated by baseline vector d.
 * The gradient tensor G_ij = dB_i/dx_j is approximated by finite difference:
 *
 * G_ij ≈ [B_i(r + d/2) - B_i(r - d/2)] / d_j
 *
 * The gradient tensor is a 3x3 traceless matrix (div B = 0) in source-free space.
 * Its invariants are used for UXO/magnetic target classification.
 *
 * @param B1, B2  Magnetic field at sensor positions [nT]
 * @param baseline Separation vector from sensor 1 to sensor 2 [m]
 * @param G        3x3 gradient tensor output (row-major) [nT/m]
 */
void compute_magnetic_gradient(const MagVector *B1, const MagVector *B2,
                                const double baseline[3], double G[9]);

/**
 * L5: Compute total field gradient magnitude (analytic signal).
 *
 * |grad|B|| = sqrt( (dB/dx)^2 + (dB/dy)^2 + (dB/dz)^2 )
 *
 * We approximate via the gradient tensor acting on the unit vector in
 * the field direction: grad(|B|) = (G^T * B_hat) where B_hat = B/|B|.
 *
 * @param G    3x3 gradient tensor [nT/m]
 * @param B    Magnetic field vector [nT]
 * @return Magnitude of the gradient of total field [nT/m]
 */
double total_field_gradient_magnitude(const double G[9], const MagVector *B);

/**
 * L5: Compute the Normalized Source Strength (NSS) from gradient tensor.
 *
 * NSS = sqrt( -lambda_2^2 - lambda_1*lambda_3 )
 *
 * where lambda_1 >= lambda_2 >= lambda_3 are eigenvalues of G.
 * NSS is nearly independent of source magnetization direction,
 * making it useful for magnetic dipole localization.
 *
 * Reference: Beiki et al., "Interpretation of magnetic gradient tensor data",
 *            Geophysics (2012)
 *
 * @param G     3x3 gradient tensor [nT/m]
 * @return NSS value [nT/m]
 */
double normalized_source_strength(const double G[9]);

/* ========================================================================
 * L6: Magnetic Anomaly Detection (MAD)
 * ======================================================================== */

/**
 * L6: Simple threshold-based MAD detector.
 *
 * Compares measured total field against expected IGRF field.
 * Anomaly = |B_measured| - |B_igrf|
 *
 * Detection rule: |anomaly| > threshold --> DETECT
 *
 * Typical thresholds:
 *   submarine detection: 0.1-5 nT (at altitude/depth)
 *   UXO detection: 1-10 nT (at 1-5 m range)
 *   geological survey: 10-1000 nT
 *
 * @param B_measured  Measured total field [nT]
 * @param B_igrf      Expected IGRF total field [nT]
 * @param threshold   Detection threshold [nT]
 * @return 1 if anomaly detected, 0 otherwise
 */
int mad_threshold_detect(double B_measured, double B_igrf, double threshold);

/**
 * L5: Orthonormalized Basis Functions (OBF) MAD algorithm.
 *
 * Advanced MAD that uses a bank of orthogonal basis functions matched
 * to dipole signatures. More sensitive than simple threshold detection.
 *
 * OBF decomposes the measurement time series onto basis functions
 * representing dipole moment components. The detection statistic is
 * the sum of squared projections.
 *
 * Reference: Ginzburg et al., "Processing of magnetic data for MAD",
 *            IEEE Trans. Geoscience and Remote Sensing (2004)
 *
 * @param signal      Time series of total field measurements [nT]
 * @param N           Length of time series
 * @param time_base   Time vector [s]
 * @param v           Relative velocity sensor/target [m/s]
 * @param cpa         Closest point of approach [m]
 * @param stat        Output detection statistic
 * @return 0 on success, -1 on error
 */
int mad_obf_detect(const double *signal, int N, const double *time_base,
                    double v, double cpa, double *stat);

/**
 * L5: Compute signal-to-noise ratio for MAD.
 *
 * SNR = (anomaly_amplitude / noise_rms)^2
 *
 * SNR_dB = 10 * log10(SNR)
 *
 * Probability of detection Pd for given Pfa:
 *   Pd = Q(Q^{-1}(Pfa) - sqrt(2*SNR))  (for coherent detection)
 *
 * @param anomaly_amplitude  Peak anomaly [nT]
 * @param noise_rms          RMS noise [nT]
 * @return SNR [dB]
 */
double mad_snr_db(double anomaly_amplitude, double noise_rms);

/**
 * L6: Magnetic dipole moment estimation from anomaly measurements.
 *
 * Given multiple measurements of magnetic anomaly around a target,
 * estimate the dipole moment vector m of the source.
 *
 * B_anomaly(r) = (mu_0 / 4*pi) * [3(m·r_hat)*r_hat - m] / r^3
 *
 * Uses linear least-squares: solve A*m = B_meas for m (3 unknown components).
 *
 * @param positions   Array of N measurement positions relative to target [m]
 * @param anomalies   Array of N anomaly vectors [nT]
 * @param N           Number of measurements (>=3)
 * @param moment      Output estimated dipole moment [A*m^2]
 * @return 0 on success, -1 if ill-conditioned
 */
int estimate_dipole_moment(const double positions[][3],
                            const double anomalies[][3],
                            int N, double moment[3]);

#endif /* GEOMAG_SENSOR_H */
