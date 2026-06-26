#ifndef INS_CALIBRATION_H
#define INS_CALIBRATION_H
#include "ins_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L5: IMU Calibration Methods
 *
 * IMU calibration determines the deterministic error parameters:
 * bias, scale factor, and misalignment (non-orthogonality) for
 * both accelerometer and gyroscope triads.
 *
 * Reference: Groves (2013), Chapter 4, "Inertial Sensor Calibration".
 * Reference: Shin & El-Sheimy (2002), "A New Calibration Method for
 *   Strapdown Inertial Navigation Systems", Z. Vermess., 127(1): 41-50.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * L1: Calibration Parameter Types
 * ------------------------------------------------------------------------- */

/** Calibration result for a single 3-axis sensor triad. */
typedef struct {
    ins_vec3_t bias;
    ins_vec3_t scale_factor;
    ins_mat3_t misalignment;
    double     residual_rms;
    int        num_positions;
} ins_calibration_result_t;

/** Full IMU calibration (accel + gyro triads). */
typedef struct {
    ins_calibration_result_t accel;
    ins_calibration_result_t gyro;
    ins_vec3_t               accel_gyro_misalign;
} ins_imu_calibration_t;

/* -------------------------------------------------------------------------
 * L5: Six-Position Static Calibration
 * ------------------------------------------------------------------------- */

/**
 * Six-position static accelerometer calibration.
 *
 * Uses Earth's gravity (g) as reference.
 * Six positions: each axis (+x,-x,+y,-y,+z,-z) pointing down.
 *
 * For each axis i:
 *   b_i = (reading_up + reading_down) / 2
 *   sf_i = (reading_down - reading_up) / (2*g) - 1
 */
int ins_calib_six_position_accel(const ins_vec3_t readings[6],
                                  double g_local,
                                  ins_calibration_result_t *result);

/* -------------------------------------------------------------------------
 * L5: Multi-Position Gyroscope Calibration
 * ------------------------------------------------------------------------- */

/**
 * Multi-position gyroscope calibration using known rotation rates.
 *
 * For navigation-grade gyros, Earth rate (15 deg/hr) can be used:
 *   gyro_north = we * cos(lat) * cos(psi)
 *   gyro_east  = we * cos(lat) * sin(psi)
 *   gyro_down  = we * sin(lat)
 */
int ins_calib_multi_position_gyro(const ins_vec3_t *rates,
                                   const ins_vec3_t *readings,
                                   size_t n_positions,
                                   ins_calibration_result_t *result);

/* -------------------------------------------------------------------------
 * L5: Cross-Axis (Misalignment) Calibration
 * ------------------------------------------------------------------------- */

/**
 * Compute misalignment correction matrix from known reference vectors.
 *
 * Correction matrix (small angles):
 *   M = [ 1       -alpha_xz  alpha_xy ]
 *       [ alpha_yz  1       -alpha_yx ]
 *       [ -alpha_zy  alpha_zx  1      ]
 */
int ins_calib_cross_axis(const ins_vec3_t *reference,
                          const ins_vec3_t *measured,
                          size_t n,
                          ins_mat3_t *M_corr);

/* -------------------------------------------------------------------------
 * L5: Thermal Calibration Model
 * ------------------------------------------------------------------------- */

/** Temperature compensation coefficients (3rd-order polynomial). */
typedef struct {
    double c0, c1, c2, c3;
    double T0;
    double T_min, T_max;
} ins_thermal_coeffs_t;

/**
 * Fit temperature compensation polynomial from calibration data.
 * Uses ordinary least squares.
 */
int ins_calib_thermal_fit(const double *temperatures, const double *biases,
                           size_t n, ins_thermal_coeffs_t *coeffs);

/**
 * Apply temperature compensation to a sensor reading.
 */
double ins_calib_thermal_apply(double reading, double T,
                                const ins_thermal_coeffs_t *coeffs);

/* -------------------------------------------------------------------------
 * L5: Full IMU Calibration Pipeline
 * ------------------------------------------------------------------------- */

/**
 * Apply full calibration correction to an IMU sample.
 * Processing: temp comp -> bias sub -> scale factor -> misalignment.
 */
void ins_calib_apply(const ins_imu_sample_t *raw,
                      const ins_imu_calibration_t *calib,
                      double temp,
                      ins_vec3_t *accel_corr,
                      ins_vec3_t *gyro_corr);

#ifdef __cplusplus
}
#endif
#endif /* INS_CALIBRATION_H */
