#ifndef INS_ERRORS_H
#define INS_ERRORS_H
#include "ins_core.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1/L4: Inertial Sensor Error Models
 *
 * IMU errors are the dominant limitation in INS performance.
 * Understanding and modeling these errors is essential for navigation
 * accuracy prediction and Kalman filter design.
 *
 * Reference: IEEE Std 952-1997, "Standard Specification Format Guide
 *   and Test Procedure for Single-Axis Interferometric Fiber Optic Gyros".
 * Reference: El-Sheimy et al. (2008), "Analysis and Modeling of Inertial
 *   Sensors Using Allan Variance", IEEE Trans. Instrum. Meas., 57(1).
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * L1: Error Model Types
 * ------------------------------------------------------------------------- */

/** IMU sensor error parameters for a single axis. */
typedef struct {
    double bias_offset;
    double bias_instability;
    double bias_random_walk;
    double scale_factor_error;
    double white_noise_std;
    double misalignment[2];
    double g_sensitivity;
} ins_axis_error_t;

/** Full 3-axis IMU error model. */
typedef struct {
    ins_axis_error_t accel_x;
    ins_axis_error_t accel_y;
    ins_axis_error_t accel_z;
    ins_axis_error_t gyro_x;
    ins_axis_error_t gyro_y;
    ins_axis_error_t gyro_z;
} ins_imu_error_model_t;

/* -------------------------------------------------------------------------
 * L1: Allan Variance Data Structures
 * ------------------------------------------------------------------------- */

/** Noise process types identified by Allan variance slopes. */
typedef enum {
    INS_NOISE_QUANTIZATION  = 0,
    INS_NOISE_WHITE         = 1,
    INS_NOISE_BIAS_INSTAB   = 2,
    INS_NOISE_RATE_RANDOM   = 3,
    INS_NOISE_RATE_RAMP     = 4
} ins_noise_type_t;

/** Allan Variance decomposition result. */
typedef struct {
    double quantization_noise;
    double angle_random_walk;
    double bias_instability;
    double rate_random_walk;
    double rate_ramp;
    double bias_instability_tau;
} ins_allan_result_t;

/* -------------------------------------------------------------------------
 * L1: IMU Error Budget Presets
 * ------------------------------------------------------------------------- */

/** Pre-defined typical error budgets for common IMU grades. */
typedef struct {
    ins_grade_t  grade;
    const char  *description;
    double       gyro_bias;
    double       gyro_arw;
    double       gyro_sf_error;
    double       accel_bias;
    double       accel_vrw;
    double       accel_sf_error;
    double       position_drift;
    const char  *typical_application;
} ins_grade_spec_t;

/* -------------------------------------------------------------------------
 * L5: Allan Variance Computation
 * ------------------------------------------------------------------------- */

/**
 * Compute the Allan variance for a sequence of sensor measurements.
 *
 * sigma^2(tau) = 1/(2*(N-2*m+1)) * sum_{k=1}^{N-2*m+1}
 *   (theta_{k+2m} - 2*theta_{k+m} + theta_k)^2 / tau^2
 *
 * where theta_k = sum_{i=1}^{k} data[i] * dt
 *
 * Reference: IEEE Std 952-1997, Annex C.
 * Complexity: O(N * M)
 */
size_t ins_allan_variance(const double *data, size_t n, double dt,
                          double *taus, double *adevs, size_t max_taus);

/**
 * Decompose Allan variance into noise component parameters.
 *
 * Fits: sigma^2(tau) = 3*Q^2/tau^2 + N^2/tau + (0.664*B)^2
 *                      + K^2*tau/3 + R^2*tau^2/2
 *
 * Uses least-squares fitting in log-log space.
 */
int ins_allan_decompose(const double *taus, const double *adevs, size_t num,
                        ins_allan_result_t *result);

/* -------------------------------------------------------------------------
 * L4: Error Propagation — Position Drift Prediction
 * ------------------------------------------------------------------------- */

/**
 * Predict free-inertial horizontal position drift from IMU error specs.
 *
 * Key contributions:
 *   - Accel bias: delta_x ~ (1/2) * accel_bias * g * t^2  (quadratic)
 *   - Gyro bias:  delta_x ~ (1/6) * gyro_bias * g * t^3   (cubic, dominant after ~10 min)
 *   - ARW/VRW:    delta_x ~ noise * sqrt(t)               (random walk)
 *
 * Reference: Groves (2013) Section 5.7, "INS Error Propagation".
 */
void ins_error_predict_drift(const ins_imu_error_model_t *model,
                              double time,
                              double *pos_drift,
                              double *vel_drift,
                              double *att_drift);

/* -------------------------------------------------------------------------
 * L5: Error Propagation — Psi-Angle Error Model (Benson, 1975)
 * ------------------------------------------------------------------------- */

/**
 * Propagate INS error state using the Psi-angle error model.
 *
 * State: x = [delta_pos^T, delta_vel^T, psi^T] (9 states)
 *
 * Continuous-time dynamics:
 *   d(delta_pos)/dt = -w_en^n x delta_pos + delta_vel
 *   d(delta_vel)/dt = -C_b^n x [f^b x psi] - (2*w_ie^n+w_en^n) x delta_vel
 *   d(psi)/dt       = -w_in^n x psi + C_b^n * delta_omega_ib^b
 *
 * Also computes 9x9 state transition matrix Phi = I + F * dt.
 */
void ins_psi_error_propagate(ins_vec3_t *delta_pos,
                              ins_vec3_t *delta_vel,
                              ins_vec3_t *psi,
                              const ins_mat3_t *C_body_ned,
                              const ins_vec3_t *f_body,
                              double lat, double alt, double dt,
                              double *Phi);

/** Get typical error budget for a given IMU grade. */
const ins_grade_spec_t *ins_grade_spec_get(ins_grade_t grade);

#ifdef __cplusplus
}
#endif
#endif /* INS_ERRORS_H */
