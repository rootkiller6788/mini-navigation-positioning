#ifndef INS_INTEGRATION_H
#define INS_INTEGRATION_H
#include "ins_core.h"
#include "ins_attitude.h"
#include "ins_mechanization.h"
#include "ins_errors.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L6/L7: GPS/INS Integration Algorithms
 *
 * Inertial navigation provides high-rate, continuous navigation but
 * errors grow unbounded. GPS provides bounded-error navigation at lower
 * rates. Integration combines both for high-rate, accurate navigation.
 *
 * Reference: Groves (2013), Chapters 14-16.
 * Reference: Farrell (2008), "Aided Navigation", Chapters 8-11.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * L1: Integration Architecture Types
 * ------------------------------------------------------------------------- */

typedef enum {
    INS_INTEG_LOOSE     = 0,
    INS_INTEG_TIGHT     = 1,
    INS_INTEG_DEEP      = 2
} ins_integration_mode_t;

/* -------------------------------------------------------------------------
 * L1: GPS Measurement Types
 * ------------------------------------------------------------------------- */

/** GPS position-velocity measurement (loosely coupled). */
typedef struct {
    ins_geodetic_t pos;
    ins_vec3_t     vel_ned;
    double         pos_std[3];
    double         vel_std[3];
    double         time;
    int            num_sats;
    double         hdop;
    double         vdop;
} ins_gps_measurement_t;

/** GPS pseudorange measurement (tightly coupled). */
typedef struct {
    int     sv_id;
    double  pseudorange;
    double  pseudorange_rate;
    double  sat_pos[3];
    double  sat_vel[3];
    double  meas_var;
    double  clock_bias;
} ins_gps_pseudorange_t;

/* -------------------------------------------------------------------------
 * L1: Kalman Filter State for INS/GPS
 * ------------------------------------------------------------------------- */

/**
 * Error-state Kalman filter (15-state).
 * x = [delta_pos^T, delta_vel^T, psi^T, accel_bias^T, gyro_bias^T]
 *   0-2:   Position error (N,E,D) [m]
 *   3-5:   Velocity error (N,E,D) [m/s]
 *   6-8:   Attitude error (psi) [rad]
 *   9-11:  Accel bias error [m/s^2]
 *   12-14: Gyro bias error [rad/s]
 */
#define INS_KF_STATE_DIM  15

typedef struct {
    double x[INS_KF_STATE_DIM];
    double P[225];
    double Q_diag[INS_KF_STATE_DIM];
    double time;
    int    initialized;
} ins_kf_state_t;

/* -------------------------------------------------------------------------
 * L5: Error-State Kalman Filter Functions
 * ------------------------------------------------------------------------- */

/** Initialize 15-state error-state Kalman filter. */
void ins_kf_init(ins_kf_state_t *kf,
                 double sigma_pos, double sigma_vel, double sigma_att);

/**
 * Predict step: propagate error state and covariance.
 * P_pred = Phi * P * Phi^T + Q_d
 */
void ins_kf_predict(ins_kf_state_t *kf,
                    const ins_mat3_t *C_body_ned,
                    const ins_vec3_t *f_body,
                    double lat, double alt, double dt,
                    double Q_accel_psd, double Q_gyro_psd);

/**
 * Update step: incorporate GPS position-velocity measurement.
 *
 * Standard Kalman update:
 *   K = P*H^T*inv(H*P*H^T + R)
 *   x = x + K*(z - H*x)
 *   P = (I - K*H)*P
 */
int ins_kf_update_gps_pos_vel(ins_kf_state_t *kf,
                               const ins_gps_measurement_t *gps,
                               const ins_nav_solution_t *ins);

/**
 * Apply estimated error corrections to INS solution (closed-loop).
 * After correction, position/velocity/attitude errors reset to zero.
 */
void ins_kf_apply_correction(ins_kf_state_t *kf, ins_nav_solution_t *ins);

/* -------------------------------------------------------------------------
 * L6: Loosely Coupled Integration Loop
 * ------------------------------------------------------------------------- */

/**
 * Run closed-loop loosely coupled GPS/INS integration.
 *
 * Timing: IMU 100 Hz, GPS 1 Hz.
 * At each IMU step: INS propagate + KF predict.
 * At each GPS step: KF update + correction apply.
 */
int ins_integrate_loose(ins_mech_state_t *mech_state,
                         ins_kf_state_t *kf,
                         const ins_imu_sample_t *imu_data,
                         size_t num_imu,
                         const ins_gps_measurement_t *gps_data,
                         size_t num_gps,
                         size_t imu_per_gps,
                         ins_nav_solution_t *output);

/* -------------------------------------------------------------------------
 * L6: Zero-Velocity Update (ZUPT)
 * ------------------------------------------------------------------------- */

/**
 * Detect zero-velocity condition using SHOE detector.
 * T(n) = (1/W) * sum [|a_k - g*u_avg|^2 / sigma_a^2 + |w_k|^2 / sigma_w^2]
 *
 * Reference: Skog et al. (2010), IEEE Trans. Biomed. Eng., 57(11).
 */
int ins_zupt_detect(const ins_vec3_t *accel, const ins_vec3_t *gyro,
                     size_t window_len,
                     double sigma_a, double sigma_w,
                     double threshold);

/** Apply zero-velocity update to Kalman filter. */
int ins_kf_update_zupt(ins_kf_state_t *kf, const ins_nav_solution_t *ins);

/* -------------------------------------------------------------------------
 * L7: GNSS-denied navigation duration estimation
 * ------------------------------------------------------------------------- */

/**
 * Estimate GNSS outage budget: max time INS can navigate within
 * specified accuracy before GNSS re-acquisition is required.
 *
 * @param model            IMU error model
 * @param max_pos_error    Max allowable horizontal position error [m]
 * @return                 Maximum outage duration [s]
 */
double ins_gnss_outage_budget(const ins_imu_error_model_t *model,
                               double max_pos_error);

#ifdef __cplusplus
}
#endif
#endif /* INS_INTEGRATION_H */
