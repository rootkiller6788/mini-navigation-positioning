/**
 * geomag_kalman.h -- Kalman Filtering for Magnetic-Aided Navigation
 *
 * L5: Kalman filter algorithms (standard, extended, unscented)
 * L6: INS/Magnetometer integration for drift correction
 *
 * The Extended Kalman Filter (EKF) fuses inertial navigation (INS)
 * with magnetic measurements to bound position error growth.
 *
 * State vector (15-DOF typical):
 *   x = [pos(3), vel(3), att(3), acc_bias(3), gyro_bias(3)]^T
 *
 * Measurement model:
 *   z_mag = h(x) + v = [F(x), D(x), I(x)]^T + v
 *
 * Where F=total field, D=declination, I=inclination are nonlinear
 * functions of position computed via IGRF model.
 *
 * Process model: standard INS strapdown equations.
 *
 * Reference:
 *   Maybeck, "Stochastic Models, Estimation, and Control" (1979)
 *   Brown & Hwang, "Introduction to Random Signals and Applied
 *     Kalman Filtering" (2012)
 *   Groves, "Principles of GNSS, Inertial, and Multisensor
 *     Integrated Navigation Systems" (2013)
 *   Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation"
 *     Proc. IEEE (2004)
 */

#ifndef GEOMAG_KALMAN_H
#define GEOMAG_KALMAN_H

#include "geomag_core.h"

/* ========================================================================
 * L2: Kalman filter state and parameter structures
 * ======================================================================== */

/**
 * L2: Standard Kalman filter state (for linear systems).
 */
typedef struct {
    int      n_states;        /* Dimension of state vector              */
    int      n_meas;          /* Dimension of measurement vector        */
    double  *x;               /* State vector [n_states]                */
    double  *P;               /* Error covariance matrix [n_states*n_states] */
    double  *F;               /* State transition matrix [n_states*n_states] */
    double  *H;               /* Measurement matrix [n_meas*n_states]    */
    double  *Q;               /* Process noise covariance [n_states*n_states] */
    double  *R;               /* Measurement noise covariance [n_meas*n_meas] */
    double  *K;               /* Kalman gain [n_states*n_meas]          */
    double  *I;               /* Identity matrix [n_states*n_states]    */
    /* Scratch buffers */
    double  *scratch_nn;      /* n_states * n_states buffer             */
    double  *scratch_nm;      /* n_states * n_meas buffer               */
    double  *scratch_mm;      /* n_meas * n_meas buffer                 */
} KalmanFilter;

/* ========================================================================
 * L2: Extended Kalman Filter state
 *
 * For nonlinear systems: x_k = f(x_{k-1}) + w_k,  z_k = h(x_k) + v_k
 * Uses Jacobians F = df/dx, H = dh/dx evaluated at current estimate.
 * ======================================================================== */
typedef struct {
    int      n_states;
    int      n_meas;
    double  *x;               /* State vector                           */
    double  *P;               /* Error covariance                       */
    double  *Q;               /* Process noise covariance               */
    double  *R;               /* Measurement noise covariance           */
    double  *K;               /* Kalman gain                            */
    double  *I;               /* Identity matrix                        */
    /* Scratch buffers */
    double  *scratch_nn;
    double  *scratch_nm;
    double  *scratch_mm;
    /* Jacobian buffers */
    double  *F_jac;           /* df/dx Jacobian                         */
    double  *H_jac;           /* dh/dx Jacobian                         */
} ExtendedKalmanFilter;

/* ========================================================================
 * L5: Standard (linear) Kalman filter API
 * ======================================================================== */

/**
 * L5: Initialize a linear Kalman filter.
 *
 * Allocates memory for all matrices. Caller must free with kalman_free().
 *
 * @param kf      Pointer to KalmanFilter struct (caller-allocated)
 * @param n_state State dimension
 * @param n_meas  Measurement dimension
 * @return 0 on success, -1 on allocation failure
 */
int kalman_init(KalmanFilter *kf, int n_state, int n_meas);

/**
 * L5: Free Kalman filter memory.
 */
void kalman_free(KalmanFilter *kf);

/**
 * L5: Set initial state and covariance.
 *
 * @param kf  Kalman filter
 * @param x0  Initial state [n_states]
 * @param P0  Initial covariance [n_states*n_states] (row-major)
 */
void kalman_set_initial(KalmanFilter *kf, const double *x0, const double *P0);

/**
 * L5: Kalman filter predict step.
 *
 * x_pred = F * x
 * P_pred = F * P * F^T + Q
 *
 * @param kf  Kalman filter with F and Q already set
 */
void kalman_predict(KalmanFilter *kf);

/**
 * L5: Kalman filter update step.
 *
 * y = z - H * x_pred                   (innovation)
 * S = H * P_pred * H^T + R             (innovation covariance)
 * K = P_pred * H^T * S^{-1}            (Kalman gain)
 * x = x_pred + K * y                    (state update)
 * P = (I - K*H) * P_pred               (covariance update)
 *
 * Joseph form is used for numerical stability:
 * P = (I - K*H)*P_pred*(I - K*H)^T + K*R*K^T
 *
 * @param kf  Kalman filter with H and R set
 * @param z   Measurement vector [n_meas]
 */
void kalman_update(KalmanFilter *kf, const double *z);

/**
 * L5: Set state transition matrix F.
 */
void kalman_set_F(KalmanFilter *kf, const double *F);

/**
 * L5: Set measurement matrix H.
 */
void kalman_set_H(KalmanFilter *kf, const double *H);

/**
 * L5: Set process noise covariance Q.
 */
void kalman_set_Q(KalmanFilter *kf, const double *Q);

/**
 * L5: Set measurement noise covariance R.
 */
void kalman_set_R(KalmanFilter *kf, const double *R);

/* ========================================================================
 * L5: Extended Kalman Filter API
 * ======================================================================== */

/**
 * L5: Initialize an Extended Kalman Filter.
 */
int ekf_init(ExtendedKalmanFilter *ekf, int n_state, int n_meas);

/**
 * L5: Free EKF memory.
 */
void ekf_free(ExtendedKalmanFilter *ekf);

/**
 * L5: Set EKF initial state and covariance.
 */
void ekf_set_initial(ExtendedKalmanFilter *ekf, const double *x0, const double *P0);

/**
 * L5: EKF predict step with user-provided nonlinear state transition.
 *
 * The caller provides a function f(x, dt) that computes the next state
 * and the Jacobian F = df/dx evaluated at current estimate.
 *
 * x_pred = f(x, dt)
 * P_pred = F * P * F^T + Q
 *
 * @param ekf       Extended Kalman filter
 * @param f         State transition function: f(x_current, dt, x_pred_out, F_jac_out)
 * @param dt        Time step
 * @param user_data User context passed to f()
 */
void ekf_predict(ExtendedKalmanFilter *ekf,
                  void (*f)(const double *x, double dt,
                            double *x_pred, double *F_jac, void *user_data),
                  double dt, void *user_data);

/**
 * L5: EKF update step with user-provided nonlinear measurement function.
 *
 * The caller provides a function h(x) that computes predicted measurement
 * and the Jacobian H = dh/dx at current estimate.
 *
 * y = z - h(x)                           (innovation)
 * S = H*P*H^T + R                        (innovation covariance)
 * K = P*H^T*inv(S)                       (Kalman gain)
 * x = x + K*y                            (state update)
 * P = (I - K*H)*P                        (covariance update, Joseph form)
 *
 * @param ekf       Extended Kalman filter
 * @param h         Measurement function: h(x, z_pred_out, H_jac_out, user_data)
 * @param z         Actual measurement vector [n_meas]
 * @param user_data User context
 */
void ekf_update(ExtendedKalmanFilter *ekf,
                 void (*h)(const double *x, double *z_pred, double *H_jac,
                           void *user_data),
                 const double *z, void *user_data);

/**
 * L5: Set EKF process and measurement noise covariance.
 */
void ekf_set_Q(ExtendedKalmanFilter *ekf, const double *Q);
void ekf_set_R(ExtendedKalmanFilter *ekf, const double *R);

/* ========================================================================
 * L6: INS/Magnetometer integration -- specialized EKF
 * ======================================================================== */

/**
 * L3: INS error state vector indices (15-state EKF).
 *
 * 0-2:   position error (NED) [m]
 * 3-5:   velocity error (NED) [m/s]
 * 6-8:   attitude error (roll, pitch, yaw) [rad]
 * 9-11:  accelerometer bias [m/s^2]
 * 12-14: gyroscope bias [rad/s]
 */
#define EKF_IDX_POS_N     0
#define EKF_IDX_POS_E     1
#define EKF_IDX_POS_D     2
#define EKF_IDX_VEL_N     3
#define EKF_IDX_VEL_E     4
#define EKF_IDX_VEL_D     5
#define EKF_IDX_ATT_R     6
#define EKF_IDX_ATT_P     7
#define EKF_IDX_ATT_Y     8
#define EKF_IDX_ACC_BIAS_N  9
#define EKF_IDX_ACC_BIAS_E  10
#define EKF_IDX_ACC_BIAS_D  11
#define EKF_IDX_GYR_BIAS_R  12
#define EKF_IDX_GYR_BIAS_P  13
#define EKF_IDX_GYR_BIAS_Y  14
#define EKF_N_STATES       15
#define EKF_N_MEAS_MAG     3    /* F, D, I or Bx, By, Bz               */

/**
 * L6: Magnetic measurement Jacobian computation.
 *
 * Computes H = dh/dx where h(x) gives predicted magnetic measurement
 * from position. The magnetic field only depends on position, so:
 *
 * H = [dF/dlat, dF/dlon, dF/dalt, 0...0]
 *
 * The spatial derivatives are computed numerically from IGRF model.
 * H is a 3 x n_states matrix (row-major).
 *
 * @param model     IGRF model for field prediction
 * @param lat       Current latitude estimate [deg]
 * @param lon       Current longitude estimate [deg]
 * @param alt       Current altitude estimate [m]
 * @param H_jac     Output Jacobian [n_meas * n_states] (row-major)
 * @param n_states  Total state dimension
 */
void mag_measurement_jacobian(const IGRFModel *model,
                               double lat, double lon, double alt,
                               double *H_jac, int n_states);

/**
 * L6: Predict magnetic measurement from state.
 *
 * z_pred = [F(lat,lon,alt), D(lat,lon,alt), I(lat,lon,alt)]^T
 *
 * @param model   IGRF model
 * @param lat,lon,alt Position estimate
 * @param z_pred  Output: [F_nT, D_deg, I_deg]
 */
void predict_magnetic_measurement(const IGRFModel *model,
                                   double lat, double lon, double alt,
                                   double z_pred[3]);

/**
 * L6: Full INS/MAG EKF measurement update.
 *
 * Convenience function that:
 *   1. Predicts magnetic measurement from current position estimate
 *   2. Computes Jacobian
 *   3. Performs EKF update
 *
 * @param ekf        Extended Kalman filter (15-state)
 * @param model      IGRF model
 * @param z_mag      Measured magnetic data [F_nT, D_deg, I_deg]
 * @param lat,lon,alt Current position from state
 */
void ins_mag_ekf_update(ExtendedKalmanFilter *ekf, const IGRFModel *model,
                         const double z_mag[3],
                         double lat, double lon, double alt);

/**
 * L6: Compute magnetic navigation solution accuracy metrics.
 *
 * Extracts position error ellipse from EKF covariance matrix.
 *
 * @param P           EKF covariance matrix [n_states*n_states]
 * @param sigma_major Output: semi-major axis of 1-sigma error ellipse [m]
 * @param sigma_minor Output: semi-minor axis [m]
 * @param orientation Output: orientation of major axis [deg from North]
 */
void ekf_position_accuracy(const double *P, int n_states,
                            double *sigma_major, double *sigma_minor,
                            double *orientation);

#endif /* GEOMAG_KALMAN_H */
