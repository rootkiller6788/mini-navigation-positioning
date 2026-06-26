/**
 * mini-uwb-localization: UWB Tracking Filters
 *
 * Implements Bayesian tracking filters: Extended Kalman Filter (EKF)
 * for nonlinear UWB measurement models with constant velocity,
 * constant acceleration, and coordinated turn motion models.
 * Includes Rauch-Tung-Striebel smoother for offline trajectory optimization.
 *
 * Reference: Bar-Shalom, Li, Kirubarajan (2001) "Estimation with
 *            Applications to Tracking and Navigation"
 * Reference: Julier & Uhlmann (2004) "Unscented Filtering and
 *            Nonlinear Estimation"
 *
 * Knowledge Coverage: L5 Algorithms (EKF)
 *                      L8 Advanced (RTS smoother)
 */

#ifndef UWB_TRACKING_H
#define UWB_TRACKING_H

#include "uwb_types.h"
#include "uwb_positioning.h"

#define EKF_STATE_DIM_2D    4
#define EKF_STATE_DIM_3D    6
#define EKF_STATE_DIM_3D_ACC 9
#define EKF_MAX_STATE_DIM   12

typedef enum {
    MOTION_STATIC              = 0,
    MOTION_CONSTANT_VELOCITY   = 1,
    MOTION_CONSTANT_ACCELERATION = 2,
    MOTION_COORDINATED_TURN    = 3
} motion_model_t;

typedef struct {
    int state_dim;
    int meas_dim;
    double state[EKF_MAX_STATE_DIM];
    double P[EKF_MAX_STATE_DIM * EKF_MAX_STATE_DIM];
    double Q[EKF_MAX_STATE_DIM * EKF_MAX_STATE_DIM];
    double R[64];
    double dt;
    motion_model_t motion_model;
    double turn_rate;
    int is_initialized;
    uint64_t last_update_ts;
} ekf_state_t;

typedef struct {
    double process_noise_pos;
    double process_noise_vel;
    double measurement_noise_range;
    double initial_pos_uncertainty;
    double initial_vel_uncertainty;
    int max_iterations;
    double convergence_threshold;
} ekf_config_t;

void ekf_init(ekf_state_t *ekf, const ekf_config_t *config, int dim,
              const uwb_pos3d_t *initial_pos, double dt);

/*
 * EKF prediction step: x_{k|k-1} = F*x_{k-1}, P_{k|k-1} = F*P*F^T + Q
 */
void ekf_predict(ekf_state_t *ekf, double dt);

/*
 * EKF update with range measurements (nonlinear measurement model).
 * Innovation: nu = z - h(x), S = H*P*H^T + R, K = P*H^T*S^-1
 * x_k = x + K*nu, P_k = (I - K*H)*P
 */
void ekf_update_range(ekf_state_t *ekf, const uwb_pos3d_t *anchors,
                      const double *ranges, int num_meas);

/*
 * EKF update with TDoA measurements.
 * h_i(x) = ||x - a_{i+1}|| - ||x - a_1|| (relative to anchor 1)
 */
void ekf_update_tdoa(ekf_state_t *ekf, const uwb_pos3d_t *anchors,
                     const double *tdoa_values, int num_tdoa);

void ekf_get_position(const ekf_state_t *ekf, uwb_pos3d_t *pos);
void ekf_get_velocity(const ekf_state_t *ekf, uwb_pos3d_t *vel);
void ekf_get_covariance(const ekf_state_t *ekf, uwb_covariance_t *cov);

/*
 * Normalized Innovation Squared (NIS) for filter consistency.
 * epsilon = nu^T*S^-1*nu, expected E[epsilon] = meas_dim.
 */
double ekf_normalized_innovation_squared(const ekf_state_t *ekf);

/*
 * State transition matrix F for given motion model.
 * CV (2D): F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]
 */
void motion_model_transition_matrix(int state_dim, double dt,
                                    motion_model_t motion, double turn_rate,
                                    double *F_out);

/*
 * Process noise covariance Q (discrete white noise acceleration).
 * Reference: Bar-Shalom et al. (2001), Section 6.3
 */
void motion_model_process_noise(int state_dim, double dt,
                                double pos_noise, double vel_noise,
                                double *Q_out);

void dead_reckon_position(uwb_pos3d_t *pos, const uwb_pos3d_t *vel, double dt);

/*
 * Rauch-Tung-Striebel fixed-interval smoother.
 * Backward recursion: C_k = P_k*F_k^T*inv(P_{k+1|k})
 * x_k^s = x_k + C_k*(x_{k+1}^s - x_{k+1|k})
 * Reference: Rauch, Tung, Striebel (1965) AIAA Journal
 */
void rts_smoother(ekf_state_t *states, int num_states, double dt);

#endif /* UWB_TRACKING_H */
