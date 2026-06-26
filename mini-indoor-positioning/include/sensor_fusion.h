/**
 * @file sensor_fusion.h
 * @brief Sensor fusion algorithms for indoor positioning
 *
 * Knowledge Coverage:
 *   L1 - Definitions: Kalman filter state/covariance/innovation,
 *        process noise, measurement noise, prediction, update
 *   L2 - Core Concepts: Bayesian filtering, recursive estimation,
 *        sensor fusion, complementary filtering
 *   L3 - Mathematical Structures: linear/nonlinear state-space models,
 *        Jacobians, information matrices, Cholesky decomposition
 *   L4 - Fundamental Laws: Kalman optimality (MMSE for linear Gaussian),
 *        Cramer-Rao lower bound, information inequality
 *   L5 - Algorithms: Kalman filter (KF), extended Kalman filter (EKF),
 *        unscented Kalman filter (UKF), particle filter (PF),
 *        complementary filter
 *   L6 - Canonical Problems: IMU+UWB fusion, foot-mounted INS+ZUPT,
 *        visual-inertial odometry (VIO) concept
 *
 * Reference: Kalman, "A New Approach to Linear Filtering and Prediction
 *            Problems," Trans. ASME, 1960.
 *            Welch & Bishop, "An Introduction to the Kalman Filter" (2006)
 *            Julier & Uhlmann, "Unscented Filtering and Nonlinear
 *            Estimation," Proc. IEEE, 2004.
 *            Arulampalam et al., "A Tutorial on Particle Filters for
 *            Online Nonlinear/Non-Gaussian Bayesian Tracking," IEEE TSP, 2002.
 *
 * Course Alignment:
 *   - MIT 6.003 (Signal Processing) — estimation theory
 *   - Stanford EE363 (Linear Dynamical Systems) — Kalman on steroids
 *   - Berkeley EE123 (DSP) — adaptive filtering
 *   - ETH 227-0427 (Signal Processing)
 *   - 清华 信号与系统 — 最优估计
 */

#ifndef SENSOR_FUSION_H
#define SENSOR_FUSION_H

#include "indoor_positioning.h"

/* ============================================================================
 * L1 - Definitions: Kalman Filter Structures
 * ============================================================================ */

#define KF_STATE_MAX_DIM     9    /**< Maximum state vector dimension */
#define UKF_SIGMA_PTS_MAX    19   /**< 2n+1 for n=9 */

/**
 * @brief Linear Kalman filter structure
 *
 * State-space model:
 *   x_k = F * x_{k-1} + B * u_k + w_k      (process)
 *   z_k = H * x_k + v_k                     (measurement)
 * where w ~ N(0,Q), v ~ N(0,R)
 *
 * L1 Definition: The Kalman filter is a recursive MMSE estimator
 * for linear Gaussian state-space models.
 */
typedef struct {
    double x[KF_STATE_MAX_DIM];            /**< State estimate vector */
    double P[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< State covariance */
    double F[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< State transition matrix */
    double H[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Measurement matrix */
    double Q[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Process noise covariance */
    double R[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Measurement noise covariance */
    int    n_state;       /**< State dimension */
    int    n_measure;     /**< Measurement dimension */
    int    initialized;   /**< 1 if filter has been initialized */
} kalman_filter_t;

/**
 * @brief Extended Kalman filter (EKF) structure
 *
 * L2 Concept: EKF linearizes nonlinear state transition f(x) and
 * measurement h(x) functions using first-order Taylor expansion
 * (Jacobian matrices).
 *
 * Nonlinear model:
 *   x_k = f(x_{k-1}, u_k) + w_k
 *   z_k = h(x_k) + v_k
 *
 * Jacobians: F = ∂f/∂x, H = ∂h/∂x
 */
typedef struct {
    double x[KF_STATE_MAX_DIM];            /**< State estimate */
    double P[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< State covariance */
    double Q[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Process noise covariance */
    double R[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Measurement noise covariance */
    int    n_state;
    int    n_measure;
    int    initialized;
} ekf_t;

/**
 * @brief Unscented Kalman filter (UKF) structure
 *
 * L8 Advanced: UKF propagates sigma points through nonlinear
 * functions without computing Jacobians, improving accuracy
 * over EKF for strongly nonlinear systems.
 */
typedef struct {
    double x[KF_STATE_MAX_DIM];            /**< State estimate */
    double P[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< State covariance */
    double Q[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Process noise covariance */
    double R[KF_STATE_MAX_DIM][KF_STATE_MAX_DIM]; /**< Measurement noise covariance */
    int    n_state;
    int    n_measure;
    double alpha;          /**< Sigma point spread parameter */
    double beta;           /**< Prior knowledge parameter (2 for Gaussian) */
    double kappa;          /**< Secondary scaling parameter */
    double lambda;         /**< Composite scaling: alpha^2*(n+kappa) - n */
    int    initialized;
} ukf_t;

/* ============================================================================
 * L5 - Kalman Filter Algorithm
 * ============================================================================ */

/**
 * @brief Initialize a Kalman filter with state dimension
 *
 * @param kf Filter to initialize
 * @param n_state State vector dimension
 * @param n_measure Measurement vector dimension
 *
 * Sets all matrices to identity/zero as appropriate.
 */
void kf_init(kalman_filter_t *kf, int n_state, int n_measure);

/**
 * @brief Set the state transition matrix F
 *
 * @param kf Filter to configure
 * @param F State transition matrix (n_state x n_state, row-major)
 */
void kf_set_F(kalman_filter_t *kf, const double *F);

/**
 * @brief Set the measurement matrix H
 *
 * @param kf Filter to configure
 * @param H Measurement matrix (n_measure x n_state, row-major)
 */
void kf_set_H(kalman_filter_t *kf, const double *H);

/**
 * @brief Set process noise covariance Q
 *
 * @param kf Filter
 * @param Q n_state x n_state covariance (row-major)
 */
void kf_set_Q(kalman_filter_t *kf, const double *Q);

/**
 * @brief Set measurement noise covariance R
 *
 * @param kf Filter
 * @param R n_measure x n_measure covariance (row-major)
 */
void kf_set_R(kalman_filter_t *kf, const double *R);

/**
 * @brief Set initial state estimate and covariance
 *
 * @param kf Filter
 * @param x0 Initial state vector (length n_state)
 * @param P0 Initial state covariance (n_state x n_state, row-major)
 */
void kf_set_initial(kalman_filter_t *kf, const double *x0, const double *P0);

/**
 * @brief Kalman filter prediction step
 *
 * x_pred = F * x
 * P_pred = F * P * F^T + Q
 *
 * @param kf Filter (state and covariance updated in place)
 *
 * L5: Kalman prediction (time update).
 * Complexity: O(n^3) for n-state, dominated by matrix multiply
 */
void kf_predict(kalman_filter_t *kf);

/**
 * @brief Kalman filter update step with measurement
 *
 * y = z - H * x           (innovation)
 * S = H * P * H^T + R     (innovation covariance)
 * K = P * H^T * S^{-1}    (Kalman gain)
 * x = x + K * y           (state update)
 * P = (I - K * H) * P     (covariance update)
 *
 * @param kf Filter (state and covariance updated)
 * @param z Measurement vector (length n_measure)
 *
 * L5: Kalman update (measurement update).
 * Complexity: O(n^2*m + m^3) for n-state, m-measure
 */
void kf_update(kalman_filter_t *kf, const double *z);

/**
 * @brief Get current state estimate
 *
 * @param kf Filter
 * @param[out] x State vector (n_state elements, pre-allocated)
 */
void kf_get_state(const kalman_filter_t *kf, double *x);

/**
 * @brief Get current state covariance diagonal (variances)
 *
 * @param kf Filter
 * @param[out] diag Diagonal elements (n_state elements, pre-allocated)
 */
void kf_get_covariance_diag(const kalman_filter_t *kf, double *diag);

/* ============================================================================
 * L5 - Extended Kalman Filter (EKF)
 * ============================================================================ */

/**
 * @brief Initialize an EKF
 */
void ekf_init(ekf_t *ekf, int n_state, int n_measure);

/**
 * @brief EKF prediction step with user-supplied state transition function
 *
 * @param ekf Filter
 * @param f Nonlinear state transition function f(x_prev, dt, user_data)
 * @param F_jacobian Jacobian of f evaluated at current state
 * @param dt Time step in seconds
 * @param user_data Opaque user data passed to f
 *
 * The caller must provide:
 * - f(x_prev, dt, user_data, x_pred): computes next state
 * - F_jacobian: n_state x n_state Jacobian of f at x_prev
 */
typedef void (*ekf_transition_fn)(const double *x, double dt,
                                  void *user_data, double *x_pred);
typedef void (*ekf_jacobian_fn)(const double *x, double dt,
                                void *user_data, double *F);

void ekf_predict(ekf_t *ekf, ekf_transition_fn f, ekf_jacobian_fn F_jacobian,
                 double dt, void *user_data);

/**
 * @brief EKF update step
 *
 * @param ekf Filter
 * @param z Measurement vector
 * @param h Nonlinear measurement function h(x, user_data, z_pred)
 * @param H Measurement Jacobian at current state
 * @param user_data Opaque user data
 */
typedef void (*ekf_measurement_fn)(const double *x, void *user_data, double *z_pred);

void ekf_update(ekf_t *ekf, const double *z,
                ekf_measurement_fn h, ekf_jacobian_fn H_jacobian,
                void *user_data);

/**
 * @brief Get EKF state
 */
void ekf_get_state(const ekf_t *ekf, double *x);

/* ============================================================================
 * L8 Advanced - Unscented Kalman Filter (UKF)
 * ============================================================================ */

/**
 * @brief Initialize UKF
 *
 * @param ukf Filter to initialize
 * @param n_state State dimension
 * @param n_measure Measurement dimension
 * @param alpha Spread parameter (1e-4 to 1, typical 1e-3)
 * @param beta Prior knowledge (2 optimal for Gaussian)
 * @param kappa Secondary scaling (0 or 3-n_state)
 */
void ukf_init(ukf_t *ukf, int n_state, int n_measure,
              double alpha, double beta, double kappa);

/**
 * @brief UKF prediction step
 *
 * Generates 2n+1 sigma points, propagates through f(x),
 * and recomputes mean and covariance.
 *
 * @param ukf Filter
 * @param f State transition function
 * @param dt Time step
 * @param user_data User data
 *
 * L8: Unscented transform for nonlinear state propagation.
 * Complexity: O(n^3) time, O(n^2) space
 */
void ukf_predict(ukf_t *ukf, ekf_transition_fn f,
                 double dt, void *user_data);

/**
 * @brief UKF update step
 *
 * @param ukf Filter
 * @param z Measurement vector
 * @param h Measurement function
 * @param user_data User data
 *
 * L8: Unscented transform for measurement update.
 */
void ukf_update(ukf_t *ukf, const double *z,
                ekf_measurement_fn h, void *user_data);

/**
 * @brief Get UKF state
 */
void ukf_get_state(const ukf_t *ukf, double *x);

/* ============================================================================
 * L5 - Complementary Filter (simple sensor fusion)
 * ============================================================================ */

/**
 * @brief Complementary filter for IMU attitude estimation
 *
 * L5: Combines high-pass filtered gyroscope integration (good short-term,
 * drifts long-term) with low-pass filtered accelerometer/magnetometer
 * (noisy short-term, stable long-term).
 *
 * angle = alpha * (angle + gyro*dt) + (1-alpha) * accel_angle
 *
 * @param current_angle Current fused angle estimate
 * @param gyro_rate Gyroscope angular rate in rad/s
 * @param accel_angle Angle from accelerometer (rad)
 * @param alpha Filter coefficient (0-1, closer to 1 trusts gyro more)
 * @param dt Sample period in seconds
 * @return Updated angle estimate
 *
 * Complexity: O(1)
 */
double complementary_filter_angle(double current_angle, double gyro_rate,
                                  double accel_angle, double alpha, double dt);

/**
 * @brief 1D complementary filter for position fusion
 *
 * Fuses high-frequency position increment (from velocity integration)
 * with low-frequency absolute position (from e.g., RSSI).
 *
 * @param fused_pos Current fused position
 * @param vel Velocity estimate in m/s
 * @param absolute_pos Position measurement from absolute source
 * @param alpha Trust coefficient (0-1)
 * @param dt Time step
 * @return Updated fused position
 */
double complementary_filter_position(double fused_pos, double vel,
                                     double absolute_pos, double alpha, double dt);

/* ============================================================================
 * L6 - Canonical Problem: IMU + Ranging Fusion
 * ============================================================================ */

/* Constant velocity model transition function (exposed for reuse) */
void cv_transition_fn(const double *x, double dt, void *user_data, double *x_pred);
void cv_jacobian_fn(const double *x, double dt, void *user_data, double *F);

/**
 * @brief Set up a constant-velocity EKF for indoor positioning
 *
 * State: [x, y, vx, vy] (4D)
 * Process: constant velocity model
 * Measurement: distance ranging (e.g., UWB/ToF)
 *
 * @param ekf Filter to configure
 * @param dt Prediction time step
 * @param process_noise_accel Process noise acceleration std (m/s^2)
 * @param measurement_noise_dist Measurement noise distance std (m)
 * @param initial_pos_x, initial_pos_y Initial position
 *
 * L6: Range-only tracking using EKF.
 */
void ekf_setup_constant_velocity_ranging(ekf_t *ekf, double dt,
                                         double process_noise_accel,
                                         double measurement_noise_dist,
                                         double initial_pos_x,
                                         double initial_pos_y);

/**
 * @brief EKF measurement prediction for range to anchor
 *
 * Computes h(x) = ||pos_anchor - pos_user|| and Jacobian H = ∂h/∂x
 */
typedef struct {
    double anchor_x;
    double anchor_y;
} ekf_range_measurement_data_t;

void ekf_range_measurement_fn(const double *x, void *user_data, double *z_pred);
void ekf_range_jacobian_fn(const double *x, double dt, void *user_data, double *F);

/* ============================================================================
 * L6 - Particle Filter
 * ============================================================================ */

#define PF_MAX_PARTICLES 500

/**
 * @brief A single particle (state hypothesis with weight)
 */
typedef struct {
    double x[KF_STATE_MAX_DIM]; /**< State vector */
    double weight;              /**< Importance weight */
} particle_t;

/**
 * @brief Particle filter for indoor positioning
 *
 * L8 Advanced: Particle filter (sequential Monte Carlo) represents
 * the posterior distribution as a set of weighted particles.
 * Handles non-Gaussian, multimodal distributions.
 *
 * Uses Sampling Importance Resampling (SIR).
 *
 * Reference: Arulampalam et al. (2002)
 */
typedef struct {
    particle_t particles[PF_MAX_PARTICLES];
    int        n_particles;
    int        n_state;
    double     effective_n_thresh; /**< Resampling threshold ratio */
    int        initialized;
} particle_filter_t;

/**
 * @brief Initialize particle filter with uniform distribution
 *
 * @param pf Particle filter
 * @param n_particles Number of particles
 * @param n_state State dimension
 * @param x_min Lower bounds for initial state
 * @param x_max Upper bounds for initial state
 */
void pf_init_uniform(particle_filter_t *pf, int n_particles, int n_state,
                     const double *x_min, const double *x_max);

/**
 * @brief Particle filter prediction (propagation) step
 *
 * @param pf Particle filter
 * @param f State transition function (applied to each particle)
 * @param dt Time step
 * @param user_data User data
 */
void pf_predict(particle_filter_t *pf, ekf_transition_fn f,
                double dt, void *user_data);

/**
 * @brief Particle filter update (weighting) step
 *
 * Updates weights based on measurement likelihood:
 *   w_i = w_i * p(z | x_i)
 *
 * @param pf Particle filter
 * @param z Measurement vector
 * @param n_measure Measurement dimension
 * @param measurement_std Measurement noise std dev (scalar, isotropic)
 */
void pf_update_gaussian(particle_filter_t *pf, const double *z,
                        int n_measure, double measurement_std);

/**
 * @brief Systematic resampling of particle filter
 *
 * Low-variance resampling algorithm.
 *
 * @param pf Particle filter (resampled in place)
 *
 * L8: Systematic resampling for particle degeneracy mitigation.
 * Complexity: O(N)
 */
void pf_resample(particle_filter_t *pf);

/**
 * @brief Get weighted mean state estimate from particle filter
 *
 * @param pf Particle filter
 * @param[out] x_mean Weighted mean (length n_state, pre-allocated)
 */
void pf_get_mean(const particle_filter_t *pf, double *x_mean);

/**
 * @brief Get weighted covariance from particle filter
 *
 * @param pf Particle filter
 * @param[out] cov Covariance matrix (n_state x n_state, row-major)
 */
void pf_get_covariance(const particle_filter_t *pf, double *cov);

/**
 * @brief Compute effective number of particles (degeneracy measure)
 *
 * N_eff = 1 / sum(w_i^2)
 *
 * @param pf Particle filter
 * @return Effective number of particles
 *
 * L8: Particle degeneracy detection.
 */
double pf_effective_particles(const particle_filter_t *pf);

/* ============================================================================
 * L3 - Matrix Utility Functions for Kalman Filtering
 * ============================================================================ */

/**
 * @brief Matrix-vector multiplication: y = A * x
 *
 * @param A Matrix (rows x cols, row-major)
 * @param x Vector (length cols)
 * @param y Output vector (length rows)
 * @param rows Number of rows
 * @param cols Number of columns
 */
void matrix_vec_mult(const double *A, const double *x, double *y,
                     int rows, int cols);

/**
 * @brief Matrix multiplication: C = A * B
 *
 * @param A Left matrix (m x k)
 * @param B Right matrix (k x n)
 * @param C Output matrix (m x n)
 * @param m Rows of A
 * @param k Cols of A / Rows of B
 * @param n Cols of B
 */
void matrix_mult(const double *A, const double *B, double *C,
                 int m, int k, int n);

/**
 * @brief Matrix transpose: B = A^T
 *
 * @param A Input matrix (rows x cols)
 * @param B Output matrix (cols x rows)
 * @param rows Rows of A
 * @param cols Cols of A
 */
void matrix_transpose(const double *A, double *B, int rows, int cols);

/**
 * @brief Cholesky decomposition: A = L * L^T
 *
 * A must be symmetric positive definite.
 *
 * @param A Input matrix (n x n, row-major), modified in-place to L
 * @param n Dimension
 * @return 0 on success, -1 if A is not positive definite
 *
 * Complexity: O(n^3/6)
 */
int cholesky_decompose(double *A, int n);

/**
 * @brief Solve L * L^T * x = b given Cholesky factor L
 *
 * @param L Lower triangular Cholesky factor (n x n)
 * @param b Right-hand side (n elements)
 * @param x Solution vector (n elements)
 * @param n Dimension
 */
void cholesky_solve(const double *L, const double *b, double *x, int n);

/**
 * @brief Compute Mahalanobis distance
 *
 * d^2 = (x - mean)^T * S^{-1} * (x - mean)
 *
 * @param x Sample vector (n elements)
 * @param mean Mean vector (n elements)
 * @param S_inv Inverse covariance matrix (n x n, row-major)
 * @param n Dimension
 * @return Squared Mahalanobis distance
 *
 * L3: Mahalanobis distance for multivariate outlier detection.
 * Complexity: O(n^2)
 */
double mahalanobis_distance(const double *x, const double *mean,
                            const double *S_inv, int n);

#endif /* SENSOR_FUSION_H */
