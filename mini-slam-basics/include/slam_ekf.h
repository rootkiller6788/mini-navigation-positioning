#ifndef SLAM_EKF_H
#define SLAM_EKF_H

/**
 * @file    slam_ekf.h
 * @brief   Extended Kalman Filter SLAM (EKF-SLAM)
 *
 * EKF-SLAM is the earliest and most studied SLAM algorithm.
 * It uses a single EKF to jointly estimate robot pose and all landmark
 * positions. The state grows with each new landmark: dim = 3 + 2×N.
 *
 * Reference:
 *   Smith, Self & Cheeseman (1990) "Estimating Uncertain Spatial
 *     Relationships in Robotics" — foundational EKF-SLAM paper.
 *   Dissanayake et al. (2001) "A Solution to the SLAM Problem" —
 *     convergence proof for EKF-SLAM.
 *   Thrun, Burgard & Fox (2005) "Probabilistic Robotics", Ch.10.
 *
 * Knowledge Coverage:
 *   L4: EKF linearization, Bayes filter for SLAM
 *   L5: EKF-SLAM prediction, update, state augmentation
 *   L6: Landmark-based SLAM with range-bearing
 */

#include "slam_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L2: EKF-SLAM Core Concepts
 * ========================================================================= */

/**
 * @brief Initialize EKF-SLAM state
 *
 * State vector: μ = [x, y, θ]^T (no landmarks initially).
 * Covariance: Σ = diag(σ_x², σ_y², σ_θ²).
 *
 * Set σ_x=σ_y=0 means known starting position (origin).
 * Set σ_x=σ_y=∞ means unknown starting position.
 *
 * @param state   allocated EKF state structure (populated by this call)
 * @param config  SLAM configuration
 * @param init_pose  initial robot pose
 * @param sigma_pos  initial position uncertainty [m]
 * @param sigma_theta initial heading uncertainty [rad]
 * @return SLAM_OK on success
 */
int slam_ekf_init(slam_ekf_state_t *state,
                  const slam_config_t *config,
                  const slam_pose2d_t *init_pose,
                  double sigma_pos,
                  double sigma_theta);

/**
 * @brief Free EKF-SLAM state memory
 */
void slam_ekf_free(slam_ekf_state_t *state);

/* =========================================================================
 * L5: EKF-SLAM Prediction (Motion Update)
 * ========================================================================= */

/**
 * @brief EKF prediction step with velocity motion model
 *
 * State propagation (nonlinear):
 *   μ̅ = g(μ, u) = μ ⊕ [v·Δt·cos(θ)·sinc(ω·Δt/2),
 *                       v·Δt·sin(θ)·sinc(ω·Δt/2),
 *                       ω·Δt]^T
 *
 * Jacobian G = ∂g/∂μ:
 *   G = [1 0 -v·Δt·sin(θ + ω·Δt/2);
 *        0 1  v·Δt·cos(θ + ω·Δt/2);
 *        0 0   1                       ]
 *
 * Covariance prediction:
 *   Σ̅ = G·Σ·G^T + F_x^T·R·F_x
 *
 * where R = motion noise covariance, F_x maps to full state.
 *
 * Complexity: O(N²) due to dense covariance update.
 *
 * @param state   EKF state (updated in-place)
 * @param vel     velocity control input
 * @return SLAM_OK on success
 */
int slam_ekf_predict_velocity(slam_ekf_state_t *state,
                               const slam_velocity_t *vel);

/**
 * @brief EKF prediction step with odometry motion model
 *
 * Uses the probabilistic odometry model (Thrun §5.4).
 * Motion Jacobian and noise depend on delta values.
 *
 * @param state  EKF state (updated in-place)
 * @param odom   odometry measurement
 * @return SLAM_OK on success
 */
int slam_ekf_predict_odometry(slam_ekf_state_t *state,
                               const slam_odometry_t *odom);

/* =========================================================================
 * L5: EKF-SLAM Update (Measurement Correction)
 * ========================================================================= */

/**
 * @brief EKF update with range-bearing observation of a known landmark
 *
 * Innovation (measurement residual):
 *   ν = z − h(μ̅, j)   where j is the landmark index
 *
 * Observation model h(μ, j) for landmark j at (m_jx, m_jy):
 *   r̂ = sqrt((m_jx − x)² + (m_jy − y)²)
 *   φ̂ = atan2(m_jy − y,  m_jx − x) − θ
 *
 * Jacobian H of h w.r.t state μ:
 *   H = [∂r̂/∂x  ∂r̂/∂y  0    0...0  ∂r̂/∂m_jx  ∂r̂/∂m_jy  0...0;
 *        ∂φ̂/∂x  ∂φ̂/∂y  −1   0...0  ∂φ̂/∂m_jx  ∂φ̂/∂m_jy  0...0 ]
 *
 * where:
 *   dx = m_jx − x,  dy = m_jy − y,  q = dx² + dy²
 *   ∂r̂/∂x = −dx/√q,  ∂r̂/∂y = −dy/√q
 *   ∂φ̂/∂x =  dy/q,   ∂φ̂/∂y = −dx/q
 *   ∂r̂/∂m_jx = dx/√q, ∂r̂/∂m_jy = dy/√q
 *   ∂φ̂/∂m_jx = −dy/q, ∂φ̂/∂m_jy = dx/q
 *
 * Kalman gain: K = Σ̅·H^T·(H·Σ̅·H^T + Q)^{-1}
 * State update:  μ = μ̅ + K·ν
 * Covariance:    Σ = (I − K·H)·Σ̅
 *
 * Complexity: O(N²) — covariance update touches all elements.
 *
 * @param state  EKF state (updated in-place)
 * @param obs    range-bearing observation
 * @param lm_idx landmark index in state vector (0-based)
 * @return SLAM_OK, SLAM_ERR_NO_ASSOC if landmark not found
 */
int slam_ekf_update_rb(slam_ekf_state_t *state,
                       const slam_obs_rb_t *obs,
                       int lm_idx);

/**
 * @brief EKF update with unknown data association
 *
 * Finds the best-matching landmark using Mahalanobis gating,
 * then delegates to slam_ekf_update_rb.
 *
 * Mahalanobis distance: d² = ν^T·S^{-1}·ν  ∼ χ²_dim(z)
 * Innovations: ν_j for each candidate landmark j
 * Innovation covariance: S_j = H_j·Σ̅·H_j^T + Q
 *
 * Accept if d² < γ (e.g., γ=5.991 for 95% confidence, 2 DOF).
 *
 * @param state     EKF state (updated in-place)
 * @param obs       range-bearing observation
 * @param matched_id output: matched landmark id (set to -1 if none)
 * @return SLAM_OK on success, SLAM_ERR_NO_ASSOC if no match found
 */
int slam_ekf_update_unknown(slam_ekf_state_t *state,
                             const slam_obs_rb_t *obs,
                             int *matched_id);

/* =========================================================================
 * L5: EKF-SLAM State Augmentation
 * ========================================================================= */

/**
 * @brief Add a new landmark to the EKF state
 *
 * For a new observation z = [r, φ]^T, the landmark position is:
 *   m_x = x + r·cos(φ + θ)
 *   m_y = y + r·sin(φ + θ)
 *
 * State augmentation:
 *   μ_new = [μ_old; m_x, m_y]^T
 *
 * Covariance augmentation — compute Jacobian of inverse observation
 * model w.r.t. state and observation noise:
 *
 *   d(m_x)/d(x,y,θ,r,φ):
 *     [1, 0, −r·sin(φ+θ), cos(φ+θ), −r·sin(φ+θ);
 *      0, 1,  r·cos(φ+θ), sin(φ+θ),  r·cos(φ+θ)]
 *
 * New covariance:
 *   Σ_new = [ Σ_old         Σ_old·G_rl^T                          ;
 *             G_rl·Σ_old     G_rl·Σ_old·G_rl^T + G_z·Q·G_z^T    ]
 *
 * Complexity: O(N) for mean, O(N²) for covariance.
 *
 * @param state  EKF state (augmented in-place)
 * @param obs    observation that triggered new landmark creation
 * @param new_id output: id assigned to new landmark
 * @return SLAM_OK on success
 */
int slam_ekf_augment(slam_ekf_state_t *state,
                     const slam_obs_rb_t *obs,
                     int *new_id);

/**
 * @brief EKF-SLAM full step: predict → associate → update or augment
 *
 * This is the main SLAM loop entry point for EKF-based systems.
 *
 * Algorithm:
 *   1. Predict using velocity/odometry
 *   2. For each observation:
 *      a. Mahalanobis gate for known landmarks
 *      b. If match found → EKF update
 *      c. If no match → state augmentation (new landmark)
 *   3. Return updated state
 *
 * @param state       EKF state
 * @param vel         control input (or NULL to skip prediction)
 * @param observations array of range-bearing observations
 * @param num_obs     number of observations
 * @return SLAM_OK on success
 */
int slam_ekf_step(slam_ekf_state_t *state,
                  const slam_velocity_t *vel,
                  const slam_obs_rb_t *observations,
                  int num_obs);

/* =========================================================================
 * L6: EKF-SLAM Consistency Measures
 * ========================================================================= */

/**
 * @brief Compute Normalized Estimation Error Squared (NEES)
 *
 * NEES = (x̂ − x_true)^T·Σ^{-1}·(x̂ − x_true)  ∼ χ²_dim(x)
 *
 * NEES tests filter consistency: the estimation error should be
 * commensurate with the filter's reported covariance.
 * Consistent filter: E[NEES] = dim(x).
 *
 * Over N Monte Carlo runs, N·NEES ∼ χ²_{N·dim(x)}.
 *
 * @param state    EKF state estimate
 * @param true_pose ground truth robot pose
 * @param nees_out output: NEES value
 * @return SLAM_OK
 */
int slam_ekf_nees(const slam_ekf_state_t *state,
                  const slam_pose2d_t *true_pose,
                  double *nees_out);

/**
 * @brief Compute Normalized Innovation Squared (NIS)
 *
 * NIS = ν^T·S^{-1}·ν  ∼ χ²_dim(z)
 *
 * NIS tests filter consistency online (no ground truth needed).
 * Consistent filter: E[NIS] = dim(z).
 *
 * @param obs      observation
 * @param state    EKF state before update
 * @param lm_idx   landmark index
 * @param nis_out  output: NIS value
 * @return SLAM_OK
 */
int slam_ekf_nis(const slam_obs_rb_t *obs,
                 const slam_ekf_state_t *state,
                 int lm_idx,
                 double *nis_out);

/* =========================================================================
 * L4: Key Theorems Implemented
 * ========================================================================= */

/**
 * @brief Verify EKF-SLAM monotonicity property
 *
 * Theorem (Dissanayake et al. 2001): In EKF-SLAM with known data
 * association, the determinant of any submatrix of the map covariance
 * matrix is monotonically non-increasing. That is, the uncertainty
 * of relative landmark positions never increases.
 *
 * This function computes det(Σ_ll) before/after an update and verifies
 * monotonicity: det(Σ_ll_new) ≤ det(Σ_ll_old).
 *
 * @param old_cov  covariance before update
 * @param new_cov  covariance after update
 * @param dim      state dimension
 * @return 1 if monotonic decrease holds, 0 otherwise
 */
int slam_ekf_verify_monotonicity(const double *old_cov,
                                  const double *new_cov,
                                  int dim);

#ifdef __cplusplus
}
#endif

#endif /* SLAM_EKF_H */
