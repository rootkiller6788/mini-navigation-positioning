#ifndef SLAM_FASTSLAM_H
#define SLAM_FASTSLAM_H

/**
 * @file    slam_fastslam.h
 * @brief   FastSLAM: Rao-Blackwellized Particle Filter for SLAM
 *
 * FastSLAM exploits the conditional independence of landmarks given
 * the robot trajectory: p(Θ, x_{1:t} | z_{1:t}, u_{1:t})
 *   = p(x_{1:t} | z_{1:t}, u_{1:t}) · Π_j p(θ_j | x_{1:t}, z_{1:t})
 *
 * This factorization allows:
 *   - Particle filter for robot trajectory (p(x_{1:t}|...))
 *   - Independent EKFs per landmark per particle (p(θ_j|...))
 *
 * Variants:
 *   FastSLAM 1.0: Uses motion model as proposal distribution
 *   FastSLAM 2.0: Uses motion+measurement as proposal (lower variance)
 *
 * Reference:
 *   Montemerlo et al. (2002) "FastSLAM: A Factored Solution to the
 *     Simultaneous Localization and Mapping Problem", AAAI.
 *   Montemerlo, Thrun et al. (2003) "FastSLAM 2.0", IJCAI.
 *   Thrun, Burgard & Fox (2005) "Probabilistic Robotics", Ch.13.
 *
 * Knowledge Coverage:
 *   L3: Rao-Blackwellization, importance sampling
 *   L4: Total probability theorem for particle filters
 *   L5: FastSLAM 1.0/2.0 particle filter
 *   L6: Particle-filter-based SLAM
 */

#include "slam_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L2: FastSLAM System
 * ========================================================================= */

/**
 * @brief Initialize FastSLAM particle set
 *
 * All particles start at init_pose with equal weight 1/M.
 * Motion noise α₁..α₄ and measurement noise σ_r, σ_b come from config.
 *
 * @param particles  output: allocated particle array
 * @param num_particles number of particles M
 * @param init_pose  initial robot pose
 * @param config     SLAM configuration
 * @return SLAM_OK
 */
int slam_fastslam_init(slam_particle_t **particles,
                       int num_particles,
                       const slam_pose2d_t *init_pose,
                       const slam_config_t *config);

/**
 * @brief Free FastSLAM particle set
 */
void slam_fastslam_free(slam_particle_t *particles,
                         int num_particles);

/* =========================================================================
 * L5: FastSLAM 1.0 — Motion Model Proposal
 * ========================================================================= */

/**
 * @brief FastSLAM 1.0: sample pose from motion model
 *
 * Proposal distribution: q = p(x_t | x_{t-1}, u_t)
 *
 * Sample from probabilistic velocity model:
 *   v̂   = v + ε_v,   ε_v   ∼ N(0, σ_v²)
 *   ω̂   = ω + ε_ω,   ε_ω   ∼ N(0, σ_ω²)
 *   γ̂   = ω̂ + ε_γ,   ε_γ   ∼ N(0, σ_γ²)   (additional heading noise)
 *
 * Then: x' = x − (v̂/ω̂)·sin(θ) + (v̂/ω̂)·sin(θ + ω̂·Δt)
 *       y' = y + (v̂/ω̂)·cos(θ) − (v̂/ω̂)·cos(θ + ω̂·Δt)
 *       θ' = θ + ω̂·Δt + γ̂·Δt
 *
 * For |ω̂| < 1e-6 (straight-line motion), use small-angle limit:
 *       x' = x + v̂·Δt·cos(θ)
 *       y' = y + v̂·Δt·sin(θ)
 *       θ' = θ
 *
 * @param particle  particle to update (pose updated in-place)
 * @param vel       velocity control
 * @param sigma_v   linear velocity noise std
 * @param sigma_omega angular velocity noise std
 * @return SLAM_OK
 */
int slam_fastslam1_sample_pose(slam_particle_t *particle,
                                const slam_velocity_t *vel,
                                double sigma_v,
                                double sigma_omega);

/**
 * @brief FastSLAM 1.0: update landmark EKFs with known association
 *
 * For each observation matched to a landmark in this particle,
 * perform a standard EKF update for that landmark only.
 * This is O(M·K) total (M particles, K landmarks per particle),
 * vs O(M·N²) for full EKF-SLAM.
 *
 * Innovation and Kalman gain are 2×2 computations.
 *
 * @param particle  particle to update
 * @param obs       range-bearing observation
 * @param lm_idx    landmark index in this particle's map
 * @return SLAM_OK
 */
int slam_fastslam1_update_landmark(slam_particle_t *particle,
                                    const slam_obs_rb_t *obs,
                                    int lm_idx);

/**
 * @brief FastSLAM 1.0: compute importance weight
 *
 * Weight update (log domain for numeric stability):
 *   w_t^(k) = w_{t-1}^(k) + log p(z_t | x_t^(k), u_t)
 *
 * For range-bearing observation with matched landmark:
 *   log w = log w_old − ½·ν^T·S^{-1}·ν − ½·log|2π·S|
 *
 * where ν = innovation, S = innovation covariance (2×2).
 *
 * @param particle   particle whose weight is updated
 * @param obs        observation
 * @param lm_idx     matched landmark index
 * @return SLAM_OK
 */
int slam_fastslam1_update_weight(slam_particle_t *particle,
                                  const slam_obs_rb_t *obs,
                                  int lm_idx);

/* =========================================================================
 * L5: FastSLAM 2.0 — Improved Proposal
 * ========================================================================= */

/**
 * @brief FastSLAM 2.0: sample from improved proposal
 *
 * Proposal: q = p(x_t | x_{t-1}, u_t, z_t, Θ_{t-1})
 *
 * FastSLAM 2.0 incorporates the current measurement into the proposal
 * to reduce sample variance. Uses EKF linearization around the motion
 * model prediction to approximate the optimal proposal.
 *
 * Steps:
 *   1. Predict: μ̅ = g(μ, u_t), Σ̅ = G·Σ·G^T + R
 *   2. Linearize observation: ẑ = h(μ̅), H = ∂h/∂μ
 *   3. Proposal covariance: Σ_prop = (Σ̅^{-1} + H^T·Q^{-1}·H)^{-1}
 *   4. Proposal mean: μ_prop = μ̅ + Σ_prop·H^T·Q^{-1}·(z_t − ẑ)
 *   5. Sample: x_t ∼ N(μ_prop, Σ_prop)
 *
 * This typically requires fewer particles than FastSLAM 1.0
 * for equivalent performance.
 *
 * @param particle     particle to update
 * @param vel          velocity control
 * @param obs          observation (used in proposal)
 * @param lm_idx       landmark index for proposal
 * @param sigma_v      velocity noise
 * @param sigma_omega  angular noise
 * @return SLAM_OK
 */
int slam_fastslam2_sample_pose(slam_particle_t *particle,
                                const slam_velocity_t *vel,
                                const slam_obs_rb_t *obs,
                                int lm_idx,
                                double sigma_v,
                                double sigma_omega);

/* =========================================================================
 * L5: Resampling
 * ========================================================================= */

/**
 * @brief Systematic resampling (Kitagawa 1996)
 *
 * Avoids particle degeneracy by resampling particles proportional to
 * their weights. Systematic resampling has lower variance than
 * multinomial resampling.
 *
 * Algorithm:
 *   Neff = 1 / Σ(w_k²)     — effective sample size
 *   If Neff < N_thresh (e.g., M/2):
 *     r ∼ U(0, 1/M)
 *     c = w_0
 *     i = 0
 *     for m = 0..M-1:
 *       U = r + m/M
 *       while U > c: i++, c += w_i
 *       new_particles[m] = copy(particles[i])
 *       new_particles[m].weight = 1/M
 *
 * @param particles   particle set (replaced with resampled set)
 * @param num_particles  number of particles M
 * @param neff_thresh  effective sample size threshold
 * @param neff_out     output: computed effective sample size
 * @return SLAM_OK
 */
int slam_fastslam_resample(slam_particle_t **particles,
                            int num_particles,
                            double neff_thresh,
                            double *neff_out);

/**
 * @brief Low-variance resampling (Probabilistic Robotics, Table 4.4)
 *
 * Deterministic resampling with O(M) complexity.
 * More efficient than multinomial and produces lower-variance estimates.
 *
 * @param particles   particle set
 * @param num_particles M
 * @return SLAM_OK
 */
int slam_fastslam_low_variance_resample(slam_particle_t *particles,
                                         int num_particles);

/* =========================================================================
 * L6: FastSLAM Full Step
 * ========================================================================= */

/**
 * @brief FastSLAM 1.0 full step: sample → associate → update → weight → resample
 *
 * Main loop for FastSLAM 1.0.
 *
 * Algorithm:
 *   For each particle k:
 *     1. Sample pose from motion model
 *     2. For each observation:
 *        a. Compute predicted observation for each landmark
 *        b. Find best match (Mahalanobis gate)
 *        c. If matched → update landmark EKF, update weight
 *        d. If new → initialize landmark in particle
 *     3. Normalize weights (convert log-weights)
 *   Resample if Neff < N_thresh
 *
 * @param particles   particle set
 * @param num_particles M
 * @param vel         control input
 * @param observations array of observations
 * @param num_obs     number of observations
 * @param config      SLAM configuration
 * @return SLAM_OK
 */
int slam_fastslam_step(slam_particle_t *particles,
                       int num_particles,
                       const slam_velocity_t *vel,
                       const slam_obs_rb_t *observations,
                       int num_obs,
                       const slam_config_t *config);

/**
 * @brief Get best particle (highest weight) — MAP estimate of robot pose
 *
 * @param particles   particle set
 * @param num_particles M
 * @param best_pose   output: pose of best particle
 * @return SLAM_OK
 */
int slam_fastslam_best_pose(const slam_particle_t *particles,
                             int num_particles,
                             slam_pose2d_t *best_pose);

/**
 * @brief Extract map from best particle
 *
 * The best particle's landmark estimates constitute the current map.
 *
 * @param particles   particle set
 * @param num_particles M
 * @param map_out     output: map populated from best particle
 * @return SLAM_OK
 */
int slam_fastslam_extract_map(const slam_particle_t *particles,
                               int num_particles,
                               slam_map2d_t *map_out);

#ifdef __cplusplus
}
#endif

#endif /* SLAM_FASTSLAM_H */
