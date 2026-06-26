#ifndef SLAM_SENSOR_H
#define SLAM_SENSOR_H

/**
 * @file    slam_sensor.h
 * @brief   Sensor models for SLAM: observation models, motion models,
 *          inverse sensor models for occupancy grids.
 *
 * Sensor models bridge raw measurements to the SLAM state representation.
 * Each model consists of:
 *   1. Forward model: h(state, landmark) → predicted measurement
 *   2. Jacobian: ∂h/∂state for linearization
 *   3. Inverse model: h^{-1}(measurement, state) → landmark position
 *
 * Reference:
 *   Thrun, Burgard & Fox (2005) "Probabilistic Robotics", Ch.5-7.
 *   Bailey & Durrant-Whyte (2006) "SLAM Part II", IEEE RAM.
 *
 * Knowledge Coverage:
 *   L1: Range-bearing, LiDAR, odometry sensor types
 *   L2: Measurement likelihood models
 *   L3: Jacobian computation for EKF linearization
 *   L6: Beam model, likelihood field for LiDAR
 */

#include "slam_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Noise sampling utilities
 * ========================================================================= */

/**
 * @brief Sample from standard normal distribution N(0,1)
 *
 * Uses Box-Muller transform:
 *   Z = √(-2·ln(U₁))·cos(2π·U₂)
 * where U₁, U₂ ∼ U(0,1) are independent.
 *
 * @return N(0,1) sample
 */
double slam_randn(void);

/**
 * @brief Sample from N(μ, σ²)
 * @param mu    mean
 * @param sigma standard deviation
 * @return N(μ,σ²) sample
 */
double slam_randn_scaled(double mu, double sigma);

/* =========================================================================
 * L2: Velocity Motion Model
 * ========================================================================= */

/**
 * @brief Velocity motion model: g(pose, velocity, dt) → new pose
 *
 * Deterministic (noise-free) prediction:
 *   For |ω| > ε:
 *     x' = x + (v/ω)·(sin(θ+ω·dt) − sin(θ))
 *     y' = y + (v/ω)·(cos(θ) − cos(θ+ω·dt))
 *     θ' = θ + ω·dt
 *
 *   For |ω| ≤ ε (straight line):
 *     x' = x + v·dt·cos(θ)
 *     y' = y + v·dt·sin(θ)
 *     θ' = θ
 *
 * @param pose  current pose
 * @param vel   velocity command
 * @param new_pose output: predicted pose
 * @return SLAM_OK
 */
int slam_motion_model_velocity(const slam_pose2d_t *pose,
                                const slam_velocity_t *vel,
                                slam_pose2d_t *new_pose);

/**
 * @brief Velocity motion model with sampled noise
 *
 * Sample v̂ = v + ε_v, ω̂ = ω + ε_ω, γ̂ = ε_γ (additional heading noise)
 * where ε_v ∼ N(0, α₁·|v| + α₂·|ω|),
 *       ε_ω ∼ N(0, α₃·|v| + α₄·|ω|),
 *       ε_γ ∼ N(0, α₅·|v| + α₆·|ω|)
 *
 * Then apply velocity model with (v̂, ω̂, γ̂).
 *
 * @param pose   current pose
 * @param vel    commanded velocity
 * @param alpha  noise coefficients α₁..α₆
 * @param new_pose output: sampled pose
 * @return SLAM_OK
 */
int slam_motion_model_velocity_noisy(const slam_pose2d_t *pose,
                                      const slam_velocity_t *vel,
                                      const double alpha[6],
                                      slam_pose2d_t *new_pose);

/**
 * @brief Velocity motion model Jacobian (state transition matrix)
 *
 * G = ∂g/∂x evaluated at current state:
 *   G = [1, 0, −v/ω·(cos(θ) − cos(θ+ω·dt))   for |ω| > ε;
 *        0, 1, −v/ω·(sin(θ+ω·dt) − sin(θ));
 *        0, 0,  1]
 *
 *   G = [1, 0, −v·dt·sin(θ)     for |ω| ≤ ε;
 *        0, 1,  v·dt·cos(θ);
 *        0, 0,  1]
 *
 * @param pose  current pose
 * @param vel   velocity
 * @param G     output: 3×3 Jacobian, row-major
 * @return SLAM_OK
 */
int slam_motion_jacobian_velocity(const slam_pose2d_t *pose,
                                   const slam_velocity_t *vel,
                                   double G[9]);

/* =========================================================================
 * L2: Odometry Motion Model
 * ========================================================================= */

/**
 * @brief Odometry motion model: apply odometry to pose
 *
 * From odometry measurement (δ_rot1, δ_trans, δ_rot2):
 *   x'  = x + δ_trans·cos(θ + δ_rot1)
 *   y'  = y + δ_trans·sin(θ + δ_rot1)
 *   θ'  = θ + δ_rot1 + δ_rot2
 *
 * @param pose  current pose
 * @param odom  odometry measurement
 * @param new_pose output
 * @return SLAM_OK
 */
int slam_motion_model_odometry(const slam_pose2d_t *pose,
                                const slam_odometry_t *odom,
                                slam_pose2d_t *new_pose);

/**
 * @brief Odometry motion model with noise
 *
 * Sampled odometry:
 *   δ̂_rot1  = δ_rot1  − ε_{α₁·|δ_rot1| + α₂·|δ_trans|}
 *   δ̂_trans = δ_trans − ε_{α₃·|δ_trans| + α₄·(|δ_rot1|+|δ_rot2|)}
 *   δ̂_rot2  = δ_rot2  − ε_{α₁·|δ_rot2| + α₂·|δ_trans|}
 *
 * @param pose  current pose
 * @param odom  measured odometry
 * @param alpha noise coefficients α₁..α₄
 * @param new_pose output
 * @return SLAM_OK
 */
int slam_motion_model_odometry_noisy(const slam_pose2d_t *pose,
                                      const slam_odometry_t *odom,
                                      const double alpha[4],
                                      slam_pose2d_t *new_pose);

/* =========================================================================
 * L3: Range-Bearing Observation Model
 * ========================================================================= */

/**
 * @brief Range-bearing observation model: predict measurement
 *
 * h(x, m_j): predict what a sensor at robot pose x would observe
 * for landmark at position m_j.
 *
 *   r̂ = √((m_jx − x)² + (m_jy − y)²)
 *   φ̂ = atan2(m_jy − y, m_jx − x) − θ
 *
 * Both r̂ and φ̂ are wrapped: r̂ ≥ 0, φ̂ ∈ [−π, π).
 *
 * @param robot_pose  robot pose x
 * @param landmark    landmark position m_j
 * @param pred        output: predicted [range, bearing]
 * @return SLAM_OK
 */
int slam_obs_model_rb(const slam_pose2d_t *robot_pose,
                      const slam_landmark2d_t *landmark,
                      double pred[2]);

/**
 * @brief Range-bearing observation Jacobian w.r.t. robot pose (2×3)
 *
 * H_pose = ∂h/∂(x,y,θ):
 *   ∂r̂/∂x = −Δx/q,   ∂r̂/∂y = −Δy/q,   ∂r̂/∂θ = 0
 *   ∂φ̂/∂x =  Δy/q²,  ∂φ̂/∂y = −Δx/q²,  ∂φ̂/∂θ = −1
 *
 * where Δx = m_jx − x, Δy = m_jy − y, q = √(Δx² + Δy²)
 *
 * @param robot_pose robot pose
 * @param landmark   landmark
 * @param H_pose     output: 2×3 Jacobian (row-major)
 * @return SLAM_OK
 */
int slam_obs_jacobian_rb_pose(const slam_pose2d_t *robot_pose,
                               const slam_landmark2d_t *landmark,
                               double H_pose[6]);

/**
 * @brief Range-bearing observation Jacobian w.r.t. landmark (2×2)
 *
 * H_lm = ∂h/∂(m_jx, m_jy):
 *   ∂r̂/∂m_jx =  Δx/q,  ∂r̂/∂m_jy =  Δy/q
 *   ∂φ̂/∂m_jx = −Δy/q², ∂φ̂/∂m_jy =  Δx/q²
 *
 * @param robot_pose robot pose
 * @param landmark   landmark
 * @param H_lm       output: 2×2 Jacobian (row-major)
 * @return SLAM_OK
 */
int slam_obs_jacobian_rb_landmark(const slam_pose2d_t *robot_pose,
                                   const slam_landmark2d_t *landmark,
                                   double H_lm[4]);

/**
 * @brief Inverse observation model: measurement → landmark position
 *
 * Given robot pose x and observation z = [r, φ]^T, compute landmark
 * position in world frame:
 *
 *   m_x = x + r·cos(φ + θ)
 *   m_y = y + r·sin(φ + θ)
 *
 * @param robot_pose robot pose
 * @param obs        range-bearing observation
 * @param landmark   output: initialized landmark position
 * @return SLAM_OK
 */
int slam_obs_inverse_rb(const slam_pose2d_t *robot_pose,
                         const slam_obs_rb_t *obs,
                         slam_landmark2d_t *landmark);

/**
 * @brief Jacobian of inverse observation model
 *
 * G_rl = ∂g/∂(x,y,θ):   G_z = ∂g/∂(r,φ):
 *   [1 0 −r·sin(φ+θ);      [cos(φ+θ)  −r·sin(φ+θ);
 *    0 1  r·cos(φ+θ)]       sin(φ+θ)   r·cos(φ+θ) ]
 *
 * Used for EKF state augmentation.
 *
 * @param pose  robot pose
 * @param obs   observation
 * @param G_rl  output: 2×3 Jacobian w.r.t. pose
 * @param G_z   output: 2×2 Jacobian w.r.t. measurement
 * @return SLAM_OK
 */
int slam_obs_inverse_jacobian_rb(const slam_pose2d_t *pose,
                                  const slam_obs_rb_t *obs,
                                  double G_rl[6],
                                  double G_z[4]);

/* =========================================================================
 * L6: LiDAR Beam Model
 * ========================================================================= */

/**
 * @brief LiDAR beam likelihood field model
 *
 * Instead of ray-casting (expensive), the likelihood field precomputes
 * the distance to the nearest obstacle for each cell. The measurement
 * probability is:
 *
 *   p(z_t^k | x_t, m) = ε·p_hit + (1−ε)·p_rand
 *
 * where:
 *   p_hit ≈ N(dist(r, θ, x_t, m); 0, σ²)
 *   dist = distance from beam endpoint to nearest occupied cell
 *
 * This is smoother and more efficient than the full beam model.
 *
 * @param scan       LiDAR scan
 * @param robot_pose proposed robot pose
 * @param grid       occupancy grid
 * @param sigma_hit  std dev of hit likelihood [m]
 * @param p_rand     uniform random measurement probability
 * @param log_likelihood output: log likelihood of scan given pose
 * @return SLAM_OK
 */
int slam_lidar_likelihood_field(const slam_lidar_scan_t *scan,
                                 const slam_pose2d_t *robot_pose,
                                 const slam_occgrid_t *grid,
                                 double sigma_hit,
                                 double p_rand,
                                 double *log_likelihood);

/**
 * @brief Occupancy grid update using inverse sensor model
 *
 * For a LiDAR beam (robot_pose, angle_k, range reading r_k):
 *   - Cells along the beam (before r_k): free
 *   - Cell at the beam endpoint (r_k): occupied
 *   - Cells beyond r_k: unknown (no update)
 *
 * Log-odds update:
 *   l_{t}(cell) = l_{t-1}(cell) + lo_occ  (if occupied)
 *   l_{t}(cell) = l_{t-1}(cell) + lo_free (if free)
 *   clamped to [lo_min, lo_max]
 *
 * Uses Bresenham's line algorithm for ray casting.
 *
 * @param grid       occupancy grid (updated in-place)
 * @param robot_pose sensor pose in world frame
 * @param scan       LiDAR scan
 * @return SLAM_OK
 */
int slam_occgrid_update_scan(slam_occgrid_t *grid,
                              const slam_pose2d_t *robot_pose,
                              const slam_lidar_scan_t *scan);

/**
 * @brief Initialize occupancy grid
 *
 * Sets all cells to prior log-odds = 0 (p=0.5, unknown).
 *
 * @param grid        grid to initialize
 * @param width       cells in x
 * @param height      cells in y
 * @param resolution  meters per cell
 * @param origin_x    world x of grid origin
 * @param origin_y    world y of grid origin
 * @return SLAM_OK
 */
int slam_occgrid_init(slam_occgrid_t *grid,
                      int width, int height,
                      double resolution,
                      double origin_x, double origin_y);

/**
 * @brief Free occupancy grid memory
 */
void slam_occgrid_free(slam_occgrid_t *grid);

#ifdef __cplusplus
}
#endif

#endif /* SLAM_SENSOR_H */
