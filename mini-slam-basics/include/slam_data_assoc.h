#ifndef SLAM_DATA_ASSOC_H
#define SLAM_DATA_ASSOC_H

/**
 * @file    slam_data_assoc.h
 * @brief   Data Association Algorithms for SLAM
 *
 * Data association is the problem of determining which landmark (if any)
 * generated a given sensor observation. This is a critical and
 * error-prone step in SLAM — wrong associations can cause filter
 * divergence.
 *
 * Key algorithms:
 *   - Nearest Neighbor (NN): pick the closest prediction
 *   - Mahalanobis gating: statistical compatibility test
 *   - Joint Compatibility Branch and Bound (JCBB): joint compatibility
 *   - ICP: Iterative Closest Point for scan matching
 *
 * Reference:
 *   Bar-Shalom & Fortmann (1988) "Tracking and Data Association"
 *   Neira & Tardos (2001) "Data Association in Stochastic Mapping",
 *     IEEE Trans. Robotics & Automation.
 *   Cox (1993) "A Review of Statistical Data Association Techniques"
 *
 * Knowledge Coverage:
 *   L3: Mahalanobis distance, chi-squared gating
 *   L5: Nearest Neighbor, JCBB, ICP
 *   L6: Feature matching for loop closure
 */

#include "slam_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L3: Mahalanobis Distance
 * ========================================================================= */

/**
 * @brief Compute Mahalanobis distance between observation and prediction
 *
 * d² = (z − ẑ)^T · S^{-1} · (z − ẑ)
 *
 * where ẑ = h(μ, m_j) is the predicted measurement and S is the
 * innovation covariance: S = H·Σ·H^T + Q
 *
 * Mahalanobis distance accounts for uncertainty correlation —
 * unlike Euclidean distance, it is scale-invariant and follows
 * a χ² distribution under Gaussian assumptions.
 *
 * Under the null hypothesis (correct association):
 *   d² ∼ χ²_{dim(z)}  (e.g., χ²₂ for range-bearing)
 *
 * Gating: accept association if d² ≤ χ²_{dim(z), α}
 *   α=0.95 → threshold=5.991 for 2 DOF
 *   α=0.99 → threshold=9.210 for 2 DOF
 *
 * @param innovation  measurement residual ν = z − ẑ (length dim_z)
 * @param S           innovation covariance (dim_z × dim_z, row-major)
 * @param dim_z       measurement dimension
 * @return Mahalanobis distance squared (≥ 0), or −1 on error
 */
double slam_mahalanobis_sq(const double *innovation,
                            const double *S,
                            int dim_z);

/**
 * @brief Compute innovation covariance for range-bearing
 *
 * S = H·Σrr·H^T + Q
 * where H is the 2×3 observation Jacobian w.r.t. robot pose,
 * Σrr is the 3×3 robot pose covariance submatrix,
 * Q = diag(σ_r², σ_φ²) is the measurement noise.
 *
 * @param H_pose   2×3 Jacobian (row-major)
 * @param Sigma_rr 3×3 robot pose covariance (row-major)
 * @param Q        2×2 measurement noise (row-major, usually diagonal)
 * @param S         output: 2×2 innovation covariance
 */
void slam_innovation_covariance_rb(const double H_pose[6],
                                    const double Sigma_rr[9],
                                    const double Q[4],
                                    double S[4]);

/* =========================================================================
 * L5: Nearest Neighbor (NN) Data Association
 * ========================================================================= */

/**
 * @brief Nearest Neighbor data association
 *
 * For a given observation, compute the Mahalanobis distance to all
 * landmarks. Accept the landmark with smallest distance if it passes
 * the chi-squared gate.
 *
 * Algorithm:
 *   d_min = ∞, best_j = −1
 *   for j = 0..N−1:
 *     ẑ_j = h(μ, m_j)    // predict measurement
 *     S_j = H·Σ·H^T + Q  // innovation covariance
 *     d² = (z − ẑ_j)^T·S_j^{-1}·(z − ẑ_j)
 *     if d² < γ and d² < d_min:
 *       d_min = d², best_j = j
 *   return best_j
 *
 * Complexity: O(N) per observation, where N = number of landmarks.
 *
 * @param robot_pose  current robot pose
 * @param Sigma_rr    robot pose covariance 3×3
 * @param obs         observation to associate
 * @param landmarks   array of candidate landmarks
 * @param num_landmarks N
 * @param gate_thresh chi-squared gate threshold (e.g., 5.991)
 * @param best_idx    output: index of best match (−1 if none)
 * @param dist_sq     output: Mahalanobis distance of best match
 * @return SLAM_OK
 */
int slam_da_nearest_neighbor(const slam_pose2d_t *robot_pose,
                              const double Sigma_rr[9],
                              const slam_obs_rb_t *obs,
                              const slam_landmark2d_t *landmarks,
                              int num_landmarks,
                              double gate_thresh,
                              int *best_idx,
                              double *dist_sq);

/* =========================================================================
 * L5: Joint Compatibility Branch and Bound (JCBB)
 * ========================================================================= */

/**
 * @brief Joint Compatibility test for a set of associations
 *
 * Joint Mahalanobis distance for K simultaneous observations:
 *   D² = ν_total^T · S_total^{-1} · ν_total
 *
 * where ν_total is the stacked innovation vector (2K × 1) and
 * S_total is the joint innovation covariance (2K × 2K).
 *
 * Individual compatibility checks can miss correlations between
 * associations — JCBB finds the largest jointly compatible set.
 *
 * @param innovations  stacked innovations [ν_0; ν_1; ...], length 2*K
 * @param S_joint      joint innovation covariance, (2K)×(2K)
 * @param K            number of simultaneous observations
 * @param gate_thresh  chi-squared threshold for 2K DOF
 * @return 1 if jointly compatible, 0 otherwise
 */
int slam_da_joint_compatible(const double *innovations,
                              const double *S_joint,
                              int K,
                              double gate_thresh);

/**
 * @brief JCBB: find the best set of jointly compatible associations
 *
 * Branch-and-bound search over the interpretation tree.
 * Each level corresponds to one observation, branches are
 * landmark associations (including "new landmark" as a branch).
 *
 * Algorithm (Neira & Tardos 2001):
 *   1. Order observations by individual compatibility
 *   2. Depth-first search over interpretation tree
 *   3. Prune branches that violate joint compatibility
 *   4. Bound: maximum possible remaining compatible pairings
 *
 * Complexity: worst-case O(N^K), but pruning makes it practical
 * for moderate N, K.
 *
 * @param robot_pose   robot pose
 * @param Sigma_rr     robot pose covariance
 * @param observations array of K observations
 * @param K            number of observations
 * @param landmarks    candidate landmarks
 * @param N            number of landmarks
 * @param gate_thresh  chi-squared threshold
 * @param associations output: association[obs_i] = lm_j (or −1 for new)
 * @param num_paired   output: number of associations found
 * @return SLAM_OK
 */
int slam_da_jcbb(const slam_pose2d_t *robot_pose,
                  const double Sigma_rr[9],
                  const slam_obs_rb_t *observations,
                  int K,
                  const slam_landmark2d_t *landmarks,
                  int N,
                  double gate_thresh,
                  int *associations,
                  int *num_paired);

/* =========================================================================
 * L5: Iterative Closest Point (ICP) for Scan Matching
 * ========================================================================= */

/**
 * @brief ICP scan matching: align two 2D point clouds
 *
 * ICP iteratively refines the transformation T = (dx, dy, dθ) that
 * best aligns source points to target points.
 *
 * Algorithm:
 *   1. For each source point, find nearest target point (NN search)
 *   2. Compute transformation that minimizes point-to-point error:
 *      E = Σ ‖R·s_i + t − m_i‖²
 *   3. Apply transformation to source
 *   4. Repeat until convergence
 *
 * Closed-form solution for step 2 (Arun et al. 1987, using SVD):
 *   H = Σ (s_i − s̄)·(m_i − m̄)^T
 *   H = U·D·V^T  (SVD)
 *   R = V·U^T
 *   t = m̄ − R·s̄
 *
 * @param src_pts    source point cloud (N×2, row-major)
 * @param tgt_pts    target point cloud (N×2)
 * @param N          number of points (must match)
 * @param max_iter   maximum ICP iterations
 * @param eps        convergence threshold on ||ΔT||
 * @param T          output: transformation [dx, dy, dθ]
 * @param rmse       output: final RMSE
 * @return SLAM_OK
 */
int slam_da_icp_2d(const double *src_pts,
                    const double *tgt_pts,
                    int N,
                    int max_iter,
                    double eps,
                    slam_pose2d_t *T,
                    double *rmse);

/* =========================================================================
 * L6: Loop Closure Detection
 * ========================================================================= */

/**
 * @brief Detect loop closure candidates by pose proximity
 *
 * A loop closure occurs when the robot revisits a previously mapped
 * area. This function finds past poses that are within a radius of
 * the current pose but separated by enough time steps.
 *
 * Criteria:
 *   - Distance between current and past pose < search_radius
 *   - Time step difference > min_steps (avoid detecting adjacent poses)
 *
 * @param poses          trajectory poses
 * @param num_poses      number of poses in trajectory
 * @param current_pose   current robot pose
 * @param search_radius  spatial search radius [m]
 * @param min_steps      minimum step difference
 * @param match_idx      output: index of best loop candidate (−1 if none)
 * @return SLAM_OK
 */
int slam_da_detect_loop(const slam_pose2d_t *poses,
                         int num_poses,
                         const slam_pose2d_t *current_pose,
                         double search_radius,
                         int min_steps,
                         int *match_idx);

#ifdef __cplusplus
}
#endif

#endif /* SLAM_DATA_ASSOC_H */
