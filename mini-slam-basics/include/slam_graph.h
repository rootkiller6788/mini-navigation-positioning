#ifndef SLAM_GRAPH_H
#define SLAM_GRAPH_H

/**
 * @file    slam_graph.h
 * @brief   Graph-based SLAM: Pose Graph Optimization
 *
 * Graph SLAM formulates the SLAM posterior as a nonlinear least-squares
 * problem on a factor graph. Vertices represent robot poses (and optionally
 * landmarks), edges represent spatial constraints between vertices.
 *
 * The optimization problem:
 *   X* = argmin_X Σ_{(i,j)∈E} e_{ij}(X)^T · Ω_{ij} · e_{ij}(X)
 *
 * where e_{ij} = z_{ij} − h(x_i, x_j) is the constraint error,
 * Ω_{ij} is the information matrix (inverse covariance), and
 * h is the measurement function mapping poses to relative constraints.
 *
 * Reference:
 *   Lu & Milios (1997) "Globally Consistent Range Scan Alignment"
 *   Gutmann & Konolige (1999) "Incremental Mapping of Large Cyclic Environments"
 *   Dellaert & Kaess (2006) "Square Root SAM", IJRR.
 *   Kümmerle et al. (2011) "g2o: A General Framework for Graph Optimization"
 *
 * Knowledge Coverage:
 *   L3: Information form, Gauss-Newton, sparse Cholesky
 *   L5: Gauss-Newton optimization, Levenberg-Marquardt
 *   L6: Pose graph optimization, loop closure
 *   L8: Sparse linear solvers, robust kernels
 */

#include "slam_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Graph SLAM System
 * ========================================================================= */

/**
 * @brief Initialize an empty pose graph
 *
 * @param graph    allocated pose graph structure
 * @param max_nodes initial capacity
 * @param max_edges initial capacity
 * @return SLAM_OK
 */
int slam_graph_init(slam_pose_graph_t *graph,
                    int max_nodes,
                    int max_edges);

/**
 * @brief Free pose graph memory
 */
void slam_graph_free(slam_pose_graph_t *graph);

/**
 * @brief Add a pose node to the graph
 *
 * Each node is a robot pose. The first node is typically fixed
 * as the origin to resolve gauge freedom.
 *
 * @param graph  pose graph
 * @param pose   node pose
 * @param fixed  whether node is anchored (not optimized)
 * @param node_id output: assigned node id
 * @return SLAM_OK
 */
int slam_graph_add_node(slam_pose_graph_t *graph,
                        const slam_pose2d_t *pose,
                        bool fixed,
                        int *node_id);

/**
 * @brief Add an edge (constraint) between two nodes
 *
 * The constraint is the measured relative pose from node_a to node_b.
 * The information matrix Ω encodes the confidence in this constraint.
 *
 * For odometry edges (sequential):
 *   constraint = odometry measurement
 *   Ω = diag(1/σ_x², 1/σ_y², 1/σ_θ²)  [diagonal approximation]
 *
 * For loop closure edges:
 *   constraint = scan matching result
 *   Ω = full 3×3 from ICP covariance estimation
 *
 * @param graph   pose graph
 * @param id_a    source node id
 * @param id_b    target node id
 * @param constraint relative pose measurement
 * @param info    3×3 information matrix (row-major)
 * @param edge_id output: assigned edge id
 * @return SLAM_OK
 */
int slam_graph_add_edge(slam_pose_graph_t *graph,
                        int id_a, int id_b,
                        const slam_pose2d_t *constraint,
                        const double info[9],
                        int *edge_id);

/* =========================================================================
 * L5: Error Functions and Jacobians for SE(2) Pose Graph
 * ========================================================================= */

/**
 * @brief Compute SE(2) relative pose error (2D)
 *
 * Given estimated poses x_i, x_j and measured relative constraint z_ij,
 * compute the error vector e = z_ij ⊖ (x_i^{-1} ∘ x_j).
 *
 * In SE(2) parameterization:
 *   Δx = cos(θ_i)·(x_j−x_i) + sin(θ_i)·(y_j−y_i)
 *   Δy = −sin(θ_i)·(x_j−x_i) + cos(θ_i)·(y_j−y_i)
 *   Δθ = normalize(θ_j − θ_i)
 *
 *   e_x = z_x − Δx
 *   e_y = z_y − Δy
 *   e_θ = normalize(z_θ − Δθ)
 *
 * @param xi     pose of node i
 * @param xj     pose of node j
 * @param z_ij   measured relative constraint
 * @param error  output: 3×1 error vector [e_x, e_y, e_θ]^T
 * @return SLAM_OK
 */
int slam_graph_error_se2(const slam_pose2d_t *xi,
                          const slam_pose2d_t *xj,
                          const slam_pose2d_t *z_ij,
                          double error[3]);

/**
 * @brief Compute Jacobians of SE(2) error w.r.t. x_i and x_j
 *
 * The error function e(x_i, x_j) = z_ij ⊖ (x_i^{-1} ∘ x_j).
 * Jacobians J_i = ∂e/∂x_i, J_j = ∂e/∂x_j, each 3×3.
 *
 * J_i = [ −R_i^T,    R_i^T·∂R(Δθ)/∂θ_i·Δp;
 *          0 0,       −1                          ]
 *
 * J_j = [  R_i^T,     0;
 *          0 0,        1                          ]
 *
 * where R_i is the 2×2 rotation of pose i,
 * Δp = [x_j−x_i, y_j−y_i]^T.
 *
 * @param xi   pose of node i
 * @param xj   pose of node j
 * @param Ji   output: 3×3 Jacobian w.r.t. x_i (row-major)
 * @param Jj   output: 3×3 Jacobian w.r.t. x_j (row-major)
 * @return SLAM_OK
 */
int slam_graph_jacobian_se2(const slam_pose2d_t *xi,
                             const slam_pose2d_t *xj,
                             double Ji[9],
                             double Jj[9]);

/* =========================================================================
 * L5: Graph Optimization
 * ========================================================================= */

/**
 * @brief Build the linear system H·Δx = −b from the pose graph
 *
 * For each edge (i,j) with error e, Jacobians J_i, J_j, and
 * information matrix Ω:
 *
 *   H_ii += J_i^T·Ω·J_i
 *   H_ij += J_i^T·Ω·J_j
 *   H_ji += J_j^T·Ω·J_i
 *   H_jj += J_j^T·Ω·J_j
 *   b_i  += J_i^T·Ω·e
 *   b_j  += J_j^T·Ω·e
 *
 * chi² = Σ e^T·Ω·e
 *
 * @param graph  pose graph (hessian/gradient populated)
 * @param chi2   output: total chi-squared error
 * @return SLAM_OK
 */
int slam_graph_build_system(slam_pose_graph_t *graph,
                             double *chi2);

/**
 * @brief Solve H·Δx = −b using Cholesky decomposition
 *
 * For small-to-medium graphs (N < ~1000), dense Cholesky is acceptable:
 *   1. H = L·L^T (Cholesky decomposition)
 *   2. L·y = −b (forward substitution)
 *   3. L^T·Δx = y (back substitution)
 *
 * For larger graphs, sparse Cholesky (CSparse/SuiteSparse) is preferred,
 * but this implementation uses a simple dense solver for clarity.
 *
 * @param graph   pose graph with H and b populated
 * @param dx      output: solution vector, length 3*num_nodes
 * @return SLAM_OK, or SLAM_ERR_SINGULAR if H is not positive definite
 */
int slam_graph_solve_cholesky(slam_pose_graph_t *graph,
                               double *dx);

/**
 * @brief Gauss-Newton optimization of the pose graph
 *
 * Algorithm:
 *   while not converged and iter < max_iter:
 *     1. Build linear system H·Δx = −b
 *     2. Solve for Δx
 *     3. Update: X ← X ⊕ Δx
 *     4. Check convergence: |Δchi²| / chi² < threshold
 *
 * Convergence is typically achieved in 5-20 iterations for pose graphs.
 *
 * @param graph    pose graph (nodes updated in-place)
 * @param config   SLAM config (max_iter, threshold)
 * @param metrics  output: chi2 history and convergence info
 * @return SLAM_OK
 */
int slam_graph_optimize_gauss_newton(slam_pose_graph_t *graph,
                                      const slam_config_t *config,
                                      slam_metrics_t *metrics);

/**
 * @brief Levenberg-Marquardt optimization
 *
 * LM adds damping to GN for better convergence from poor initial guesses:
 *   (H + λ·diag(H))·Δx = −b
 *
 * λ adjustment (Nielsen's strategy):
 *   If chi²_new < chi²:    λ ← λ·max(1/3, 1−(2ρ−1)³)   (accept step)
 *   If chi²_new ≥ chi²:    λ ← λ·ν, ν ← 2·ν             (reject, increase damping)
 *
 * where ρ = (chi² − chi²_new) / (Δx^T·(λ·diag(H)·Δx − b))
 * is the gain ratio (predicted vs actual improvement).
 *
 * @param graph    pose graph
 * @param config   config with lm_lambda_init and convergence params
 * @param metrics  output metrics
 * @return SLAM_OK
 */
int slam_graph_optimize_lm(slam_pose_graph_t *graph,
                            const slam_config_t *config,
                            slam_metrics_t *metrics);

/* =========================================================================
 * L8: Robust Kernels
 * ========================================================================= */

/**
 * @brief Apply Huber robust kernel to an edge error
 *
 * Huber loss reduces the influence of outlier constraints (e.g., false
 * loop closures):
 *
 *   ρ(e) = { ½·e²,               if |e| ≤ δ
 *          { δ·(|e| − δ/2),     if |e| > δ
 *
 * Weight for information matrix: w = { 1,      if |e| ≤ δ
 *                                      { δ/|e|, if |e| > δ
 *
 * The edge information is scaled by w: Ω_robust = w·Ω.
 *
 * @param error     3×1 error vector
 * @param info      3×3 information matrix (modified in-place)
 * @param delta     Huber kernel threshold
 * @param weight    output: computed weight
 * @return SLAM_OK
 */
int slam_graph_huber_kernel(const double error[3],
                             double info[9],
                             double delta,
                             double *weight);

#ifdef __cplusplus
}
#endif

#endif /* SLAM_GRAPH_H */
