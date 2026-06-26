/**
 * @file    slam_graph.c
 * @brief   Graph-based SLAM: Pose Graph Optimization
 *
 * Implements Gauss-Newton and Levenberg-Marquardt optimization for
 * pose graphs in SE(2). The pose graph formulation is:
 *
 *   X* = argmin_X Σ_{ij} e_{ij}(X)^T · Ω_{ij} · e_{ij}(X)
 *
 * where e_{ij} = z_{ij} ⊖ (x_i^{-1} ∘ x_j) is the SE(2) constraint error.
 *
 * Solves H·Δx = −b using dense Cholesky decomposition.
 * Supports Huber robust kernel for outlier loop closures.
 *
 * Reference:
 *   Lu & Milios (1997)
 *   Dellaert & Kaess (2006) "Square Root SAM", IJRR
 *   Kümmerle et al. (2011) "g2o: A General Framework for Graph Optimization"
 */

#include "slam_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

/* External declarations */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern void slam_pose_compose(const slam_pose2d_t *a,
                               const slam_pose2d_t *b,
                               slam_pose2d_t *c);
extern void slam_pose_inverse(const slam_pose2d_t *a, slam_pose2d_t *inv);
extern int slam_cholesky_3x3(const double A[9], double L[9]);
extern void slam_forward_sub_3x3(const double L[9], const double b[3],
                                  double y[3]);
extern void slam_back_sub_3x3(const double L[9], const double y[3],
                               double x[3]);
extern void slam_matvec_mul(const double *A, const double *x,
                             int m, int n, double *y);
extern void slam_eye(double *M, int n);
extern void slam_mat_copy(const double *src, double *dst, int n);
extern int slam_solve_cholesky_3x3(const double A[9], const double b[3],
                                    double x[3]);

/* =========================================================================
 * L1: Graph Management
 * ========================================================================= */

int slam_graph_init(slam_pose_graph_t *graph,
                    int max_nodes, int max_edges) {
    if (!graph || max_nodes <= 0 || max_edges <= 0)
        return SLAM_ERR_INVALID_PARAM;

    graph->num_nodes     = 0;
    graph->node_capacity = max_nodes;
    graph->num_edges     = 0;
    graph->edge_capacity = max_edges;

    graph->nodes = (slam_graph_node_t *)calloc(max_nodes,
                                                sizeof(slam_graph_node_t));
    graph->edges = (slam_graph_edge_t *)calloc(max_edges,
                                                sizeof(slam_graph_edge_t));
    graph->hessian  = (double *)calloc(3 * max_nodes * 3 * max_nodes,
                                        sizeof(double));
    graph->gradient = (double *)calloc(3 * max_nodes, sizeof(double));
    graph->chi2_history = (double *)calloc(100, sizeof(double));
    graph->chi2_len = 0;

    if (!graph->nodes || !graph->edges || !graph->hessian
        || !graph->gradient || !graph->chi2_history) {
        slam_graph_free(graph);
        return SLAM_ERR_MEMORY;
    }

    return SLAM_OK;
}

void slam_graph_free(slam_pose_graph_t *graph) {
    if (graph) {
        free(graph->nodes);
        free(graph->edges);
        free(graph->hessian);
        free(graph->gradient);
        free(graph->chi2_history);
        memset(graph, 0, sizeof(slam_pose_graph_t));
    }
}

int slam_graph_add_node(slam_pose_graph_t *graph,
                        const slam_pose2d_t *pose,
                        bool fixed, int *node_id) {
    if (!graph || !pose || !node_id) return SLAM_ERR_NULL_PTR;
    if (graph->num_nodes >= graph->node_capacity) return SLAM_ERR_MEMORY;

    int id = graph->num_nodes;
    graph->nodes[id].id    = id;
    graph->nodes[id].pose  = *pose;
    graph->nodes[id].is_fixed = fixed;
    graph->num_nodes++;
    *node_id = id;
    return SLAM_OK;
}

int slam_graph_add_edge(slam_pose_graph_t *graph,
                        int id_a, int id_b,
                        const slam_pose2d_t *constraint,
                        const double info[9], int *edge_id) {
    if (!graph || !constraint || !info || !edge_id)
        return SLAM_ERR_NULL_PTR;
    if (id_a < 0 || id_a >= graph->num_nodes) return SLAM_ERR_INVALID_PARAM;
    if (id_b < 0 || id_b >= graph->num_nodes) return SLAM_ERR_INVALID_PARAM;
    if (graph->num_edges >= graph->edge_capacity) return SLAM_ERR_MEMORY;

    int eid = graph->num_edges;
    graph->edges[eid].id_a       = id_a;
    graph->edges[eid].id_b       = id_b;
    graph->edges[eid].constraint = *constraint;
    memcpy(graph->edges[eid].info_matrix, info, 9 * sizeof(double));
    graph->edges[eid].is_loop_closure = false;
    graph->edges[eid].inlier_count = 0;
    graph->num_edges++;
    *edge_id = eid;
    return SLAM_OK;
}

/* =========================================================================
 * L5: SE(2) Error and Jacobian
 * ========================================================================= */

int slam_graph_error_se2(const slam_pose2d_t *xi,
                          const slam_pose2d_t *xj,
                          const slam_pose2d_t *z_ij,
                          double error[3]) {
    if (!xi || !xj || !z_ij || !error) return SLAM_ERR_NULL_PTR;

    /* Δ = inv(xi) ⊕ xj  (relative pose in xi's frame) */
    double dx = xj->x - xi->x;
    double dy = xj->y - xi->y;
    double ci = cos(xi->theta), si = sin(xi->theta);

    double delta_x  =  dx * ci + dy * si;
    double delta_y  = -dx * si + dy * ci;
    double delta_th = slam_normalize_angle(xj->theta - xi->theta);

    /* error = z_ij ⊖ delta */
    error[0] = z_ij->x - delta_x;
    error[1] = z_ij->y - delta_y;
    error[2] = slam_normalize_angle(z_ij->theta - delta_th);

    return SLAM_OK;
}

int slam_graph_jacobian_se2(const slam_pose2d_t *xi,
                             const slam_pose2d_t *xj,
                             double Ji[9], double Jj[9]) {
    if (!xi || !xj || !Ji || !Jj) return SLAM_ERR_NULL_PTR;

    double dx = xj->x - xi->x;
    double dy = xj->y - xi->y;
    double ci = cos(xi->theta), si = sin(xi->theta);

    /* Rotation matrix of xi: R_i = [ci, -si; si, ci] */
    /* R_i^T = [ci, si; -si, ci] */

    /* J_i = [-R_i^T,  R_i^T·∂R(Δθ)/∂θ_i·Δp;
     *          0 0,  -1                          ]
     * where ∂R(Δθ)/∂θ_i ≈ [−sin(Δθ), cos(Δθ); −cos(Δθ), −sin(Δθ)]
     *
     * Simplified: J_i(1:2, 1:2) = -R_i^T
     * J_i(1:2, 3) = [−Δy; Δx]   (in xi's frame after rotation, ≈ [dy; -dx])
     * with some detail...
     *
     * Standard g2o-style Jacobians for SE(2):
     * Let (x_j - x_i)_i = R_i^T · (p_j - p_i).
     * Then:
     * J_i = [-R_i^T,  [delta_y; -delta_x];
     *         0 0  ,  -1                ]
     *
     * J_j = [ R_i^T,  0 0;
     *          0 0 ,   1  ]
     */

    double delta_x =  dx * ci + dy * si;   /* cos component */
    double delta_y = -dx * si + dy * ci;   /* sin component */

    /* J_i */
    Ji[0] = -ci;  Ji[1] = -si;  Ji[2] =  delta_y;
    Ji[3] =  si;  Ji[4] = -ci;  Ji[5] = -delta_x;
    Ji[6] =  0;   Ji[7] =  0;   Ji[8] = -1.0;

    /* J_j */
    Jj[0] =  ci;  Jj[1] =  si;  Jj[2] = 0;
    Jj[3] = -si;  Jj[4] =  ci;  Jj[5] = 0;
    Jj[6] =  0;   Jj[7] =  0;   Jj[8] = 1.0;

    return SLAM_OK;
}

/* =========================================================================
 * L5: Build Linear System
 * ========================================================================= */

int slam_graph_build_system(slam_pose_graph_t *graph, double *chi2) {
    if (!graph || !chi2) return SLAM_ERR_NULL_PTR;

    int N = graph->num_nodes;
    int dim = 3 * N;

    /* Zero out H and b */
    memset(graph->hessian,  0, dim * dim * sizeof(double));
    memset(graph->gradient, 0, dim * sizeof(double));
    *chi2 = 0.0;

    for (int e = 0; e < graph->num_edges; e++) {
        int ia = graph->edges[e].id_a;
        int ib = graph->edges[e].id_b;
        const double *Omega = graph->edges[e].info_matrix;
        const slam_pose2d_t *z = &graph->edges[e].constraint;

        slam_pose2d_t *xa = &graph->nodes[ia].pose;
        slam_pose2d_t *xb = &graph->nodes[ib].pose;

        /* Compute error */
        double error[3];
        slam_graph_error_se2(xa, xb, z, error);

        /* Compute Jacobians */
        double Ja[9], Jb[9];
        slam_graph_jacobian_se2(xa, xb, Ja, Jb);

        /* chi2 += error^T · Omega · error */
        double Oe[3];
        slam_matvec_mul(Omega, error, 3, 3, Oe);
        *chi2 += error[0]*Oe[0] + error[1]*Oe[1] + error[2]*Oe[2];

        /* b_a += J_a^T · Omega · error */
        /* b_b += J_b^T · Omega · error */
        for (int i = 0; i < 3; i++) {
            int bi = 3*ia + i;
            for (int j = 0; j < 3; j++) {
                graph->gradient[bi] += Ja[j*3 + i] * Oe[j];
            }
        }
        for (int i = 0; i < 3; i++) {
            int bi = 3*ib + i;
            for (int j = 0; j < 3; j++) {
                graph->gradient[bi] += Jb[j*3 + i] * Oe[j];
            }
        }

        /* H_aa += J_a^T · Omega · J_a */
        /* H_ab += J_a^T · Omega · J_b */
        /* H_ba += J_b^T · Omega · J_a  (symmetric with H_ab) */
        /* H_bb += J_b^T · Omega · J_b */
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                /* J_a^T · Omega · J_a, Omega diagonal simplified */

                /* H_aa */
                double haa = 0.0;
                for (int k = 0; k < 3; k++)
                    haa += Ja[k*3 + r] * Omega[k*3 + k] * Ja[k*3 + c];
                graph->hessian[(3*ia+r) * dim + (3*ia+c)] += haa;

                /* H_ab */
                double hab = 0.0;
                for (int k = 0; k < 3; k++)
                    hab += Ja[k*3 + r] * Omega[k*3 + k] * Jb[k*3 + c];
                graph->hessian[(3*ia+r) * dim + (3*ib+c)] += hab;

                /* H_ba */
                double hba = 0.0;
                for (int k = 0; k < 3; k++)
                    hba += Jb[k*3 + r] * Omega[k*3 + k] * Ja[k*3 + c];
                graph->hessian[(3*ib+r) * dim + (3*ia+c)] += hba;

                /* H_bb */
                double hbb = 0.0;
                for (int k = 0; k < 3; k++)
                    hbb += Jb[k*3 + r] * Omega[k*3 + k] * Jb[k*3 + c];
                graph->hessian[(3*ib+r) * dim + (3*ib+c)] += hbb;
            }
        }
    }

    /* Fix anchor nodes: set diagonal to large value */
    for (int n = 0; n < N; n++) {
        if (graph->nodes[n].is_fixed) {
            for (int i = 0; i < 3; i++) {
                int idx = 3*n + i;
                graph->hessian[idx * dim + idx] += 1e10;
                graph->gradient[idx] = 0.0;
            }
        }
    }

    return SLAM_OK;
}

/* =========================================================================
 * L5: Cholesky Solver for Pose Graph
 * ========================================================================= */

int slam_graph_solve_cholesky(slam_pose_graph_t *graph, double *dx) {
    if (!graph || !dx) return SLAM_ERR_NULL_PTR;

    int dim = 3 * graph->num_nodes;

    /* Solve H·dx = -b using dense Cholesky */
    /* H is SPD by construction (information matrices are SPD) */
    double *L = (double *)calloc(dim * dim, sizeof(double));
    if (!L) return SLAM_ERR_MEMORY;

    /* Cholesky: H = L·L^T */
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i*dim + k] * L[j*dim + k];
            }
            if (i == j) {
                double diag = graph->hessian[i*dim + i] - sum;
                if (diag <= 0) {
                    /* Add damping for numerical stability */
                    diag = 1e-12;
                }
                L[i*dim + i] = sqrt(diag);
            } else {
                L[i*dim + j] = (graph->hessian[i*dim + j] - sum) / L[j*dim + j];
            }
        }
    }

    /* Forward substitution: L·y = -b */
    double *y = (double *)calloc(dim, sizeof(double));
    if (!y) { free(L); return SLAM_ERR_MEMORY; }

    for (int i = 0; i < dim; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i*dim + j] * y[j];
        }
        y[i] = (-graph->gradient[i] - sum) / L[i*dim + i];
    }

    /* Back substitution: L^T·dx = y */
    for (int i = dim - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < dim; j++) {
            sum += L[j*dim + i] * dx[j];
        }
        dx[i] = (y[i] - sum) / L[i*dim + i];
    }

    free(L); free(y);
    return SLAM_OK;
}

/* =========================================================================
 * L5: Gauss-Newton Optimization
 * ========================================================================= */

int slam_graph_optimize_gauss_newton(slam_pose_graph_t *graph,
                                      const slam_config_t *config,
                                      slam_metrics_t *metrics) {
    if (!graph || !config) return SLAM_ERR_NULL_PTR;

    int N = graph->num_nodes;
    int dim = 3 * N;
    double prev_chi2 = DBL_MAX;
    graph->chi2_len = 0;

    double *dx = (double *)calloc(dim, sizeof(double));
    if (!dx) return SLAM_ERR_MEMORY;

    for (int iter = 0; iter < config->max_graph_iterations; iter++) {
        double chi2;
        slam_graph_build_system(graph, &chi2);

        /* Convergence check */
        if (graph->chi2_len < 100) {
            graph->chi2_history[graph->chi2_len++] = chi2;
        }

        if (fabs(prev_chi2 - chi2) / fmax(chi2, 1e-12)
            < config->convergence_thresh) {
            break;
        }
        prev_chi2 = chi2;

        /* Solve */
        slam_graph_solve_cholesky(graph, dx);

        /* Update poses: x ← x ⊕ Δx */
        for (int n = 0; n < N; n++) {
            if (graph->nodes[n].is_fixed) continue;
            double dxi = dx[3*n + 0];
            double dyi = dx[3*n + 1];
            double dti = dx[3*n + 2];
            /* compose change in local frame */
            double ci = cos(graph->nodes[n].pose.theta);
            double si = sin(graph->nodes[n].pose.theta);
            graph->nodes[n].pose.x += dxi * ci - dyi * si;
            graph->nodes[n].pose.y += dxi * si + dyi * ci;
            graph->nodes[n].pose.theta = slam_normalize_angle(
                graph->nodes[n].pose.theta + dti);
        }
    }

    if (metrics) {
        metrics->cpu_time_per_step = 0.0; /* not profiled */
    }

    free(dx);
    return SLAM_OK;
}

/* =========================================================================
 * L8: Levenberg-Marquardt Optimization
 * ========================================================================= */

int slam_graph_optimize_lm(slam_pose_graph_t *graph,
                            const slam_config_t *config,
                            slam_metrics_t *metrics) {
    (void)metrics;  /* reserved for future metric collection */
    if (!graph || !config) return SLAM_ERR_NULL_PTR;

    int N = graph->num_nodes;
    int dim = 3 * N;
    double lambda = config->lm_lambda_init;
    double nu = 2.0;
    double chi2 = 0.0;

    double *dx_backup = (double *)calloc(dim, sizeof(double));
    slam_pose2d_t *pose_backup = (slam_pose2d_t *)calloc(
        N, sizeof(slam_pose2d_t));
    double *H_augmented = (double *)calloc(dim * dim, sizeof(double));

    if (!dx_backup || !pose_backup || !H_augmented) {
        free(dx_backup); free(pose_backup); free(H_augmented);
        return SLAM_ERR_MEMORY;
    }

    graph->chi2_len = 0;

    for (int iter = 0; iter < config->max_graph_iterations; iter++) {
        /* Build system */
        slam_graph_build_system(graph, &chi2);

        if (graph->chi2_len < 100) {
            graph->chi2_history[graph->chi2_len++] = chi2;
        }

        /* Save current state */
        for (int n = 0; n < N; n++) pose_backup[n] = graph->nodes[n].pose;

        /* Augment Hessian: H_aug = H + λ·diag(H) */
        memcpy(H_augmented, graph->hessian, dim*dim * sizeof(double));
        for (int i = 0; i < dim; i++) {
            H_augmented[i*dim + i] += lambda * graph->hessian[i*dim + i];
        }

        /* Swap hessian temporarily */
        double *temp_hessian = graph->hessian;
        graph->hessian = H_augmented;

        /* Solve */
        double *dx = (double *)calloc(dim, sizeof(double));
        if (!dx) { graph->hessian = temp_hessian; goto cleanup; }
        slam_graph_solve_cholesky(graph, dx);
        graph->hessian = temp_hessian;

        /* Update */
        for (int n = 0; n < N; n++) {
            if (graph->nodes[n].is_fixed) continue;
            double ci = cos(graph->nodes[n].pose.theta);
            double si = sin(graph->nodes[n].pose.theta);
            graph->nodes[n].pose.x += dx[3*n]*ci - dx[3*n+1]*si;
            graph->nodes[n].pose.y += dx[3*n]*si + dx[3*n+1]*ci;
            graph->nodes[n].pose.theta = slam_normalize_angle(
                graph->nodes[n].pose.theta + dx[3*n+2]);
        }

        /* Evaluate new chi2 */
        double chi2_new;
        slam_graph_build_system(graph, &chi2_new);

        /* Gain ratio */
        double num = chi2 - chi2_new;
        double den = 0.0;
        for (int i = 0; i < dim; i++) {
            den += dx[i] * (lambda * graph->hessian[i*dim + i] * dx[i]
                            - graph->gradient[i]);
        }
        double rho = (den > 1e-12) ? num / den : 0.0;

        if (rho > 0) {
            /* Accept step */
            lambda = lambda * fmax(1.0/3.0, 1.0 - pow(2.0*rho - 1.0, 3));
            nu = 2.0;
        } else {
            /* Reject step, restore poses */
            for (int n = 0; n < N; n++)
                graph->nodes[n].pose = pose_backup[n];
            lambda *= nu;
            nu *= 2.0;
        }

        free(dx);

        if (fabs(chi2 - chi2_new) / fmax(chi2, 1e-12)
            < config->convergence_thresh && rho > 0) {
            break;
        }
        if (rho > 0) chi2 = chi2_new;
    }

    cleanup:
    free(dx_backup); free(pose_backup); free(H_augmented);
    return SLAM_OK;
}

/* =========================================================================
 * L8: Huber Robust Kernel
 * ========================================================================= */

int slam_graph_huber_kernel(const double error[3],
                             double info[9],
                             double delta,
                             double *weight) {
    if (!error || !info || !weight) return SLAM_ERR_NULL_PTR;

    /* Compute weighted error norm: ||e||_Ω² = e^T·Ω·e */
    double sq_error = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            sq_error += error[i] * info[i*3 + j] * error[j];

    double abs_error = sqrt(fmax(sq_error, 0.0));

    /* Huber weight */
    if (abs_error <= delta) {
        *weight = 1.0;
    } else {
        *weight = delta / abs_error;
    }

    /* Scale information matrix by weight */
    /* Ω_robust = w·Ω  (for the first derivative; the second
     * derivative has an additional Jacobian-of-weight term) */
    for (int i = 0; i < 9; i++) {
        info[i] *= (*weight);
    }

    return SLAM_OK;
}
