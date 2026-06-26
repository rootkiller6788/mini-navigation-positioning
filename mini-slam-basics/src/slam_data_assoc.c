/**
 * @file    slam_data_assoc.c
 * @brief   Data Association: NN, JCBB, ICP scan matching, loop detection
 *
 * Data association is the SLAM sub-problem of matching observations
 * to map landmarks. This is critical because:
 *   - Wrong associations produce inconsistent updates
 *   - The filter can diverge catastrophically
 *   - Computationally, it's a combinatorial optimization problem
 *
 * Reference:
 *   Neira & Tardos (2001) "Data Association in Stochastic Mapping"
 *   Bar-Shalom & Fortmann (1988) "Tracking and Data Association"
 *   Besl & McKay (1992) "A Method for Registration of 3-D Shapes" (ICP)
 */

#include "slam_data_assoc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

/* External */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern int slam_inv2x2(const double A[4], double Ainv[4]);
extern void slam_matvec_mul(const double *A, const double *x,
                             int m, int n, double *y);
extern double slam_det2x2(const double A[4]);

/* =========================================================================
 * L3: Mahalanobis Distance
 * ========================================================================= */

double slam_mahalanobis_sq(const double *innovation,
                            const double *S,
                            int dim_z) {
    if (!innovation || !S || dim_z <= 0) return -1.0;

    if (dim_z == 2) {
        /* Specialized 2×2 for efficiency */
        double S_inv[4];
        if (!slam_inv2x2(S, S_inv)) return 1e99;

        return innovation[0]*S_inv[0]*innovation[0]
             + innovation[0]*S_inv[1]*innovation[1]
             + innovation[1]*S_inv[2]*innovation[0]
             + innovation[1]*S_inv[3]*innovation[1];
    }

    /* General case: S^{-1} via LU, then ν^T·S^{-1}·ν */
    /* For dim_z=1: */
    if (dim_z == 1) {
        if (fabs(S[0]) < 1e-12) return 1e99;
        return innovation[0] * innovation[0] / S[0];
    }

    /* For larger dim_z, would need full matrix inverse */
    return -1.0;
}

void slam_innovation_covariance_rb(const double H_pose[6],
                                    const double Sigma_rr[9],
                                    const double Q[4],
                                    double S[4]) {
    if (!H_pose || !Sigma_rr || !Q || !S) return;

    /* S = H·Σ_rr·H^T + Q */
    /* H is 2×3, Σ_rr is 3×3 */
    double H_Sigma[6] = {0}; /* 2×3 */
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += H_pose[i*3 + k] * Sigma_rr[k*3 + j];
            }
            H_Sigma[i*3 + j] = sum;
        }
    }

    /* S = H_Sigma · H^T */
    S[0] = H_Sigma[0]*H_pose[0] + H_Sigma[1]*H_pose[1] + H_Sigma[2]*H_pose[2]
         + Q[0];
    S[1] = H_Sigma[0]*H_pose[3] + H_Sigma[1]*H_pose[4] + H_Sigma[2]*H_pose[5]
         + Q[1];
    S[2] = H_Sigma[3]*H_pose[0] + H_Sigma[4]*H_pose[1] + H_Sigma[5]*H_pose[2]
         + Q[2];
    S[3] = H_Sigma[3]*H_pose[3] + H_Sigma[4]*H_pose[4] + H_Sigma[5]*H_pose[5]
         + Q[3];
}

/* =========================================================================
 * L5: Nearest Neighbor Data Association
 * ========================================================================= */

int slam_da_nearest_neighbor(const slam_pose2d_t *robot_pose,
                              const double Sigma_rr[9],
                              const slam_obs_rb_t *obs,
                              const slam_landmark2d_t *landmarks,
                              int num_landmarks,
                              double gate_thresh,
                              int *best_idx,
                              double *dist_sq) {
    if (!robot_pose || !obs || !landmarks || !best_idx || !dist_sq)
        return SLAM_ERR_NULL_PTR;

    *best_idx = -1;
    *dist_sq = 1e99;

    double rx = robot_pose->x, ry = robot_pose->y, rt = robot_pose->theta;
    double Q[4] = {obs->range_std * obs->range_std, 0,
                   0, obs->bearing_std * obs->bearing_std};

    for (int j = 0; j < num_landmarks; j++) {
        const slam_landmark2d_t *lm = &landmarks[j];
        if (!lm->is_active) continue;

        double dx = lm->x - rx, dy = lm->y - ry;
        double q = dx*dx + dy*dy;
        if (q < 1e-12) continue;
        double sqrt_q = sqrt(q);

        /* Predicted measurement */
        double z_hat_r = sqrt_q;
        double z_hat_b = slam_normalize_angle(atan2(dy, dx) - rt);

        /* Innovation */
        double nu[2] = {
            obs->range - z_hat_r,
            slam_angle_diff(obs->bearing, z_hat_b)
        };

        /* Observation Jacobian w.r.t. robot pose (2×3) */
        double H[6] = {
            -dx/sqrt_q, -dy/sqrt_q, 0,
             dy/q,      -dx/q,     -1
        };

        /* Innovation covariance */
        double S[4] = {0};
        /* S = H·Σ_rr·H^T + Q */
        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 2; b++) {
                double sum = 0.0;
                for (int p = 0; p < 3; p++)
                    for (int qq = 0; qq < 3; qq++)
                        sum += H[a*3+p] * Sigma_rr[p*3+qq] * H[b*3+qq];
                S[a*2+b] = sum;
            }
        }
        S[0] += Q[0]; S[3] += Q[3];

        /* Mahalanobis distance */
        double d2 = slam_mahalanobis_sq(nu, S, 2);
        if (d2 < 0) continue;

        if (d2 < gate_thresh && d2 < *dist_sq) {
            *dist_sq = d2;
            *best_idx = j;
        }
    }

    return SLAM_OK;
}

/* =========================================================================
 * L5: Joint Compatibility Branch and Bound (JCBB)
 * ========================================================================= */

int slam_da_joint_compatible(const double *innovations,
                              const double *S_joint,
                              int K,
                              double gate_thresh) {
    if (!innovations || !S_joint || K <= 0) return 0;

    if (K == 1) {
        double d2 = slam_mahalanobis_sq(innovations, S_joint, 2);
        return (d2 >= 0 && d2 < gate_thresh) ? 1 : 0;
    }

    /* For K > 1, compute D² = ν^T·S^{-1}·ν */
    /* This requires full matrix inverse — for now, check compatibility
     * of each observation individually (permissive) */
    int dim = 2 * K;
    for (int k = 0; k < K; k++) {
        double d2 = slam_mahalanobis_sq(&innovations[2*k],
                                         &S_joint[2*k * dim + 2*k], 2);
        /* Note: this uses diagonal blocks only, not the full joint test */
        if (d2 < 0 || d2 >= gate_thresh) return 0;
    }
    return 1;
}

/* JCBB: recursive branch-and-bound */
typedef struct {
    const slam_pose2d_t *robot_pose;
    const double *Sigma_rr;
    const slam_obs_rb_t *observations;
    const slam_landmark2d_t *landmarks;
    int N;         /* number of landmarks */
    double gate;
    int *best_assoc;
    int best_count;
    double best_quality;
} jcbb_ctx_t;

static void jcbb_search(jcbb_ctx_t *ctx, int *current_assoc,
                         int obs_idx, int K, int num_paired) {
    /* Base: all observations processed */
    if (obs_idx >= K) {
        if (num_paired > ctx->best_count) {
            ctx->best_count = num_paired;
            memcpy(ctx->best_assoc, current_assoc, K * sizeof(int));
        }
        return;
    }

    /* Branch: try each landmark (including -1 for "new") */
    /* Bound: can we beat best_count? */
    int remaining = K - obs_idx;
    if (num_paired + remaining < ctx->best_count) return;

    const slam_obs_rb_t *obs = &ctx->observations[obs_idx];
    double rx = ctx->robot_pose->x;
    double ry = ctx->robot_pose->y;
    double rt = ctx->robot_pose->theta;

    /* Option 1: associate to each landmark */
    for (int j = 0; j < ctx->N; j++) {
        /* Quick gate check */
        double dx = ctx->landmarks[j].x - rx;
        double dy = ctx->landmarks[j].y - ry;
        double q = dx*dx + dy*dy;
        if (q < 1e-12) continue;
        double sqrt_q = sqrt(q);

        double zr = sqrt_q;
        double zb = slam_normalize_angle(atan2(dy, dx) - rt);
        double nu[2] = {obs->range - zr,
                        slam_angle_diff(obs->bearing, zb)};
        double d2 = nu[0]*nu[0]/(obs->range_std*obs->range_std)
                  + nu[1]*nu[1]/(obs->bearing_std*obs->bearing_std);

        if (d2 < ctx->gate) {
            current_assoc[obs_idx] = j;
            jcbb_search(ctx, current_assoc, obs_idx + 1, K, num_paired + 1);
        }
    }

    /* Option 2: new landmark */
    current_assoc[obs_idx] = -1;
    jcbb_search(ctx, current_assoc, obs_idx + 1, K, num_paired);
}

int slam_da_jcbb(const slam_pose2d_t *robot_pose,
                  const double Sigma_rr[9],
                  const slam_obs_rb_t *observations,
                  int K,
                  const slam_landmark2d_t *landmarks,
                  int N,
                  double gate_thresh,
                  int *associations,
                  int *num_paired) {
    if (!robot_pose || !observations || !landmarks || !associations
        || !num_paired)
        return SLAM_ERR_NULL_PTR;

    jcbb_ctx_t ctx;
    ctx.robot_pose   = robot_pose;
    ctx.Sigma_rr     = Sigma_rr;
    ctx.observations = observations;
    ctx.landmarks    = landmarks;
    ctx.N            = N;
    ctx.gate         = gate_thresh;
    ctx.best_assoc   = associations;
    ctx.best_count   = 0;

    int *current = (int *)malloc(K * sizeof(int));
    if (!current) return SLAM_ERR_MEMORY;
    for (int i = 0; i < K; i++) current[i] = -1;

    jcbb_search(&ctx, current, 0, K, 0);

    *num_paired = ctx.best_count;
    free(current);
    return SLAM_OK;
}

/* =========================================================================
 * L5: ICP for 2D Point Cloud Alignment
 * ========================================================================= */

int slam_da_icp_2d(const double *src_pts,
                    const double *tgt_pts,
                    int N,
                    int max_iter,
                    double eps,
                    slam_pose2d_t *T,
                    double *rmse) {
    if (!src_pts || !tgt_pts || !T || !rmse || N <= 0)
        return SLAM_ERR_INVALID_PARAM;

    /* Initialize transformation to zero */
    T->x = 0; T->y = 0; T->theta = 0;

    double *src_trans = (double *)malloc(2 * N * sizeof(double));
    int *correspondences = (int *)malloc(N * sizeof(int));
    if (!src_trans || !correspondences) {
        free(src_trans); free(correspondences);
        return SLAM_ERR_MEMORY;
    }

    double prev_error = DBL_MAX;

    for (int iter = 0; iter < max_iter; iter++) {
        /* 1. Transform source points by current T */
        double ct = cos(T->theta), st = sin(T->theta);
        for (int i = 0; i < N; i++) {
            double sx = src_pts[2*i], sy = src_pts[2*i + 1];
            src_trans[2*i]     = ct * sx - st * sy + T->x;
            src_trans[2*i + 1] = st * sx + ct * sy + T->y;
        }

        /* 2. Find nearest neighbors */
        for (int i = 0; i < N; i++) {
            double min_dist = DBL_MAX;
            int best_j = 0;
            for (int j = 0; j < N; j++) {
                double dx = src_trans[2*i] - tgt_pts[2*j];
                double dy = src_trans[2*i + 1] - tgt_pts[2*j + 1];
                double d2 = dx*dx + dy*dy;
                if (d2 < min_dist) { min_dist = d2; best_j = j; }
            }
            correspondences[i] = best_j;
        }

        /* 3. Compute centroids */
        double src_cx = 0, src_cy = 0, tgt_cx = 0, tgt_cy = 0;
        for (int i = 0; i < N; i++) {
            int j = correspondences[i];
            src_cx += src_pts[2*i]; src_cy += src_pts[2*i + 1];
            tgt_cx += tgt_pts[2*j]; tgt_cy += tgt_pts[2*j + 1];
        }
        src_cx /= N; src_cy /= N;
        tgt_cx /= N; tgt_cy /= N;

        /* 4. Compute cross-covariance H */
        double H11 = 0, H12 = 0, H21 = 0, H22 = 0;
        for (int i = 0; i < N; i++) {
            int j = correspondences[i];
            double sx = src_pts[2*i] - src_cx;
            double sy = src_pts[2*i + 1] - src_cy;
            double tx = tgt_pts[2*j] - tgt_cx;
            double ty = tgt_pts[2*j + 1] - tgt_cy;
            H11 += sx * tx; H12 += sx * ty;
            H21 += sy * tx; H22 += sy * ty;
        }

        /* 5. SVD in 2D: closed-form rotation
         * Using atan2-based solution from the cross-covariance matrix */
        (void)H11; (void)H21; /* used above via H_norm computation which is
                                 folded into the atan2 solution */
        double cos_th = (H11 + H22) / (sqrt((H11+H22)*(H11+H22) + (H12-H21)*(H12-H21)) + 1e-12);
        double sin_th = (H12 - H21) / (sqrt((H11+H22)*(H11+H22) + (H12-H21)*(H12-H21)) + 1e-12);
        double dtheta = atan2(sin_th, cos_th);

        /* 6. Translation */
        double dx = tgt_cx - (cos_th * src_cx - sin_th * src_cy);
        double dy = tgt_cy - (sin_th * src_cx + cos_th * src_cy);

        /* 7. Update cumulative transformation */
        T->x += T->x * cos(dtheta) - T->y * sin(dtheta) + dx;
        T->y += T->x * sin(dtheta) + T->y * cos(dtheta) + dy;
        T->theta = slam_normalize_angle(T->theta + dtheta);

        /* Recompute: T_new = Δ ⊕ T_old */
        /* Better: apply increment to src points directly */
        double ct2 = cos(T->theta), st2 = sin(T->theta);
        double total_error = 0.0;
        for (int i = 0; i < N; i++) {
            double sx = src_pts[2*i], sy = src_pts[2*i + 1];
            double wx = ct2 * sx - st2 * sy + T->x;
            double wy = st2 * sx + ct2 * sy + T->y;
            int j = correspondences[i];
            double err = sqrt((wx - tgt_pts[2*j])*(wx - tgt_pts[2*j])
                            + (wy - tgt_pts[2*j+1])*(wy - tgt_pts[2*j+1]));
            total_error += err;
        }
        total_error /= N;

        if (fabs(prev_error - total_error) < eps) break;
        prev_error = total_error;
    }

    *rmse = prev_error;
    free(src_trans); free(correspondences);
    return SLAM_OK;
}

/* =========================================================================
 * L6: Loop Closure Detection
 * ========================================================================= */

int slam_da_detect_loop(const slam_pose2d_t *poses,
                         int num_poses,
                         const slam_pose2d_t *current_pose,
                         double search_radius,
                         int min_steps,
                         int *match_idx) {
    if (!poses || !current_pose || !match_idx) return SLAM_ERR_NULL_PTR;

    *match_idx = -1;

    /* Search backwards from num_poses-1, skipping recent poses */
    int search_end = num_poses - min_steps;
    if (search_end < 0) search_end = 0;

    double best_dist = search_radius + 1e-6;
    int best_i = -1;

    for (int i = 0; i < search_end; i++) {
        double dx = current_pose->x - poses[i].x;
        double dy = current_pose->y - poses[i].y;
        double dist = sqrt(dx*dx + dy*dy);

        /* Also check orientation similarity */
        double dtheta = fabs(slam_angle_diff(current_pose->theta,
                                              poses[i].theta));

        /* Loop closure: close in position AND rough orientation */
        if (dist < best_dist && dtheta < M_PI/2) {
            best_dist = dist;
            best_i = i;
        }
    }

    if (best_i >= 0 && best_dist <= search_radius) {
        *match_idx = best_i;
    }

    return SLAM_OK;
}
