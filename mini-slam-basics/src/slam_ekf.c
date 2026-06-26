/**
 * @file    slam_ekf.c
 * @brief   EKF-SLAM implementation: prediction, update, state augmentation
 *
 * Implements the full EKF-SLAM pipeline with velocity motion model
 * and range-bearing observations.
 *
 * EKF-SLAM Theorem (Dissanayake et al. 2001):
 *   - The determinant of any submatrix of the map covariance matrix
 *     decreases monotonically as observations are made.
 *   - In the limit of infinite observations, the relative positions
 *     of all landmarks become fully correlated.
 *   - The absolute accuracy of the map is bounded below by the
 *     initial robot position uncertainty.
 *
 * Reference:
 *   Smith, Self & Cheeseman (1990)
 *   Dissanayake, Newman, Clark et al. (2001) IEEE T-RO
 *   Thrun, Burgard & Fox (2005), Chapter 10.
 */

#include "slam_ekf.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* Declaration of functions from slam_core.c */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern void slam_pose_compose(const slam_pose2d_t *a,
                               const slam_pose2d_t *b,
                               slam_pose2d_t *c);
extern void slam_matvec_mul(const double *A, const double *x,
                             int m, int n, double *y);
extern void slam_matmul(const double *A, const double *B,
                         int m, int k, int n, double *C);
extern void slam_matmul_AT_B(const double *A, const double *B,
                              int k, int m, int n, double *C);
extern void slam_matmul_A_BT(const double *A, const double *B,
                              int m, int k, int n, double *C);
extern void slam_covariance_propagate(const double *Sigma,
                                       const double *J, const double *Q,
                                       int m, double *S);
extern double slam_det2x2(const double A[4]);
extern int slam_inv2x2(const double A[4], double Ainv[4]);
extern void slam_eye(double *M, int n);
extern void slam_mat_copy(const double *src, double *dst, int n);
extern void slam_mat_add(const double *A, const double *B,
                          double *C, int n);
extern double slam_det3x3(const double M[9]);
extern double slam_subdet(const double *M, int n, int k);
extern int slam_map2d_add_landmark(slam_map2d_t *map,
                                    const slam_landmark2d_t *lm, int *idx);

/* =========================================================================
 * L2: EKF-SLAM Initialization
 * ========================================================================= */

int slam_ekf_init(slam_ekf_state_t *state,
                  const slam_config_t *config,
                  const slam_pose2d_t *init_pose,
                  double sigma_pos,
                  double sigma_theta) {
    if (!state || !config || !init_pose) return SLAM_ERR_NULL_PTR;
    if (sigma_pos < 0 || sigma_theta < 0) return SLAM_ERR_INVALID_PARAM;

    state->robot_pose    = *init_pose;
    state->num_landmarks = 0;
    state->state_dim     = 3;  /* initially just pose */
    state->cov_stride    = 3 + 2 * config->max_landmarks;
    state->status        = SLAM_STATUS_INITIALIZING;

    /* Allocate mean vector and covariance */
    int stride = state->cov_stride;
    state->state_mean = (double *)calloc(stride, sizeof(double));
    state->covariance = (double *)calloc(stride * stride, sizeof(double));

    if (!state->state_mean || !state->covariance) {
        free(state->state_mean);
        free(state->covariance);
        return SLAM_ERR_MEMORY;
    }

    /* Set initial mean */
    state->state_mean[0] = init_pose->x;
    state->state_mean[1] = init_pose->y;
    state->state_mean[2] = init_pose->theta;

    /* Set initial covariance: diag(σ_x², σ_y², σ_θ²) */
    state->covariance[0 * stride + 0] = sigma_pos * sigma_pos;
    state->covariance[1 * stride + 1] = sigma_pos * sigma_pos;
    state->covariance[2 * stride + 2] = sigma_theta * sigma_theta;

    state->timestamp = 0;
    return SLAM_OK;
}

void slam_ekf_free(slam_ekf_state_t *state) {
    if (state) {
        free(state->state_mean);
        free(state->covariance);
        state->state_mean = NULL;
        state->covariance = NULL;
    }
}

/* =========================================================================
 * L5: EKF Prediction with Velocity Model
 * ========================================================================= */

int slam_ekf_predict_velocity(slam_ekf_state_t *state,
                               const slam_velocity_t *vel) {
    if (!state || !vel) return SLAM_ERR_NULL_PTR;

    double v = vel->v, w = vel->omega, dt = vel->dt;
    double theta = state->robot_pose.theta;
    int s = state->cov_stride;  /* stride for all matrix indexing */
    int d = state->state_dim;   /* active state dimension */

    /* 1. State prediction (nonlinear) */
    if (fabs(w) < 1e-8) {
        state->robot_pose.x += v * dt * cos(theta);
        state->robot_pose.y += v * dt * sin(theta);
    } else {
        double v_over_w = v / w;
        state->robot_pose.x += v_over_w * (sin(theta + w * dt) - sin(theta));
        state->robot_pose.y += v_over_w * (cos(theta) - cos(theta + w * dt));
        state->robot_pose.theta += w * dt;
    }
    state->robot_pose.theta = slam_normalize_angle(state->robot_pose.theta);
    state->state_mean[0] = state->robot_pose.x;
    state->state_mean[1] = state->robot_pose.y;
    state->state_mean[2] = state->robot_pose.theta;

    /* 2. Jacobian G (active state only, packed to d×d) */
    double *G = (double *)calloc(d * d, sizeof(double));
    if (!G) return SLAM_ERR_MEMORY;
    slam_eye(G, d);
    if (fabs(w) < 1e-8) {
        G[0*d + 2] = -v * dt * sin(theta);
        G[1*d + 2] =  v * dt * cos(theta);
    } else {
        double v_over_w = v / w;
        double theta_new = theta + w * dt;
        G[0*d + 2] = -v_over_w * (cos(theta) - cos(theta_new));
        G[1*d + 2] = -v_over_w * (sin(theta_new) - sin(theta));
    }

    /* 3-5. Covariance prediction using only the active d×d sub-block.
     * Copy active sub-block from full covariance (stride s) to packed (stride d),
     * operate on packed matrices, then copy back. */
    double *cov_packed = (double *)calloc(d * d, sizeof(double));
    double *G_Sigma = (double *)calloc(d * d, sizeof(double));
    double *new_cov = (double *)calloc(d * d, sizeof(double));
    if (!cov_packed || !G_Sigma || !new_cov) {
        free(G); free(cov_packed); free(G_Sigma); free(new_cov);
        return SLAM_ERR_MEMORY;
    }

    /* Copy active d×d from full covariance */
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            cov_packed[i*d + j] = state->covariance[i*s + j];

    /* Temporaries for R (motion noise, only top-left 3×3 non-zero) */
    double sigma_v = 0.1, sigma_omega = 0.05, sigma_gamma = 0.02;
    double Q_ctrl[9] = {0};
    Q_ctrl[0] = sigma_v * sigma_v;
    Q_ctrl[4] = sigma_omega * sigma_omega;
    Q_ctrl[8] = sigma_gamma * sigma_gamma;

    /* R = F_x · Q_ctrl · F_x^T, F_x is (d×3), non-zero only in first 3 rows */
    double *Fx = (double *)calloc(d * 3, sizeof(double));
    if (!Fx) { free(G); free(cov_packed); free(G_Sigma); free(new_cov); return SLAM_ERR_MEMORY; }
    Fx[0*3+0] = cos(theta) * dt;
    Fx[1*3+0] = sin(theta) * dt;
    Fx[2*3+1] = dt;
    Fx[2*3+2] = dt;

    double R[81] = {0}; /* d×d max, only top-left 3×3 used */
    double Fx_Q[9] = {0};
    slam_matmul(Fx, Q_ctrl, d, 3, 3, Fx_Q);
    slam_matmul_A_BT(Fx_Q, Fx, d, 3, d, R);

    /* Σ = G·Σ·G^T + R, all packed d×d */
    slam_matmul(G, cov_packed, d, d, d, G_Sigma);
    slam_matmul_A_BT(G_Sigma, G, d, d, d, new_cov);
    for (int i = 0; i < d*d; i++) new_cov[i] += R[i];

    /* Copy back to full covariance */
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            state->covariance[i*s + j] = new_cov[i*d + j];

    free(G); free(Fx); free(cov_packed); free(G_Sigma); free(new_cov);
    state->timestamp++;
    return SLAM_OK;
}

/* =========================================================================
 * L5: EKF Update with Range-Bearing Observation
 * ========================================================================= */

int slam_ekf_update_rb(slam_ekf_state_t *state,
                       const slam_obs_rb_t *obs,
                       int lm_idx) {
    if (!state || !obs) return SLAM_ERR_NULL_PTR;
    if (lm_idx < 0 || lm_idx >= state->num_landmarks)
        return SLAM_ERR_NO_ASSOC;

    int s = state->cov_stride;
    int d = state->state_dim;

    double rx = state->robot_pose.x;
    double ry = state->robot_pose.y;
    double rt = state->robot_pose.theta;

    int lm_base = 3 + 2 * lm_idx;
    double mx = state->state_mean[lm_base];
    double my = state->state_mean[lm_base + 1];

    /* 1. Predicted measurement & innovation */
    double dx = mx - rx, dy = my - ry;
    double q = dx*dx + dy*dy;
    if (q < 1e-12) return SLAM_ERR_SINGULAR;
    double sqrt_q = sqrt(q);

    double nu[2] = {
        obs->range - sqrt_q,
        slam_angle_diff(obs->bearing,
                         slam_normalize_angle(atan2(dy, dx) - rt))
    };

    /* 2. Observation Jacobians */
    double H_pose[6] = {-dx/sqrt_q, -dy/sqrt_q, 0,
                         dy/q,      -dx/q,     -1};
    double H_lm[4]   = {dx/sqrt_q, dy/sqrt_q,
                        -dy/q,      dx/q};
    double Qm[4] = {obs->range_std * obs->range_std, 0,
                    0, obs->bearing_std * obs->bearing_std};

    /* 3. Build S = H·Σ·H^T + Q using packed covariance */
    /* Copy active d×d submatrix from full covariance */
    double *cov_p = (double *)calloc(d * d, sizeof(double));
    if (!cov_p) return SLAM_ERR_MEMORY;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            cov_p[i*d + j] = state->covariance[i*s + j];

    /* S = H_pose·Σ_rr·H_pose^T + H_lm·Σ_ll·H_lm^T + cross + Q */
    double S[4] = {0};
    for (int a = 0; a < 2; a++) {
        for (int b = 0; b < 2; b++) {
            /* Robot pose contribution */
            for (int p = 0; p < 3; p++)
                for (int qq = 0; qq < 3; qq++)
                    S[a*2+b] += H_pose[a*3+p] * cov_p[p*d + qq] * H_pose[b*3+qq];
            /* Cross terms */
            for (int p = 0; p < 3; p++)
                for (int qq = 0; qq < 2; qq++)
                    S[a*2+b] += H_pose[a*3+p] * cov_p[p*d + lm_base + qq] * H_lm[b*2+qq]
                              + H_lm[a*2+qq] * cov_p[(lm_base+qq)*d + p] * H_pose[b*3+p];
            /* Landmark contribution */
            for (int p = 0; p < 2; p++)
                for (int qq = 0; qq < 2; qq++)
                    S[a*2+b] += H_lm[a*2+p] * cov_p[(lm_base+p)*d + lm_base + qq] * H_lm[b*2+qq];
        }
    }
    S[0] += Qm[0]; S[3] += Qm[3];

    double S_inv[4];
    if (!slam_inv2x2(S, S_inv)) { free(cov_p); return SLAM_ERR_SINGULAR; }

    /* 4. Kalman gain K = Σ·H^T·S^{-1} */
    double *K = (double *)calloc(d * 2, sizeof(double));
    if (!K) { free(cov_p); return SLAM_ERR_MEMORY; }
    for (int row = 0; row < d; row++) {
        double kr = 0, kb = 0;
        for (int c = 0; c < 3; c++) {
            kr += cov_p[row*d + c] * H_pose[c];
            kb += cov_p[row*d + c] * H_pose[3+c];
        }
        kr += cov_p[row*d + lm_base]     * H_lm[0]
            + cov_p[row*d + lm_base + 1] * H_lm[1];
        kb += cov_p[row*d + lm_base]     * H_lm[2]
            + cov_p[row*d + lm_base + 1] * H_lm[3];
        K[row*2+0] = kr * S_inv[0] + kb * S_inv[2];
        K[row*2+1] = kr * S_inv[1] + kb * S_inv[3];
    }

    /* 5. State update: μ = μ + K·ν */
    for (int row = 0; row < d; row++) {
        state->state_mean[row] += K[row*2+0] * nu[0] + K[row*2+1] * nu[1];
    }
    state->robot_pose.x     = state->state_mean[0];
    state->robot_pose.y     = state->state_mean[1];
    state->robot_pose.theta = slam_normalize_angle(state->state_mean[2]);
    state->state_mean[2]    = state->robot_pose.theta;

    /* 6. Covariance update (Joseph form): Σ = Σ − K·H·Σ */
    /* Build H_full: 2×d packed */
    double *H_full = (double *)calloc(2 * d, sizeof(double));
    if (!H_full) { free(K); free(cov_p); return SLAM_ERR_MEMORY; }
    H_full[0*d+0] = H_pose[0]; H_full[0*d+1] = H_pose[1]; H_full[0*d+2] = H_pose[2];
    H_full[0*d+lm_base]     = H_lm[0]; H_full[0*d+lm_base+1] = H_lm[1];
    H_full[1*d+0] = H_pose[3]; H_full[1*d+1] = H_pose[4]; H_full[1*d+2] = H_pose[5];
    H_full[1*d+lm_base]     = H_lm[2]; H_full[1*d+lm_base+1] = H_lm[3];

    double *KH = (double *)calloc(d * d, sizeof(double));
    if (!KH) { free(K); free(cov_p); free(H_full); return SLAM_ERR_MEMORY; }
    slam_matmul(K, H_full, d, 2, d, KH);
    double *KH_S = (double *)calloc(d * d, sizeof(double));
    if (!KH_S) { free(K); free(cov_p); free(H_full); free(KH); return SLAM_ERR_MEMORY; }
    slam_matmul(KH, cov_p, d, d, d, KH_S);
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            state->covariance[i*s + j] = cov_p[i*d + j] - KH_S[i*d + j];

    free(K); free(cov_p); free(H_full); free(KH); free(KH_S);
    state->timestamp++;
    return SLAM_OK;
}

/* =========================================================================
 * L5: EKF Update with Unknown Association
 * ========================================================================= */

int slam_ekf_update_unknown(slam_ekf_state_t *state,
                             const slam_obs_rb_t *obs,
                             int *matched_id) {
    if (!state || !obs || !matched_id) return SLAM_ERR_NULL_PTR;

    *matched_id = -1;
    double best_d2 = 1e99;
    int best_j = -1;

    double rx = state->robot_pose.x;
    double ry = state->robot_pose.y;
    double rt = state->robot_pose.theta;
    int s = state->cov_stride;

    /* Extract robot pose covariance */
    double Sigma_rr[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Sigma_rr[i*3+j] = state->covariance[i * s + j];

    double Q[4] = {obs->range_std * obs->range_std, 0,
                   0, obs->bearing_std * obs->bearing_std};

    for (int j = 0; j < state->num_landmarks; j++) {
        int bj = 3 + 2 * j;
        double mx = state->state_mean[bj];
        double my = state->state_mean[bj + 1];

        double dx = mx - rx, dy = my - ry;
        double q = dx*dx + dy*dy;
        if (q < 1e-12) continue;
        double sqrt_q = sqrt(q);

        double z_hat_r = sqrt_q;
        double z_hat_b = slam_normalize_angle(atan2(dy, dx) - rt);

        double nu[2] = {obs->range - z_hat_r,
                        slam_angle_diff(obs->bearing, z_hat_b)};

        /* Build H for this landmark */
        double H_p[6] = {-dx/sqrt_q, -dy/sqrt_q, 0,
                          dy/q,      -dx/q,     -1};

        /* Innovation covariance: H_pose·Σ_rr·H_pose^T + Q */
        double S[4] = {0};
        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 2; b++) {
                double sum = 0.0;
                for (int p = 0; p < 3; p++)
                    for (int qq = 0; qq < 3; qq++)
                        sum += H_p[a*3+p] * Sigma_rr[p*3+qq] * H_p[b*3+qq];
                S[a*2+b] = sum;
            }
        }
        S[0] += Q[0]; S[3] += Q[3];

        double S_inv[4];
        if (!slam_inv2x2(S, S_inv)) continue;

        /* Mahalanobis distance */
        double d2 = nu[0]*S_inv[0]*nu[0] + nu[0]*S_inv[1]*nu[1]
                  + nu[1]*S_inv[2]*nu[0] + nu[1]*S_inv[3]*nu[1];

        if (d2 < 5.991 && d2 < best_d2) {
            best_d2 = d2;
            best_j = j;
        }
    }

    if (best_j >= 0) {
        *matched_id = best_j;
        return slam_ekf_update_rb(state, obs, best_j);
    }
    return SLAM_ERR_NO_ASSOC;
}

/* =========================================================================
 * L5: State Augmentation — adding new landmarks
 * ========================================================================= */

int slam_ekf_augment(slam_ekf_state_t *state,
                     const slam_obs_rb_t *obs,
                     int *new_id) {
    if (!state || !obs || !new_id) return SLAM_ERR_NULL_PTR;

    int N = state->num_landmarks;
    int d_old = state->state_dim;   /* = 3 + 2*N */
    int d_new = d_old + 2;
    int s = state->cov_stride;      /* full allocation stride */

    double rx = state->robot_pose.x;
    double ry = state->robot_pose.y;
    double rt = state->robot_pose.theta;

    /* 1. New landmark position (inverse observation) */
    double mx = rx + obs->range * cos(obs->bearing + rt);
    double my = ry + obs->range * sin(obs->bearing + rt);
    state->state_mean[d_old]     = mx;
    state->state_mean[d_old + 1] = my;

    /* 2. Jacobians */
    double r = obs->range, phi = obs->bearing;
    double sp = sin(phi + rt), cp = cos(phi + rt);
    double G_rl[6] = {1, 0, -r * sp,  0, 1, r * cp};
    double G_z[4]  = {cp, -r * sp, sp, r * cp};

    /* 3. Copy active covariance to packed */
    double *cov_p = (double *)calloc(d_old * d_old, sizeof(double));
    if (!cov_p) return SLAM_ERR_MEMORY;
    for (int i = 0; i < d_old; i++)
        for (int j = 0; j < d_old; j++)
            cov_p[i*d_old + j] = state->covariance[i*s + j];

    /* 4. Cross-covariance: Σ_rl = Σ_old · G_rl^T (d_old × 2) */
    /* G_rl_full: 2 × d_old, sparse (non-zero only in first 3 cols) */
    double *G_rl_full = (double *)calloc(2 * d_old, sizeof(double));
    if (!G_rl_full) { free(cov_p); return SLAM_ERR_MEMORY; }
    G_rl_full[0*d_old+0] = G_rl[0]; G_rl_full[0*d_old+1] = G_rl[1];
    G_rl_full[0*d_old+2] = G_rl[2];
    G_rl_full[1*d_old+0] = G_rl[3]; G_rl_full[1*d_old+1] = G_rl[4];
    G_rl_full[1*d_old+2] = G_rl[5];

    double *Sigma_rl = (double *)calloc(d_old * 2, sizeof(double));
    if (!Sigma_rl) { free(cov_p); free(G_rl_full); return SLAM_ERR_MEMORY; }
    slam_matmul_A_BT(cov_p, G_rl_full, d_old, d_old, 2, Sigma_rl);

    /* 5. New landmark covariance: Σ_ll = G_rl·Σ·G_rl^T + G_z·Q·G_z^T */
    double Q[4] = {obs->range_std*obs->range_std, 0,
                   0, obs->bearing_std*obs->bearing_std};
    double G_Sigma[4];
    slam_matmul(G_rl_full, cov_p, 2, d_old, d_old, G_Sigma);
    double G_Sigma_GT[4];
    slam_matmul_A_BT(G_Sigma, G_rl_full, 2, d_old, 2, G_Sigma_GT);
    double tmp[4], GzQGzT[4];
    slam_matmul(G_z, Q, 2, 2, 2, tmp);
    slam_matmul_A_BT(tmp, G_z, 2, 2, 2, GzQGzT);
    for (int i = 0; i < 4; i++) G_Sigma_GT[i] += GzQGzT[i];

    /* 6. Write back to full covariance */
    /* Old block stays */
    for (int i = 0; i < d_old; i++)
        for (int j = 0; j < d_old; j++)
            state->covariance[i*s + j] = cov_p[i*d_old + j];
    /* Cross blocks */
    for (int i = 0; i < d_old; i++) {
        state->covariance[i*s + d_old]     = Sigma_rl[i*2 + 0];
        state->covariance[i*s + d_old + 1] = Sigma_rl[i*2 + 1];
        state->covariance[d_old*s + i]     = Sigma_rl[i*2 + 0];
        state->covariance[(d_old+1)*s + i] = Sigma_rl[i*2 + 1];
    }
    /* New landmark block */
    state->covariance[d_old*s + d_old]         = G_Sigma_GT[0];
    state->covariance[d_old*s + d_old + 1]     = G_Sigma_GT[1];
    state->covariance[(d_old+1)*s + d_old]     = G_Sigma_GT[2];
    state->covariance[(d_old+1)*s + d_old + 1] = G_Sigma_GT[3];

    state->num_landmarks++;
    state->state_dim = d_new;
    *new_id = N;

    free(cov_p); free(G_rl_full); free(Sigma_rl);
    state->timestamp++;
    return SLAM_OK;
}

/* =========================================================================
 * L6: EKF-SLAM Full Step
 * ========================================================================= */

int slam_ekf_step(slam_ekf_state_t *state,
                  const slam_velocity_t *vel,
                  const slam_obs_rb_t *observations,
                  int num_obs) {
    if (!state) return SLAM_ERR_NULL_PTR;

    /* 1. Prediction */
    if (vel) {
        int rc = slam_ekf_predict_velocity(state, vel);
        if (rc != SLAM_OK) return rc;
    }

    /* 2. Process each observation */
    if (observations && num_obs > 0) {
        for (int i = 0; i < num_obs; i++) {
            int matched;
            int rc = slam_ekf_update_unknown(state, &observations[i], &matched);
            if (rc == SLAM_ERR_NO_ASSOC) {
                /* New landmark */
                int new_id;
                slam_ekf_augment(state, &observations[i], &new_id);
            } else if (rc != SLAM_OK) {
                /* Silently skip problematic observations in this implementation */
            }
        }
    }

    state->status = SLAM_STATUS_RUNNING;
    return SLAM_OK;
}

/* =========================================================================
 * L6: Consistency Measures
 * ========================================================================= */

int slam_ekf_nees(const slam_ekf_state_t *state,
                  const slam_pose2d_t *true_pose,
                  double *nees_out) {
    if (!state || !true_pose || !nees_out) return SLAM_ERR_NULL_PTR;

    /* Error for robot pose only (first 3 states) */
    double err[3] = {
        state->state_mean[0] - true_pose->x,
        state->state_mean[1] - true_pose->y,
        slam_angle_diff(state->state_mean[2], true_pose->theta)
    };

    int s = state->cov_stride;
    double Sigma_rr[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Sigma_rr[i*3+j] = state->covariance[i * s + j];

    /* NEES = err^T · Σ^{-1} · err */
    double Sigma_inv[9];
    int ok = 1;
    double det = slam_det3x3(Sigma_rr);
    if (fabs(det) < 1e-12) ok = 0;
    if (ok) {
        /* Use cofactor inverse for 3×3 */
        double invdet = 1.0 / det;
        Sigma_inv[0] = (Sigma_rr[4]*Sigma_rr[8] - Sigma_rr[5]*Sigma_rr[7]) * invdet;
        Sigma_inv[1] = (Sigma_rr[2]*Sigma_rr[7] - Sigma_rr[1]*Sigma_rr[8]) * invdet;
        Sigma_inv[2] = (Sigma_rr[1]*Sigma_rr[5] - Sigma_rr[2]*Sigma_rr[4]) * invdet;
        Sigma_inv[3] = (Sigma_rr[5]*Sigma_rr[6] - Sigma_rr[3]*Sigma_rr[8]) * invdet;
        Sigma_inv[4] = (Sigma_rr[0]*Sigma_rr[8] - Sigma_rr[2]*Sigma_rr[6]) * invdet;
        Sigma_inv[5] = (Sigma_rr[3]*Sigma_rr[2] - Sigma_rr[0]*Sigma_rr[5]) * invdet;
        Sigma_inv[6] = (Sigma_rr[3]*Sigma_rr[7] - Sigma_rr[4]*Sigma_rr[6]) * invdet;
        Sigma_inv[7] = (Sigma_rr[6]*Sigma_rr[1] - Sigma_rr[0]*Sigma_rr[7]) * invdet;
        Sigma_inv[8] = (Sigma_rr[0]*Sigma_rr[4] - Sigma_rr[1]*Sigma_rr[3]) * invdet;

        *nees_out = 0.0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                *nees_out += err[i] * Sigma_inv[i*3+j] * err[j];
    } else {
        *nees_out = -1.0;
    }

    return ok ? SLAM_OK : SLAM_ERR_SINGULAR;
}

int slam_ekf_nis(const slam_obs_rb_t *obs,
                 const slam_ekf_state_t *state,
                 int lm_idx,
                 double *nis_out) {
    if (!obs || !state || !nis_out) return SLAM_ERR_NULL_PTR;
    if (lm_idx < 0 || lm_idx >= state->num_landmarks)
        return SLAM_ERR_NO_ASSOC;

    double rx = state->robot_pose.x;
    double ry = state->robot_pose.y;
    double rt = state->robot_pose.theta;
    int bj = 3 + 2 * lm_idx;
    double mx = state->state_mean[bj];
    double my = state->state_mean[bj + 1];
    int s = state->cov_stride;

    double dx = mx - rx, dy = my - ry;
    double q = dx*dx + dy*dy;
    if (q < 1e-12) { *nis_out = -1; return SLAM_ERR_SINGULAR; }
    double sqrt_q = sqrt(q);

    double nu[2] = {
        obs->range - sqrt_q,
        slam_angle_diff(obs->bearing,
                         slam_normalize_angle(atan2(dy, dx) - rt))
    };

    double Sigma_rr[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Sigma_rr[i*3+j] = state->covariance[i*s + j];

    double H_p[6] = {-dx/sqrt_q, -dy/sqrt_q, 0,
                      dy/q,      -dx/q,     -1};
    double Qm[4] = {obs->range_std*obs->range_std, 0,
                    0, obs->bearing_std*obs->bearing_std};

    double S[4] = {0};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int p = 0; p < 3; p++)
                for (int qq = 0; qq < 3; qq++)
                    S[a*2+b] += H_p[a*3+p] * Sigma_rr[p*3+qq] * H_p[b*3+qq];
    S[0] += Qm[0]; S[3] += Qm[3];

    double S_inv[4];
    if (!slam_inv2x2(S, S_inv)) { *nis_out = -1; return SLAM_ERR_SINGULAR; }

    *nis_out = nu[0]*S_inv[0]*nu[0] + nu[0]*S_inv[1]*nu[1]
             + nu[1]*S_inv[2]*nu[0] + nu[1]*S_inv[3]*nu[1];
    return SLAM_OK;
}

/* =========================================================================
 * L4: Monotonicity Verification (Dissanayake et al. 2001)
 * ========================================================================= */

int slam_ekf_verify_monotonicity(const double *old_cov,
                                  const double *new_cov,
                                  int dim) {
    if (!old_cov || !new_cov || dim < 3) return 0;

    /* Verify det(Σ_map) decreases for each 2×2 landmark sub-block */
    /* Iterate over landmarks (starting from index 3) */
    for (int lm_base = 3; lm_base + 1 < dim; lm_base += 2) {
        /* 2×2 sub-block for this landmark */
        double old_block[4] = {
            old_cov[lm_base * dim + lm_base],
            old_cov[lm_base * dim + lm_base + 1],
            old_cov[(lm_base+1) * dim + lm_base],
            old_cov[(lm_base+1) * dim + lm_base + 1]
        };
        double new_block[4] = {
            new_cov[lm_base * dim + lm_base],
            new_cov[lm_base * dim + lm_base + 1],
            new_cov[(lm_base+1) * dim + lm_base],
            new_cov[(lm_base+1) * dim + lm_base + 1]
        };

        double det_old = slam_det2x2(old_block);
        double det_new = slam_det2x2(new_block);

        /* Monotonicity: det_new ≤ det_old (within numerical tolerance) */
        if (det_new > det_old + 1e-10) return 0;
    }
    return 1;
}
