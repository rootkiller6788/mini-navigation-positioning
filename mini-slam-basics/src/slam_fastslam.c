/**
 * @file    slam_fastslam.c
 * @brief   FastSLAM: Rao-Blackwellized Particle Filter SLAM
 *
 * FastSLAM splits the SLAM posterior:
 *   p(x_{1:t}, Θ | z_{1:t}, u_{1:t})
 *   = p(x_{1:t} | z_{1:t}, u_{1:t}) × Π_j p(θ_j | x_{1:t}, z_{1:t})
 *
 * The robot trajectory is estimated by a particle filter (non-Gaussian,
 * multi-modal capable), while each landmark is estimated by a
 * separate low-dimensional EKF conditional on the trajectory.
 *
 * This factorization gives O(M·K) complexity vs O(N²) for EKF-SLAM,
 * where M = #particles, K = #landmarks.
 *
 * Reference:
 *   Montemerlo, Thrun, Koller & Wegbreit (2002), AAAI
 *   Montemerlo, Thrun (2003) "FastSLAM 2.0", IJCAI
 *   Thrun, Burgard & Fox (2005), Chapter 13
 */

#include "slam_fastslam.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* External declarations */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern void slam_pose_compose(const slam_pose2d_t *a,
                               const slam_pose2d_t *b,
                               slam_pose2d_t *c);
extern void slam_matmul(const double *A, const double *B,
                         int m, int k, int n, double *C);
extern void slam_matmul_A_BT(const double *A, const double *B,
                              int m, int k, int n, double *C);
extern int slam_inv2x2(const double A[4], double Ainv[4]);
extern double slam_randn(void);
extern double slam_randn_scaled(double mu, double sigma);
extern int slam_map2d_add_landmark(slam_map2d_t *map,
                                    const slam_landmark2d_t *lm, int *idx);

/* =========================================================================
 * Simple PRNG for reproducibility
 * ========================================================================= */

/* Marsaglia's xorshift128+ */
static uint64_t rng_state[2] = {123456789ULL, 987654321ULL};

static uint64_t xorshift128plus(void) {
    uint64_t s1 = rng_state[0];
    uint64_t s0 = rng_state[1];
    rng_state[0] = s0;
    s1 ^= s1 << 23;
    rng_state[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng_state[1] + s0;
}

static double uniform_random(void) {
    /* Generate uniform in [0, 1) */
    uint64_t x = xorshift128plus();
    return (x >> 11) * 0x1.0p-53;
}

double slam_randn(void) {
    /* Box-Muller: two uniforms → two normals */
    double u1 = uniform_random();
    double u2 = uniform_random();
    /* Avoid log(0) */
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

double slam_randn_scaled(double mu, double sigma) {
    return mu + sigma * slam_randn();
}

/* =========================================================================
 * L2: Initialization
 * ========================================================================= */

int slam_fastslam_init(slam_particle_t **particles_out,
                       int num_particles,
                       const slam_pose2d_t *init_pose,
                       const slam_config_t *config) {
    if (!particles_out || !init_pose || !config || num_particles <= 0)
        return SLAM_ERR_INVALID_PARAM;

    slam_particle_t *p = (slam_particle_t *)calloc(num_particles,
                                                    sizeof(slam_particle_t));
    if (!p) return SLAM_ERR_MEMORY;

    double init_weight = 1.0 / num_particles;
    double init_logw = log(init_weight);
    int lm_cap = config->max_landmarks;

    for (int i = 0; i < num_particles; i++) {
        p[i].pose       = *init_pose;
        p[i].weight     = init_weight;
        p[i].log_weight = init_logw;
        p[i].num_landmarks = 0;
        p[i].lm_capacity   = lm_cap;
        p[i].landmarks = (slam_landmark2d_t *)calloc(lm_cap,
                                                      sizeof(slam_landmark2d_t));
        p[i].lm_covariances = (double *)calloc(lm_cap * 4, sizeof(double));
        if (!p[i].landmarks || !p[i].lm_covariances) {
            /* Cleanup on failure */
            for (int j = 0; j <= i; j++) {
                free(p[j].landmarks);
                free(p[j].lm_covariances);
            }
            free(p);
            return SLAM_ERR_MEMORY;
        }
        p[i].traj_len = 0;
    }

    *particles_out = p;
    return SLAM_OK;
}

void slam_fastslam_free(slam_particle_t *particles, int num_particles) {
    if (particles) {
        for (int i = 0; i < num_particles; i++) {
            free(particles[i].landmarks);
            free(particles[i].lm_covariances);
        }
        free(particles);
    }
}

/* =========================================================================
 * L5: FastSLAM 1.0 — Motion Model Proposal
 * ========================================================================= */

int slam_fastslam1_sample_pose(slam_particle_t *particle,
                                const slam_velocity_t *vel,
                                double sigma_v,
                                double sigma_omega) {
    if (!particle || !vel) return SLAM_ERR_NULL_PTR;

    /* Sample noisy velocities */
    double v_hat = vel->v + slam_randn_scaled(0, sigma_v);
    double w_hat = vel->omega + slam_randn_scaled(0, sigma_omega);
    double gamma_hat = slam_randn_scaled(0, sigma_omega * 0.5);
    double dt = vel->dt;

    double x = particle->pose.x;
    double y = particle->pose.y;
    double theta = particle->pose.theta;

    if (fabs(w_hat) < 1e-8) {
        /* Straight-line motion */
        x += v_hat * dt * cos(theta);
        y += v_hat * dt * sin(theta);
        theta += gamma_hat * dt;
    } else {
        /* Arc motion */
        double vw = v_hat / w_hat;
        x += vw * (sin(theta + w_hat * dt) - sin(theta));
        y += vw * (cos(theta) - cos(theta + w_hat * dt));
        theta += w_hat * dt + gamma_hat * dt;
    }

    particle->pose.x     = x;
    particle->pose.y     = y;
    particle->pose.theta = slam_normalize_angle(theta);

    /* Record trajectory */
    if (particle->traj_len < 1024) {
        particle->trajectory_x[particle->traj_len] = x;
        particle->trajectory_y[particle->traj_len] = y;
        particle->traj_len++;
    }

    return SLAM_OK;
}

/* =========================================================================
 * L5: FastSLAM 1.0 — Landmark Update
 * ========================================================================= */

/**
 * @brief EKF update for a single landmark (2D)
 *
 * This is the key efficiency of FastSLAM: each landmark is updated
 * independently. The innovation is 2D, covariance is 2×2.
 *
 * Complexity: O(1) per landmark, vs O(N²) for full EKF update.
 */
int slam_fastslam1_update_landmark(slam_particle_t *particle,
                                    const slam_obs_rb_t *obs,
                                    int lm_idx) {
    if (!particle || !obs) return SLAM_ERR_NULL_PTR;
    if (lm_idx < 0 || lm_idx >= particle->num_landmarks)
        return SLAM_ERR_NO_ASSOC;

    slam_landmark2d_t *lm = &particle->landmarks[lm_idx];
    double *cov = &particle->lm_covariances[lm_idx * 4];

    double rx = particle->pose.x;
    double ry = particle->pose.y;
    double rt = particle->pose.theta;

    /* 1. Predicted measurement */
    double dx = lm->x - rx, dy = lm->y - ry;
    double q = dx*dx + dy*dy;
    if (q < 1e-12) return SLAM_ERR_SINGULAR;
    double sqrt_q = sqrt(q);

    double z_hat_r = sqrt_q;
    double z_hat_b = slam_normalize_angle(atan2(dy, dx) - rt);

    double nu[2] = {
        obs->range - z_hat_r,
        slam_angle_diff(obs->bearing, z_hat_b)
    };

    /* 2. Observation Jacobian w.r.t. landmark position (2×2) */
    double H[4] = {
         dx / sqrt_q,  dy / sqrt_q,
        -dy / q,       dx / q
    };

    /* 3. Innovation covariance S = H·Σ_lm·H^T + Q */
    double Q[4] = {obs->range_std * obs->range_std, 0,
                   0, obs->bearing_std * obs->bearing_std};

    double H_Sigma[4];
    slam_matmul(H, cov, 2, 2, 2, H_Sigma);
    double S[4];
    slam_matmul_A_BT(H_Sigma, H, 2, 2, 2, S);
    S[0] += Q[0]; S[3] += Q[3];

    /* 4. Kalman gain K = Σ·H^T·S^{-1} (2×2) */
    double Sigma_HT[4] = {
        cov[0]*H[0] + cov[1]*H[1],
        cov[2]*H[0] + cov[3]*H[1],
        cov[0]*H[2] + cov[1]*H[3],
        cov[2]*H[2] + cov[3]*H[3]
    };

    double S_inv[4];
    if (!slam_inv2x2(S, S_inv)) return SLAM_ERR_SINGULAR;

    double K[4];
    slam_matmul(Sigma_HT, S_inv, 2, 2, 2, K);

    /* 5. Update landmark state */
    lm->x += K[0]*nu[0] + K[2]*nu[1];
    lm->y += K[1]*nu[0] + K[3]*nu[1];

    /* 6. Joseph-form covariance update: Σ = (I − K·H)·Σ */
    double KH[4];
    slam_matmul(K, H, 2, 2, 2, KH);
    double new_cov[4] = {
        cov[0] - (KH[0]*cov[0] + KH[1]*cov[2]),  /* KH[0,0]·Σ[0,0] + KH[0,1]·Σ[0,1] */
        cov[1] - (KH[0]*cov[1] + KH[1]*cov[3]),
        cov[2] - (KH[2]*cov[0] + KH[3]*cov[2]),
        cov[3] - (KH[2]*cov[1] + KH[3]*cov[3])
    };
    memcpy(cov, new_cov, 4 * sizeof(double));

    return SLAM_OK;
}

/* =========================================================================
 * L5: Importance Weight Update
 * ========================================================================= */

int slam_fastslam1_update_weight(slam_particle_t *particle,
                                  const slam_obs_rb_t *obs,
                                  int lm_idx) {
    if (!particle || !obs) return SLAM_ERR_NULL_PTR;
    if (lm_idx < 0 || lm_idx >= particle->num_landmarks)
        return SLAM_ERR_NO_ASSOC;

    slam_landmark2d_t *lm = &particle->landmarks[lm_idx];
    double *cov = &particle->lm_covariances[lm_idx * 4];

    double rx = particle->pose.x;
    double ry = particle->pose.y;
    double rt = particle->pose.theta;

    double dx = lm->x - rx, dy = lm->y - ry;
    double q = dx*dx + dy*dy;
    if (q < 1e-12) return SLAM_ERR_SINGULAR;
    double sqrt_q = sqrt(q);

    double z_hat_r = sqrt_q;
    double z_hat_b = slam_normalize_angle(atan2(dy, dx) - rt);

    double nu[2] = {
        obs->range - z_hat_r,
        slam_angle_diff(obs->bearing, z_hat_b)
    };

    double H[4] = {
         dx / sqrt_q,  dy / sqrt_q,
        -dy / q,       dx / q
    };

    double Q[4] = {obs->range_std * obs->range_std, 0,
                   0, obs->bearing_std * obs->bearing_std};

    /* S = H·Σ_lm·H^T + Q */
    double H_Sigma[4];
    slam_matmul(H, cov, 2, 2, 2, H_Sigma);
    double S[4];
    slam_matmul_A_BT(H_Sigma, H, 2, 2, 2, S);
    S[0] += Q[0]; S[3] += Q[3];

    double S_inv[4];
    if (!slam_inv2x2(S, S_inv)) return SLAM_ERR_SINGULAR;

    /* log w += −0.5·ν^T·S^{-1}·ν − 0.5·log|2π·S| */
    double d2 = nu[0]*S_inv[0]*nu[0] + nu[0]*S_inv[1]*nu[1]
              + nu[1]*S_inv[2]*nu[0] + nu[1]*S_inv[3]*nu[1];

    double det_S = S[0]*S[3] - S[1]*S[2];
    if (det_S <= 0) det_S = 1e-12;

    particle->log_weight += -0.5 * d2 - 0.5 * log(2.0 * M_PI * 2.0 * M_PI * det_S);

    return SLAM_OK;
}

/* =========================================================================
 * L5: FastSLAM 2.0 — Improved Proposal
 * ========================================================================= */

int slam_fastslam2_sample_pose(slam_particle_t *particle,
                                const slam_velocity_t *vel,
                                const slam_obs_rb_t *obs,
                                int lm_idx,
                                double sigma_v,
                                double sigma_omega) {
    if (!particle || !vel || !obs) return SLAM_ERR_NULL_PTR;

    /* If no landmark reference, fall back to motion model */
    if (lm_idx < 0 || lm_idx >= particle->num_landmarks) {
        return slam_fastslam1_sample_pose(particle, vel, sigma_v, sigma_omega);
    }

    slam_landmark2d_t *lm = &particle->landmarks[lm_idx];

    /* 1. Predict using motion model (mean) */
    double v = vel->v, w = vel->omega, dt = vel->dt;
    double x0 = particle->pose.x;
    double y0 = particle->pose.y;
    double t0 = particle->pose.theta;

    double x_pred, y_pred, t_pred;
    if (fabs(w) < 1e-8) {
        x_pred = x0 + v * dt * cos(t0);
        y_pred = y0 + v * dt * sin(t0);
        t_pred = t0;
    } else {
        double vw = v / w;
        x_pred = x0 + vw * (sin(t0 + w*dt) - sin(t0));
        y_pred = y0 + vw * (cos(t0) - cos(t0 + w*dt));
        t_pred = t0 + w * dt;
    }
    t_pred = slam_normalize_angle(t_pred);

    /* 2. Motion covariance (3×3) */
    double Sigma_pred[9] = {0};
    Sigma_pred[0] = sigma_v * sigma_v;
    Sigma_pred[4] = sigma_v * sigma_v;
    Sigma_pred[8] = sigma_omega * sigma_omega;

    /* 3. Observation Jacobian at predicted pose (2×3) */
    double dx = lm->x - x_pred, dy = lm->y - y_pred;
    double q = dx*dx + dy*dy;
    if (q < 1e-10) {
        /* Too close, fall back to motion model */
        particle->pose.x = x_pred;
        particle->pose.y = y_pred;
        particle->pose.theta = t_pred;
        return SLAM_OK;
    }
    double sqrt_q = sqrt(q);

    double H[6] = {-dx/sqrt_q, -dy/sqrt_q, 0,
                    dy/q,      -dx/q,     -1};

    /* 4. Measurement noise */
    double Q[4] = {obs->range_std * obs->range_std, 0,
                   0, obs->bearing_std * obs->bearing_std};

    /* 5. Proposal covariance: (Sigma_pred^{-1} + H^T Q^{-1} H)^{-1}
     * For this simplified implementation, we use a diagonal approximation.
     * Sigma_inv_prop = H^T Q^{-1} H (diagonal approx) */
    double Sigma_inv[9] = {0};
    Sigma_inv[0] = H[0]*H[0]/Q[0] + H[3]*H[3]/Q[3] + 1.0/Sigma_pred[0];
    Sigma_inv[4] = H[1]*H[1]/Q[0] + H[4]*H[4]/Q[3] + 1.0/Sigma_pred[4];
    Sigma_inv[8] = H[2]*H[2]/Q[0] + H[5]*H[5]/Q[3] + 1.0/Sigma_pred[8];
    Sigma_inv[1] = H[0]*H[1]/Q[0] + H[3]*H[4]/Q[3];
    Sigma_inv[3] = Sigma_inv[1];
    Sigma_inv[2] = H[0]*H[2]/Q[0] + H[3]*H[5]/Q[3];
    Sigma_inv[6] = Sigma_inv[2];
    Sigma_inv[5] = H[1]*H[2]/Q[0] + H[4]*H[5]/Q[3];
    Sigma_inv[7] = Sigma_inv[5];

    /* 6. Predicted measurement */
    double z_hat_r = sqrt_q;
    double z_hat_b = slam_normalize_angle(atan2(dy, dx) - t_pred);
    double nu_z[2] = {obs->range - z_hat_r,
                      slam_angle_diff(obs->bearing, z_hat_b)};

    /* 7. Proposal mean shift = Σ_prop·H^T·Q^{-1}·ν */
    double dx_shift = (H[0]*nu_z[0]/Q[0] + H[3]*nu_z[1]/Q[3]) / Sigma_inv[0];
    double dy_shift = (H[1]*nu_z[0]/Q[0] + H[4]*nu_z[1]/Q[3]) / Sigma_inv[4];
    double dt_shift = (H[2]*nu_z[0]/Q[0] + H[5]*nu_z[1]/Q[3]) / Sigma_inv[8];

    /* 8. Sample from proposal N(μ_prop, Σ_prop) */
    double sigma_prop_x = sqrt(1.0 / fmax(Sigma_inv[0], 1e-6));
    double sigma_prop_y = sqrt(1.0 / fmax(Sigma_inv[4], 1e-6));
    double sigma_prop_t = sqrt(1.0 / fmax(Sigma_inv[8], 1e-6));

    particle->pose.x     = x_pred + dx_shift + slam_randn_scaled(0, sigma_prop_x);
    particle->pose.y     = y_pred + dy_shift + slam_randn_scaled(0, sigma_prop_y);
    particle->pose.theta = slam_normalize_angle(
        t_pred + dt_shift + slam_randn_scaled(0, sigma_prop_t));

    /* Record trajectory */
    if (particle->traj_len < 1024) {
        particle->trajectory_x[particle->traj_len] = particle->pose.x;
        particle->trajectory_y[particle->traj_len] = particle->pose.y;
        particle->traj_len++;
    }

    return SLAM_OK;
}

/* =========================================================================
 * L5: Resampling
 * ========================================================================= */

int slam_fastslam_resample(slam_particle_t **particles,
                            int num_particles,
                            double neff_thresh,
                            double *neff_out) {
    if (!particles || !*particles || num_particles <= 0)
        return SLAM_ERR_INVALID_PARAM;

    slam_particle_t *p = *particles;

    /* 1. Normalize weights (convert from log-space) */
    double max_logw = p[0].log_weight;
    for (int i = 1; i < num_particles; i++) {
        if (p[i].log_weight > max_logw) max_logw = p[i].log_weight;
    }

    double sum_w = 0.0;
    for (int i = 0; i < num_particles; i++) {
        p[i].weight = exp(p[i].log_weight - max_logw);
        sum_w += p[i].weight;
    }
    for (int i = 0; i < num_particles; i++) {
        p[i].weight /= sum_w;
    }

    /* 2. Effective sample size */
    double neff = 0.0;
    for (int i = 0; i < num_particles; i++) {
        neff += p[i].weight * p[i].weight;
    }
    neff = 1.0 / neff;
    if (neff_out) *neff_out = neff;

    /* 3. Resample only if needed */
    if (neff >= neff_thresh) return SLAM_OK;

    /* 4. Systematic resampling */
    slam_particle_t *new_p = (slam_particle_t *)calloc(num_particles,
                                                        sizeof(slam_particle_t));
    if (!new_p) return SLAM_ERR_MEMORY;

    /* Allocate landmark arrays for new particles */
    for (int i = 0; i < num_particles; i++) {
        new_p[i].landmarks = (slam_landmark2d_t *)calloc(
            p[0].lm_capacity, sizeof(slam_landmark2d_t));
        new_p[i].lm_covariances = (double *)calloc(
            p[0].lm_capacity * 4, sizeof(double));
        new_p[i].lm_capacity = p[0].lm_capacity;
        if (!new_p[i].landmarks || !new_p[i].lm_covariances) {
            for (int j = 0; j <= i; j++) {
                free(new_p[j].landmarks);
                free(new_p[j].lm_covariances);
            }
            free(new_p);
            return SLAM_ERR_MEMORY;
        }
    }

    /* Systematic resample */
    double inv_M = 1.0 / num_particles;
    double r = uniform_random() * inv_M;
    double c = p[0].weight;
    int idx = 0;

    for (int m = 0; m < num_particles; m++) {
        double U = r + m * inv_M;
        while (U > c && idx < num_particles - 1) {
            idx++;
            c += p[idx].weight;
        }
        /* Copy particle idx */
        new_p[m].pose          = p[idx].pose;
        new_p[m].weight        = inv_M;
        new_p[m].log_weight    = log(inv_M);
        new_p[m].num_landmarks = p[idx].num_landmarks;
        memcpy(new_p[m].landmarks, p[idx].landmarks,
               p[idx].num_landmarks * sizeof(slam_landmark2d_t));
        memcpy(new_p[m].lm_covariances, p[idx].lm_covariances,
               p[idx].num_landmarks * 4 * sizeof(double));
        memcpy(new_p[m].trajectory_x, p[idx].trajectory_x,
               p[idx].traj_len * sizeof(double));
        memcpy(new_p[m].trajectory_y, p[idx].trajectory_y,
               p[idx].traj_len * sizeof(double));
        new_p[m].traj_len = p[idx].traj_len;
    }

    /* Free old particles */
    slam_fastslam_free(p, num_particles);

    *particles = new_p;
    return SLAM_OK;
}

int slam_fastslam_low_variance_resample(slam_particle_t *particles,
                                         int num_particles) {
    if (!particles || num_particles <= 0) return SLAM_ERR_INVALID_PARAM;

    /* Low-variance resampling (Probabilistic Robotics, Table 4.4) */
    /* This is a simpler O(M) resampler that modifies in-place */
    /* First convert log-weights to weights */
    double max_lw = particles[0].log_weight;
    for (int i = 1; i < num_particles; i++)
        if (particles[i].log_weight > max_lw)
            max_lw = particles[i].log_weight;

    double *cw = (double *)malloc(num_particles * sizeof(double));
    if (!cw) return SLAM_ERR_MEMORY;

    cw[0] = exp(particles[0].log_weight - max_lw);
    for (int i = 1; i < num_particles; i++) {
        cw[i] = cw[i-1] + exp(particles[i].log_weight - max_lw);
    }
    /* Normalize */
    double total = cw[num_particles-1];
    for (int i = 0; i < num_particles; i++) cw[i] /= total;

    /* Systematic resampling, in-place */
    double *new_logw = (double *)malloc(num_particles * sizeof(double));
    slam_pose2d_t *new_poses = (slam_pose2d_t *)malloc(
        num_particles * sizeof(slam_pose2d_t));
    if (!new_logw || !new_poses) {
        free(cw); free(new_logw); free(new_poses);
        return SLAM_ERR_MEMORY;
    }

    double r = uniform_random() / num_particles;
    int j = 0;
    for (int i = 0; i < num_particles; i++) {
        double U = r + (double)i / num_particles;
        while (U > cw[j] && j < num_particles - 1) j++;
        new_poses[i] = particles[j].pose;
        new_logw[i] = log(1.0 / num_particles);
    }

    for (int i = 0; i < num_particles; i++) {
        particles[i].pose       = new_poses[i];
        particles[i].log_weight = new_logw[i];
        particles[i].weight     = 1.0 / num_particles;
    }

    free(cw); free(new_logw); free(new_poses);
    return SLAM_OK;
}

/* =========================================================================
 * L6: FastSLAM Full Step
 * ========================================================================= */

int slam_fastslam_step(slam_particle_t *particles,
                       int num_particles,
                       const slam_velocity_t *vel,
                       const slam_obs_rb_t *observations,
                       int num_obs,
                       const slam_config_t *config) {
    if (!particles || !config) return SLAM_ERR_NULL_PTR;

    for (int k = 0; k < num_particles; k++) {
        /* 1. Sample pose from motion model (FastSLAM 1.0 proposal) */
        slam_fastslam1_sample_pose(&particles[k], vel,
                                    config->sigma_v, config->sigma_omega);

        /* 2. Process each observation */
        if (observations && num_obs > 0) {
            for (int i = 0; i < num_obs; i++) {
                const slam_obs_rb_t *obs = &observations[i];

                /* Find best matching landmark (NN) */
                int best_j = -1;
                double best_d2 = 1e99;

                for (int j = 0; j < particles[k].num_landmarks; j++) {
                    slam_landmark2d_t *lm = &particles[k].landmarks[j];
                    double dx = lm->x - particles[k].pose.x;
                    double dy = lm->y - particles[k].pose.y;
                    double qq = dx*dx + dy*dy;
                    if (qq < 1e-12) continue;
                    double sq = sqrt(qq);

                    double zr = sq;
                    double zb = slam_normalize_angle(
                        atan2(dy, dx) - particles[k].pose.theta);

                    double nu_r = obs->range - zr;
                    double nu_b = slam_angle_diff(obs->bearing, zb);

                    double d2 = (nu_r*nu_r)/(obs->range_std*obs->range_std)
                              + (nu_b*nu_b)/(obs->bearing_std*obs->bearing_std);

                    if (d2 < config->mahalanobis_gate && d2 < best_d2) {
                        best_d2 = d2;
                        best_j = j;
                    }
                }

                if (best_j >= 0) {
                    /* Update matched landmark */
                    slam_fastslam1_update_landmark(&particles[k], obs, best_j);
                    slam_fastslam1_update_weight(&particles[k], obs, best_j);
                } else if (particles[k].num_landmarks < particles[k].lm_capacity) {
                    /* Initialize new landmark */
                    int idx = particles[k].num_landmarks;
                    double rx = particles[k].pose.x;
                    double ry = particles[k].pose.y;
                    double rt = particles[k].pose.theta;

                    particles[k].landmarks[idx].x =
                        rx + obs->range * cos(obs->bearing + rt);
                    particles[k].landmarks[idx].y =
                        ry + obs->range * sin(obs->bearing + rt);
                    particles[k].landmarks[idx].id = idx;
                    particles[k].landmarks[idx].is_active = true;

                    /* Initialize covariance */
                    double *cov = &particles[k].lm_covariances[idx * 4];
                    cov[0] = obs->range_std * obs->range_std;
                    cov[3] = obs->bearing_std * obs->bearing_std;
                    cov[1] = cov[2] = 0.0;

                    particles[k].num_landmarks++;
                }
            }
        }
    }

    /* 3. Resample */
    double neff;
    slam_fastslam_resample(&particles, num_particles,
                            num_particles / 2.0, &neff);

    return SLAM_OK;
}

/* =========================================================================
 * L6: Query Best Estimate
 * ========================================================================= */

int slam_fastslam_best_pose(const slam_particle_t *particles,
                             int num_particles,
                             slam_pose2d_t *best_pose) {
    if (!particles || !best_pose) return SLAM_ERR_NULL_PTR;

    int best_i = 0;
    double best_w = particles[0].log_weight;
    for (int i = 1; i < num_particles; i++) {
        if (particles[i].log_weight > best_w) {
            best_w = particles[i].log_weight;
            best_i = i;
        }
    }

    *best_pose = particles[best_i].pose;
    return SLAM_OK;
}

int slam_fastslam_extract_map(const slam_particle_t *particles,
                               int num_particles,
                               slam_map2d_t *map_out) {
    if (!particles || !map_out) return SLAM_ERR_NULL_PTR;

    /* Find best particle */
    int best_i = 0;
    double best_w = particles[0].log_weight;
    for (int i = 1; i < num_particles; i++) {
        if (particles[i].log_weight > best_w) {
            best_w = particles[i].log_weight;
            best_i = i;
        }
    }

    /* Copy landmarks from best particle */
    const slam_particle_t *best = &particles[best_i];
    for (int j = 0; j < best->num_landmarks; j++) {
        int idx;
        slam_map2d_add_landmark(map_out, &best->landmarks[j], &idx);
    }

    return SLAM_OK;
}
