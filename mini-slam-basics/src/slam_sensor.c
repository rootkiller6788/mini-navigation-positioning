/**
 * @file    slam_sensor.c
 * @brief   Sensor models: motion models, observation models, LiDAR
 *
 * Implements all sensor models used in SLAM.
 *
 * Reference:
 *   Thrun, Burgard & Fox (2005), Chapters 5, 6, 7, 9.
 *   Bailey & Durrant-Whyte (2006) "SLAM Part II"
 */

#include "slam_sensor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* External */
extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);
extern double slam_randn(void);
extern double slam_randn_scaled(double mu, double sigma);

/* =========================================================================
 * L2: Velocity Motion Model
 * ========================================================================= */

int slam_motion_model_velocity(const slam_pose2d_t *pose,
                                const slam_velocity_t *vel,
                                slam_pose2d_t *new_pose) {
    if (!pose || !vel || !new_pose) return SLAM_ERR_NULL_PTR;

    double v = vel->v, w = vel->omega, dt = vel->dt;
    double x0 = pose->x, y0 = pose->y, t0 = pose->theta;

    if (fabs(w) < 1e-8) {
        new_pose->x     = x0 + v * dt * cos(t0);
        new_pose->y     = y0 + v * dt * sin(t0);
        new_pose->theta = t0;
    } else {
        double vw = v / w;
        new_pose->x     = x0 + vw * (sin(t0 + w*dt) - sin(t0));
        new_pose->y     = y0 + vw * (cos(t0) - cos(t0 + w*dt));
        new_pose->theta = slam_normalize_angle(t0 + w * dt);
    }

    return SLAM_OK;
}

int slam_motion_model_velocity_noisy(const slam_pose2d_t *pose,
                                      const slam_velocity_t *vel,
                                      const double alpha[6],
                                      slam_pose2d_t *new_pose) {
    if (!pose || !vel || !alpha || !new_pose) return SLAM_ERR_NULL_PTR;

    double v = vel->v, w = vel->omega, dt = vel->dt;

    /* Sample noise */
    double v_hat  = v + slam_randn_scaled(0, alpha[0]*fabs(v) + alpha[1]*fabs(w));
    double w_hat  = w + slam_randn_scaled(0, alpha[2]*fabs(v) + alpha[3]*fabs(w));
    double gamma  = slam_randn_scaled(0, alpha[4]*fabs(v) + alpha[5]*fabs(w));

    double x0 = pose->x, y0 = pose->y, t0 = pose->theta;

    if (fabs(w_hat) < 1e-8) {
        new_pose->x     = x0 + v_hat * dt * cos(t0);
        new_pose->y     = y0 + v_hat * dt * sin(t0);
        new_pose->theta = slam_normalize_angle(t0 + gamma * dt);
    } else {
        double vw = v_hat / w_hat;
        new_pose->x     = x0 + vw * (sin(t0 + w_hat*dt) - sin(t0));
        new_pose->y     = y0 + vw * (cos(t0) - cos(t0 + w_hat*dt));
        new_pose->theta = slam_normalize_angle(t0 + w_hat*dt + gamma*dt);
    }

    return SLAM_OK;
}

int slam_motion_jacobian_velocity(const slam_pose2d_t *pose,
                                   const slam_velocity_t *vel,
                                   double G[9]) {
    if (!pose || !vel || !G) return SLAM_ERR_NULL_PTR;

    double v = vel->v, w = vel->omega, dt = vel->dt;
    double t0 = pose->theta;

    /* Identity */
    memset(G, 0, 9 * sizeof(double));
    G[0] = 1.0; G[4] = 1.0; G[8] = 1.0;

    if (fabs(w) < 1e-8) {
        G[2] = -v * dt * sin(t0);
        G[5] =  v * dt * cos(t0);
    } else {
        G[2] = -v/w * (cos(t0) - cos(t0 + w*dt));
        G[5] = -v/w * (sin(t0 + w*dt) - sin(t0));
    }

    return SLAM_OK;
}

/* =========================================================================
 * L2: Odometry Motion Model
 * ========================================================================= */

int slam_motion_model_odometry(const slam_pose2d_t *pose,
                                const slam_odometry_t *odom,
                                slam_pose2d_t *new_pose) {
    if (!pose || !odom || !new_pose) return SLAM_ERR_NULL_PTR;

    double dr1 = odom->delta_rot1;
    double dt  = odom->delta_trans;
    double dr2 = odom->delta_rot2;
    double t0  = pose->theta;

    new_pose->x     = pose->x + dt * cos(t0 + dr1);
    new_pose->y     = pose->y + dt * sin(t0 + dr1);
    new_pose->theta = slam_normalize_angle(t0 + dr1 + dr2);

    return SLAM_OK;
}

int slam_motion_model_odometry_noisy(const slam_pose2d_t *pose,
                                      const slam_odometry_t *odom,
                                      const double alpha[4],
                                      slam_pose2d_t *new_pose) {
    if (!pose || !odom || !alpha || !new_pose) return SLAM_ERR_NULL_PTR;

    double dr1 = odom->delta_rot1;
    double dt  = odom->delta_trans;
    double dr2 = odom->delta_rot2;

    /* Sample noisy odometry */
    double dr1_hat = dr1 - slam_randn_scaled(0,
        alpha[0]*fabs(dr1) + alpha[1]*fabs(dt));
    double dt_hat  = dt - slam_randn_scaled(0,
        alpha[2]*fabs(dt) + alpha[3]*(fabs(dr1)+fabs(dr2)));
    double dr2_hat = dr2 - slam_randn_scaled(0,
        alpha[0]*fabs(dr2) + alpha[1]*fabs(dt));

    double t0 = pose->theta;

    new_pose->x     = pose->x + dt_hat * cos(t0 + dr1_hat);
    new_pose->y     = pose->y + dt_hat * sin(t0 + dr1_hat);
    new_pose->theta = slam_normalize_angle(t0 + dr1_hat + dr2_hat);

    return SLAM_OK;
}

/* =========================================================================
 * L3: Range-Bearing Observation Model
 * ========================================================================= */

int slam_obs_model_rb(const slam_pose2d_t *robot_pose,
                      const slam_landmark2d_t *landmark,
                      double pred[2]) {
    if (!robot_pose || !landmark || !pred) return SLAM_ERR_NULL_PTR;

    double dx = landmark->x - robot_pose->x;
    double dy = landmark->y - robot_pose->y;
    double q = dx*dx + dy*dy;

    pred[0] = sqrt(q);
    pred[1] = slam_normalize_angle(atan2(dy, dx) - robot_pose->theta);

    return SLAM_OK;
}

int slam_obs_jacobian_rb_pose(const slam_pose2d_t *robot_pose,
                               const slam_landmark2d_t *landmark,
                               double H_pose[6]) {
    if (!robot_pose || !landmark || !H_pose) return SLAM_ERR_NULL_PTR;

    double dx = landmark->x - robot_pose->x;
    double dy = landmark->y - robot_pose->y;
    double q = dx*dx + dy*dy;

    if (q < 1e-12) {
        memset(H_pose, 0, 6 * sizeof(double));
        return SLAM_OK;
    }

    double sqrt_q = sqrt(q);

    H_pose[0] = -dx / sqrt_q;  /* ∂r̂/∂x */
    H_pose[1] = -dy / sqrt_q;  /* ∂r̂/∂y */
    H_pose[2] =  0.0;          /* ∂r̂/∂θ */
    H_pose[3] =  dy / q;       /* ∂φ̂/∂x */
    H_pose[4] = -dx / q;       /* ∂φ̂/∂y */
    H_pose[5] = -1.0;          /* ∂φ̂/∂θ */

    return SLAM_OK;
}

int slam_obs_jacobian_rb_landmark(const slam_pose2d_t *robot_pose,
                                   const slam_landmark2d_t *landmark,
                                   double H_lm[4]) {
    if (!robot_pose || !landmark || !H_lm) return SLAM_ERR_NULL_PTR;

    double dx = landmark->x - robot_pose->x;
    double dy = landmark->y - robot_pose->y;
    double q = dx*dx + dy*dy;

    if (q < 1e-12) {
        memset(H_lm, 0, 4 * sizeof(double));
        return SLAM_OK;
    }

    double sqrt_q = sqrt(q);

    H_lm[0] =  dx / sqrt_q;  /* ∂r̂/∂m_x */
    H_lm[1] =  dy / sqrt_q;  /* ∂r̂/∂m_y */
    H_lm[2] = -dy / q;       /* ∂φ̂/∂m_x */
    H_lm[3] =  dx / q;       /* ∂φ̂/∂m_y */

    return SLAM_OK;
}

int slam_obs_inverse_rb(const slam_pose2d_t *robot_pose,
                         const slam_obs_rb_t *obs,
                         slam_landmark2d_t *landmark) {
    if (!robot_pose || !obs || !landmark) return SLAM_ERR_NULL_PTR;

    double r = obs->range;
    double phi = obs->bearing;
    double theta = robot_pose->theta;

    double cos_angle = cos(phi + theta);
    double sin_angle = sin(phi + theta);

    landmark->x = robot_pose->x + r * cos_angle;
    landmark->y = robot_pose->y + r * sin_angle;

    return SLAM_OK;
}

int slam_obs_inverse_jacobian_rb(const slam_pose2d_t *pose,
                                  const slam_obs_rb_t *obs,
                                  double G_rl[6],
                                  double G_z[4]) {
    if (!pose || !obs || !G_rl || !G_z) return SLAM_ERR_NULL_PTR;

    double r = obs->range;
    double phi = obs->bearing;
    double theta = pose->theta;
    double cp = cos(phi + theta);
    double sp = sin(phi + theta);

    /* G_rl = ∂g/∂(x,y,θ) */
    G_rl[0] = 1.0;  G_rl[1] = 0.0;  G_rl[2] = -r * sp;
    G_rl[3] = 0.0;  G_rl[4] = 1.0;  G_rl[5] =  r * cp;

    /* G_z = ∂g/∂(r,φ) */
    G_z[0] = cp;  G_z[1] = -r * sp;
    G_z[2] = sp;  G_z[3] =  r * cp;

    return SLAM_OK;
}

/* =========================================================================
 * L6: LiDAR Beam Model and Occupancy Grid
 * ========================================================================= */

int slam_lidar_likelihood_field(const slam_lidar_scan_t *scan,
                                 const slam_pose2d_t *robot_pose,
                                 const slam_occgrid_t *grid,
                                 double sigma_hit,
                                 double p_rand,
                                 double *log_likelihood) {
    if (!scan || !robot_pose || !grid || !log_likelihood)
        return SLAM_ERR_NULL_PTR;

    double log_lik = 0.0;
    double var = sigma_hit * sigma_hit;
    double norm_factor = 1.0 / sqrt(2.0 * M_PI * var);

    for (int k = 0; k < scan->num_ranges; k++) {
        double r = scan->ranges[k];
        if (r < scan->range_min || r > scan->range_max) continue;

        double angle = scan->angle_min + k * scan->angle_inc;
        double world_angle = angle + robot_pose->theta;

        /* Beam endpoint in world frame */
        double ex = robot_pose->x + r * cos(world_angle);
        double ey = robot_pose->y + r * sin(world_angle);

        /* Find distance to nearest obstacle in grid */
        /* Grid cell containing endpoint */
        int cx = (int)((ex - grid->origin_x) / grid->resolution);
        int cy = (int)((ey - grid->origin_y) / grid->resolution);

        /* Check neighbors for nearest occupied cell */
        double min_dist = 1e99;
        bool found = false;
        int search_radius = (int)ceil(3.0 * sigma_hit / grid->resolution) + 1;
        for (int dx = -search_radius; dx <= search_radius; dx++) {
            for (int dy = -search_radius; dy <= search_radius; dy++) {
                int nx = cx + dx, ny = cy + dy;
                if (nx < 0 || nx >= grid->width || ny < 0 || ny >= grid->height)
                    continue;
                double lo = grid->log_odds[ny * grid->width + nx];
                /* Occupied if log-odds > 0 (p > 0.5) */
                if (lo > 0.5) {
                    double wx = grid->origin_x + (nx + 0.5) * grid->resolution;
                    double wy = grid->origin_y + (ny + 0.5) * grid->resolution;
                    double dist = sqrt((ex - wx)*(ex - wx) + (ey - wy)*(ey - wy));
                    if (dist < min_dist) {
                        min_dist = dist;
                        found = true;
                    }
                }
            }
        }

        double p_hit;
        if (found) {
            p_hit = norm_factor * exp(-min_dist*min_dist / (2.0*var));
        } else {
            p_hit = norm_factor * exp(-(3.0*sigma_hit)*(3.0*sigma_hit)/(2.0*var))
                    * 0.1;
        }

        double p = 0.5 * p_hit + 0.5 * p_rand;
        if (p < 1e-100) p = 1e-100;
        log_lik += log(p);
    }

    *log_likelihood = log_lik;
    return SLAM_OK;
}

/**
 * @brief Bresenham line rasterization for ray casting
 *
 * Returns grid cells along the ray from start to end.
 * This is used for occupancy grid updates from LiDAR scans.
 */
static int bresenham_line(int x0, int y0, int x1, int y1,
                          int *out_x, int *out_y, int max_points) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2;
    int count = 0;

    while (count < max_points) {
        out_x[count] = x0;
        out_y[count] = y0;
        count++;
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
    return count;
}

int slam_occgrid_update_scan(slam_occgrid_t *grid,
                              const slam_pose2d_t *robot_pose,
                              const slam_lidar_scan_t *scan) {
    if (!grid || !robot_pose || !scan) return SLAM_ERR_NULL_PTR;

    int *bx = (int *)malloc(grid->width * grid->height * sizeof(int));
    int *by = (int *)malloc(grid->width * grid->height * sizeof(int));
    if (!bx || !by) { free(bx); free(by); return SLAM_ERR_MEMORY; }

    /* Robot pose in grid coordinates */
    int rx = (int)((robot_pose->x - grid->origin_x) / grid->resolution);
    int ry = (int)((robot_pose->y - grid->origin_y) / grid->resolution);

    for (int k = 0; k < scan->num_ranges; k++) {
        double r = scan->ranges[k];
        if (r < scan->range_min || r > scan->range_max) continue;

        double angle = scan->angle_min + k * scan->angle_inc;
        double world_angle = angle + robot_pose->theta;

        /* Beam endpoint */
        double ex = robot_pose->x + r * cos(world_angle);
        double ey = robot_pose->y + r * sin(world_angle);
        int end_x = (int)((ex - grid->origin_x) / grid->resolution);
        int end_y = (int)((ey - grid->origin_y) / grid->resolution);

        /* Clip to grid */
        if (end_x < 0 || end_x >= grid->width ||
            end_y < 0 || end_y >= grid->height) continue;
        if (rx < 0 || rx >= grid->width || ry < 0 || ry >= grid->height)
            continue;

        /* Mark cells along beam as free */
        int n_cells = bresenham_line(rx, ry, end_x, end_y, bx, by,
                                      grid->width * grid->height);
        for (int c = 0; c < n_cells - 1; c++) {
            if (bx[c] >= 0 && bx[c] < grid->width &&
                by[c] >= 0 && by[c] < grid->height) {
                double *lo = &grid->log_odds[by[c] * grid->width + bx[c]];
                *lo -= grid->lo_free;
                if (*lo < grid->lo_min) *lo = grid->lo_min;
            }
        }

        /* Mark endpoint as occupied */
        if (end_x >= 0 && end_x < grid->width &&
            end_y >= 0 && end_y < grid->height) {
            double *lo = &grid->log_odds[end_y * grid->width + end_x];
            *lo += grid->lo_occ;
            if (*lo > grid->lo_max) *lo = grid->lo_max;
        }
    }

    free(bx); free(by);
    return SLAM_OK;
}

int slam_occgrid_init(slam_occgrid_t *grid,
                      int width, int height,
                      double resolution,
                      double origin_x, double origin_y) {
    if (!grid || width <= 0 || height <= 0) return SLAM_ERR_INVALID_PARAM;

    grid->width      = width;
    grid->height     = height;
    grid->resolution = resolution;
    grid->origin_x   = origin_x;
    grid->origin_y   = origin_y;
    grid->lo_occ     = 0.85;   /* p(occ)=0.7 */
    grid->lo_free    = 0.40;   /* p(free)=0.4 */
    grid->lo_min     = -4.0;
    grid->lo_max     =  4.0;

    int n_cells = width * height;
    grid->log_odds = (double *)calloc(n_cells, sizeof(double));
    if (!grid->log_odds) return SLAM_ERR_MEMORY;

    /* All cells start unknown (log-odds = 0, p = 0.5) */
    return SLAM_OK;
}

void slam_occgrid_free(slam_occgrid_t *grid) {
    if (grid) {
        free(grid->log_odds);
        grid->log_odds = NULL;
    }
}
