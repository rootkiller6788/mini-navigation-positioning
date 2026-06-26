/**
 * @file    example_fastslam.c
 * @brief   FastSLAM Particle Filter SLAM Demo
 *
 * Demonstrates FastSLAM 1.0 with a robot navigating and mapping
 * landmarks using a Rao-Blackwellized particle filter.
 *
 * Scenario:
 *   - Robot starts at origin, drives a figure-8 pattern
 *   - 10 point landmarks scattered in the environment
 *   - 30 particles with independent EKFs per landmark
 *   - Systematic resampling when Neff < N/2
 *
 * Reference: Montemerlo et al. (2002), Thrun Ch.13
 */

#include "slam_types.h"
#include "slam_fastslam.h"
#include "slam_sensor.h"

extern void slam_config_default(slam_config_t *cfg);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

extern double slam_normalize_angle(double theta);

#define NUM_STEPS       100
#define NUM_LANDMARKS   10
#define NUM_PARTICLES   30

int main(void) {
    printf("=== FastSLAM Particle Filter Demo ===\n");
    printf("Robot drives figure-8, mapping %d landmarks with %d particles.\n\n",
           NUM_LANDMARKS, NUM_PARTICLES);

    srand((unsigned)time(NULL));

    /* 1. Setup ground truth landmarks */
    double lm_gt_x[NUM_LANDMARKS], lm_gt_y[NUM_LANDMARKS];
    for (int i = 0; i < NUM_LANDMARKS; i++) {
        double angle = 2.0 * M_PI * i / NUM_LANDMARKS;
        double radius = 5.0 + 2.0 * ((double)rand() / RAND_MAX);
        lm_gt_x[i] = radius * cos(angle);
        lm_gt_y[i] = radius * sin(angle);
    }

    /* 2. Initialize FastSLAM */
    slam_config_t config;
    slam_config_default(&config);
    config.sigma_r = 0.15;
    config.sigma_b = 0.03;
    config.sigma_v = 0.1;
    config.sigma_omega = 0.05;
    config.max_landmarks = NUM_LANDMARKS;

    slam_particle_t *particles;
    slam_pose2d_t init_pose = {0, 0, 0};
    slam_fastslam_init(&particles, NUM_PARTICLES, &init_pose, &config);

    /* 3. True trajectory */
    double true_x = 0, true_y = 0, true_th = 0;

    printf("%6s %10s %10s %10s %8s %6s\n",
           "Step", "Est X", "Est Y", "Err[m]", "Neff", "#LM");
    printf("------ ---------- ---------- ---------- -------- ------\n");

    for (int t = 0; t < NUM_STEPS; t++) {
        /* Figure-8 motion */
        double v = 0.8;
        double w = 0.3 * sin(2.0 * M_PI * t / (double)NUM_STEPS);
        double dt = 0.1;

        /* True motion */
        if (fabs(w) < 1e-8) {
            true_x += v * dt * cos(true_th);
            true_y += v * dt * sin(true_th);
        } else {
            true_x += v/w * (sin(true_th + w*dt) - sin(true_th));
            true_y += v/w * (cos(true_th) - cos(true_th + w*dt));
        }
        true_th = slam_normalize_angle(true_th + w * dt);

        /* Generate observations */
        int num_obs = 0;
        slam_obs_rb_t observations[20];

        for (int j = 0; j < NUM_LANDMARKS; j++) {
            double dx = lm_gt_x[j] - true_x;
            double dy = lm_gt_y[j] - true_y;
            double range = sqrt(dx*dx + dy*dy);

            if (range < config.max_range && range > 0.3 && num_obs < 20) {
                double bearing = slam_normalize_angle(atan2(dy, dx) - true_th);
                double noisy_r = range + ((double)rand()/RAND_MAX - 0.5) * 2.0 * config.sigma_r;
                double noisy_b = bearing + ((double)rand()/RAND_MAX - 0.5) * 2.0 * config.sigma_b;

                observations[num_obs].range       = noisy_r;
                observations[num_obs].bearing     = noisy_b;
                observations[num_obs].range_std   = config.sigma_r;
                observations[num_obs].bearing_std = config.sigma_b;
                observations[num_obs].landmark_id = -1;
                observations[num_obs].timestamp   = (uint64_t)t;
                num_obs++;
            }
        }

        /* FastSLAM step */
        slam_velocity_t vel = {v, w, dt};
        slam_fastslam_step(particles, NUM_PARTICLES, &vel,
                            observations, num_obs, &config);

        /* Get best estimate */
        slam_pose2d_t best_pose;
        slam_fastslam_best_pose(particles, NUM_PARTICLES, &best_pose);

        double err = sqrt((best_pose.x - true_x)*(best_pose.x - true_x)
                        + (best_pose.y - true_y)*(best_pose.y - true_y));

        /* Compute effective sample size */
        double sum_w = 0, sum_w2 = 0;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            sum_w += particles[i].weight;
            sum_w2 += particles[i].weight * particles[i].weight;
        }
        double neff = sum_w2 > 0 ? (sum_w * sum_w) / sum_w2 : 0;

        /* Count landmarks in best particle */
        int best_i = 0;
        double best_w = particles[0].log_weight;
        for (int i = 1; i < NUM_PARTICLES; i++) {
            if (particles[i].log_weight > best_w) {
                best_w = particles[i].log_weight;
                best_i = i;
            }
        }

        if (t % 10 == 0) {
            printf("%6d %10.3f %10.3f %10.3f %8.1f %6d\n",
                   t, best_pose.x, best_pose.y, err, neff,
                   particles[best_i].num_landmarks);
        }
    }

    /* Final results */
    slam_pose2d_t final_pose;
    slam_fastslam_best_pose(particles, NUM_PARTICLES, &final_pose);

    /* Get map from best particle */
    int best_i = 0;
    double best_w = particles[0].log_weight;
    for (int i = 1; i < NUM_PARTICLES; i++) {
        if (particles[i].log_weight > best_w) {
            best_w = particles[i].log_weight;
            best_i = i;
        }
    }

    printf("------ ---------- ---------- ---------- -------- ------\n");
    printf("\nResults:\n");
    printf("  Final estimated pose: (%.3f, %.3f, %.2f°)\n",
           final_pose.x, final_pose.y, final_pose.theta * 180.0 / M_PI);
    printf("  True final pose:      (%.3f, %.3f, %.2f°)\n",
           true_x, true_y, true_th * 180.0 / M_PI);
    printf("  Landmarks mapped:      %d / %d\n",
           particles[best_i].num_landmarks, NUM_LANDMARKS);

    printf("  Map estimate (best particle):\n");
    for (int j = 0; j < particles[best_i].num_landmarks; j++) {
        double err_j = sqrt(
            (particles[best_i].landmarks[j].x - lm_gt_x[j]) *
            (particles[best_i].landmarks[j].x - lm_gt_x[j])
          + (particles[best_i].landmarks[j].y - lm_gt_y[j]) *
            (particles[best_i].landmarks[j].y - lm_gt_y[j]));
        printf("    LM%d: est=(%.2f, %.2f)  true=(%.2f, %.2f)  err=%.3f m\n",
               j,
               particles[best_i].landmarks[j].x,
               particles[best_i].landmarks[j].y,
               lm_gt_x[j], lm_gt_y[j], err_j);
    }

    slam_fastslam_free(particles, NUM_PARTICLES);
    printf("\nFastSLAM demo complete.\n");
    return 0;
}
