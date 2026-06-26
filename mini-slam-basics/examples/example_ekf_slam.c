/**
 * @file    example_ekf_slam.c
 * @brief   End-to-end EKF-SLAM simulation with range-bearing landmarks
 *
 * Simulates a differential-drive robot moving in a 2D environment with
 * point landmarks. The robot uses EKF-SLAM to simultaneously estimate
 * its pose and build a map of landmarks.
 *
 * Scenario: Robot drives a circular path with 8 landmarks in its field
 * of view. Demonstrates prediction, update, data association, and
 * state augmentation.
 *
 * Reference: Dissanayake et al. (2001), Thrun Ch.10
 */

#include "slam_types.h"
#include "slam_ekf.h"
#include "slam_sensor.h"

extern void slam_config_default(slam_config_t *cfg);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

extern double slam_normalize_angle(double theta);
extern double slam_angle_diff(double a, double b);

#define NUM_STEPS       200
#define NUM_LANDMARKS   8

int main(void) {
    printf("=== EKF-SLAM Simulation ===\n");
    printf("Robot drives a circular path, mapping %d landmarks via range-bearing.\n\n",
           NUM_LANDMARKS);

    /* 1. Setup ground truth landmarks (placed in a circle around origin) */
    double lm_gt_x[NUM_LANDMARKS], lm_gt_y[NUM_LANDMARKS];
    for (int i = 0; i < NUM_LANDMARKS; i++) {
        double angle = 2.0 * M_PI * i / NUM_LANDMARKS;
        lm_gt_x[i] = 8.0 * cos(angle);
        lm_gt_y[i] = 8.0 * sin(angle);
    }

    /* 2. Initialize EKF-SLAM */
    slam_config_t config;
    slam_config_default(&config);
    config.sigma_r = 0.15;
    config.sigma_b = 0.03;

    slam_ekf_state_t ekf;
    slam_pose2d_t init = {0, 0, 0};
    slam_ekf_init(&ekf, &config, &init, 0.01, 0.01);

    /* 3. True robot trajectory */
    double true_x = 0, true_y = 0, true_th = 0;
    double total_err = 0.0;

    printf("%6s %10s %10s %10s %6s %6s\n",
           "Step", "Est X", "Est Y", "Err[m]", "#LM", "NEES");
    printf("------ ---------- ---------- ---------- ------ ------\n");

    for (int t = 0; t < NUM_STEPS; t++) {
        /* Motion: circular path */
        double v = 0.5;
        double w = 0.15;
        double dt = 0.1;

        /* True motion (no noise for ground truth) */
        if (fabs(w) < 1e-8) {
            true_x += v * dt * cos(true_th);
            true_y += v * dt * sin(true_th);
        } else {
            true_x += v/w * (sin(true_th + w*dt) - sin(true_th));
            true_y += v/w * (cos(true_th) - cos(true_th + w*dt));
        }
        true_th = slam_normalize_angle(true_th + w * dt);

        /* EKF predict */
        slam_velocity_t vel = {v, w, dt};
        slam_ekf_predict_velocity(&ekf, &vel);

        /* Generate observations of visible landmarks */
        for (int j = 0; j < NUM_LANDMARKS; j++) {
            double dx = lm_gt_x[j] - true_x;
            double dy = lm_gt_y[j] - true_y;
            double range = sqrt(dx*dx + dy*dy);

            if (range < config.max_range && range > 0.5) {
                /* Generate noisy observation */
                double bearing = slam_normalize_angle(atan2(dy, dx) - true_th);
                double noisy_r = range + ((double)rand()/RAND_MAX - 0.5) * 2.0 * config.sigma_r;
                double noisy_b = bearing + ((double)rand()/RAND_MAX - 0.5) * 2.0 * config.sigma_b;

                slam_obs_rb_t obs = {noisy_r, noisy_b, -1,
                                      config.sigma_r, config.sigma_b, (uint64_t)t};

                /* Try to associate or augment */
                int matched_id;
                int rc = slam_ekf_update_unknown(&ekf, &obs, &matched_id);
                if (rc == SLAM_ERR_NO_ASSOC) {
                    slam_ekf_augment(&ekf, &obs, &matched_id);
                }
            }
        }

        /* Compute estimation error */
        double err = sqrt((ekf.robot_pose.x - true_x)*(ekf.robot_pose.x - true_x)
                        + (ekf.robot_pose.y - true_y)*(ekf.robot_pose.y - true_y));
        total_err += err;

        /* NEES check */
        slam_pose2d_t true_pose = {true_x, true_y, true_th};
        double nees;
        slam_ekf_nees(&ekf, &true_pose, &nees);

        /* Print status every 20 steps */
        if (t % 20 == 0) {
            printf("%6d %10.3f %10.3f %10.3f %6d %6.2f\n",
                   t, ekf.robot_pose.x, ekf.robot_pose.y,
                   err, ekf.num_landmarks, nees);
        }
    }

    double avg_err = total_err / NUM_STEPS;
    printf("------ ---------- ---------- ---------- ------ ------\n");
    printf("\nResults:\n");
    printf("  Final pose:       (%.3f, %.3f, %.2f°)\n",
           ekf.robot_pose.x, ekf.robot_pose.y,
           ekf.robot_pose.theta * 180.0 / M_PI);
    printf("  True final pose:  (%.3f, %.3f, %.2f°)\n",
           true_x, true_y, true_th * 180.0 / M_PI);
    printf("  Landmarks mapped: %d / %d\n",
           ekf.num_landmarks, NUM_LANDMARKS);
    printf("  Average error:    %.3f m\n", avg_err);
    printf("  Map estimate:\n");

    for (int j = 0; j < ekf.num_landmarks && j < NUM_LANDMARKS; j++) {
        int bj = 3 + 2*j;
        printf("    LM%d: est=(%.2f, %.2f)  true=(%.2f, %.2f)  err=%.3f m\n",
               j, ekf.state_mean[bj], ekf.state_mean[bj+1],
               lm_gt_x[j], lm_gt_y[j],
               sqrt((ekf.state_mean[bj]-lm_gt_x[j])*(ekf.state_mean[bj]-lm_gt_x[j])
                  + (ekf.state_mean[bj+1]-lm_gt_y[j])*(ekf.state_mean[bj+1]-lm_gt_y[j])));
    }

    slam_ekf_free(&ekf);
    printf("\nEKF-SLAM demo complete.\n");
    return 0;
}
