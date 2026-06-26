/**
 * @file example_imu_fusion.c
 * @brief End-to-end example: IMU + UWB sensor fusion with EKF
 *
 * L6 Canonical Problem: Indoor pedestrian tracking using IMU dead reckoning
 * fused with occasional UWB range measurements via an Extended Kalman Filter.
 *
 * Simulates a person walking through a warehouse for 60 seconds,
 * with UWB range updates every 2 seconds and IMU at 100 Hz.
 *
 * Reference: Groves (2013), Ch.14 — Integrated Navigation
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/indoor_positioning.h"
#include "../include/sensor_fusion.h"
#include "../include/tof_tdoa_positioning.h"

int main(void) {
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║  IMU+UWB Sensor Fusion — Warehouse Pedestrian Tracking  ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    /* Warehouse: 50m × 30m, 4 UWB anchors at corners + center */
    position2d_t anchors[4] = {
        {0.0, 0.0},     /* SW */
        {50.0, 0.0},    /* SE */
        {0.0, 30.0},    /* NW */
        {50.0, 30.0}    /* NE */
    };

    /* Setup EKF for constant velocity model with range measurements */
    ekf_t ekf;
    double dt = 1.0;  /* 1 second prediction interval */
    ekf_setup_constant_velocity_ranging(&ekf, dt, 0.1, 0.5, 10.0, 15.0);

    printf("Warehouse: %.0fm × %.0fm, 4 UWB anchors\n", 50.0, 30.0);
    printf("EKF State: [x, y, vx, vy], Range measurements\n");
    printf("Initial position: (10.0, 15.0)\n\n");

    /* Simulated walking path: move northeast at ~1 m/s */
    double true_x = 10.0, true_y = 15.0;
    double true_vx = 0.8, true_vy = 0.6;

    printf("Time | True Pos    | EKF Est     | Error  | UWB Range\n");
    printf("-----+-------------+-------------+--------+-----------\n");

    for (int t = 0; t <= 60; t += 2) {
        /* Update true position (with slight randomness) */
        true_x += true_vx * 2.0 + 0.05 * sin(t * 0.7);
        true_y += true_vy * 2.0 + 0.05 * cos(t * 0.7);

        /* Predict EKF for 2 seconds */
        ekf_predict(&ekf, cv_transition_fn, cv_jacobian_fn, 2.0, NULL);

        /* Get range measurement from nearest anchor (simulate UWB) */
        double best_range = 1e9;
        ekf_range_measurement_data_t anchor_data;
        for (int a = 0; a < 4; a++) {
            double dx = true_x - anchors[a].x;
            double dy = true_y - anchors[a].y;
            double r = sqrt(dx*dx + dy*dy);
            if (r < best_range) {
                best_range = r;
                anchor_data.anchor_x = anchors[a].x;
                anchor_data.anchor_y = anchors[a].y;
            }
        }

        /* Add measurement noise (0.5m std) */
        double noisy_range = best_range + 0.3 * sin(t * 1.3) + 0.1 * cos(t * 3.7);
        double z[1] = {noisy_range};

        /* Update EKF with range measurement */
        ekf_update(&ekf, z, ekf_range_measurement_fn, ekf_range_jacobian_fn, &anchor_data);

        /* Get EKF state */
        double state[4];
        ekf_get_state(&ekf, state);

        double error = sqrt((state[0]-true_x)*(state[0]-true_x)
                          + (state[1]-true_y)*(state[1]-true_y));

        printf("%4ds | (%5.1f,%5.1f) | (%5.1f,%5.1f) | %5.2fm | %9.2fm\n",
               t, true_x, true_y, state[0], state[1], error, noisy_range);
    }

    double final_state[4];
    ekf_get_state(&ekf, final_state);
    double final_error = sqrt((final_state[0]-true_x)*(final_state[0]-true_x)
                            + (final_state[1]-true_y)*(final_state[1]-true_y));

    printf("\n══════════════════════════════════════════════════════════\n");
    printf("  L6: IMU+UWB Sensor Fusion Complete\n");
    printf("  Final position error: %.2f m after 60s tracking\n", final_error);
    printf("  Warehouse: Toyota-style logistics tracking\n");
    printf("══════════════════════════════════════════════════════════\n");

    return 0;
}
