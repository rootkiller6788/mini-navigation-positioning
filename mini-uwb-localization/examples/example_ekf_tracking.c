#include <stdio.h>
#include <math.h>
#include "../include/uwb_types.h"
#include "../include/uwb_tracking.h"

int main(void) {
    /* Setup 4 anchors in corners of a 20m x 15m warehouse */
    uwb_pos3d_t anchors[4] = {
        {0.0, 0.0, 3.0},
        {20.0, 0.0, 3.0},
        {20.0, 15.0, 3.0},
        {0.0, 15.0, 3.0}
    };

    /* EKF config */
    ekf_config_t cfg = {
        .process_noise_pos = 0.01,
        .process_noise_vel = 0.1,
        .measurement_noise_range = 0.01,
        .initial_pos_uncertainty = 1.0,
        .initial_vel_uncertainty = 0.1,
        .max_iterations = 10,
        .convergence_threshold = 1e-6
    };

    /* Initialize EKF at known start position */
    ekf_state_t ekf;
    uwb_pos3d_t start = {2.0, 2.0, 1.0};
    double dt = 0.1; /* 100 ms update rate */
    ekf_init(&ekf, &cfg, EKF_STATE_DIM_3D, &start, dt);

    printf("EKF Tracking Demo: tag walks through 20m x 15m warehouse\n");
    printf("Start: (2.0, 2.0), velocity: (1.0, 0.5) m/s\n\n");

    /* Simulate 50 steps of tag walking with velocity (1.0, 0.5) m/s */
    double true_x = 2.0, true_y = 2.0, true_z = 1.0;
    double vx = 1.0, vy = 0.5;
    int step;

    for (step = 0; step < 50; step++) {
        /* True motion */
        true_x += vx * dt;
        true_y += vy * dt;

        /* Generate noisy range measurements */
        double ranges[4];
        int i;
        for (i = 0; i < 4; i++) {
            double dx = true_x - anchors[i].x;
            double dy = true_y - anchors[i].y;
            double dz = true_z - anchors[i].z;
            ranges[i] = sqrt(dx*dx + dy*dy + dz*dz);
            /* Add measurement noise (sigma = 0.1 m) */
            ranges[i] += ((double)((i*17+step*31)%200)/100.0 - 1.0) * 0.1;
        }

        /* EKF predict + update */
        ekf_predict(&ekf, dt);
        ekf_update_range(&ekf, anchors, ranges, 4);

        /* Get estimate */
        uwb_pos3d_t est;
        ekf_get_position(&ekf, &est);

        if (step % 10 == 0) {
            double err = sqrt((est.x-true_x)*(est.x-true_x) +
                              (est.y-true_y)*(est.y-true_y));
            printf("t=%.1fs: est=(%.2f, %.2f) true=(%.2f, %.2f) err=%.3fm\n",
                   step*dt, est.x, est.y, true_x, true_y, err);
        }
    }

    printf("\nFinal estimate: (%.3f, %.3f)\n", ekf.state[0], ekf.state[1]);
    printf("True position:  (%.3f, %.3f)\n", true_x, true_y);

    return 0;
}
