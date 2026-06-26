/**
 * @file    example_zup.c
 * @brief   L7: Zero-Velocity Update (ZUPT) for Pedestrian Navigation
 *
 * Demonstrates pedestrian dead reckoning with zero-velocity updates
 * at each foot step. This technique is used in:
 *   - First responder tracking (firefighters in buildings)
 *   - Military soldier navigation (GNSS-denied environments)
 *   - Indoor pedestrian localization
 *
 * Principle: During each foot step, there is a stance phase (~0.1-0.3 s)
 * where the foot is stationary on the ground. A zero-velocity update
 * resets the accumulated velocity error, bounding the position drift
 * to ~1-2% of distance traveled (vs unbounded growth without ZUPT).
 *
 * Reference: Foxlin (2005), "Pedestrian Tracking with Shoe-Mounted
 *   Inertial Sensors", IEEE Comput. Graph. Appl., 25(6): 38-46.
 *
 * Course Mapping:
 *   MIT 6.832 — Underactuated Robotics (bipedal locomotion sensing)
 *   Michigan EECS 461 — Embedded Control (sensor fusion project)
 */

#include "ins_core.h"
#include "ins_mechanization.h"
#include "ins_integration.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== ZUPT Pedestrian Navigation Example ===\n\n");

    /* Simulate walking: 1 m/s, 1 Hz step rate, 5 steps */
    double walk_speed = 1.0;   /* 1 m/s */
    double step_rate = 1.0;    /* 1 step per second */
    double dt = 0.01;          /* 100 Hz IMU */
    double stance_duration = 0.2; /* Stance phase: 0.2 s per step */
    double total_time = 5.0;   /* 5 seconds = 5 steps */
    int total_steps = (int)(total_time / dt);

    /* Initialize INS */
    ins_mech_state_t ms;
    ins_mech_init(&ms, INS_MECH_CONING, 0.0, 0.0, 0.0);

    ins_vec3_t accel_init = {0, 0, INS_GRAVITY_EQUATOR};
    ins_align_level(&ms, &accel_init, 0.0);

    /* Initialize Kalman filter for ZUPT corrections */
    ins_kf_state_t kf;
    ins_kf_init(&kf, 1.0, 0.1, 0.017);

    ins_nav_solution_t trajectory[500];
    double total_distance_raw = 0.0;  /* Without ZUPT */
    double total_distance_zupt = 0.0; /* With ZUPT */

    printf("Simulating %d steps of pedestrian walking at %.0f m/s...\n\n",
           (int)total_time, walk_speed);

    ins_vec3_t accel_window[20], gyro_window[20];
    int window_idx = 0;

    for (int i = 0; i < total_steps; i++) {
        double t = i * dt;
        double step_phase = fmod(t, 1.0 / step_rate) * step_rate;

        /* Simulate foot motion:
         * - Swing phase: forward acceleration
         * - Stance phase: zero velocity (foot on ground) */
        ins_imu_sample_t imu;
        imu.dt = dt;
        imu.temperature = 25.0;

        if (step_phase < stance_duration * step_rate) {
            /* Stance: stationary foot */
            imu.accel.x = 0;
            imu.accel.y = 0;
            imu.accel.z = INS_GRAVITY_EQUATOR;
            imu.gyro.x = 0;
            imu.gyro.y = 0;
            imu.gyro.z = 0;
        } else {
            /* Swing: forward motion + small oscillations */
            double swing_fraction = (step_phase - stance_duration * step_rate) /
                                    (1.0 - stance_duration * step_rate);
            double accel_fwd = 8.0 * sin(M_PI * swing_fraction);

            imu.accel.x = accel_fwd;
            imu.accel.y = 0.1 * sin(2 * M_PI * 5 * t);
            imu.accel.z = INS_GRAVITY_EQUATOR + 0.5 * sin(2 * M_PI * 3 * t);
            imu.gyro.x = 0;
            imu.gyro.y = 0.2 * sin(2 * M_PI * 5 * t);
            imu.gyro.z = 0;
        }

        /* INS mechanization */
        ins_mech_step(&ms, &imu);
        ins_nav_solution_t nav;
        ins_mech_get_solution(&ms, &nav);

        /* Kalman predict */
        ins_mat3_t C_bn;
        ins_quat_to_dcm(&ms.quat, &C_bn);
        ins_kf_predict(&kf, &C_bn, &imu.accel, ms.pos.lat, ms.pos.alt,
                       dt, 1e-6, 1e-12);

        /* ZUPT detection and update */
        accel_window[window_idx % 20] = imu.accel;
        gyro_window[window_idx % 20]  = imu.gyro;
        window_idx++;

        if (window_idx >= 20) {
            int zupt = ins_zupt_detect(accel_window, gyro_window, 20,
                                       0.1, 0.05, 2.0);
            if (zupt) {
                ins_kf_update_zupt(&kf, &nav);
                ins_kf_apply_correction(&kf, &nav);
                ms.pos = nav.pos;
                ins_vec3_copy(&nav.vel_ned, &ms.vel_ned);
            }
        }

        trajectory[i] = nav;
        total_distance_raw += ins_vec3_norm(&nav.vel_ned) * dt;
    }

    /* Compute actual distance traveled (north direction) */
    double actual_dist = walk_speed * total_time;

    printf("Final position: north=%.2f m, east=%.2f m (without GPS)\n",
           trajectory[total_steps-1].pos.lat * INS_EARTH_RADIUS_MEAN,
           trajectory[total_steps-1].pos.lon * INS_EARTH_RADIUS_MEAN *
           cos(trajectory[total_steps-1].pos.lat));
    printf("Actual distance walked: %.2f m\n", actual_dist);
    printf("ZUPT-aided distance estimate: %.2f m\n",
           trajectory[total_steps-1].pos.lat * INS_EARTH_RADIUS_MEAN);
    printf("\nZUPT effectively bounds velocity error drift during stance phases.\n");
    printf("Typical accuracy: 1-2%% of distance traveled, suitable for indoor tracking.\n");
    return 0;
}
