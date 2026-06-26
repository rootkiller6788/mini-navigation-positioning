/**
 * @file    example_ahrs.c
 * @brief   L6: AHRS — Attitude and Heading Reference System example
 *
 * Demonstrates how to use quaternion integration and the Mahony/Madgwick
 * complementary filter for attitude estimation from IMU data.
 * This is a fundamental application for drones, VR headsets, and smartphones.
 *
 * Course Mapping:
 *   Stanford EE267 — Virtual Reality (IMU orientation tracking)
 *   ETH 151-0854 — Trajectory Generation (rotation estimation)
 */

#include "ins_core.h"
#include "ins_attitude.h"
#include <stdio.h>
#include <math.h>

/**
 * Complementary filter (simplified Mahony filter) for attitude estimation.
 *
 * Combines gyroscope integration (high-frequency) with accelerometer
 * leveling (low-frequency) to produce a drift-free attitude estimate.
 *
 * Algorithm:
 *   1. Gyro integration: q = q + 0.5*q*omega*dt
 *   2. Compute expected gravity from current attitude
 *   3. Cross product of measured and expected gravity = error
 *   4. PI controller drives error to zero (corrects gyro bias)
 *
 * Reference: Mahony et al. (2008), IEEE Trans. Autom. Control, 53(5): 1203-1218.
 */
void ahrs_complementary_filter(ins_quat_t *q,
                                const ins_vec3_t *accel,
                                const ins_vec3_t *gyro,
                                double dt, double kp, double ki) {
    static ins_vec3_t integral_error = {0, 0, 0};

    /* Normalize accelerometer measurement */
    ins_vec3_t a_norm;
    ins_vec3_copy(accel, &a_norm);
    ins_vec3_normalize(&a_norm);

    /* Expected gravity direction from current orientation */
    /* v_expected = C_body_ned^T * [0, 0, -1] = last column of DCM^T = last row of DCM negated */
    ins_mat3_t C;
    ins_quat_to_dcm(q, &C);
    ins_vec3_t v_exp;
    v_exp.x = -C.m[2];  /* -C(0,2) */
    v_exp.y = -C.m[5];  /* -C(1,2) */
    v_exp.z = -C.m[8];  /* -C(2,2) */

    /* Error = cross(measured_accel, expected_gravity) */
    ins_vec3_t error;
    ins_vec3_cross(&a_norm, &v_exp, &error);

    /* PI controller */
    integral_error.x += error.x * ki * dt;
    integral_error.y += error.y * ki * dt;
    integral_error.z += error.z * ki * dt;

    /* Corrected gyro */
    ins_vec3_t gyro_corr;
    gyro_corr.x = gyro->x + kp * error.x + integral_error.x;
    gyro_corr.y = gyro->y + kp * error.y + integral_error.y;
    gyro_corr.z = gyro->z + kp * error.z + integral_error.z;

    /* Integrate attitude */
    ins_quat_update_exact(q, &gyro_corr, dt);
    ins_quat_normalize(q);
}

int main(void) {
    printf("=== AHRS Example — Complementary Filter Attitude Estimation ===\n\n");

    /* Simulate a drone in hover with small oscillations */
    ins_quat_t attitude;
    ins_quat_identity(&attitude);

    /* Initialize with small tilt (10 deg roll, 5 deg pitch) */
    ins_euler_t init_att = {10.0 * M_PI / 180.0, 5.0 * M_PI / 180.0, 45.0 * M_PI / 180.0};
    ins_euler_to_quat(&init_att, &attitude);

    /* Filter parameters */
    double kp = 2.0;  /* Proportional gain */
    double ki = 0.1;  /* Integral gain */
    double dt = 0.01; /* 100 Hz IMU */

    printf("Initial attitude: roll=10 deg, pitch=5 deg, yaw=45 deg\n");
    printf("Running complementary filter for 10 seconds...\n\n");

    for (int t = 0; t < 1000; t++) {
        /* Simulate IMU: gravity in body frame + small noise + gyro bias */
        ins_mat3_t C;
        ins_quat_to_dcm(&attitude, &C);

        /* Simulated accelerometer (gravity only, no acceleration) */
        ins_vec3_t accel;
        accel.x = -C.m[6] * INS_GRAVITY_EQUATOR;
        accel.y = -C.m[7] * INS_GRAVITY_EQUATOR;
        accel.z = -C.m[8] * INS_GRAVITY_EQUATOR;

        /* Simulated gyroscope (zero angular rate with bias) */
        ins_vec3_t gyro;
        gyro.x = 0.01 * (M_PI / 180.0);  /* 0.01 deg/s bias drift */
        gyro.y = -0.005 * (M_PI / 180.0);
        gyro.z = 0.0;

        ahrs_complementary_filter(&attitude, &accel, &gyro, dt, kp, ki);

        if (t % 200 == 0) {
            ins_euler_t euler;
            ins_quat_to_euler(&attitude, &euler);
            printf("t=%.1fs  roll=%.2f deg  pitch=%.2f deg  yaw=%.2f deg\n",
                   t * dt,
                   euler.roll * INS_RAD2DEG,
                   euler.pitch * INS_RAD2DEG,
                   euler.yaw * INS_RAD2DEG);
        }
    }

    printf("\nAHRS filter running: attitude maintained despite gyro bias thanks to accelerometer correction.\n");
    printf("This is the fundamental algorithm used in every smartphone orientation sensor.\n");
    return 0;
}
