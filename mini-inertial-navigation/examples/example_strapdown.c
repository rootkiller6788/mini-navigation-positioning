/**
 * @file    example_strapdown.c
 * @brief   L6: Complete Strapdown INS Navigation example
 *
 * Demonstrates a full strapdown inertial navigation solution from
 * simulated IMU data. A vehicle moves in a circular trajectory and
 * the INS computes its position, velocity, and attitude.
 *
 * This is the canonical problem for inertial navigation:
 *   Given: Initial position, velocity, attitude + IMU data
 *   Find: Vehicle trajectory over time
 *
 * Course Mapping:
 *   TU Munich — Navigation (inertial systems lab exercise)
 *   Tsinghua — Inertial Navigation (strapdown simulation)
 */

#include "ins_core.h"
#include "ins_attitude.h"
#include "ins_mechanization.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Strapdown INS Navigation Example ===\n\n");

    /* Initialize mechanization: starting at Beijing (approx) */
    ins_mech_state_t ms;
    double lat0 = 39.9 * INS_DEG2RAD;
    double lon0 = 116.4 * INS_DEG2RAD;
    double alt0 = 50.0;

    ins_mech_init(&ms, INS_MECH_CONING, lat0, lon0, alt0);

    /* Level alignment (assuming stationary start) */
    ins_vec3_t accel_init = {0, 0, INS_GRAVITY_EQUATOR};
    ins_align_level(&ms, &accel_init, 0.0);

    printf("Initial position: lat=%.4f deg, lon=%.4f deg, alt=%.1f m\n",
           lat0 * INS_RAD2DEG, lon0 * INS_RAD2DEG, alt0);

    /* Simulate circular trajectory (constant speed, constant turn rate) */
    double vt = 30.0;        /* 30 m/s (~108 km/h) */
    double turn_rate = 5.0 * INS_DEG2RAD;  /* 5 deg/s turn */
    double radius = vt / turn_rate;
    double dt = 0.01;
    int total_steps = 1500;  /* 15 seconds */

    printf("Simulating circular trajectory: speed=%.0f m/s, turn=%.0f deg/s\n",
           vt, turn_rate * INS_RAD2DEG);
    printf("Turn radius: %.1f m\n\n", radius);

    ins_nav_solution_t trajectory[1500];

    for (int i = 0; i < total_steps; i++) {
        double t = i * dt;

        /* IMU truth for circular motion:
         * Specific force in body frame: f_b = a_body - C_n_b * g
         * In a coordinated turn: a_body = [v_turn^2/r, 0, 0] (centripetal, in turn plane)
         * Actually for coordinated turn, specific force is along body z-axis
         * with centripetal lateral force applied via roll angle
         */

        /* Gyro: turn rate about body z-axis */
        ins_imu_sample_t imu;
        imu.gyro.x = 0.0;
        imu.gyro.y = 0.0;
        imu.gyro.z = turn_rate;

        /* Accel: gravity in body frame + centripetal
         * For a banked turn at roll angle phi:
         * Lateral accel = v^2/r = g*tan(phi) */
        double roll = atan(vt * vt / (radius * INS_GRAVITY_EQUATOR));
        ins_euler_t bank_att = {roll, 0.0, turn_rate * t};
        ins_mat3_t C_bn;
        ins_euler_to_dcm(&bank_att, &C_bn);

        /* f_body = C_n_b * a_n - C_n_b * g_n = C_n_b * [0,0,0]^T - C_n_b * g_n */
        /* Actually for coordinated turn, the accelerometer measures reaction to
           combined gravity and centripetal force */
        ins_vec3_t accel_ned = {0, 0, -INS_GRAVITY_EQUATOR};
        ins_mat3_t C_nb;
        ins_mat3_transpose(&C_bn, &C_nb);
        ins_mat3_mul_vec(&C_nb, &accel_ned, &imu.accel);

        /* Add centripetal in body frame: a_lateral = v^2/r along body y-axis */
        imu.accel.y += vt * vt / radius;

        imu.dt = dt;
        imu.temperature = 25.0;

        ins_mech_step(&ms, &imu);
        ins_mech_get_solution(&ms, &trajectory[i]);

        if (i % 300 == 0) {
            printf("t=%.1fs  lat=%.6f deg  lon=%.6f deg  hdg=%.1f deg  v=%.1f m/s\n",
                   t,
                   trajectory[i].pos.lat * INS_RAD2DEG,
                   trajectory[i].pos.lon * INS_RAD2DEG,
                   trajectory[i].yaw * INS_RAD2DEG,
                   ins_vec3_norm(&trajectory[i].vel_ned));
        }
    }

    printf("\n... trajectory computed for %.0f seconds.\n", total_steps * dt);
    printf("INS successfully navigated the vehicle through a circular path.\n");
    printf("In a real system, GPS updates would bound the drift error.\n");
    return 0;
}
