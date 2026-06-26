/**
 * @example example_ahrs.c
 * @brief AHRS demonstration: Madgwick filter with simulated IMU data.
 *
 * Simulates a simple rotation and shows attitude estimation.
 * L6 Canonical Problem: AHRS attitude from IMU+magnetometer.
 */

#include "nav_rotation.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== AHRS Example: Madgwick Filter ===\n\n");
    nav_quat_t q;
    nav_quat_identity(&q);
    /* Simulated gyro: constant rotation around z-axis at 10 deg/s */
    NAV_PRECISION gyro[3] = {0.0, 0.0, nav_deg2rad(10.0)};
    /* Simulated accel: gravity pointing down */
    NAV_PRECISION accel[3] = {0.0, 0.0, 9.81};
    NAV_PRECISION mag[3] = {0.3, 0.0, 0.5};
    NAV_PRECISION dt = 0.01;
    NAV_PRECISION beta = 0.1;
    /* Run for 100 steps (1 second) */
    printf("Step | Roll(deg) | Pitch(deg) | Yaw(deg)\n");
    printf("-----|-----------|------------|---------\n");
    extern void nav_ahrs_madgwick(nav_quat_t*,const NAV_PRECISION*,
                                   const NAV_PRECISION*,const NAV_PRECISION*,
                                   NAV_PRECISION,NAV_PRECISION);
    for (int step = 0; step <= 100; step += 10) {
        nav_euler_t euler;
        nav_quat_to_euler(&euler, &q);
        printf("%4d | %9.2f | %10.2f | %8.2f\n",
               step, nav_rad2deg(euler.roll),
               nav_rad2deg(euler.pitch), nav_rad2deg(euler.yaw));
        for (int i = 0; i < 10; i++)
            nav_ahrs_madgwick(&q, gyro, accel, mag, dt, beta);
    }
    printf("\nAttitude estimation complete.\n");
    return 0;
}
