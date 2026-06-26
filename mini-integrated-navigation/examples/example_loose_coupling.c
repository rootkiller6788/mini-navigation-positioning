/**
 * @example example_loose_coupling.c
 * @brief Loosely Coupled INS/GNSS Integration Demonstration
 *
 * Simulates IMU data with GNSS position updates every second.
 * L6 Canonical Problem: Loosely coupled INS/GNSS integration.
 */

#include "nav_integration.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

extern void nav_ahrs_madgwick(nav_quat_t*,const NAV_PRECISION*,
                               const NAV_PRECISION*,const NAV_PRECISION*,
                               NAV_PRECISION,NAV_PRECISION);

int main(void) {
    printf("=== Loosely Coupled INS/GNSS Example ===\n\n");
    nav_loose_integration_t integ;
    nav_loose_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.type = NAV_INTEG_LOOSE;
    cfg.use_gnss_velocity = 1;
    cfg.pos_noise_n = cfg.pos_noise_e = 1.0;
    cfg.pos_noise_d = 4.0;
    cfg.vel_noise_n = cfg.vel_noise_e = 0.01;
    cfg.vel_noise_d = 0.04;
    if (nav_loose_init(&integ, &cfg) != 0) {
        printf("Failed to initialize integration.\n");
        return 1;
    }
    /* Simulate 10 seconds of IMU data at 100Hz with GNSS at 1Hz */
    printf("Time(s) | Lat(deg) | Lon(deg) | Alt(m) | Vn(m/s) | Ve(m/s)\n");
    printf("--------|----------|----------|--------|----------|---------\n");
    for (int t = 0; t <= 10; t++) {
        /* Simulate stationary IMU with noise */
        nav_imu_meas_t imu;
        imu.gyro_x = 0.0;   imu.gyro_y = 0.0;   imu.gyro_z = 7.29e-5;
        imu.accel_x = 0.0;  imu.accel_y = 0.0;  imu.accel_z = 9.81;
        imu.dt = 0.01;
        for (int i = 0; i < 100; i++) {
            nav_loose_predict(&integ, &imu);
        }
        /* GNSS update every second */
        nav_gnss_solution_t gnss;
        memset(&gnss, 0, sizeof(gnss));
        gnss.pos.latitude = 0.001 * t;  /* moving north */
        gnss.pos.longitude = 0.0;
        gnss.pos.altitude = 100.0;
        gnss.vel_enu.x = 0.1;
        gnss.vel_enu.y = 0.0;
        gnss.timestamp = t * 1000000;
        nav_loose_update(&integ, &gnss);
        nav_solution_t sol;
        nav_loose_get_solution(&sol, &integ);
        printf("%6d  | %8.4f | %8.4f | %6.1f | %8.4f | %8.4f\n",
               t, nav_rad2deg(sol.pos.latitude),
               nav_rad2deg(sol.pos.longitude), sol.pos.altitude,
               sol.vel_ned.x, sol.vel_ned.y);
    }
    printf("\nIntegration simulation complete.\n");
    nav_ekf_free(integ.ekf);
    return 0;
}
