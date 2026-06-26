/**
 * @file nav_ins_gnss_loose.c
 * @brief Loosely Coupled INS/GNSS Integration
 *
 * L6 Canonical Problem: Loosely coupled INS/GNSS using 15-state EKF.
 * GNSS provides position and velocity updates; INS propagates between.
 *
 * Reference: Groves, "Principles of GNSS, Inertial, and Multisensor
 * Navigation", Chapter 14.
 */

#include "nav_integration.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

int nav_loose_init(nav_loose_integration_t *integ,
                    const nav_loose_config_t *config) {
    if (!integ || !config) return -1;
    memset(integ, 0, sizeof(nav_loose_integration_t));
    integ->config = *config;
    /* Allocate 15-state EKF (6 measurements: pos NED + vel NED) */
    int m_dim = 3; /* position only by default */
    if (config->use_gnss_velocity) m_dim = 6;
    integ->ekf = nav_ekf_alloc(15, m_dim, 0);
    if (!integ->ekf) return -1;
    /* Initialize EKF with position/velocity measurement model */
    /* H maps INS errors to position and velocity errors */
    memset(integ->ekf->base.H, 0, m_dim * 15 * sizeof(NAV_PRECISION));
    /* Position measurement: z = [lat_gnss-lat_ins, lon_gnss-lon_ins, h_gnss-h_ins] */
    /* H_pos = [I_3x3, 0_3x12] */
    for (int i = 0; i < 3; i++)
        integ->ekf->base.H[i*15+i] = 1.0;
    /* Velocity measurement (if used) */
    if (config->use_gnss_velocity) {
        for (int i = 0; i < 3; i++)
            integ->ekf->base.H[(3+i)*15+(3+i)] = 1.0;
    }
    /* Set measurement noise */
    NAV_PRECISION R[36] = {0};
    R[0] = config->pos_noise_n; R[7] = config->pos_noise_e; R[14] = config->pos_noise_d;
    if (config->use_gnss_velocity) {
        R[21] = config->vel_noise_n; R[28] = config->vel_noise_e; R[35] = config->vel_noise_d;
    }
    nav_kf_set_R(&integ->ekf->base, R);
    /* Initialize INS with default state */
    integ->ins.pos.latitude = 0.0;
    integ->ins.pos.longitude = 0.0;
    integ->ins.pos.altitude = 0.0;
    integ->ins.vel_ned.x = 0.0;
    integ->ins.vel_ned.y = 0.0;
    integ->ins.vel_ned.z = 0.0;
    nav_quat_identity(&integ->ins.quat);
    nav_quat_to_dcm(&integ->ins.dcm, &integ->ins.quat);
    return 0;
}

void nav_loose_predict(nav_loose_integration_t *integ,
                        const nav_imu_meas_t *imu) {
    if (!integ || !imu) return;
    /* Correct IMU measurements with current bias estimates */
    NAV_PRECISION gyro_c[3], accel_c[3];
    NAV_PRECISION gyro_r[3] = {imu->gyro_x, imu->gyro_y, imu->gyro_z};
    NAV_PRECISION accel_r[3] = {imu->accel_x, imu->accel_y, imu->accel_z};
    nav_imu_correct_bias(gyro_c, accel_c,
                          gyro_r, accel_r,
                          integ->ins.gyro_bias, integ->ins.accel_bias);
    /* INS mechanization */
    nav_ins_mechanize(&integ->ins, gyro_c, accel_c, imu->dt);
    /* Error state covariance propagation */
    if (integ->ekf) {
        NAV_PRECISION Phi[225], Qd[225];
        nav_ins_error_transition(Phi, &integ->ins, imu->dt);
        nav_ins_process_noise(Qd, &integ->ins, 1e-6, 1e-4, 1e-10, 1e-6, imu->dt);
        nav_kf_set_F(&integ->ekf->base, Phi);
        nav_kf_set_Q(&integ->ekf->base, Qd);
        nav_kf_predict(&integ->ekf->base);
    }
    integ->dt_accum += imu->dt;
}

int nav_loose_update(nav_loose_integration_t *integ,
                      const nav_gnss_solution_t *gnss) {
    if (!integ || !gnss || !integ->ekf) return -1;
    /* Build measurement: position innovation (and velocity if used) */
    int m = integ->ekf->base.m;
    NAV_PRECISION *z = (NAV_PRECISION*)malloc(m * sizeof(NAV_PRECISION));
    if (!z) return -1;
    /* Position innovation: convert to local NED frame */
    NAV_PRECISION dlat = gnss->pos.latitude - integ->ins.pos.latitude;
    NAV_PRECISION dlon = gnss->pos.longitude - integ->ins.pos.longitude;
    NAV_PRECISION dalt = gnss->pos.altitude - integ->ins.pos.altitude;
    /* Convert to meters using radii of curvature */
    NAV_PRECISION Rm = nav_meridian_radius(integ->ins.pos.latitude);
    NAV_PRECISION Rn = nav_transverse_radius(integ->ins.pos.latitude);
    NAV_PRECISION h_Rm = Rm + integ->ins.pos.altitude;
    NAV_PRECISION h_Rn = Rn + integ->ins.pos.altitude;
    NAV_PRECISION cl = cos(integ->ins.pos.latitude);
    z[0] = dlat * h_Rm; /* N position error */
    z[1] = dlon * h_Rn * cl; /* E position error */
    z[2] = -dalt; /* D position error (positive down) */
    if (integ->config.use_gnss_velocity) {
        z[3] = gnss->vel_enu.y - integ->ins.vel_ned.y; /* v_east */
        z[4] = gnss->vel_enu.x - integ->ins.vel_ned.x; /* v_north */
        z[5] = -(-integ->ins.vel_ned.z); /* v_down */
    }
    int ret = nav_kf_update(&integ->ekf->base, z);
    free(z);
    integ->gnss_valid = (ret == 0);
    if (ret == 0) {
        integ->first_fix = 1;
        integ->last_gnss_time = gnss->timestamp;
        /* Apply error corrections to INS (closed-loop) */
        nav_loose_correct(integ);
    }
    integ->dt_accum = 0.0;
    return ret;
}

void nav_loose_correct(nav_loose_integration_t *integ) {
    if (!integ) return;
    const NAV_PRECISION *x = nav_kf_get_x(&integ->ekf->base);
    if (!x) return;
    /* Correct position */
    NAV_PRECISION Rm = nav_meridian_radius(integ->ins.pos.latitude);
    NAV_PRECISION Rn = nav_transverse_radius(integ->ins.pos.latitude);
    NAV_PRECISION cl = cos(integ->ins.pos.latitude);
    integ->ins.pos.latitude  -= x[0] / (Rm + integ->ins.pos.altitude);
    integ->ins.pos.longitude -= x[1] / ((Rn + integ->ins.pos.altitude) * cl);
    integ->ins.pos.altitude  += x[2];
    /* Correct velocity */
    integ->ins.vel_ned.x -= x[3];
    integ->ins.vel_ned.y -= x[4];
    integ->ins.vel_ned.z -= x[5];
    /* Correct attitude */
    nav_vector3_t psi;
    psi.x = x[6]; psi.y = x[7]; psi.z = x[8];
    nav_quat_apply_correction(&integ->ins.quat, &psi);
    nav_quat_to_dcm(&integ->ins.dcm, &integ->ins.quat);
    /* Update bias estimates */
    integ->ins.gyro_bias[0] += x[9];
    integ->ins.gyro_bias[1] += x[10];
    integ->ins.gyro_bias[2] += x[11];
    integ->ins.accel_bias[0] += x[12];
    integ->ins.accel_bias[1] += x[13];
    integ->ins.accel_bias[2] += x[14];
    /* Reset error state to zero */
    memset(integ->ekf->base.x, 0, 15 * sizeof(NAV_PRECISION));
}

void nav_loose_get_solution(nav_solution_t *sol,
                             const nav_loose_integration_t *integ) {
    if (!sol || !integ) return;
    memset(sol, 0, sizeof(nav_solution_t));
    sol->pos = integ->ins.pos;
    sol->vel_ned = integ->ins.vel_ned;
    { nav_euler_t eul; nav_quat_to_euler(&eul, &integ->ins.quat);
      sol->roll = eul.roll; sol->pitch = eul.pitch; sol->yaw = eul.yaw; }
    memcpy(sol->quat, &integ->ins.quat, 4*sizeof(NAV_PRECISION));
    memcpy(sol->gyro_bias, integ->ins.gyro_bias, 3*sizeof(NAV_PRECISION));
    memcpy(sol->accel_bias, integ->ins.accel_bias, 3*sizeof(NAV_PRECISION));
    if (integ->ekf) {
        const NAV_PRECISION *P = nav_kf_get_P(&integ->ekf->base);
        if (P) {
            for (int i = 0; i < 3; i++) {
                sol->pos_cov[i*4] = P[i*16];
                sol->vel_cov[i*4] = P[(3+i)*16];
                sol->att_cov[i*4] = P[(6+i)*16];
            }
        }
    }
}
