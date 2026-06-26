/**
 * @file nav_imu.h
 * @brief Inertial Measurement Unit Models and Mechanization
 *
 * L2 Core Concepts: Strapdown inertial navigation, IMU error modeling.
 * L5 Algorithms: IMU mechanization, coning/sculling compensation.
 *
 * Reference: Titterton & Weston, "Strapdown Inertial Navigation Technology"
 *            Savage (1998), "Strapdown Analytics"
 */

#ifndef NAV_IMU_H
#define NAV_IMU_H

#include "nav_common.h"
#include "nav_rotation.h"

typedef struct {
    nav_imu_error_t gyro;
    nav_imu_error_t accel;
    NAV_PRECISION   sample_rate_hz;
    NAV_PRECISION   gyro_range;
    NAV_PRECISION   accel_range;
    int             axes_configured;
} nav_imu_config_t;

typedef struct {
    nav_geodetic_t  pos;
    nav_vector3_t   vel_ned;
    nav_quat_t      quat;
    nav_dcm_t       dcm;
    NAV_PRECISION   gyro_bias[3];
    NAV_PRECISION   accel_bias[3];
    NAV_PRECISION   P[225];
    uint64_t        timestamp;
} nav_ins_state_t;

void nav_imu_correct_errors(NAV_PRECISION gyro_out[3],
                             NAV_PRECISION accel_out[3],
                             const NAV_PRECISION gyro_raw[3],
                             const NAV_PRECISION accel_raw[3],
                             const NAV_PRECISION gyro_bias[3],
                             const NAV_PRECISION accel_bias[3],
                             const nav_imu_error_t *gyro_err,
                             const nav_imu_error_t *accel_err);

void nav_imu_correct_bias(NAV_PRECISION gyro_out[3],
                           NAV_PRECISION accel_out[3],
                           const NAV_PRECISION gyro_raw[3],
                           const NAV_PRECISION accel_raw[3],
                           const NAV_PRECISION gyro_bias[3],
                           const NAV_PRECISION accel_bias[3]);

void nav_ins_mechanize(nav_ins_state_t *state,
                        const NAV_PRECISION gyro_c[3],
                        const NAV_PRECISION accel_c[3],
                        NAV_PRECISION dt);

void nav_ins_mechanize_2sample(nav_ins_state_t *state,
                                const NAV_PRECISION gyro_prev[3],
                                const NAV_PRECISION gyro_curr[3],
                                const NAV_PRECISION accel_prev[3],
                                const NAV_PRECISION accel_curr[3],
                                NAV_PRECISION dt);

void nav_ins_update_position(nav_ins_state_t *state, NAV_PRECISION dt);
void nav_ins_update_velocity(nav_ins_state_t *state,
                              const NAV_PRECISION accel_ned[3],
                              NAV_PRECISION dt);
void nav_ins_update_attitude(nav_ins_state_t *state,
                              const NAV_PRECISION gyro_b[3],
                              NAV_PRECISION dt);

int nav_ins_coarse_alignment(nav_ins_state_t *state,
                              const NAV_PRECISION accel_mean[3],
                              const NAV_PRECISION gyro_mean[3],
                              NAV_PRECISION lat_rad);

int nav_ins_zupt(nav_ins_state_t *state,
                 const NAV_PRECISION gyro[3],
                 const NAV_PRECISION accel[3],
                 NAV_PRECISION dt);

void nav_ins_nhc(nav_ins_state_t *state, NAV_PRECISION R[9]);

void nav_ins_error_transition(NAV_PRECISION Phi[225],
                               const nav_ins_state_t *state,
                               NAV_PRECISION dt);

void nav_ins_process_noise(NAV_PRECISION Qd[225],
                            const nav_ins_state_t *state,
                            NAV_PRECISION gyro_arw,
                            NAV_PRECISION accel_vrw,
                            NAV_PRECISION gyro_bi,
                            NAV_PRECISION accel_bi,
                            NAV_PRECISION dt);

#endif /* NAV_IMU_H */
