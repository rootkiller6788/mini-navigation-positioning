/**
 * @file nav_integration.h
 * @brief INS/GNSS Integration Architectures
 *
 * L2 Core Concepts: Loose, tight, deep coupling architectures.
 * L4 Fundamental Laws: Error-state observability, complementary filtering.
 * L6 Canonical Problems: Loosely coupled INS/GNSS, TC INS/GNSS.
 *
 * Reference: Groves, "Principles of GNSS, Inertial, and Multisensor Navigation"
 *            Grewal et al., "Global Positioning Systems, Inertial Navigation, and Integration"
 */

#ifndef NAV_INTEGRATION_H
#define NAV_INTEGRATION_H

#include "nav_common.h"
#include "nav_rotation.h"
#include "nav_kalman.h"
#include "nav_imu.h"
#include "nav_gnss.h"

typedef enum {
    NAV_INTEG_LOOSE = 0,
    NAV_INTEG_TIGHT,
    NAV_INTEG_DEEP,
    NAV_INTEG_ULTRA_TIGHT
} nav_integration_type_t;

typedef struct {
    nav_integration_type_t type;
    int     use_gnss_velocity;
    int     use_heading_update;
    int     use_nhc;
    int     use_zupt;
    NAV_PRECISION pos_noise_n, pos_noise_e, pos_noise_d;
    NAV_PRECISION vel_noise_n, vel_noise_e, vel_noise_d;
    NAV_PRECISION lever_arm[3];
    int     rate_ratio;
} nav_loose_config_t;

typedef struct {
    nav_integration_type_t type;
    NAV_PRECISION pr_noise;
    NAV_PRECISION prr_noise;
    NAV_PRECISION lever_arm[3];
    int     max_svs;
    int     use_carrier_phase;
    int     use_doppler;
} nav_tight_config_t;

typedef struct {
    nav_ins_state_t    ins;
    nav_ekf_t         *ekf;
    nav_loose_config_t config;
    NAV_PRECISION      P[225];
    NAV_PRECISION      x_err[15];
    NAV_PRECISION      dt_accum;
    int                gnss_valid;
    int                first_fix;
    uint64_t           last_gnss_time;
} nav_loose_integration_t;

int nav_loose_init(nav_loose_integration_t *integ, const nav_loose_config_t *config);
void nav_loose_predict(nav_loose_integration_t *integ, const nav_imu_meas_t *imu);
int nav_loose_update(nav_loose_integration_t *integ, const nav_gnss_solution_t *gnss);
void nav_loose_correct(nav_loose_integration_t *integ);
void nav_loose_get_solution(nav_solution_t *sol, const nav_loose_integration_t *integ);

typedef struct {
    nav_ins_state_t    ins;
    nav_ekf_t         *ekf;
    nav_tight_config_t config;
    NAV_PRECISION      P[289];
    NAV_PRECISION      x_err[17];
    NAV_PRECISION      dt_accum;
    int                n_active;
    uint64_t           last_gnss_time;
} nav_tight_integration_t;

int nav_tight_init(nav_tight_integration_t *integ, const nav_tight_config_t *config);
void nav_tight_predict(nav_tight_integration_t *integ, const nav_imu_meas_t *imu);
int nav_tight_update(nav_tight_integration_t *integ, const nav_gnss_sv_t *svs, int n_svs, const nav_gnss_ephemeris_t *eph, int n_eph);
void nav_tight_correct(nav_tight_integration_t *integ);
void nav_tight_get_solution(nav_solution_t *sol, const nav_tight_integration_t *integ);

void nav_tight_pr_innovation(NAV_PRECISION *innov, const nav_gnss_sv_t *sv, const NAV_PRECISION rx_ecef[3], NAV_PRECISION clk_bias);
void nav_tight_prr_innovation(NAV_PRECISION *innov, const nav_gnss_sv_t *sv, const NAV_PRECISION rx_ecef[3], const NAV_PRECISION rx_vel[3], NAV_PRECISION clk_drift);
void nav_los_vector(NAV_PRECISION los[3], const NAV_PRECISION rx_ecef[3], const NAV_PRECISION sat_ecef[3]);
void nav_lever_arm_compensate(NAV_PRECISION gnss_ecef[3], const NAV_PRECISION ins_ecef[3], const NAV_PRECISION R_ned_body[9], const NAV_PRECISION lever_arm[3], NAV_PRECISION lat, NAV_PRECISION lon);

#endif /* NAV_INTEGRATION_H */
