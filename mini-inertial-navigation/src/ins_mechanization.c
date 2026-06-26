/**
 * @file    ins_mechanization.c
 * @brief   Strapdown INS mechanization: attitude, velocity, position propagation
 *
 * Knowledge Coverage:
 *   L5 (Algorithms): Strapdown mechanization, velocity/position integration
 *   L6 (Canonical Problems): Complete INS navigation loop
 *   L4 (Fundamental Laws): Schuler pendulum, Newton's laws in rotating frame
 *
 * The strapdown INS equations in the NED frame:
 *
 * Velocity: dv^n/dt = C_b^n * f^b - (2*w_ie^n + w_en^n) x v^n + g^n
 *   where f^b = specific force measured by accelerometers
 *         C_b^n = body-to-NED transformation from attitude
 *         w_ie^n = Earth rate in NED
 *         w_en^n = transport rate
 *         g^n = gravity vector in NED
 *
 * Position: d(lat)/dt = v_n / (M + h)
 *           d(lon)/dt = v_e / ((N + h) * cos(lat))
 *           d(alt)/dt = -v_d
 *
 * Reference:
 *   Titterton & Weston (2004), Chapters 11-12.
 *   Savage (2000), "Strapdown Analytics", Strapdown Associates.
 *   Groves (2013), Chapter 5, "Inertial Navigation".
 *
 * Course Mapping:
 *   MIT 6.832 - Underactuated Robotics (state estimation)
 *   Stanford AA272 - Global Positioning Systems (INS mechanization)
 *   Georgia Tech ECE 6601 - Communications (navigation)
 *   TU Munich - Navigation (inertial systems)
 *   Tsinghua - Inertial Navigation (precision instrument)
 */

#include "ins_mechanization.h"
#include "ins_core.h"
#include "ins_attitude.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

/* =========================================================================
 * L5: Mechanization Initialization
 * ========================================================================= */

void ins_mech_init(ins_mech_state_t *state, ins_mech_method_t method,
                   double init_lat, double init_lon, double init_alt) {
    if (!state) return;

    memset(state, 0, sizeof(*state));
    state->method = method;
    state->pos.lat = init_lat;
    state->pos.lon = init_lon;
    state->pos.alt = init_alt;
    ins_quat_identity(&state->quat);
    ins_vec3_zero(&state->vel_ned);
    ins_vec3_zero(&state->omega_prev);
    ins_vec3_zero(&state->accel_prev);
    state->dt_prev = 0.0;
    state->initialized = 0;
    state->iteration = 0;
}

/* =========================================================================
 * L5: Level Alignment
 *
 * For a stationary vehicle, accelerometers measure:
 *   a_body = C_body_ned * a_ned
 * where a_ned = -g^n + disturbance (nominally [0, 0, g] in NED
 * because gravity is downward and accelerometer measures reaction)
 *
 * Actually: accelerometer measures specific force = a_kinematic - g
 * For stationary: a_kinematic = 0, so f_body = -g_body = C_body_ned^T * [0, 0, -g]
 * but sign conventions vary. We assume:
 *   accel_reading_body ~ C_body_ned^T * [0, 0, -g]
 *
 * With acceleration frame conventions:
 *   a_x_body ~ g*sin(pitch) ~ 0 for level
 *   a_y_body ~ -g*sin(roll)*cos(pitch)
 *   a_z_body ~ g*cos(roll)*cos(pitch) ~ -g
 *
 * Standard level formulas (assuming a measures reaction to gravity):
 *   roll  = atan2(-a_y, -a_z)
 *   pitch = atan2(a_x, sqrt(a_y^2 + a_z^2))
 *   yaw   = heading_input (not observable from accelerometer alone)
 * ========================================================================= */

int ins_align_level(ins_mech_state_t *state, const ins_vec3_t *accel_body,
                     double heading) {
    if (!state || !accel_body) return -1;

    double ax = accel_body->x;
    double ay = accel_body->y;
    double az = accel_body->z;
    double mag = sqrt(ax * ax + ay * ay + az * az);

    /* Sanity check: acceleration magnitude should be ~g */
    double g_nominal = INS_GRAVITY_EQUATOR;
    if (mag < 0.1 * g_nominal || mag > 10.0 * g_nominal) return -1;

    ins_euler_t att;
    att.roll  = atan2(-ay, -az);
    att.pitch = atan2(ax, sqrt(ay * ay + az * az));
    att.yaw   = heading;

    ins_euler_to_quat(&att, &state->quat);
    ins_quat_normalize(&state->quat);

    return 0;
}

/* =========================================================================
 * L5: Main Strapdown Mechanization Step
 *
 * Full navigation update from one IMU sample.
 *
 * The velocity update in the rotating NED frame:
 *   v^n(t+dt) = v^n(t) + dv
 * where dv = (C_b^n * f^b + g^n - (2*w_ie^n + w_en^n) x v^n) * dt
 *
 * The position update integrates NED velocity to geodetic coordinates.
 * ========================================================================= */

int ins_mech_step(ins_mech_state_t *state, const ins_imu_sample_t *imu) {
    if (!state || !imu) return -1;
    if (imu->dt <= 0) return -1;

    double dt = imu->dt;
    ins_vec3_t omega, accel;
    ins_vec3_copy(&imu->gyro, &omega);
    ins_vec3_copy(&imu->accel, &accel);

    /* Step 1: Attitude Update */
    if (state->initialized && state->method == INS_MECH_CONING) {
        ins_quat_update_coning(&state->quat, &state->omega_prev, &omega, dt);
    } else if (state->method == INS_MECH_EULER) {
        ins_quat_update_euler(&state->quat, &omega, dt);
    } else {
        /* Default to exact update (works for RK2, RK4 with midpoint averaging) */
        ins_quat_update_exact(&state->quat, &omega, dt);
    }

    /* Step 2: Get DCM from current attitude */
    ins_mat3_t C_bn;
    ins_quat_to_dcm(&state->quat, &C_bn);

    /* Step 3: Transform specific force to NED frame */
    ins_vec3_t f_ned;
    ins_mat3_mul_vec(&C_bn, &accel, &f_ned);

    /* Step 4: Velocity update */
    /* Get gravity vector in NED */
    ins_vec3_t g_ned;
    ins_gravity_ned(state->pos.lat, state->pos.alt, &g_ned);

    /* Get Coriolis + transport rate acceleration */
    ins_vec3_t a_coriolis;
    ins_coriolis_accel(&state->vel_ned, state->pos.lat, state->pos.alt,
                       &a_coriolis);

    /* v_new = v_old + (f_ned + g_ned + a_coriolis) * dt */
    /* Note: a_coriolis already contains the negative of the cross product,
     * so adding it accounts for -(2*w_ie+w_en) x v */
    state->vel_ned.x += (f_ned.x + g_ned.x + a_coriolis.x) * dt;
    state->vel_ned.y += (f_ned.y + g_ned.y + a_coriolis.y) * dt;
    state->vel_ned.z += (f_ned.z + g_ned.z + a_coriolis.z) * dt;

    /* Step 5: Position update */
    ins_pos_update(&state->pos, &state->vel_ned, dt);

    /* Store for next iteration */
    ins_vec3_copy(&omega, &state->omega_prev);
    ins_vec3_copy(&accel, &state->accel_prev);
    state->dt_prev = dt;
    state->initialized = 1;
    state->iteration++;

    return 0;
}

/* =========================================================================
 * L5: Navigation Solution Extraction
 * ========================================================================= */

void ins_mech_get_solution(const ins_mech_state_t *state,
                            ins_nav_solution_t *nav) {
    if (!state || !nav) return;

    ins_vec3_copy(&state->vel_ned, &nav->vel_ned);
    nav->pos = state->pos;
    nav->q[0] = state->quat.w;
    nav->q[1] = state->quat.x;
    nav->q[2] = state->quat.y;
    nav->q[3] = state->quat.z;
    nav->time = 0.0;

    /* Compute Euler angles */
    ins_euler_t att;
    ins_quat_to_euler(&state->quat, &att);
    nav->roll  = att.roll;
    nav->pitch = att.pitch;
    nav->yaw   = att.yaw;

    /* Compute ECEF velocity (approximate) */
    ins_ecef_t ecef_pos;
    ins_geodetic_to_ecef(&state->pos, &ecef_pos);
    nav->vel_ecef.x = state->vel_ned.x;  /* Simplified — full transform needs NED->ECEF */
    nav->vel_ecef.y = state->vel_ned.y;
    nav->vel_ecef.z = state->vel_ned.z;

    /* Body-frame velocity */
    ins_mat3_t C_bn;
    ins_quat_to_dcm(&state->quat, &C_bn);
    ins_mat3_t C_nb;
    ins_mat3_transpose(&C_bn, &C_nb);
    ins_mat3_mul_vec(&C_nb, &state->vel_ned, &nav->vel_body);
}

/* =========================================================================
 * L5: Standalone Velocity Update
 * ========================================================================= */

void ins_vel_update(ins_vec3_t *vel_ned, const ins_mat3_t *C_body_ned,
                    const ins_vec3_t *f_body, double lat, double alt,
                    double dt) {
    if (!vel_ned || !C_body_ned || !f_body || dt <= 0) return;

    /* Transform specific force to NED */
    ins_vec3_t f_ned;
    ins_mat3_mul_vec(C_body_ned, f_body, &f_ned);

    /* Gravity */
    ins_vec3_t g_ned;
    ins_gravity_ned(lat, alt, &g_ned);

    /* Coriolis */
    ins_vec3_t a_coriolis;
    ins_coriolis_accel(vel_ned, lat, alt, &a_coriolis);

    /* Integrate */
    vel_ned->x += (f_ned.x + g_ned.x + a_coriolis.x) * dt;
    vel_ned->y += (f_ned.y + g_ned.y + a_coriolis.y) * dt;
    vel_ned->z += (f_ned.z + g_ned.z + a_coriolis.z) * dt;
}

void ins_pos_update(ins_geodetic_t *pos, const ins_vec3_t *vel_ned, double dt) {
    if (!pos || !vel_ned || dt <= 0) return;

    double M = ins_meridian_radius(pos->lat);
    double N = ins_prime_vertical_radius(pos->lat);
    double cos_lat = cos(pos->lat);

    if (fabs(cos_lat) < 1e-12) {
        /* Near pole: longitude update is singular */
        pos->lat += vel_ned->x * dt / (M + pos->alt);
        pos->alt += -vel_ned->z * dt;
        return;
    }

    pos->lat += vel_ned->x * dt / (M + pos->alt);
    pos->lon += vel_ned->y * dt / ((N + pos->alt) * cos_lat);
    pos->alt += -vel_ned->z * dt;

    /* Wrap longitude to [-pi, pi] */
    pos->lon = ins_angle_wrap(pos->lon);
}

/* =========================================================================
 * L6: Complete Strapdown Navigation Loop
 * ========================================================================= */

int ins_navigate(const ins_mech_state_t *state,
                 const ins_imu_sample_t *imu_data,
                 size_t num_samples,
                 ins_nav_solution_t *trajectory) {
    if (!state || !imu_data || !trajectory || num_samples == 0) return -1;

    /* Make a mutable copy of the initial state */
    ins_mech_state_t ms;
    memcpy(&ms, state, sizeof(ms));

    for (size_t i = 0; i < num_samples; i++) {
        int ret = ins_mech_step(&ms, &imu_data[i]);
        if (ret != 0) return -1;

        ins_mech_get_solution(&ms, &trajectory[i]);
        trajectory[i].time = (double)(i + 1) * imu_data[0].dt;
    }

    return 0;
}
