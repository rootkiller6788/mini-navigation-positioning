/**
 * @file nav_imu.c
 * @brief IMU Error Compensation, Strapdown Mechanization, and Alignment
 *
 * L2: Strapdown INS concept
 * L5: IMU mechanization, coning/sculling, coarse alignment, ZUPT, NHC
 * L4: INS error dynamics propagation
 */

#include "nav_imu.h"
#include <math.h>
#include <string.h>

void nav_imu_correct_bias(NAV_PRECISION gyro_out[3], NAV_PRECISION accel_out[3],
                           const NAV_PRECISION gyro_raw[3],
                           const NAV_PRECISION accel_raw[3],
                           const NAV_PRECISION gyro_bias[3],
                           const NAV_PRECISION accel_bias[3]) {
    if (!gyro_out || !accel_out || !gyro_raw || !accel_raw ||
        !gyro_bias || !accel_bias) return;
    for (int i = 0; i < 3; i++) {
        gyro_out[i]  = gyro_raw[i]  - gyro_bias[i];
        accel_out[i] = accel_raw[i] - accel_bias[i];
    }
}

void nav_imu_correct_errors(NAV_PRECISION gyro_out[3], NAV_PRECISION accel_out[3],
                             const NAV_PRECISION gyro_raw[3],
                             const NAV_PRECISION accel_raw[3],
                             const NAV_PRECISION gyro_bias[3],
                             const NAV_PRECISION accel_bias[3],
                             const nav_imu_error_t *gyro_err,
                             const nav_imu_error_t *accel_err) {
    if (!gyro_out || !accel_out || !gyro_raw || !accel_raw ||
        !gyro_bias || !accel_bias) return;
    for (int i = 0; i < 3; i++) {
        NAV_PRECISION sf = (gyro_err ? gyro_err->scale_factor_ppm * 1e-6 : 0.0);
        gyro_out[i] = (gyro_raw[i] - gyro_bias[i]) / (1.0 + sf);
        NAV_PRECISION asf = (accel_err ? accel_err->scale_factor_ppm * 1e-6 : 0.0);
        accel_out[i] = (accel_raw[i] - accel_bias[i]) / (1.0 + asf);
    }
    /* Cross-axis compensation (simplified first-order) */
    if (gyro_err) {
        gyro_out[0] -= gyro_err->misalignment[0] * gyro_raw[1] +
                        gyro_err->misalignment[1] * gyro_raw[2];
        gyro_out[1] -= gyro_err->misalignment[2] * gyro_raw[0];
    }
    if (accel_err) {
        accel_out[0] -= accel_err->misalignment[0] * accel_raw[1] +
                         accel_err->misalignment[1] * accel_raw[2];
        accel_out[1] -= accel_err->misalignment[2] * accel_raw[0];
    }
}

void nav_ins_update_position(nav_ins_state_t *state, NAV_PRECISION dt) {
    (void)dt;
    if (!state) return;
    NAV_PRECISION lat = state->pos.latitude;
    NAV_PRECISION alt = state->pos.altitude;
    NAV_PRECISION Rm = nav_meridian_radius(lat);
    NAV_PRECISION Rn = nav_transverse_radius(lat);
    NAV_PRECISION vn = state->vel_ned.x;
    NAV_PRECISION ve = state->vel_ned.y;
    NAV_PRECISION vd = state->vel_ned.z;
    state->pos.latitude  += vn / (Rm + alt) * dt;
    state->pos.longitude += ve / ((Rn + alt) * cos(lat)) * dt;
    state->pos.altitude  -= vd * dt;
}

void nav_ins_update_velocity(nav_ins_state_t *state,
                              const NAV_PRECISION accel_ned[3],
                              NAV_PRECISION dt) {
    (void)dt;
    if (!state || !accel_ned) return;
    NAV_PRECISION lat = state->pos.latitude;
    NAV_PRECISION alt = state->pos.altitude;
    NAV_PRECISION vn = state->vel_ned.x;
    NAV_PRECISION ve = state->vel_ned.y;
    NAV_PRECISION vd = state->vel_ned.z;
    /* Coriolis and transport rate */
    NAV_PRECISION omega_ie[3], omega_en[3];
    nav_earth_rotation_ned(lat, omega_ie);
    nav_transport_rate_ned(lat, alt, vn, ve, omega_en);
    /* Gravity */
    NAV_PRECISION g = nav_normal_gravity(lat, alt);
    /* Velocity update: dv/dt = f_ned - (2*omega_ie + omega_en) x v + g */
    NAV_PRECISION coriolis_n = -(2.0*omega_ie[2] + omega_en[2])*ve +
                                 (2.0*omega_ie[1] + omega_en[1])*vd;
    NAV_PRECISION coriolis_e =  (2.0*omega_ie[2] + omega_en[2])*vn -
                                 (2.0*omega_ie[0] + omega_en[0])*vd;
    NAV_PRECISION coriolis_d = -(2.0*omega_ie[1] + omega_en[1])*vn +
                                 (2.0*omega_ie[0] + omega_en[0])*ve;
    state->vel_ned.x += (accel_ned[0] + coriolis_n) * dt;
    state->vel_ned.y += (accel_ned[1] + coriolis_e) * dt;
    state->vel_ned.z += (accel_ned[2] + coriolis_d + g) * dt;
}

void nav_ins_update_attitude(nav_ins_state_t *state,
                              const NAV_PRECISION gyro_b[3],
                              NAV_PRECISION dt) {
    (void)dt;
    if (!state || !gyro_b) return;
    nav_vector3_t omega;
    omega.x = gyro_b[0]; omega.y = gyro_b[1]; omega.z = gyro_b[2];
    nav_quat_kinematics(&state->quat, &omega, dt);
    nav_quat_to_dcm(&state->dcm, &state->quat);
}

void nav_ins_mechanize(nav_ins_state_t *state,
                        const NAV_PRECISION gyro_c[3],
                        const NAV_PRECISION accel_c[3],
                        NAV_PRECISION dt) {
    (void)dt;
    if (!state || !gyro_c || !accel_c) return;
    /* Transform accel to NED */
    nav_vector3_t accel_b, accel_ned;
    accel_b.x = accel_c[0]; accel_b.y = accel_c[1]; accel_b.z = accel_c[2];
    nav_dcm_rotate_vector(&accel_ned, &state->dcm, &accel_b);
    NAV_PRECISION accel[3] = {accel_ned.x, accel_ned.y, accel_ned.z};
    nav_ins_update_velocity(state, accel, dt);
    nav_ins_update_position(state, dt);
    nav_ins_update_attitude(state, gyro_c, dt);
}

void nav_ins_mechanize_2sample(nav_ins_state_t *state,
                                const NAV_PRECISION gyro_prev[3],
                                const NAV_PRECISION gyro_curr[3],
                                const NAV_PRECISION accel_prev[3],
                                const NAV_PRECISION accel_curr[3],
                                NAV_PRECISION dt) {
    (void)dt;
    if (!state || !gyro_prev || !gyro_curr || !accel_prev || !accel_curr) return;
    /* Coning compensation: delta_theta = (gyro_prev + gyro_curr) / 2 */
    NAV_PRECISION gyro_avg[3], accel_avg[3];
    for (int i = 0; i < 3; i++) {
        gyro_avg[i]  = 0.5 * (gyro_prev[i] + gyro_curr[i]);
        accel_avg[i] = 0.5 * (accel_prev[i] + accel_curr[i]);
    }
    /* Coning correction: 1/12 * (delta_theta_prev x delta_theta_curr) */
    NAV_PRECISION coning[3];
    coning[0] = (gyro_prev[1]*gyro_curr[2] - gyro_prev[2]*gyro_curr[1]) / 12.0;
    coning[1] = (gyro_prev[2]*gyro_curr[0] - gyro_prev[0]*gyro_curr[2]) / 12.0;
    coning[2] = (gyro_prev[0]*gyro_curr[1] - gyro_prev[1]*gyro_curr[0]) / 12.0;
    for (int i = 0; i < 3; i++)
        gyro_avg[i] += coning[i];
    /* Sculling compensation */
    NAV_PRECISION scul[3];
    scul[0] = (gyro_prev[1]*accel_curr[2] - gyro_prev[2]*accel_curr[1] +
               accel_prev[1]*gyro_curr[2] - accel_prev[2]*gyro_curr[1]) / 12.0;
    scul[1] = (gyro_prev[2]*accel_curr[0] - gyro_prev[0]*accel_curr[2] +
               accel_prev[2]*gyro_curr[0] - accel_prev[0]*gyro_curr[2]) / 12.0;
    scul[2] = (gyro_prev[0]*accel_curr[1] - gyro_prev[1]*accel_curr[0] +
               accel_prev[0]*gyro_curr[1] - accel_prev[1]*gyro_curr[0]) / 12.0;
    for (int i = 0; i < 3; i++)
        accel_avg[i] += scul[i];
    nav_ins_mechanize(state, gyro_avg, accel_avg, dt);
}

int nav_ins_coarse_alignment(nav_ins_state_t *state,
                              const NAV_PRECISION accel_mean[3],
                              const NAV_PRECISION gyro_mean[3],
                              NAV_PRECISION lat_rad) {
    if (!state || !accel_mean || !gyro_mean) return -1;
    /* Leveling: pitch = asin(fx/g), roll = -asin(fy/(g*cos(pitch))) */
    NAV_PRECISION g = nav_normal_gravity(lat_rad, 0.0);
    NAV_PRECISION fx = accel_mean[0], fy = accel_mean[1];
    (void)accel_mean[2];
    NAV_PRECISION pitch = asin(fx / g);
    NAV_PRECISION roll  = -asin(fy / (g * cos(pitch)));
    /* Gyrocompassing: yaw from Earth rate */
    NAV_PRECISION omega_ie[3];
    nav_earth_rotation_ned(lat_rad, omega_ie);
    NAV_PRECISION wy = gyro_mean[1];
    NAV_PRECISION wz = gyro_mean[2];
    NAV_PRECISION yaw = atan2(-wy, wz) + lat_rad;
    /* Set initial attitude */
    nav_euler_t euler;
    euler.roll = roll; euler.pitch = pitch; euler.yaw = yaw;
    nav_euler_to_quat(&state->quat, &euler);
    nav_quat_to_dcm(&state->dcm, &state->quat);
    /* Initialize position */
    state->pos.latitude = lat_rad;
    state->pos.longitude = 0.0;
    state->pos.altitude = 0.0;
    /* Zero velocity */
    state->vel_ned.x = 0.0;
    state->vel_ned.y = 0.0;
    state->vel_ned.z = 0.0;
    /* Zero bias estimates initially */
    memset(state->gyro_bias, 0, 3 * sizeof(NAV_PRECISION));
    memset(state->accel_bias, 0, 3 * sizeof(NAV_PRECISION));
    /* Init covariance: position(10m)^2, velocity(0.1m/s)^2, attitude(1deg)^2 */
    memset(state->P, 0, 225 * sizeof(NAV_PRECISION));
    NAV_PRECISION deg1 = nav_deg2rad(1.0);
    for (int i = 0; i < 3; i++) {
        state->P[i*16] = 100.0;       /* pos: 100 m^2 */
        state->P[(3+i)*16] = 0.01;    /* vel: 0.01 (m/s)^2 */
        state->P[(6+i)*16] = deg1*deg1; /* att: rad^2 */
        state->P[(9+i)*16] = 1e-8;    /* gyro bias */
        state->P[(12+i)*16] = 1e-4;   /* accel bias */
    }
    return 0;
}

int nav_ins_zupt(nav_ins_state_t *state,
                 const NAV_PRECISION gyro[3],
                 const NAV_PRECISION accel[3],
                 NAV_PRECISION dt) {
    (void)dt;
    if (!state || !gyro || !accel) return -1;
    /* Static detector: accel magnitude near g, gyro magnitude near 0 */
    NAV_PRECISION acc_mag = sqrt(accel[0]*accel[0] + accel[1]*accel[1] +
                                  accel[2]*accel[2]);
    NAV_PRECISION gyro_mag = sqrt(gyro[0]*gyro[0] + gyro[1]*gyro[1] +
                                   gyro[2]*gyro[2]);
    NAV_PRECISION g_local = nav_normal_gravity(state->pos.latitude,
                                                state->pos.altitude);
    /* Thresholds: accel within 2% of g, gyro < 0.05 rad/s (~3 deg/s) */
    if (fabs(acc_mag - g_local) < 0.02 * g_local && gyro_mag < 0.05) {
        /* Zero velocity update: set velocity to zero */
        state->vel_ned.x = 0.0;
        state->vel_ned.y = 0.0;
        state->vel_ned.z = 0.0;
        /* Estimate gyro bias from static measurement */
        for (int i = 0; i < 3; i++)
            state->gyro_bias[i] = 0.999 * state->gyro_bias[i] + 0.001 * gyro[i];
        return 1; /* stationary detected */
    }
    return 0; /* moving */
}

void nav_ins_nhc(nav_ins_state_t *state, NAV_PRECISION R[9]) {
    if (!state || !R) return;
    /* Non-holonomic constraint: lateral and vertical velocity in body frame = 0.
     * v_body_y = 0, v_body_z = 0 (for land vehicles, x-forward).
     * In NED: R_body_to_ned * [vx, 0, 0]^T = v_ned_expected.
     * We measure: [0, 1, 0]*R_ned_to_body * v_ned = 0, similarly for z. */
    /* R is body-to-NED DCM (3x3 row-major). Compute NED-to-body = R^T. */
    NAV_PRECISION vy_body = R[3]*state->vel_ned.x + R[4]*state->vel_ned.y +
                             R[5]*state->vel_ned.z;
    NAV_PRECISION vz_body = R[6]*state->vel_ned.x + R[7]*state->vel_ned.y +
                             R[8]*state->vel_ned.z;
    /* Project velocity onto forward direction only */
    NAV_PRECISION vx_body = R[0]*state->vel_ned.x + R[1]*state->vel_ned.y +
                             R[2]*state->vel_ned.z;
    /* Reconstruct NED velocity from forward component only */
    state->vel_ned.x = R[0] * vx_body;
    state->vel_ned.y = R[1] * vx_body;
    state->vel_ned.z = R[2] * vx_body;
    (void)vy_body; (void)vz_body;
}

void nav_ins_error_transition(NAV_PRECISION Phi[225],
                               const nav_ins_state_t *state,
                               NAV_PRECISION dt) {
    (void)dt;
    if (!Phi || !state) return;
    /* Build 15x15 INS error state transition matrix (psi-angle formulation).
     * Phi = I + F * dt (first-order approximation for small dt).
     * F blocks:
     *   F_rr, F_rv, F_rpsi, 0, 0
     *   F_vr, F_vv, F_vpsi, 0, F_vba
     *   F_psir, F_psiv, F_psipsi, F_psibg, 0
     *   0, 0, 0, F_bgbg, 0
     *   0, 0, 0, 0, F_baba
     */
    NAV_PRECISION lat = state->pos.latitude;
    NAV_PRECISION alt = state->pos.altitude;
    NAV_PRECISION vn = state->vel_ned.x;
    NAV_PRECISION ve = state->vel_ned.y;
    NAV_PRECISION Rm = nav_meridian_radius(lat);
    NAV_PRECISION Rn = nav_transverse_radius(lat);
    NAV_PRECISION h_Rm = Rm + alt;
    NAV_PRECISION h_Rn = Rn + alt;
    (void)NAV_EARTH_ROTATION_RATE;
    NAV_PRECISION cl = cos(lat), tl = tan(lat);
    (void)sin(lat);
    NAV_PRECISION omega_ie[3];
    nav_earth_rotation_ned(lat, omega_ie);
    NAV_PRECISION omega_en[3];
    nav_transport_rate_ned(lat, alt, vn, ve, omega_en);
    /* Initialize to identity */
    memset(Phi, 0, 225 * sizeof(NAV_PRECISION));
    for (int i = 0; i < 15; i++) Phi[i*16] = 1.0;
    /* Position error rows */
    Phi[0*15+4] = 1.0/h_Rm * dt;
    Phi[1*15+3] = 1.0/(h_Rn*cl) * dt;
    Phi[1*15+0] = ve*tl/(h_Rn*cl) * dt;
    Phi[2*15+5] = -dt;
    /* Velocity error rows (simplified) */
    NAV_PRECISION f_n = 0.0, f_e = 0.0, f_d = -nav_normal_gravity(lat, alt);
    Phi[3*15+6] = -f_d * dt; Phi[3*15+7] =  f_e * dt;
    Phi[3*15+13] = dt;
    Phi[4*15+6] =  f_d * dt; Phi[4*15+8] = -f_n * dt;
    Phi[4*15+14] = dt;
    Phi[5*15+6] = -f_e * dt; Phi[5*15+7] =  f_n * dt;
    Phi[5*15+12] = dt;
    /* Attitude error: psi_dot = -omega_in^n x psi */
    Phi[6*15+9] = -dt;
    Phi[7*15+10] = -dt;
    Phi[8*15+11] = -dt;
    /* Bias models: first-order Gauss-Markov */
    NAV_PRECISION tau_g = 3600.0; /* 1 hour correlation time */
    NAV_PRECISION tau_a = 3600.0;
    Phi[9*15+9]   = 1.0 - dt/tau_g;
    Phi[10*15+10] = 1.0 - dt/tau_g;
    Phi[11*15+11] = 1.0 - dt/tau_g;
    Phi[12*15+12] = 1.0 - dt/tau_a;
    Phi[13*15+13] = 1.0 - dt/tau_a;
    Phi[14*15+14] = 1.0 - dt/tau_a;
}

void nav_ins_process_noise(NAV_PRECISION Qd[225],
                            const nav_ins_state_t *state,
                            NAV_PRECISION gyro_arw,
                            NAV_PRECISION accel_vrw,
                            NAV_PRECISION gyro_bi,
                            NAV_PRECISION accel_bi,
                            NAV_PRECISION dt) {
    (void)dt;
    if (!Qd || !state) return;
    memset(Qd, 0, 225 * sizeof(NAV_PRECISION));
    /* Simplified discrete-time process noise.
     * For ARW/VRW: Q = sigma^2 * dt * I (for velocity/attitude states)
     * For bias instability: Q = 2*sigma^2/tau * dt * I */
    NAV_PRECISION q_att = gyro_arw * gyro_arw * dt;
    NAV_PRECISION q_vel = accel_vrw * accel_vrw * dt;
    NAV_PRECISION q_gb  = 2.0 * gyro_bi * gyro_bi / 3600.0 * dt;
    NAV_PRECISION q_ab  = 2.0 * accel_bi * accel_bi / 3600.0 * dt;
    for (int i = 0; i < 3; i++) {
        Qd[(6+i)*16] = q_att;
        Qd[(3+i)*16] = q_vel;
        Qd[(9+i)*16] = q_gb;
        Qd[(12+i)*16] = q_ab;
    }
}
