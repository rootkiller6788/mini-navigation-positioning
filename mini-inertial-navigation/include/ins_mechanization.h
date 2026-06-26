#ifndef INS_MECHANIZATION_H
#define INS_MECHANIZATION_H
#include "ins_core.h"
#include "ins_attitude.h"
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L5/L6: Strapdown Inertial Navigation System Mechanization
 *
 * The strapdown INS mechanization propagates the navigation state
 * (attitude, velocity, position) by integrating inertial sensor
 * measurements (specific force and angular rate) over time.
 *
 * Core Processing Chain:
 *   1. Attitude Update: integrate gyro measurements to update quaternion
 *   2. Specific Force Transform: rotate accelerometer readings to NED
 *   3. Velocity Update: integrate NED-frame specific force (accounting for
 *      gravity, Coriolis, and transport rate)
 *   4. Position Update: integrate velocity to update lat/lon/alt
 *
 * Reference: Titterton & Weston (2004), Chapters 11-12.
 * Reference: Savage (2000), "Strapdown Analytics", Strapdown Associates.
 * ========================================================================= */

typedef enum {
    INS_MECH_EULER    = 0,
    INS_MECH_RK2      = 1,
    INS_MECH_RK4      = 2,
    INS_MECH_CONING   = 3
} ins_mech_method_t;

typedef struct {
    ins_mech_method_t method;
    ins_quat_t        quat;
    ins_vec3_t        vel_ned;
    ins_geodetic_t    pos;
    ins_vec3_t        omega_prev;
    ins_vec3_t        accel_prev;
    double            dt_prev;
    int               initialized;
    unsigned long     iteration;
} ins_mech_state_t;

/* -------------------------------------------------------------------------
 * L5: Strapdown Mechanization Functions
 * ------------------------------------------------------------------------- */

/** Initialize mechanization state. */
void ins_mech_init(ins_mech_state_t *state, ins_mech_method_t method,
                   double init_lat, double init_lon, double init_alt);

/**
 * Align INS by setting initial attitude from level accelerometer readings
 * and known heading. Assumes stationary: accelerometer measures gravity.
 *
 * roll  = atan2(-a_y, -a_z)
 * pitch = atan2(a_x, sqrt(a_y^2 + a_z^2))
 * yaw   = heading_input
 *
 * @return 0 on success, -1 if accel near zero
 */
int ins_align_level(ins_mech_state_t *state, const ins_vec3_t *accel_body,
                     double heading);

/**
 * Main strapdown mechanization step.
 *
 * Velocity Update Equation:
 *   dv^n/dt = C_b^n * f^b - (2*w_ie^n + w_en^n) x v^n + g^n
 *
 * Position Update Equations:
 *   d(lat)/dt = v_north / (M + h)
 *   d(lon)/dt = v_east / ((N + h) * cos(lat))
 *   d(alt)/dt = -v_down
 *
 * @return 0 on success, -1 on numerical error
 */
int ins_mech_step(ins_mech_state_t *state, const ins_imu_sample_t *imu);

/** Extract navigation solution from mechanization state. */
void ins_mech_get_solution(const ins_mech_state_t *state,
                            ins_nav_solution_t *nav);

/* -------------------------------------------------------------------------
 * L5: Standalone Velocity/Position Update Functions
 * ------------------------------------------------------------------------- */

/**
 * Single velocity update step.
 * v(t+dt) = v(t) + (C_b^n * f^b + g^n - (2*w_ie^n + w_en^n) x v^n) * dt
 */
void ins_vel_update(ins_vec3_t *vel_ned, const ins_mat3_t *C_body_ned,
                    const ins_vec3_t *f_body, double lat, double alt,
                    double dt);

/**
 * Single position update step from NED velocity.
 * lat(t+dt) = lat(t) + v_north * dt / (M + h)
 * lon(t+dt) = lon(t) + v_east * dt / ((N + h) * cos(lat))
 * alt(t+dt) = alt(t) - v_down * dt
 */
void ins_pos_update(ins_geodetic_t *pos, const ins_vec3_t *vel_ned, double dt);

/* -------------------------------------------------------------------------
 * L6: Canonical Problem — Complete INS Navigation Loop
 * ------------------------------------------------------------------------- */

/**
 * Run complete INS navigation: process IMU samples to produce trajectory.
 * This is the canonical endpoint problem for strapdown INS.
 */
int ins_navigate(const ins_mech_state_t *state,
                 const ins_imu_sample_t *imu_data,
                 size_t num_samples,
                 ins_nav_solution_t *trajectory);

#ifdef __cplusplus
}
#endif
#endif /* INS_MECHANIZATION_H */
