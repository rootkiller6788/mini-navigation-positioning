/**
 * @file nav_common.h
 * @brief Integrated Navigation - Common Definitions & Types
 *
 * Covers L1: Core definitions for navigation states, frames, and errors.
 * Reference: Jekeli, "Inertial Navigation Systems with Geodetic Applications"
 *            Groves, "Principles of GNSS, Inertial, and Multisensor Navigation"
 */

#ifndef NAV_COMMON_H
#define NAV_COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef NAV_PRECISION
#define NAV_PRECISION           double
#endif

#define NAV_PI                  3.14159265358979323846
#define NAV_DEG2RAD             (NAV_PI / 180.0)
#define NAV_RAD2DEG             (180.0 / NAV_PI)
#define NAV_EARTH_RADIUS_M      6378137.0
#define NAV_EARTH_FLATTENING    (1.0 / 298.257223563)
#define NAV_EARTH_ECCENTRICITY2 (6.69437999014e-3)
#define NAV_EARTH_ROTATION_RATE 7.2921151467e-5
#define NAV_GRAVITY_EQUATOR     9.7803267715
#define NAV_GRAVITY_POLE        9.8321863685
#define NAV_C_LIGHT             299792458.0

#define NAV_STATE_DIM           15
#define NAV_IMU_MEAS_DIM        6
#define NAV_GNSS_MEAS_DIM       4
#define NAV_QUAT_DIM            4
#define NAV_ROT_MATRIX_DIM      9

typedef enum {
    NAV_FRAME_ECI = 0,
    NAV_FRAME_ECEF,
    NAV_FRAME_NED,
    NAV_FRAME_ENU,
    NAV_FRAME_BODY,
    NAV_FRAME_WANDER,
    NAV_FRAME_COUNT
} nav_frame_t;

typedef struct {
    NAV_PRECISION latitude;
    NAV_PRECISION longitude;
    NAV_PRECISION altitude;
} nav_geodetic_t;

typedef struct {
    NAV_PRECISION x;
    NAV_PRECISION y;
    NAV_PRECISION z;
} nav_vector3_t;

typedef struct {
    NAV_PRECISION dr_n;    NAV_PRECISION dr_e;    NAV_PRECISION dr_d;
    NAV_PRECISION dv_n;    NAV_PRECISION dv_e;    NAV_PRECISION dv_d;
    NAV_PRECISION psi_n;   NAV_PRECISION psi_e;   NAV_PRECISION psi_d;
    NAV_PRECISION bg_x;    NAV_PRECISION bg_y;    NAV_PRECISION bg_z;
    NAV_PRECISION ba_x;    NAV_PRECISION ba_y;    NAV_PRECISION ba_z;
} nav_ins_error_t;

typedef struct {
    NAV_PRECISION gyro_x;
    NAV_PRECISION gyro_y;
    NAV_PRECISION gyro_z;
    NAV_PRECISION accel_x;
    NAV_PRECISION accel_y;
    NAV_PRECISION accel_z;
    NAV_PRECISION dt;
    uint64_t      timestamp;
} nav_imu_meas_t;

typedef struct {
    NAV_PRECISION bias;
    NAV_PRECISION scale_factor_ppm;
    NAV_PRECISION noise_psd;
    NAV_PRECISION bias_instability;
    NAV_PRECISION misalignment[3];
    NAV_PRECISION g_sensitivity[3];
} nav_imu_error_t;

typedef struct {
    uint8_t       sv_id;
    NAV_PRECISION pseudorange;
    NAV_PRECISION pseudorange_rate;
    NAV_PRECISION carrier_phase;
    NAV_PRECISION cn0;
    NAV_PRECISION sat_x, sat_y, sat_z;
    NAV_PRECISION sat_vx, sat_vy, sat_vz;
    NAV_PRECISION iono_delay;
    NAV_PRECISION tropo_delay;
    NAV_PRECISION sat_clk_bias;
    uint64_t      tow;
} nav_gnss_sv_t;

typedef struct {
    nav_geodetic_t pos;
    nav_vector3_t  vel_enu;
    NAV_PRECISION  clock_bias;
    NAV_PRECISION  clock_drift;
    NAV_PRECISION  gdop, pdop, hdop, vdop, tdop;
    NAV_PRECISION  pos_sigma_e, pos_sigma_n, pos_sigma_u;
    uint8_t        num_svs;
    uint64_t       timestamp;
} nav_gnss_solution_t;

typedef struct {
    nav_geodetic_t pos;
    nav_vector3_t  vel_ned;
    nav_vector3_t  vel_enu;
    NAV_PRECISION  roll, pitch, yaw;
    NAV_PRECISION  quat[4];
    NAV_PRECISION  gyro_bias[3];
    NAV_PRECISION  accel_bias[3];
    NAV_PRECISION  pos_cov[9];
    NAV_PRECISION  vel_cov[9];
    NAV_PRECISION  att_cov[9];
    uint64_t       timestamp;
} nav_solution_t;

static inline NAV_PRECISION nav_deg2rad(NAV_PRECISION deg) {
    return deg * NAV_DEG2RAD;
}

static inline NAV_PRECISION nav_rad2deg(NAV_PRECISION rad) {
    return rad * NAV_RAD2DEG;
}

static inline NAV_PRECISION nav_wrap_pi(NAV_PRECISION angle) {
    while (angle > NAV_PI)  angle -= 2.0 * NAV_PI;
    while (angle < -NAV_PI) angle += 2.0 * NAV_PI;
    return angle;
}

static inline NAV_PRECISION nav_wrap_2pi(NAV_PRECISION angle) {
    angle = fmod(angle, 2.0 * NAV_PI);
    if (angle < 0.0) angle += 2.0 * NAV_PI;
    return angle;
}

static inline NAV_PRECISION nav_meridian_radius(NAV_PRECISION lat_rad) {
    NAV_PRECISION s = sin(lat_rad);
    NAV_PRECISION e2_s2 = NAV_EARTH_ECCENTRICITY2 * s * s;
    return NAV_EARTH_RADIUS_M * (1.0 - NAV_EARTH_ECCENTRICITY2) /
           pow(1.0 - e2_s2, 1.5);
}

static inline NAV_PRECISION nav_transverse_radius(NAV_PRECISION lat_rad) {
    NAV_PRECISION s = sin(lat_rad);
    return NAV_EARTH_RADIUS_M / sqrt(1.0 - NAV_EARTH_ECCENTRICITY2 * s * s);
}

static inline NAV_PRECISION nav_normal_gravity(NAV_PRECISION lat_rad,
                                                NAV_PRECISION alt_m) {
    NAV_PRECISION s2 = sin(lat_rad) * sin(lat_rad);
    NAV_PRECISION k = 0.00193185265241;
    NAV_PRECISION e2 = NAV_EARTH_ECCENTRICITY2;
    NAV_PRECISION g0 = NAV_GRAVITY_EQUATOR * (1.0 + k * s2) /
                       sqrt(1.0 - e2 * s2);
    NAV_PRECISION fac = -3.086e-6;
    return g0 + fac * alt_m;
}

static inline void nav_earth_rotation_ned(NAV_PRECISION lat_rad,
                                           NAV_PRECISION omega[3]) {
    omega[0] = NAV_EARTH_ROTATION_RATE * cos(lat_rad);
    omega[1] = 0.0;
    omega[2] = -NAV_EARTH_ROTATION_RATE * sin(lat_rad);
}

static inline void nav_transport_rate_ned(NAV_PRECISION lat_rad,
                                           NAV_PRECISION alt_m,
                                           NAV_PRECISION vn, NAV_PRECISION ve,
                                           NAV_PRECISION rho[3]) {
    NAV_PRECISION Rm = nav_meridian_radius(lat_rad);
    NAV_PRECISION Rn = nav_transverse_radius(lat_rad);
    rho[0] = ve / (Rn + alt_m);
    rho[1] = -vn / (Rm + alt_m);
    rho[2] = -ve * tan(lat_rad) / (Rn + alt_m);
}

#endif /* NAV_COMMON_H */
