#ifndef INS_CORE_H
#define INS_CORE_H
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Core Definitions — Inertial Navigation System Types
 *
 * Reference: Titterton & Weston (2004), "Strapdown Inertial Navigation
 *   Technology", 2nd ed., IET Radar Sonar Navig.
 * Reference: Groves (2013), "Principles of GNSS, Inertial, and
 *   Multisensor Integrated Navigation Systems", 2nd ed., Artech House.
 * Reference: Farrell (2008), "Aided Navigation: GPS with High Rate Sensors",
 *   McGraw-Hill.
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * L1: Coordinate Frames
 * ------------------------------------------------------------------------- */

typedef enum {
    INS_FRAME_ECI   = 0,
    INS_FRAME_ECEF  = 1,
    INS_FRAME_NED   = 2,
    INS_FRAME_ENU   = 3,
    INS_FRAME_BODY  = 4,
    INS_FRAME_PLATFORM = 5
} ins_frame_t;

/* -------------------------------------------------------------------------
 * L1: Earth Model Parameters (WGS84 Ellipsoid)
 * ------------------------------------------------------------------------- */

#define INS_WGS84_A    6378137.0
#define INS_WGS84_F    (1.0 / 298.257223563)
#define INS_WGS84_B    (INS_WGS84_A * (1.0 - INS_WGS84_F))
#define INS_WGS84_E2   (2.0 * INS_WGS84_F - INS_WGS84_F * INS_WGS84_F)
#define INS_WGS84_EP2  (INS_WGS84_E2 / (1.0 - INS_WGS84_E2))
#define INS_WGS84_WE   (7.2921151467e-5)
#define INS_WGS84_GM   (3.986004418e14)

#define INS_GRAVITY_EQUATOR  9.7803253359
#define INS_GRAVITY_POLE     9.8321849378
#define INS_K_SOMIGLIANA     ((INS_WGS84_B * INS_GRAVITY_POLE) / (INS_WGS84_A * INS_GRAVITY_EQUATOR) - 1.0)

/* Fallback for M_PI if not provided by math.h */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* -------------------------------------------------------------------------
 * L1: 3D Vector and Matrix Types
 * ------------------------------------------------------------------------- */

typedef struct {
    double x, y, z;
} ins_vec3_t;

typedef struct {
    double m[9];
} ins_mat3_t;

/* -------------------------------------------------------------------------
 * L1: Geodetic Position
 * ------------------------------------------------------------------------- */

typedef struct {
    double lat;
    double lon;
    double alt;
} ins_geodetic_t;

typedef struct {
    double x, y, z;
} ins_ecef_t;

typedef struct {
    double north;
    double east;
    double down;
} ins_ned_t;

/* -------------------------------------------------------------------------
 * L1: Navigation State Storage
 * ------------------------------------------------------------------------- */

typedef struct {
    ins_geodetic_t pos;
    ins_vec3_t     vel_ned;
    ins_vec3_t     vel_ecef;
    ins_vec3_t     vel_body;
    double         roll;
    double         pitch;
    double         yaw;
    double         q[4];
    double         time;
} ins_nav_solution_t;

/* -------------------------------------------------------------------------
 * L1: IMU Sensor Configuration
 * ------------------------------------------------------------------------- */

typedef enum {
    INS_GRADE_CONSUMER   = 0,
    INS_GRADE_INDUSTRIAL = 1,
    INS_GRADE_TACTICAL   = 2,
    INS_GRADE_NAVIGATION = 3,
    INS_GRADE_STRATEGIC  = 4
} ins_grade_t;

typedef struct {
    ins_vec3_t accel;
    ins_vec3_t gyro;
    double     temperature;
    double     dt;
} ins_imu_sample_t;

/* -------------------------------------------------------------------------
 * L2: Earth/Gravity Model Functions
 * ------------------------------------------------------------------------- */

/**
 * Compute normal gravity at given geodetic latitude using Somigliana formula.
 * g(phi) = g_e * (1 + k*sin^2(phi)) / sqrt(1 - e^2*sin^2(phi))
 * Reference: Hofmann-Wellenhof & Moritz (2006), "Physical Geodesy", 2nd ed.
 */
double ins_gravity_wgs84(double lat, double alt);

/**
 * Compute gravity vector in NED frame accounting for gravitation and
 * centrifugal acceleration due to Earth rotation.
 */
void ins_gravity_ned(double lat, double alt, ins_vec3_t *g_ned);

/* -------------------------------------------------------------------------
 * L2: Coordinate Transform Functions
 * ------------------------------------------------------------------------- */

/** Convert geodetic to ECEF Cartesian using standard closed-form. */
void ins_geodetic_to_ecef(const ins_geodetic_t *geo, ins_ecef_t *ecef);

/** Convert ECEF Cartesian to geodetic using Bowring/Zhu iterative method. */
void ins_ecef_to_geodetic(const ins_ecef_t *ecef, ins_geodetic_t *geo);

/** Compute ECEF-to-NED rotation matrix: C_ned^ecef = R_y(-lat-pi/2) * R_z(lon) */
void ins_ecef_to_ned_dcm(double lat, double lon, ins_mat3_t *C);

/**
 * Transport rate: angular rate of NED frame relative to ECEF.
 * omega_en^n = [lon_dot*cos(lat), -lat_dot, -lon_dot*sin(lat)]^T
 */
void ins_transport_rate(const ins_vec3_t *v_ned, double lat, double alt,
                        ins_vec3_t *omega);

/** Earth rotation rate in NED frame: omega_ie^n = [we*cos(lat), 0, -we*sin(lat)]^T */
void ins_earth_rate_ned(double lat, ins_vec3_t *omega);

/* -------------------------------------------------------------------------
 * L2: Vector/Matrix Utilities
 * ------------------------------------------------------------------------- */

void ins_vec3_zero(ins_vec3_t *v);
void ins_vec3_set(ins_vec3_t *v, double x, double y, double z);
void ins_vec3_copy(const ins_vec3_t *src, ins_vec3_t *dst);
void ins_vec3_add(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c);
void ins_vec3_sub(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c);
void ins_vec3_scale(const ins_vec3_t *a, double s, ins_vec3_t *b);
double ins_vec3_dot(const ins_vec3_t *a, const ins_vec3_t *b);
void ins_vec3_cross(const ins_vec3_t *a, const ins_vec3_t *b, ins_vec3_t *c);
double ins_vec3_norm(const ins_vec3_t *v);
int ins_vec3_normalize(ins_vec3_t *v);
void ins_vec3_skew(const ins_vec3_t *v, ins_mat3_t *S);
void ins_mat3_mul_vec(const ins_mat3_t *A, const ins_vec3_t *x, ins_vec3_t *y);
void ins_mat3_mul(const ins_mat3_t *A, const ins_mat3_t *B, ins_mat3_t *C);
void ins_mat3_transpose(const ins_mat3_t *A, ins_mat3_t *AT);
void ins_mat3_identity(ins_mat3_t *M);
double ins_vec3_angle(const ins_vec3_t *a, const ins_vec3_t *b);

/* -------------------------------------------------------------------------
 * L2: Meridian Radius and Related Quantities
 * ------------------------------------------------------------------------- */

/** Meridian radius of curvature: M = a*(1-e^2)/(1-e^2*sin^2(lat))^(3/2) */
double ins_meridian_radius(double lat);

/** Prime vertical radius of curvature: N = a/sqrt(1-e^2*sin^2(lat)) */
double ins_prime_vertical_radius(double lat);

/* -------------------------------------------------------------------------
 * L3: Coriolis Compensation
 * ------------------------------------------------------------------------- */

/**
 * Coriolis acceleration in NED frame.
 * a_coriolis = -(2*omega_ie^n + omega_en^n) x v^n
 * Essential for accurate strapdown velocity integration.
 */
void ins_coriolis_accel(const ins_vec3_t *v_ned, double lat, double alt,
                        ins_vec3_t *a_coriolis);

/* -------------------------------------------------------------------------
 * Utility Constants
 * ------------------------------------------------------------------------- */

#define INS_DEG2RAD   (M_PI / 180.0)
#define INS_RAD2DEG   (180.0 / M_PI)
#define INS_EARTH_RADIUS_MEAN  6371000.0
#define INS_SCHULER_PERIOD     5067.0
#define INS_SCHULER_FREQ       (2.0 * M_PI / INS_SCHULER_PERIOD)

#ifdef __cplusplus
}
#endif
#endif /* INS_CORE_H */
