#ifndef GNSS_COMMON_H
#define GNSS_COMMON_H
#include <stddef.h>
#include <stdint.h>
#define _USE_MATH_DEFINES
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * L1: Core GNSS Constants — defined per ICD-GPS-200 & WGS84
 * ========================================================================= */

#define GNSS_C_LIGHT          299792458.0       /* Speed of light [m/s] */
#define GNSS_L1_FREQ          1575.42e6         /* GPS L1 carrier [Hz] */
#define GNSS_L2_FREQ          1227.60e6         /* GPS L2 carrier [Hz] */
#define GNSS_L5_FREQ          1176.45e6         /* GPS L5 carrier [Hz] */
#define GNSS_L1_WAVELENGTH    (GNSS_C_LIGHT / GNSS_L1_FREQ)  /* ~0.1903 m */
#define GNSS_L2_WAVELENGTH    (GNSS_C_LIGHT / GNSS_L2_FREQ)  /* ~0.2442 m */
#define GNSS_L5_WAVELENGTH    (GNSS_C_LIGHT / GNSS_L5_FREQ)  /* ~0.2548 m */
#define GNSS_CA_CHIP_RATE     1.023e6           /* C/A code chip rate [cps] */
#define GNSS_CA_CODE_LENGTH   1023              /* C/A code chips */
#define GNSS_CA_PERIOD        0.001             /* C/A code period [s] */

/* WGS84 Ellipsoid Parameters (NIMA TR8350.2) */
#define GNSS_WGS84_A          6378137.0         /* Semi-major axis [m] */
#define GNSS_WGS84_F          0.0033528106647474805  /* 1 / 298.257223563 */
#define GNSS_WGS84_E2         (2.0 * GNSS_WGS84_F - GNSS_WGS84_F * GNSS_WGS84_F)
#define GNSS_WGS84_EP2        (GNSS_WGS84_E2 / (1.0 - GNSS_WGS84_E2))
#define GNSS_WGS84_B          (GNSS_WGS84_A * (1.0 - GNSS_WGS84_F)) /* polar semi-axis */

/* Earth rotation */
#define GNSS_OMEGA_E          7.2921151467e-5   /* Earth rotation rate [rad/s] */
#define GNSS_GM_EARTH         3.986005e14       /* WGS84 GM [m³/s²] */

/* GPS time */
#define GNSS_SECONDS_PER_WEEK 604800
#define GNSS_GPS_EPOCH_JD     2444244.5         /* GPS epoch: 1980-01-06 */

/* -------------------------------------------------------------------------
 * L1: Coordinate Systems — three fundamental frames
 * ------------------------------------------------------------------------- */

/**
 * @brief Earth-Centered Earth-Fixed (ECEF) Cartesian coordinates
 *
 * X-axis: intersection of prime meridian and equator
 * Y-axis: 90° east in equator plane
 * Z-axis: north along rotation axis (right-handed)
 *
 * WGS84 reference frame, units in meters.
 */
typedef struct {
    double x, y, z;     /* Cartesian coordinates [m] */
} gnss_ecef_t;

/**
 * @brief Geodetic coordinates (Latitude, Longitude, Altitude)
 *
 * Latitude φ ∈ [-π/2, π/2], Longitude λ ∈ [-π, π) in radians.
 * Altitude h: height above WGS84 ellipsoid [m].
 */
typedef struct {
    double lat;         /* Geodetic latitude [rad] */
    double lon;         /* Geodetic longitude [rad] */
    double alt;         /* Geodetic height above ellipsoid [m] */
} gnss_lla_t;

/**
 * @brief East-North-Up (ENU) local tangent plane coordinates
 *
 * Origin at a reference LLA point. East-x, North-y, Up-z form a
 * right-handed local Cartesian system.
 */
typedef struct {
    double e, n, u;     /* East, North, Up [m] */
} gnss_enu_t;

/* -------------------------------------------------------------------------
 * L1: GPS Time Structures
 * ------------------------------------------------------------------------- */

/** GPS week number and seconds-of-week */
typedef struct {
    int32_t  week;          /* GPS week number (modulo 1024 or 8192) */
    double   sow;           /* Seconds of week [s] (0 .. 604800) */
} gnss_gpstime_t;

/** UTC time with leap-second awareness */
typedef struct {
    int16_t  year, month, day;
    int16_t  hour, minute;
    double   second;
    int16_t  leap_seconds;  /* ΔT = GPS - UTC (currently +18 s) */
} gnss_utctime_t;

/* -------------------------------------------------------------------------
 * L1: Satellite ID and constellation enums
 * ------------------------------------------------------------------------- */

/** GNSS constellation type */
typedef enum {
    GNSS_CONST_GPS       = 0,
    GNSS_CONST_GLONASS   = 1,
    GNSS_CONST_GALILEO   = 2,
    GNSS_CONST_BEIDOU    = 3,
    GNSS_CONST_QZSS      = 4,
    GNSS_CONST_NAVIC     = 5,
    GNSS_CONST_SBAS      = 6
} gnss_constellation_t;

/** Per-satellite identification within a constellation */
typedef struct {
    gnss_constellation_t constellation;
    int32_t              prn;         /* PRN / slot number */
    int32_t              svn;         /* Space Vehicle Number */
} gnss_satid_t;

/* -------------------------------------------------------------------------
 * L2: WGS84 Ellipsoid descriptor (captures the datum)
 * ------------------------------------------------------------------------- */

typedef struct {
    double a;           /* Semi-major axis [m] */
    double f;           /* Flattening */
    double e2;          /* First eccentricity squared */
    double ep2;         /* Second eccentricity squared = e²/(1-e²) */
    double b;           /* Semi-minor axis [m] */
} gnss_ellipsoid_t;

/* -------------------------------------------------------------------------
 * L3: 3×3 matrix — used for coordinate rotations and DOP
 * ------------------------------------------------------------------------- */

typedef struct {
    double m[3][3];
} gnss_mat33_t;

/* -------------------------------------------------------------------------
 * L3: 4×4 matrix — used for augmented least squares (position + clock)
 * ------------------------------------------------------------------------- */

typedef struct {
    double m[4][4];
} gnss_mat44_t;

/* -------------------------------------------------------------------------
 * API: Coordinate conversions
 * ------------------------------------------------------------------------- */

/**
 * @brief LLA → ECEF conversion (closed-form)
 *
 * Uses standard geodetic formula:
 *   N(φ) = a / √(1 - e²·sin²φ)
 *   x = (N + h)·cos φ·cos λ
 *   y = (N + h)·cos φ·sin λ
 *   z = (N·(1-e²) + h)·sin φ
 */
gnss_ecef_t gnss_lla_to_ecef(gnss_lla_t lla);
gnss_ecef_t gnss_lla_to_ecef_ell(gnss_lla_t lla, const gnss_ellipsoid_t *ell);

/**
 * @brief ECEF → LLA conversion (Bowring 1985 iterative method)
 *
 * Iteratively refines latitude until convergence (typically 2-3 iterations).
 * Handles polar singularities gracefully.
 *
 * Reference: Bowring, B.R. (1985). "The accuracy of geodetic latitude
 * and height equations". Survey Review, 28(218), 202-206.
 */
gnss_lla_t gnss_ecef_to_lla(gnss_ecef_t ecef);
gnss_lla_t gnss_ecef_to_lla_ell(gnss_ecef_t ecef, const gnss_ellipsoid_t *ell);

/**
 * @brief ECEF → ENU at reference point (closed-form rotation)
 *
 * ENU = R·Δecef where R rotates from ECEF axes to local ENU.
 */
gnss_enu_t gnss_ecef_to_enu(gnss_ecef_t ecef, gnss_lla_t ref);
gnss_ecef_t gnss_enu_to_ecef(gnss_enu_t enu,   gnss_lla_t ref);

/** @brief Compute prime vertical radius of curvature N(φ) */
double gnss_rad_curvature_prime(double lat_rad);
double gnss_rad_curvature_prime_ell(double lat_rad, const gnss_ellipsoid_t *ell);

/** @brief Compute meridian radius of curvature M(φ) */
double gnss_rad_curvature_meridian(double lat_rad);
double gnss_rad_curvature_meridian_ell(double lat_rad, const gnss_ellipsoid_t *ell);

/**
 * @brief Sagnac-effect Earth rotation correction for pseudorange
 *
 * During signal travel time τ, the Earth rotates by ω_e·τ.
 * The satellite position must be rotated back:
 *   correction = (ω_e / c) · (x_s·y_r - y_s·x_r)
 */
double gnss_sagnac_correction(gnss_ecef_t sat_pos, gnss_ecef_t rx_pos);
double gnss_sagnac_rotate(gnss_ecef_t sat_pos, double travel_time);

/* -------------------------------------------------------------------------
 * API: Time conversions
 * ------------------------------------------------------------------------- */

/** @brief Convert GPS week + seconds-of-week to UTC */
gnss_utctime_t gnss_gpstime_to_utc(gnss_gpstime_t gps, int16_t leap_sec);

/** @brief Convert calendar UTC to GPS time */
gnss_gpstime_t gnss_utc_to_gpstime(gnss_utctime_t utc);

/** @brief Compute Julian Date from UTC */
double gnss_utc_to_jd(gnss_utctime_t utc);

/** @brief Compute geometric range ||sat - rx|| */
double gnss_geometric_range(gnss_ecef_t sat, gnss_ecef_t rx);

/** @brief Compute satellite elevation angle at receiver [rad] */
double gnss_sat_elevation(gnss_ecef_t sat, gnss_ecef_t rx);

/** @brief Compute satellite azimuth angle at receiver [rad] */
double gnss_sat_azimuth(gnss_ecef_t sat, gnss_ecef_t rx);

/** @brief Get WGS84 ellipsoid descriptor */
gnss_ellipsoid_t gnss_wgs84_ellipsoid(void);

/* =========================================================================
 * L3: Vector/matrix operations (3-vector, 3×3, 4×4)
 * ========================================================================= */

/** 3-vector operations */
void   gnss_vec3_sub(const double a[3], const double b[3], double result[3]);
void   gnss_vec3_add(const double a[3], const double b[3], double result[3]);
void   gnss_vec3_scale(const double v[3], double s, double result[3]);
double gnss_vec3_dot(const double a[3], const double b[3]);
void   gnss_vec3_cross(const double a[3], const double b[3], double result[3]);
double gnss_vec3_norm(const double v[3]);
void   gnss_vec3_normalize(const double v[3], double result[3]);

/** 3×3 matrix operations */
void gnss_mat33_identity(gnss_mat33_t *m);
void gnss_mat33_multiply(const gnss_mat33_t *a, const gnss_mat33_t *b, gnss_mat33_t *c);
void gnss_mat33_transpose(const gnss_mat33_t *a, gnss_mat33_t *at);
void gnss_mat33_vec_multiply(const gnss_mat33_t *a, const double v[3], double result[3]);
int  gnss_mat33_inverse(const gnss_mat33_t *a, gnss_mat33_t *ainv);

/** 4×4 matrix operations */
void gnss_mat44_identity(gnss_mat44_t *m);
void gnss_mat44_multiply(const gnss_mat44_t *a, const gnss_mat44_t *b, gnss_mat44_t *c);
void gnss_mat44_transpose(const gnss_mat44_t *a, gnss_mat44_t *at);
void gnss_mat44_vec_multiply(const gnss_mat44_t *a, const double v[4], double result[4]);
int  gnss_mat44_inverse(gnss_mat44_t *a, gnss_mat44_t *inv);
int  gnss_mat44_cholesky(const gnss_mat44_t *a, gnss_mat44_t *L);

#ifdef __cplusplus
}
#endif
#endif /* GNSS_COMMON_H */
