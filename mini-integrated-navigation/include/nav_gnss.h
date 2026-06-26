/**
 * @file nav_gnss.h
 * @brief GNSS Measurement Models and Position Solution
 *
 * L2 Core Concepts: GNSS positioning, pseudorange measurement model.
 * L5 Algorithms: Weighted Least Squares PVT, satellite position computation.
 * L6 Canonical Problem: GNSS standalone positioning.
 *
 * Reference: Misra & Enge, "Global Positioning System: Signals,
 * Measurements, and Performance"
 *            Tsui, "Fundamentals of Global Positioning System Receivers"
 */

#ifndef NAV_GNSS_H
#define NAV_GNSS_H

#include "nav_common.h"

/* ---------- GNSS Configuration --------------------------------------- */

typedef enum {
    NAV_GNSS_GPS = 0,
    NAV_GNSS_GLONASS,
    NAV_GNSS_GALILEO,
    NAV_GNSS_BEIDOU,
    NAV_GNSS_COUNT
} nav_gnss_constellation_t;

typedef struct {
    nav_gnss_constellation_t constellation;
    NAV_PRECISION iono_alpha[4];  /* Klobuchar model parameters */
    NAV_PRECISION iono_beta[4];
    NAV_PRECISION tropo_refractivity;
    NAV_PRECISION elevation_mask; /* rad, min elevation */
    NAV_PRECISION cn0_mask;       /* dB-Hz, min signal strength */
    int           use_iono_model;
    int           use_tropo_model;
} nav_gnss_config_t;

/* ---------- Satellite Ephemeris (simplified Keplerian) --------------- */

typedef struct {
    uint8_t       sv_id;
    uint32_t      week;
    NAV_PRECISION toe;             /* time of ephemeris, s */
    NAV_PRECISION sqrt_a;         /* sqrt(semi-major axis), m^0.5 */
    NAV_PRECISION ecc;            /* eccentricity */
    NAV_PRECISION i0;             /* inclination at reference, rad */
    NAV_PRECISION omega0;         /* longitude of ascending node, rad */
    NAV_PRECISION w;              /* argument of perigee, rad */
    NAV_PRECISION M0;             /* mean anomaly at reference, rad */
    NAV_PRECISION delta_n;        /* mean motion correction, rad/s */
    NAV_PRECISION idot;           /* rate of inclination, rad/s */
    NAV_PRECISION omegadot;       /* rate of LAN, rad/s */
    NAV_PRECISION cuc, cus;       /* latitude argument corrections, rad */
    NAV_PRECISION crc, crs;       /* radius corrections, m */
    NAV_PRECISION cic, cis;       /* inclination corrections, rad */
    NAV_PRECISION af0, af1, af2;  /* clock correction polynomial */
    NAV_PRECISION tgd;            /* group delay, s */
    int           healthy;
} nav_gnss_ephemeris_t;

/* ---------- GNSS PVT Solution Functions ------------------------------ */

/**
 * @brief Weighted Least Squares GNSS position solution.
 *
 *  Solves for [x, y, z, c*dt] using pseudorange measurements.
 *  Iterative Gauss-Newton method for linearization around approximate position.
 *
 *  d_rho = sqrt((x_sv - x_r)^2 + (y_sv - y_r)^2 + (z_sv - z_r)^2) + c*dt
 *
 * @param solution [out] computed GNSS solution
 * @param svs [in] array of satellite measurements
 * @param n_svs [in] number of satellites (>= 4)
 * @param config [in] GNSS configuration
 * @return 0 on success, -1 if insufficient satellites or no convergence
 */
int nav_gnss_wls_pvt(nav_gnss_solution_t *solution,
                     const nav_gnss_sv_t *svs,
                     int n_svs,
                     const nav_gnss_config_t *config);

/**
 * @brief Compute satellite ECEF position from ephemeris at given time.
 *
 *  Implements the full GPS ICD-200 ephemeris model.
 *
 * @param pos [out] ECEF position {x, y, z}, m
 * @param vel [out] ECEF velocity, m/s (NULL to skip)
 * @param clk [out] satellite clock correction, m (NULL to skip)
 * @param eph [in] satellite ephemeris data
 * @param t_tx [in] transmission time, seconds of week
 */
void nav_gnss_sat_position(NAV_PRECISION pos[3],
                            NAV_PRECISION vel[3],
                            NAV_PRECISION *clk,
                            const nav_gnss_ephemeris_t *eph,
                            NAV_PRECISION t_tx);

/**
 * @brief Compute azimuth and elevation of satellite from receiver position.
 * @param az [out] azimuth, rad [0, 2*PI]
 * @param el [out] elevation, rad [-PI/2, PI/2]
 * @param rx_pos [in] receiver ECEF position, m
 * @param sat_pos [in] satellite ECEF position, m
 */
void nav_gnss_azel(NAV_PRECISION *az, NAV_PRECISION *el,
                   const NAV_PRECISION rx_pos[3],
                   const NAV_PRECISION sat_pos[3]);

/**
 * @brief Klobuchar ionospheric delay model.
 *
 *  I_z = 5e-9 + A * cos(2*PI*(t - phi)/P)  (vertical delay, seconds)
 *  where A and P are computed from alpha/beta parameters.
 *
 * @param delay_m [out] ionospheric delay, m
 * @param lat_rad [in] receiver latitude, rad
 * @param lon_rad [in] receiver longitude, rad
 * @param az_rad [in] satellite azimuth, rad
 * @param el_rad [in] satellite elevation, rad
 * @param tow [in] GPS time of week, s
 * @param alpha [in] Klobuchar alpha parameters
 * @param beta [in] Klobuchar beta parameters
 */
void nav_iono_klobuchar(NAV_PRECISION *delay_m,
                         NAV_PRECISION lat_rad, NAV_PRECISION lon_rad,
                         NAV_PRECISION az_rad, NAV_PRECISION el_rad,
                         NAV_PRECISION tow,
                         const NAV_PRECISION alpha[4],
                         const NAV_PRECISION beta[4]);

/**
 * @brief Simple tropospheric delay model (Saastamoinen).
 * @param delay_m [out] tropospheric delay, m
 * @param el_rad [in] satellite elevation, rad
 * @param alt_m [in] receiver altitude, m
 */
void nav_tropo_saastamoinen(NAV_PRECISION *delay_m,
                             NAV_PRECISION el_rad, NAV_PRECISION alt_m);

/**
 * @brief Compute DOP values from geometry matrix G (n x 4).
 *  H = inv(G^T * G)
 *  GDOP = sqrt(trace(H))
 *  PDOP = sqrt(H_11 + H_22 + H_33)
 *  HDOP = sqrt(H_11 + H_22)
 *  VDOP = sqrt(H_33)
 *  TDOP = sqrt(H_44)
 */
void nav_gnss_compute_dop(NAV_PRECISION *gdop, NAV_PRECISION *pdop,
                           NAV_PRECISION *hdop, NAV_PRECISION *vdop,
                           NAV_PRECISION *tdop,
                           const NAV_PRECISION *G, int n_svs);

/**
 * @brief Correct pseudorange for satellite clock, iono, and tropo errors.
 */
NAV_PRECISION nav_gnss_correct_pseudorange(NAV_PRECISION raw_pr,
                                            NAV_PRECISION sat_clk,
                                            NAV_PRECISION iono,
                                            NAV_PRECISION tropo);

/**
 * @brief Convert ECEF position to geodetic (WGS84).
 *  Iterative Bowring or Zhu method.
 */
void nav_ecef_to_geodetic(nav_geodetic_t *geo, const NAV_PRECISION ecef[3]);

/**
 * @brief Convert geodetic to ECEF position.
 */
void nav_geodetic_to_ecef(NAV_PRECISION ecef[3], const nav_geodetic_t *geo);

/**
 * @brief Compute ENU-to-ECEF rotation matrix at given reference point.
 */
void nav_enu_to_ecef_matrix(NAV_PRECISION R[9],
                             NAV_PRECISION lat_rad, NAV_PRECISION lon_rad);

#endif /* NAV_GNSS_H */
