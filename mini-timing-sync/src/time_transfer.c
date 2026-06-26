/**
 * @file time_transfer.c
 * @brief Time transfer methods implementation
 *
 * Implements: Two-way time transfer with asymmetry,
 * common-view and all-in-view GPS time transfer,
 * TWSTFT delay modeling, Sagnac correction,
 * White Rabbit link model and phase tracking.
 */

#include "time_transfer.h"
#include <math.h>
#include <string.h>

/* Physical constants */
#define EARTH_ROTATION_RATE  7.2921159e-5  /* rad/s */
#define SPEED_OF_LIGHT       299792458.0    /* m/s */
#define EARTH_RADIUS_M       6371000.0      /* Earth mean radius [m] */

/* ===================================================================
 * L4: Two-Way Time Transfer with Asymmetric Delay
 * =================================================================== */

int twtt_with_asymmetry(const PtpTimestamps *ts, double asymmetry_ns,
                        double *offset, double *delay)
{
    if (!ts || !offset || !delay) return -1;
    if (!timing_timestamp_valid(&ts->t1) ||
        !timing_timestamp_valid(&ts->t2) ||
        !timing_timestamp_valid(&ts->t3) ||
        !timing_timestamp_valid(&ts->t4)) {
        return -1;
    }

    double t2_m_t1 = timing_timestamp_diff_ns(&ts->t2, &ts->t1);
    double t4_m_t3 = timing_timestamp_diff_ns(&ts->t4, &ts->t3);

    /* With asymmetry compensation:
     *   offset = ((t2 - t1) - (t4 - t3)) / 2 + asymmetry
     *   delay  = ((t2 - t1) + (t4 - t3)) / 2
     *
     * Where asymmetry = (forward_delay - reverse_delay) / 2
     *   forward_delay  = (t2 - t1) - offset_master
     *   reverse_delay  = (t4 - t3) + offset_master
     */
    *delay  = (t2_m_t1 + t4_m_t3) / 2.0;
    *offset = (t2_m_t1 - t4_m_t3) / 2.0 + asymmetry_ns;

    return 0;
}

/* ===================================================================
 * Helper: Geodetic to ECEF coordinate conversion
 * =================================================================== */

static void geodetic_to_ecef(double lat_deg, double lon_deg, double alt_m,
                             double *x, double *y, double *z)
{
    /* WGS-84 ellipsoid parameters */
    double a = 6378137.0;          /* Semi-major axis [m] */
    double f = 1.0 / 298.257223563; /* Flattening */
    double e2 = 2.0 * f - f * f;    /* First eccentricity squared */

    double lat_rad = lat_deg * M_PI / 180.0;
    double lon_rad = lon_deg * M_PI / 180.0;

    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    /* Radius of curvature in prime vertical */
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);

    *x = (N + alt_m) * cos_lat * cos(lon_rad);
    *y = (N + alt_m) * cos_lat * sin(lon_rad);
    *z = (N * (1.0 - e2) + alt_m) * sin_lat;
}

/* ===================================================================
 * L5: Common-View GPS Time Transfer
 * =================================================================== */

TimeTransferResult common_view_transfer(const GpsCommonViewObs *obs_A,
                                        const GpsCommonViewObs *obs_B,
                                        double sta_A_pos_x, double sta_A_pos_y,
                                        double sta_A_pos_z,
                                        double sta_B_pos_x, double sta_B_pos_y,
                                        double sta_B_pos_z)
{
    TimeTransferResult result;
    memset(&result, 0, sizeof(result));

    if (!obs_A || !obs_B) return result;
    (void)sta_A_pos_x; (void)sta_A_pos_y; (void)sta_A_pos_z;
    (void)sta_B_pos_x; (void)sta_B_pos_y; (void)sta_B_pos_z;

    /* Common-view: same satellite observed by both stations.
     *
     * offset_AB = (obs_A - obs_B) + corrections
     *
     * The key insight: satellite clock error cancels out
     * when both stations observe the same satellite at
     * the same GPS time.
     */

    /* Satellite clock cancels, leaving:
     *   offset_AB = (local_A - GPS_time) - (local_B - GPS_time)
     *             = local_A - local_B
     *
     * Plus differential path corrections.
     */

    /* Differential ionospheric and tropospheric delays */
    double diff_iono = obs_A->iono_delay_ns - obs_B->iono_delay_ns;
    double diff_tropo = obs_A->tropo_delay_ns - obs_B->tropo_delay_ns;

    /* Simple common-view: offset is difference of clock offsets
     * corrected for differential atmospheric delays */
    result.offset_ns = (obs_A->local_clock_offset_ns
                       - obs_B->local_clock_offset_ns)
                       + diff_iono + diff_tropo;

    /* Uncertainty estimate (simplified) */
    result.uncertainty_ns = 5.0; /* ~5 ns for single-satellite common-view */

    result.method = TT_COMMON_VIEW;
    result.valid = 1;
    result.measurement_time = obs_A->gps_time;

    return result;
}

/* ===================================================================
 * L5: All-in-View GPS Time Transfer
 * =================================================================== */

TimeTransferResult all_in_view_transfer(const GpsCommonViewObs *obs_A,
                                        int num_obs_A,
                                        const GpsCommonViewObs *obs_B,
                                        int num_obs_B, int method)
{
    TimeTransferResult result;
    memset(&result, 0, sizeof(result));

    if (!obs_A || !obs_B || num_obs_A < 1 || num_obs_B < 1) return result;

    /* All-in-view: stations observe different satellites.
     * Requires precise satellite clock/orbit products.
     *
     * For this simplified implementation:
     *   offset_AB = mean(obs_A) - mean(obs_B)
     *
     * Method 0: simple mean
     * Method 1: weighted by elevation (higher elevation = lower noise)
     */

    double sum_A = 0.0, sum_B = 0.0, weight_sum_A = 0.0, weight_sum_B = 0.0;
    int count_A = 0, count_B = 0;

    for (int i = 0; i < num_obs_A; i++) {
        double weight;

        if (method == 1) {
            /* Weight by sin(elevation), rejecting low-elevation (< 10 deg) */
            double elev_rad = obs_A[i].elevation_deg * M_PI / 180.0;
            if (elev_rad < 10.0 * M_PI / 180.0) continue; /* Skip low elevations */
            weight = sin(elev_rad);
        } else {
            weight = 1.0;
        }

        sum_A += obs_A[i].local_clock_offset_ns * weight;
        weight_sum_A += weight;
        count_A++;
    }

    for (int i = 0; i < num_obs_B; i++) {
        double weight;

        if (method == 1) {
            double elev_rad = obs_B[i].elevation_deg * M_PI / 180.0;
            if (elev_rad < 10.0 * M_PI / 180.0) continue;
            weight = sin(elev_rad);
        } else {
            weight = 1.0;
        }

        sum_B += obs_B[i].local_clock_offset_ns * weight;
        weight_sum_B += weight;
        count_B++;
    }

    if (count_A == 0 || count_B == 0 || weight_sum_A == 0.0 || weight_sum_B == 0.0) {
        return result;
    }

    double mean_A = sum_A / weight_sum_A;
    double mean_B = sum_B / weight_sum_B;

    result.offset_ns = mean_A - mean_B;
    result.uncertainty_ns = 3.0; /* ~3 ns typical for all-in-view with IGS products */
    result.method = TT_ALL_IN_VIEW;
    result.valid = 1;

    return result;
}

/* ===================================================================
 * L5: Time Transfer Uncertainty Budget
 * =================================================================== */

void time_transfer_uncertainty(const TimeTransferLink *link,
                               const AsymmetryModel *asymmetry,
                               double temp_c,
                               TimeTransferResult *result)
{
    if (!link || !asymmetry || !result) return;

    /* Uncertainty budget components (all in ns):
     *
     * 1. Timestamp resolution: typically 0.1-10 ns
     * 2. Asymmetry uncertainty: from calibration and dynamics
     * 3. Temperature sensitivity: temp_coeff * delta_T
     * 4. Multipath: 0.5-5 ns depending on environment
     * 5. Instrumentation noise: 0.1-2 ns
     */

    double U_ts_resolution = 1.0; /* 1 ns timestamp resolution */

    double U_asymmetry = asymmetry->asymmetry_uncertainty_ns
                       + asymmetry->dynamic_asymmetry_ns;

    /* Temperature effect on delay (referenced to 25 C) */
    double delta_T = temp_c - 25.0;
    double U_temp = fabs(link->temperature_coeff_ps_per_c)
                  * fabs(delta_T) * 1.0e-3; /* ps -> ns */

    double U_multipath = 2.0; /* 2 ns typical multipath */

    double U_instr_noise = 0.5; /* 0.5 ns instrumentation */

    /* Combined uncertainty (RSS) */
    result->uncertainty_ns = sqrt(
        U_ts_resolution * U_ts_resolution +
        U_asymmetry * U_asymmetry +
        U_temp * U_temp +
        U_multipath * U_multipath +
        U_instr_noise * U_instr_noise
    );
}

/* ===================================================================
 * L6: TWSTFT Delay Model
 * =================================================================== */

int twstft_delay_model(double sat_longitude_deg,
                       double sta_A_lat_deg, double sta_A_lon_deg,
                       double sta_A_alt_m,
                       double sta_B_lat_deg, double sta_B_lon_deg,
                       double sta_B_alt_m,
                       double *sagnac_correction, double *total_delay)
{
    if (!sagnac_correction || !total_delay) return -1;

    /* Geostationary satellite at altitude ~35786000 m */
    double sat_alt_m = 35786000.0;

    /* Convert station coordinates to ECEF */
    double Ax, Ay, Az, Bx, By, Bz;
    geodetic_to_ecef(sta_A_lat_deg, sta_A_lon_deg, sta_A_alt_m, &Ax, &Ay, &Az);
    geodetic_to_ecef(sta_B_lat_deg, sta_B_lon_deg, sta_B_alt_m, &Bx, &By, &Bz);

    /* Satellite position in ECEF */
    double sat_lon_rad = sat_longitude_deg * M_PI / 180.0;
    double sat_r = EARTH_RADIUS_M + sat_alt_m;
    double Sx = sat_r * cos(sat_lon_rad);
    double Sy = sat_r * sin(sat_lon_rad);
    double Sz = 0.0;

    /* Compute slant ranges */
    double dx_A = Sx - Ax, dy_A = Sy - Ay, dz_A = Sz - Az;
    double range_A = sqrt(dx_A*dx_A + dy_A*dy_A + dz_A*dz_A);

    double dx_B = Sx - Bx, dy_B = Sy - By, dz_B = Sz - Bz;
    double range_B = sqrt(dx_B*dx_B + dy_B*dy_B + dz_B*dz_B);

    /* Total path delay (A -> satellite -> B) */
    double path_length_m = range_A + range_B;
    *total_delay = path_length_m / SPEED_OF_LIGHT * 1.0e9; /* ns */

    /* Sagnac correction */
    *sagnac_correction = sagnac_correction_ns(sta_A_lat_deg, sta_A_lon_deg,
                                              sta_B_lat_deg, sta_B_lon_deg);

    return 0;
}

/* ===================================================================
 * L6: Sagnac Effect Correction
 * =================================================================== */

double sagnac_correction_ns(double sta_A_lat_deg, double sta_A_lon_deg,
                            double sta_B_lat_deg, double sta_B_lon_deg)
{
    /* Sagnac effect for Earth rotation:
     *
     *   Delta_t = (2 * omega_E * A_E) / c^2
     *
     * Where A_E is the projected area on the equatorial plane
     * of the polygon formed by the signal path and the Earth's axis.
     *
     * For two stations on the surface:
     *   Delta_t = (omega_E * R^2 * cos(lat_A) * cos(lat_B)
     *             * sin(lon_B - lon_A)) / c^2
     */

    double lat_A_rad = sta_A_lat_deg * M_PI / 180.0;
    double lat_B_rad = sta_B_lat_deg * M_PI / 180.0;
    double lon_A_rad = sta_A_lon_deg * M_PI / 180.0;
    double lon_B_rad = sta_B_lon_deg * M_PI / 180.0;

    double area = EARTH_RADIUS_M * EARTH_RADIUS_M
                * cos(lat_A_rad) * cos(lat_B_rad)
                * sin(lon_B_rad - lon_A_rad) / 2.0;

    double sagnac_s = 2.0 * EARTH_ROTATION_RATE * area
                    / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    /* Convert to nanoseconds */
    return sagnac_s * 1.0e9;
}

/* ===================================================================
 * L8: White Rabbit Link Model
 * =================================================================== */

void white_rabbit_link_model(double fiber_length_m,
                             double fiber_refractive_index,
                             double tx_delay_ns, double rx_delay_ns,
                             double *forward_delay,
                             double *reverse_delay)
{
    if (!forward_delay || !reverse_delay) return;

    /* Fiber propagation delay:
     *   t_fiber = L * n / c
     *
     * Where L = fiber length, n = refractive index (~1.467 for SMF-28),
     * c = speed of light in vacuum.
     */
    if (fiber_refractive_index <= 0.0) {
        fiber_refractive_index = 1.467; /* SMF-28 at 1550nm */
    }

    double fiber_delay_ns = fiber_length_m * fiber_refractive_index
                          / SPEED_OF_LIGHT * 1.0e9;

    /* Total one-way delay = fiber + TX electronics + RX electronics */
    *forward_delay = fiber_delay_ns + tx_delay_ns + rx_delay_ns;

    /* WR model: symmetry assumption
     *   reverse_delay = forward_delay (same fiber in both directions)
     *
     * In practice, bidirectional transmission on single fiber
     * ensures excellent symmetry. Remaining asymmetry:
     *   - Different wavelengths (1490/1310 nm -> chromatic dispersion)
     *   - TX/RX delay differences
     */
    *reverse_delay = *forward_delay;
}

/* ===================================================================
 * L8: White Rabbit Phase Tracking
 * =================================================================== */

double white_rabbit_phase_track(double main_phase_ps,
                                double helper_phase_ps,
                                double ref_freq_hz)
{
    /* White Rabbit uses DDMTD (Digital Dual Mixer Time Difference)
     * for sub-picosecond phase measurement.
     *
     * Architecture:
     * - Main PLL: locks recovered clock to local reference
     * - Helper PLL: offset by ~10 kHz for phase magnification
     * - Phase detector measures main vs helper
     *
     * Phase tracking for alignment:
     *   delta_phi = (main_phase - helper_phase)
     *
     * Convert phase difference to time:
     *   delta_t = delta_phi / (2*pi*f_ref)
     *
     * For 125 MHz reference: 1 radian = 1273 ps
     */
    if (ref_freq_hz <= 0.0) return 0.0;

    double delta_phase_ps = main_phase_ps - helper_phase_ps;

    /* DDMTD magnification: phase difference is amplified by
     * the frequency offset ratio.
     *   phase_magnification = f_ref / (f_helper - f_main)
     *
     * Typical: 125 MHz / 10 kHz = 12500x magnification
     */

    /* Phase wrapping to [-T/2, T/2] where T = 1/f_ref */
    double period_ps = 1.0e12 / ref_freq_hz;

    while (delta_phase_ps > period_ps / 2.0) {
        delta_phase_ps -= period_ps;
    }
    while (delta_phase_ps < -period_ps / 2.0) {
        delta_phase_ps += period_ps;
    }

    return delta_phase_ps;
}
