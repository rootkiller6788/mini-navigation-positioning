/* =========================================================================
 * gnss_common.c — Coordinate transforms, WGS84, time, vector/matrix ops
 *
 * Covers L1 (coordinates/Ellipsoid), L2 (ECEF↔LLA↔ENU),
 * L3 (3×3/4×4 matrix algebra), L4 (Sagnac effect).
 *
 * References:
 * - NIMA TR8350.2, WGS84 technical report
 * - Bowring, B.R. (1985). Survey Review 28(218), 202-206
 * - IS-GPS-200 §20.3.3.3.3.3 (Sagnac correction)
 * ========================================================================= */

#include "gnss_common.h"
#include <string.h>

/* -------------------------------------------------------------------------
 * L1: WGS84 ellipsoid singleton
 * ------------------------------------------------------------------------- */

gnss_ellipsoid_t gnss_wgs84_ellipsoid(void) {
    gnss_ellipsoid_t ell;
    ell.a   = GNSS_WGS84_A;
    ell.f   = GNSS_WGS84_F;
    ell.e2  = GNSS_WGS84_E2;
    ell.ep2 = GNSS_WGS84_EP2;
    ell.b   = GNSS_WGS84_B;
    return ell;
}

/* -------------------------------------------------------------------------
 * L2: Prime vertical and meridian radii of curvature
 * ------------------------------------------------------------------------- */

double gnss_rad_curvature_prime(double lat_rad) {
    double s = sin(lat_rad);
    return GNSS_WGS84_A / sqrt(1.0 - GNSS_WGS84_E2 * s * s);
}

double gnss_rad_curvature_prime_ell(double lat_rad, const gnss_ellipsoid_t *ell) {
    double s = sin(lat_rad);
    return ell->a / sqrt(1.0 - ell->e2 * s * s);
}

double gnss_rad_curvature_meridian(double lat_rad) {
    double s = sin(lat_rad);
    double d = 1.0 - GNSS_WGS84_E2 * s * s;
    return GNSS_WGS84_A * (1.0 - GNSS_WGS84_E2) / (d * sqrt(d));
}

double gnss_rad_curvature_meridian_ell(double lat_rad, const gnss_ellipsoid_t *ell) {
    double s = sin(lat_rad);
    double d = 1.0 - ell->e2 * s * s;
    return ell->a * (1.0 - ell->e2) / (d * sqrt(d));
}

/* -------------------------------------------------------------------------
 * L2: LLA → ECEF (closed-form)
 *
 * N(φ) = a / √(1 - e²·sin²φ)
 * x = (N + h)·cos φ·cos λ
 * y = (N + h)·cos φ·sin λ
 * z = (N·(1-e²) + h)·sin φ
 * ------------------------------------------------------------------------- */

gnss_ecef_t gnss_lla_to_ecef(gnss_lla_t lla) {
    gnss_ellipsoid_t ell = gnss_wgs84_ellipsoid();
    return gnss_lla_to_ecef_ell(lla, &ell);
}

gnss_ecef_t gnss_lla_to_ecef_ell(gnss_lla_t lla, const gnss_ellipsoid_t *ell) {
    gnss_ecef_t ecef;
    double s_lat = sin(lla.lat);
    double c_lat = cos(lla.lat);
    double s_lon = sin(lla.lon);
    double c_lon = cos(lla.lon);
    double N = ell->a / sqrt(1.0 - ell->e2 * s_lat * s_lat);
    ecef.x = (N + lla.alt) * c_lat * c_lon;
    ecef.y = (N + lla.alt) * c_lat * s_lon;
    ecef.z = (N * (1.0 - ell->e2) + lla.alt) * s_lat;
    return ecef;
}

/* -------------------------------------------------------------------------
 * L2: ECEF → LLA (Bowring 1985 iterative method)
 *
 * Uses Bowring's approximation for rapid convergence:
 *   p = √(x² + y²)
 *   θ = atan(z·a / (p·b))
 *   φ = atan((z + ep²·b·sin³θ) / (p - e²·a·cos³θ))
 *
 * Iteratively refines φ until change < 1e-12 rad (typically 2-3 iterations).
 * ------------------------------------------------------------------------- */

gnss_lla_t gnss_ecef_to_lla(gnss_ecef_t ecef) {
    gnss_ellipsoid_t ell = gnss_wgs84_ellipsoid();
    return gnss_ecef_to_lla_ell(ecef, &ell);
}

gnss_lla_t gnss_ecef_to_lla_ell(gnss_ecef_t ecef, const gnss_ellipsoid_t *ell) {
    gnss_lla_t lla;
    double x = ecef.x, y = ecef.y, z = ecef.z;
    double p = sqrt(x*x + y*y);

    /* Handle polar special case (p → 0) */
    if (p < 1e-12) {
        lla.lon = 0.0;
        lla.lat = (z >= 0.0) ? (M_PI / 2.0) : (-M_PI / 2.0);
        lla.alt = fabs(z) - ell->b;
        return lla;
    }

    /* Longitude (exact) */
    lla.lon = atan2(y, x);

    /* Iterative latitude refinement:
     * From: z = (N(1-e²)+h)·sin(φ), p = (N+h)·cos(φ)
     * ⇒ tan(φ) = (z + e²·N·sin(φ)) / p
     * Iteration: φ_{k+1} = atan2(z + e²·N_k·sin(φ_k), p)
     * where N_k = a / sqrt(1 - e²·sin²(φ_k))
     * Initial guess: φ₀ = atan2(z, p·(1-e²)) */
    double lat = atan2(z, p * (1.0 - ell->e2));
    double prev_lat;
    int iter;
    for (iter = 0; iter < 10; iter++) {
        prev_lat = lat;
        double s = sin(lat);
        double N = ell->a / sqrt(1.0 - ell->e2 * s * s);
        lat = atan2(z + ell->e2 * N * s, p);
        if (fabs(lat - prev_lat) < 1e-14) break;
    }
    lla.lat = lat;

    /* Height from final latitude */
    double s = sin(lat);
    double N = ell->a / sqrt(1.0 - ell->e2 * s * s);
    if (fabs(lat) < 1.5533) {
        lla.alt = p / cos(lat) - N;
    } else {
        lla.alt = z / s - N * (1.0 - ell->e2);
        if (fabs(s) < 1e-12) lla.alt = fabs(z) - ell->b;
    }

    return lla;
}

/* -------------------------------------------------------------------------
 * L2: ECEF ↔ ENU at reference LLA
 *
 * Rotation from ECEF to ENU:
 *   [e]   [   -sin(λ)         cos(λ)         0  ] [Δx]
 *   [n] = [ -sin(φ)cos(λ)  -sin(φ)sin(λ)   cos(φ) ] [Δy]
 *   [u]   [  cos(φ)cos(λ)   cos(φ)sin(λ)   sin(φ) ] [Δz]
 *
 * This rotation is orthonormal: R⁻¹ = Rᵀ.
 * ------------------------------------------------------------------------- */

gnss_enu_t gnss_ecef_to_enu(gnss_ecef_t ecef, gnss_lla_t ref) {
    gnss_ecef_t ref_ecef = gnss_lla_to_ecef(ref);
    double dx = ecef.x - ref_ecef.x;
    double dy = ecef.y - ref_ecef.y;
    double dz = ecef.z - ref_ecef.z;

    double s_lat = sin(ref.lat), c_lat = cos(ref.lat);
    double s_lon = sin(ref.lon), c_lon = cos(ref.lon);

    gnss_enu_t enu;
    enu.e = -s_lon       * dx + c_lon       * dy;
    enu.n = -s_lat*c_lon * dx - s_lat*s_lon * dy + c_lat * dz;
    enu.u =  c_lat*c_lon * dx + c_lat*s_lon * dy + s_lat * dz;
    return enu;
}

gnss_ecef_t gnss_enu_to_ecef(gnss_enu_t enu, gnss_lla_t ref) {
    gnss_ecef_t ref_ecef = gnss_lla_to_ecef(ref);
    double s_lat = sin(ref.lat), c_lat = cos(ref.lat);
    double s_lon = sin(ref.lon), c_lon = cos(ref.lon);

    gnss_ecef_t ecef;
    ecef.x = -s_lon*enu.e - s_lat*c_lon*enu.n + c_lat*c_lon*enu.u + ref_ecef.x;
    ecef.y =  c_lon*enu.e - s_lat*s_lon*enu.n + c_lat*s_lon*enu.u + ref_ecef.y;
    ecef.z =                c_lat       *enu.n + s_lat       *enu.u + ref_ecef.z;
    return ecef;
}

/* -------------------------------------------------------------------------
 * L4: Sagnac effect (Earth rotation correction)
 *
 * During signal travel time τ = ||sat - rx|| / c, the Earth rotates by
 * ω_e·τ radians. The satellite position must be corrected for this rotation.
 *
 * Correction (in meters):
 *   Δρ_Sagnac = (ω_e / c) · (x_sat · y_rx - y_sat · x_rx)
 *
 * Satellite ECEF position rotated backward by ω_e·τ:
 *   x_s'  =  x_s·cos(θ) + y_s·sin(θ)
 *   y_s'  = -x_s·sin(θ) + y_s·cos(θ)
 *   z_s'  =  z_s
 * where θ = ω_e · τ = ω_e · ||sat - rx|| / c
 * ------------------------------------------------------------------------- */

double gnss_sagnac_correction(gnss_ecef_t sat_pos, gnss_ecef_t rx_pos) {
    return (GNSS_OMEGA_E / GNSS_C_LIGHT) *
           (sat_pos.x * rx_pos.y - sat_pos.y * rx_pos.x);
}

double gnss_sagnac_rotate(gnss_ecef_t sat_pos, double travel_time) {
    double theta = GNSS_OMEGA_E * travel_time;
    double c_th = cos(theta);
    double s_th = sin(theta);
    /* Returns rotated X coordinate; caller must also adjust Y */
    return sat_pos.x * c_th + sat_pos.y * s_th;
    /* Y' = -sat_pos.x * s_th + sat_pos.y * c_th */
}

/* -------------------------------------------------------------------------
 * L1: Time conversions
 * ------------------------------------------------------------------------- */

gnss_utctime_t gnss_gpstime_to_utc(gnss_gpstime_t gps, int16_t leap_sec) {
    gnss_utctime_t utc;
    double jd = GNSS_GPS_EPOCH_JD + gps.week * 7.0 + gps.sow / 86400.0;
    jd -= leap_sec / 86400.0;

    /* Convert JD to calendar date (Fliegel & Van Flandern algorithm) */
    double jd_ut = jd + 0.5;
    int32_t Z = (int32_t)jd_ut;
    double F = jd_ut - Z;
    int32_t A = Z;
    if (Z >= 2299161) {
        int32_t alpha = (int32_t)((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - alpha / 4;
    }
    int32_t B = A + 1524;
    int32_t C = (int32_t)((B - 122.1) / 365.25);
    int32_t D = (int32_t)(365.25 * C);
    int32_t E = (int32_t)((B - D) / 30.6001);

    utc.day   = B - D - (int32_t)(30.6001 * E);
    utc.month = (E < 14) ? (E - 1) : (E - 13);
    utc.year  = (utc.month > 2) ? (C - 4716) : (C - 4715);

    double day_frac = F + ((double)(B - D - (int32_t)(30.6001 * E))) - utc.day;
    if (day_frac >= 1.0) { day_frac -= 1.0; utc.day++; }
    if (day_frac < 0.0)  { day_frac += 1.0; utc.day--; }

    double total_seconds = day_frac * 86400.0;
    utc.hour   = (int16_t)(total_seconds / 3600.0);
    utc.minute = (int16_t)((total_seconds - utc.hour * 3600.0) / 60.0);
    utc.second = total_seconds - utc.hour * 3600.0 - utc.minute * 60.0;
    utc.leap_seconds = leap_sec;

    return utc;
}

gnss_gpstime_t gnss_utc_to_gpstime(gnss_utctime_t utc) {
    gnss_gpstime_t gps;
    /* Compute Julian date */
    int32_t y = utc.year, m = utc.month;
    if (m <= 2) { y--; m += 12; }
    int32_t A = y / 100;
    int32_t B = 2 - A + A / 4;
    double jd = (int32_t)(365.25 * (y + 4716))
              + (int32_t)(30.6001 * (m + 1))
              + utc.day + B - 1524.5;
    jd += utc.hour / 24.0 + utc.minute / 1440.0
       + utc.second / 86400.0;
    jd += utc.leap_seconds / 86400.0;

    double total_days = jd - GNSS_GPS_EPOCH_JD;
    gps.week = (int32_t)(total_days / 7.0);
    gps.sow  = (total_days - gps.week * 7.0) * 86400.0;
    if (gps.sow < 0.0) { gps.week--; gps.sow += 604800.0; }
    if (gps.sow >= 604800.0) { gps.week++; gps.sow -= 604800.0; }
    return gps;
}

double gnss_utc_to_jd(gnss_utctime_t utc) {
    int32_t y = utc.year, m = utc.month;
    if (m <= 2) { y--; m += 12; }
    int32_t A = y / 100;
    int32_t B = 2 - A + A / 4;
    double jd = (int32_t)(365.25 * (y + 4716))
              + (int32_t)(30.6001 * (m + 1))
              + utc.day + B - 1524.5;
    jd += utc.hour / 24.0 + utc.minute / 1440.0
       + utc.second / 86400.0;
    return jd;
}

/* -------------------------------------------------------------------------
 * L1: Geometric range and satellite angles
 * ------------------------------------------------------------------------- */

double gnss_geometric_range(gnss_ecef_t sat, gnss_ecef_t rx) {
    double dx = sat.x - rx.x;
    double dy = sat.y - rx.y;
    double dz = sat.z - rx.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double gnss_sat_elevation(gnss_ecef_t sat, gnss_ecef_t rx) {
    gnss_lla_t rx_lla = gnss_ecef_to_lla(rx);
    gnss_enu_t enu = gnss_ecef_to_enu(sat, rx_lla);
    double h_len = sqrt(enu.e*enu.e + enu.n*enu.n);
    return atan2(enu.u, h_len); /* elevation = angle above horizon */
}

double gnss_sat_azimuth(gnss_ecef_t sat, gnss_ecef_t rx) {
    gnss_lla_t rx_lla = gnss_ecef_to_lla(rx);
    gnss_enu_t enu = gnss_ecef_to_enu(sat, rx_lla);
    return atan2(enu.e, enu.n); /* azimuth from North, clockwise = positive */
}

/* =========================================================================
 * L3: 3-vector operations
 * ========================================================================= */

void gnss_vec3_sub(const double a[3], const double b[3], double result[3]) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

void gnss_vec3_add(const double a[3], const double b[3], double result[3]) {
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}

void gnss_vec3_scale(const double v[3], double s, double result[3]) {
    result[0] = v[0] * s;
    result[1] = v[1] * s;
    result[2] = v[2] * s;
}

double gnss_vec3_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void gnss_vec3_cross(const double a[3], const double b[3], double result[3]) {
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

double gnss_vec3_norm(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void gnss_vec3_normalize(const double v[3], double result[3]) {
    double n = gnss_vec3_norm(v);
    if (n > 1e-30) {
        result[0] = v[0] / n;
        result[1] = v[1] / n;
        result[2] = v[2] / n;
    } else {
        result[0] = result[1] = result[2] = 0.0;
    }
}

/* =========================================================================
 * L3: 3×3 matrix operations
 * ========================================================================= */

void gnss_mat33_identity(gnss_mat33_t *m) {
    m->m[0][0] = 1.0; m->m[0][1] = 0.0; m->m[0][2] = 0.0;
    m->m[1][0] = 0.0; m->m[1][1] = 1.0; m->m[1][2] = 0.0;
    m->m[2][0] = 0.0; m->m[2][1] = 0.0; m->m[2][2] = 1.0;
}

void gnss_mat33_multiply(const gnss_mat33_t *a, const gnss_mat33_t *b,
                          gnss_mat33_t *c) {
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            double sum = 0.0;
            for (k = 0; k < 3; k++) {
                sum += a->m[i][k] * b->m[k][j];
            }
            c->m[i][j] = sum;
        }
    }
}

void gnss_mat33_transpose(const gnss_mat33_t *a, gnss_mat33_t *at) {
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            at->m[j][i] = a->m[i][j];
}

void gnss_mat33_vec_multiply(const gnss_mat33_t *a, const double v[3],
                              double result[3]) {
    result[0] = a->m[0][0]*v[0] + a->m[0][1]*v[1] + a->m[0][2]*v[2];
    result[1] = a->m[1][0]*v[0] + a->m[1][1]*v[1] + a->m[1][2]*v[2];
    result[2] = a->m[2][0]*v[0] + a->m[2][1]*v[1] + a->m[2][2]*v[2];
}

int gnss_mat33_inverse(const gnss_mat33_t *a, gnss_mat33_t *ainv) {
    /* Cofactor expansion (3×3 is small enough for explicit formula) */
    const double (*m)[3] = a->m;
    double det = m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
               - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
               + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

    if (fabs(det) < 1e-30) return -1;

    double inv_det = 1.0 / det;
    ainv->m[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * inv_det;
    ainv->m[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) * inv_det;
    ainv->m[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * inv_det;
    ainv->m[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) * inv_det;
    ainv->m[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * inv_det;
    ainv->m[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) * inv_det;
    ainv->m[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * inv_det;
    ainv->m[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) * inv_det;
    ainv->m[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * inv_det;
    return 0;
}

/* =========================================================================
 * L3: 4×4 matrix operations
 * ========================================================================= */

void gnss_mat44_identity(gnss_mat44_t *m) {
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            m->m[i][j] = (i == j) ? 1.0 : 0.0;
}

void gnss_mat44_multiply(const gnss_mat44_t *a, const gnss_mat44_t *b,
                          gnss_mat44_t *c) {
    int i, j, k;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            double sum = 0.0;
            for (k = 0; k < 4; k++)
                sum += a->m[i][k] * b->m[k][j];
            c->m[i][j] = sum;
        }
    }
}

void gnss_mat44_transpose(const gnss_mat44_t *a, gnss_mat44_t *at) {
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            at->m[j][i] = a->m[i][j];
}

void gnss_mat44_vec_multiply(const gnss_mat44_t *a, const double v[4],
                              double result[4]) {
    int i, j;
    for (i = 0; i < 4; i++) {
        result[i] = 0.0;
        for (j = 0; j < 4; j++)
            result[i] += a->m[i][j] * v[j];
    }
}

int gnss_mat44_inverse(gnss_mat44_t *a, gnss_mat44_t *inv) {
    /* Gauss-Jordan elimination with partial (row) pivoting on augmented [A|I] */
    int i, j, k;
    double aug[4][8]; /* rows: 4, cols: A(4) + I(4) = 8 */

    /* Initialize augmented matrix [A | I] */
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            aug[i][j] = a->m[i][j];
            aug[i][j+4] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Forward elimination with partial pivoting */
    for (k = 0; k < 4; k++) {
        /* Find pivot row */
        int pivot = k;
        double maxv = fabs(aug[k][k]);
        for (i = k+1; i < 4; i++) {
            if (fabs(aug[i][k]) > maxv) {
                maxv = fabs(aug[i][k]);
                pivot = i;
            }
        }
        if (maxv < 1e-30) return -1; /* singular */

        /* Swap rows */
        if (pivot != k) {
            for (j = 0; j < 8; j++) {
                double tmp = aug[k][j];
                aug[k][j] = aug[pivot][j];
                aug[pivot][j] = tmp;
            }
        }

        /* Scale pivot row */
        double piv_val = aug[k][k];
        for (j = 0; j < 8; j++) aug[k][j] /= piv_val;

        /* Eliminate other rows */
        for (i = 0; i < 4; i++) {
            if (i == k) continue;
            double factor = aug[i][k];
            for (j = 0; j < 8; j++)
                aug[i][j] -= factor * aug[k][j];
        }
    }

    /* Extract inverse from right half */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            inv->m[i][j] = aug[i][j+4];

    return 0;
}

int gnss_mat44_cholesky(const gnss_mat44_t *a, gnss_mat44_t *L) {
    /* Cholesky decomposition: A = L·Lᵀ, L is lower triangular
     * Only valid for symmetric positive-definite A (e.g., HᵀH matrix) */
    int i, j, k;
    memset(L->m, 0, sizeof(L->m));

    for (i = 0; i < 4; i++) {
        for (j = 0; j <= i; j++) {
            double sum = a->m[i][j];
            for (k = 0; k < j; k++)
                sum -= L->m[i][k] * L->m[j][k];

            if (i == j) {
                if (sum <= 1e-30) return -1; /* Not positive definite */
                L->m[i][j] = sqrt(sum);
            } else {
                L->m[i][j] = sum / L->m[j][j];
            }
        }
    }
    return 0;
}
