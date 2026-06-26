/**
 * geomag_model.c -- IGRF/WMM Geomagnetic Field Model Implementation
 *
 * L4: Spherical harmonic synthesis of geomagnetic potential
 * L5: Schmidt semi-normalized Legendre recurrence
 * L5: Field computation from Gauss coefficients
 * L5: Secular variation prediction, dipole approximation
 *
 * Reference: Alken et al., "IGRF-13", Earth Planets Space (2021)
 *            Langel, "The Main Field", Geomagnetism Vol.1 (1987)
 *            Lowes, "Mean-square values...", JGR (1966)
 */

#include "geomag_model.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * L5: Flat array index for Gauss coefficients
 * index(n,m) = n*(n+1)/2 + m   (n=1..Nmax, m=0..n)
 * Count = Nmax*(Nmax+3)/2
 * ======================================================================== */
int schmidt_index(int n, int m) {
    return n * (n + 1) / 2 + m;
}

/* Precompute sqrt((2m+1)/(2m)) for diagonal Legendre recurrence */
static void precompute_schmidt_diag_factors(int nmax, double *factors) {
    factors[0] = 1.0;
    for (int m = 1; m <= nmax; m++)
        factors[m] = sqrt((2.0 * m + 1.0) / (2.0 * m));
}

/* ========================================================================
 * L5: Compute Schmidt semi-normalized Legendre functions P_n^m(cos theta)
 * and dP_n^m/dtheta via forward column recurrence.
 *
 * Recurrence (Holmes & Featherstone, 2002):
 *   P_0^0 = 1
 *   P_m^m = diag_factor[m] * sin(theta) * P_{m-1}^{m-1}
 *   P_n^m = [(2n-1)*cos(theta)*P_{n-1}^m - (n+m-1)*P_{n-2}^m] / (n-m)
 *
 * dP_n^m/dtheta = [n*cos(theta)*P_n^m - (n+m)*P_{n-1}^m] / sin(theta)
 *
 * Schmidt normalization integral: int_0^pi [P_n^m]^2 sin(theta) dtheta = 1/(2n+1)
 *
 * Complexity: O(nmax^2) time, O(nmax^2) space.
 * ======================================================================== */
void compute_schmidt_legendre(int nmax, double colatitude, LegendreState *state) {
    if (nmax < 1 || !state || !state->values) return;

    double ct = cos(colatitude);
    double st = sin(colatitude);
    state->sin_theta = st;
    state->cos_theta = ct;
    state->nmax = nmax;

    int total = (nmax + 1) * (nmax + 2) / 2;
    memset(state->values, 0, total * sizeof(double));
    if (state->derivatives)
        memset(state->derivatives, 0, total * sizeof(double));

    double *diag_factor = (double *)malloc((nmax + 1) * sizeof(double));
    precompute_schmidt_diag_factors(nmax, diag_factor);

    /* P_0^0 = 1 */
    state->values[0] = 1.0;
    if (state->derivatives) state->derivatives[0] = 0.0;

    for (int m = 1; m <= nmax; m++) {
        int idx_mm = schmidt_index(m, m);
        int idx_prev = schmidt_index(m - 1, m - 1);
        state->values[idx_mm] = diag_factor[m] * st * state->values[idx_prev];

        if (state->derivatives) {
            if (st > 1e-15)
                state->derivatives[idx_mm] = m * ct / st * state->values[idx_mm];
            else
                state->derivatives[idx_mm] = 0.0;
        }

        for (int n = m + 1; n <= nmax; n++) {
            int idx_nm = schmidt_index(n, m);
            int idx_n1m = schmidt_index(n - 1, m);
            int idx_n2m = schmidt_index(n - 2, m);
            double a_coef = (2.0 * n - 1.0) / (n - m);
            double b_coef = (n + m - 1.0) / (n - m);
            state->values[idx_nm] = a_coef * ct * state->values[idx_n1m]
                                  - b_coef * state->values[idx_n2m];
            if (state->derivatives) {
                if (st > 1e-15)
                    state->derivatives[idx_nm] =
                        (n * ct * state->values[idx_nm]
                         - (n + m) * state->values[idx_n1m]) / st;
                else
                    state->derivatives[idx_nm] = 0.0;
            }
        }
    }
    free(diag_factor);
}

LegendreState *alloc_legendre_state(int nmax) {
    LegendreState *state = (LegendreState *)malloc(sizeof(LegendreState));
    if (!state) return NULL;
    int total = (nmax + 1) * (nmax + 2) / 2;
    state->nmax = nmax;
    state->values = (double *)malloc(total * sizeof(double));
    state->derivatives = (double *)malloc(total * sizeof(double));
    state->sin_theta = 0.0;
    state->cos_theta = 1.0;
    if (!state->values || !state->derivatives) {
        free(state->values); free(state->derivatives); free(state);
        return NULL;
    }
    return state;
}

void free_legendre_state(LegendreState *state) {
    if (state) { free(state->values); free(state->derivatives); free(state); }
}

/* Convert geodetic to geocentric spherical (r, theta=colatitude, phi) */
static void geodetic_to_geocentric_spherical(const GeodeticCoord *loc,
                                              double *r, double *theta, double *phi) {
    ECEFCoord ecef;
    geodetic_to_ecef(loc, &ecef);
    *r = sqrt(ecef.x * ecef.x + ecef.y * ecef.y + ecef.z * ecef.z);
    *theta = (*r > 1e-10) ? acos(ecef.z / *r) : 0.0;
    *phi = atan2(ecef.y, ecef.x);
}

/* Rotate field from geocentric spherical to geodetic NED */
static void rotate_geo_to_ned(double B_r, double B_theta, double B_phi,
                               double geodetic_lat, double geocentric_lat,
                               MagVector *B_ned) {
    double delta = (geodetic_lat - geocentric_lat) * DEG2RAD;
    double cd = cos(delta), sd = sin(delta);
    double B_south = -B_theta;
    B_ned->bx = -B_south * cd + B_r * sd;  /* North */
    B_ned->by = B_phi;                       /* East  */
    B_ned->bz = -B_south * sd - B_r * cd;  /* Down  */
}

/* ========================================================================
 * L4+L5: Full IGRF field computation
 *
 * B_r = sum_{n=1}^{nmax} (a/r)^{n+2} (n+1) sum_{m=0}^{n}
 *       [g_n^m cos(m phi) + h_n^m sin(m phi)] P_n^m(cos theta)
 * B_theta = - sum (a/r)^{n+2} sum [...] dP_n^m/dtheta
 * B_phi = -1/sin(theta) sum (a/r)^{n+2} sum m [-g sin(m phi) + h cos(m phi)] P_n^m
 *
 * Complexity: O(nmax^2) time, O(nmax^2) space.
 * Reference: Langel (1987), Eq. 24-27.
 * ======================================================================== */
int igrf_compute_field(const IGRFModel *model, const GeodeticCoord *loc,
                        MagVector *B) {
    if (!model || !loc || !B || !model->coeffs) return -1;
    int nmax = model->nmax;
    if (nmax < 1) return -1;

    double r, theta, phi;
    geodetic_to_geocentric_spherical(loc, &r, &theta, &phi);
    double sin_theta = sin(theta);

    LegendreState *leg = alloc_legendre_state(nmax);
    if (!leg) return -1;
    compute_schmidt_legendre(nmax, theta, leg);

    double geoc_lat = (M_PI / 2.0 - theta) * RAD2DEG;
    double B_r = 0.0, B_theta = 0.0, B_phi = 0.0;

    double *ar_pow = (double *)malloc((nmax + 1) * sizeof(double));
    if (!ar_pow) { free_legendre_state(leg); return -1; }
    double ar = GEOMAG_REF_RADIUS / r;
    ar_pow[1] = ar * ar * ar;
    for (int n = 2; n <= nmax; n++) ar_pow[n] = ar_pow[n - 1] * ar;

    for (int n = 1; n <= nmax; n++) {
        double ar_np2 = ar_pow[n];
        double sum_r = 0.0, sum_theta = 0.0, sum_phi = 0.0;
        for (int m = 0; m <= n; m++) {
            int idx = schmidt_index(n, m) - 1;
            GaussCoeff *coeff = &model->coeffs[idx];
            double P_nm = leg->values[idx];
            double dP_nm = leg->derivatives[idx];
            double cos_mphi = cos(m * phi), sin_mphi = sin(m * phi);
            double g_term = coeff->g_nm * cos_mphi + coeff->h_nm * sin_mphi;
            sum_r += g_term * P_nm;
            sum_theta += g_term * dP_nm;
            if (m > 0) {
                double h_term = -coeff->g_nm * sin_mphi + coeff->h_nm * cos_mphi;
                sum_phi += m * h_term * P_nm;
            }
        }
        B_r += (n + 1) * ar_np2 * sum_r;
        B_theta += ar_np2 * sum_theta;
        B_phi += ar_np2 * sum_phi;
    }
    B_theta = -B_theta;
    if (fabs(sin_theta) > 1e-15) B_phi = -B_phi / sin_theta;
    else B_phi = 0.0;

    rotate_geo_to_ned(B_r, B_theta, B_phi, loc->lat, geoc_lat, B);
    free(ar_pow);
    free_legendre_state(leg);
    return 0;
}
/* ========================================================================
 * L5: Secular variation computation
 * Same algorithm as igrf_compute_field but with dg/dh coefficients.
 * ======================================================================== */
int igrf_compute_secular_variation(const IGRFModel *model,
                                    const GeodeticCoord *loc, MagVector *dB_dt) {
    if (!model || !loc || !dB_dt || !model->coeffs) return -1;
    int nmax = model->nmax;

    double r, theta, phi;
    geodetic_to_geocentric_spherical(loc, &r, &theta, &phi);
    double sin_theta = sin(theta);

    LegendreState *leg = alloc_legendre_state(nmax);
    if (!leg) return -1;
    compute_schmidt_legendre(nmax, theta, leg);
    double geoc_lat = (M_PI / 2.0 - theta) * RAD2DEG;

    double B_r = 0.0, B_theta = 0.0, B_phi = 0.0;
    double *ar_pow = (double *)malloc((nmax + 1) * sizeof(double));
    if (!ar_pow) { free_legendre_state(leg); return -1; }
    double ar = GEOMAG_REF_RADIUS / r;
    ar_pow[1] = ar * ar * ar;
    for (int n = 2; n <= nmax; n++) ar_pow[n] = ar_pow[n - 1] * ar;

    for (int n = 1; n <= nmax; n++) {
        double ar_np2 = ar_pow[n];
        double sum_r = 0.0, sum_theta = 0.0, sum_phi = 0.0;
        for (int m = 0; m <= n; m++) {
            int idx = schmidt_index(n, m) - 1;
            GaussCoeff *coeff = &model->coeffs[idx];
            double P_nm = leg->values[idx];
            double dP_nm = leg->derivatives[idx];
            double cos_mphi = cos(m * phi), sin_mphi = sin(m * phi);
            double g_term = coeff->dg_nm * cos_mphi + coeff->dh_nm * sin_mphi;
            sum_r += g_term * P_nm;
            sum_theta += g_term * dP_nm;
            if (m > 0) {
                double h_term = -coeff->dg_nm * sin_mphi + coeff->dh_nm * cos_mphi;
                sum_phi += m * h_term * P_nm;
            }
        }
        B_r += (n + 1) * ar_np2 * sum_r;
        B_theta += ar_np2 * sum_theta;
        B_phi += ar_np2 * sum_phi;
    }
    B_theta = -B_theta;
    if (fabs(sin_theta) > 1e-15) B_phi = -B_phi / sin_theta;
    else B_phi = 0.0;
    rotate_geo_to_ned(B_r, B_theta, B_phi, loc->lat, geoc_lat, dB_dt);
    free(ar_pow);
    free_legendre_state(leg);
    return 0;
}

/* L5: Predict field at future year by linear SV extrapolation */
int igrf_predict_field(const IGRFModel *model, const GeodeticCoord *loc,
                        double year, MagVector *B) {
    if (!model || !loc || !B) return -1;
    double dt = year - model->epoch;
    if (fabs(dt) > 5.0) return -1;
    MagVector B_epoch;
    if (igrf_compute_field(model, loc, &B_epoch) != 0) return -1;
    MagVector dB;
    if (igrf_compute_secular_variation(model, loc, &dB) != 0) return -1;
    B->bx = B_epoch.bx + dB.bx * dt;
    B->by = B_epoch.by + dB.by * dt;
    B->bz = B_epoch.bz + dB.bz * dt;
    return 0;
}

/* ========================================================================
 * L4: Geocentric tilted dipole approximation
 *
 * B_r = 2*(a/r)^3 * [g10*cos(theta) + (g11*cos(phi)+h11*sin(phi))*sin(theta)]
 * B_theta = -(a/r)^3 * [g10*sin(theta) - (g11*cos(phi)+h11*sin(phi))*cos(theta)]
 * B_phi = (a/r)^3 * [g11*sin(phi) - h11*cos(phi)]
 *
 * Reference: Merrill et al., "The Magnetic Field of the Earth" (1996), Ch.1.
 * ======================================================================== */
int compute_dipole_field(double g10, double g11, double h11,
                          const GeodeticCoord *loc, MagVector *B) {
    if (!loc || !B) return -1;
    double r, theta, phi;
    geodetic_to_geocentric_spherical(loc, &r, &theta, &phi);
    double st = sin(theta), ct = cos(theta);
    double cp = cos(phi), sp = sin(phi);
    double ar3 = pow(GEOMAG_REF_RADIUS / r, 3.0);
    double g_term = g11 * cp + h11 * sp;

    double B_r = 2.0 * ar3 * (g10 * ct + g_term * st);
    double B_theta = ar3 * (g10 * st - g_term * ct);
    B_theta = -B_theta;
    double B_phi = -ar3 * (g11 * sp - h11 * cp);

    double geoc_lat = (M_PI / 2.0 - theta) * RAD2DEG;
    rotate_geo_to_ned(B_r, B_theta, B_phi, loc->lat, geoc_lat, B);
    return 0;
}

/* ========================================================================
 * L5: Find magnetic dip poles by iterative gradient-ascent on inclination.
 *
 * Starting from dipole pole positions, iteratively adjust lat/lon to
 * maximize |inclination| (toward +/- 90 degrees).
 *
 * Complexity: ~50 IGRF evaluations per pole.
 * ======================================================================== */
int compute_magnetic_poles(const IGRFModel *model,
                            GeodeticCoord *north_pole, GeodeticCoord *south_pole) {
    if (!model || !north_pole || !south_pole) return -1;

    double g10 = model->coeffs[schmidt_index(1, 0) - 1].g_nm;
    double g11 = model->coeffs[schmidt_index(1, 1) - 1].g_nm;
    double h11 = model->coeffs[schmidt_index(1, 1) - 1].h_nm;

    double dipole_lat_n = atan2(g10, sqrt(g11*g11 + h11*h11)) * RAD2DEG;
    double dipole_lon_n = atan2(h11, g11) * RAD2DEG;

    GeodeticCoord targets[2];
    targets[0].lat = dipole_lat_n;  targets[0].lon = dipole_lon_n; targets[0].alt = 0.0;
    targets[1].lat = -dipole_lat_n; targets[1].lon = dipole_lon_n + 180.0;
    if (targets[1].lon > 180.0) targets[1].lon -= 360.0;
    targets[1].alt = 0.0;

    GeodeticCoord *outs[2] = { north_pole, south_pole };

    for (int pole = 0; pole < 2; pole++) {
        double lat = targets[pole].lat, lon = targets[pole].lon;
        MagVector B;
        MagneticElements elem;
        double target_incl = (pole == 0) ? 90.0 : -90.0;

        for (int iter = 0; iter < 50; iter++) {
            GeodeticCoord trial = { lat, lon, 0.0 };
            if (igrf_compute_field(model, &trial, &B) != 0) break;
            compute_magnetic_elements(&B, &elem);
            if (fabs(elem.inclination - target_incl) < 0.01) break;

            double eps = 0.1;
            GeodeticCoord t2; MagVector B2; MagneticElements elem2;
            t2 = trial; t2.lat += eps;
            igrf_compute_field(model, &t2, &B2);
            compute_magnetic_elements(&B2, &elem2);
            double grad_lat = (elem2.inclination - elem.inclination) / eps;

            t2 = trial; t2.lon += eps;
            igrf_compute_field(model, &t2, &B2);
            compute_magnetic_elements(&B2, &elem2);
            double grad_lon = (elem2.inclination - elem.inclination) / eps;

            double step = (target_incl - elem.inclination) * 0.5;
            double gn = grad_lat*grad_lat + grad_lon*grad_lon;
            if (gn < 1e-10) break;
            lat += step * grad_lat / gn;
            lon += step * grad_lon / gn;
            if (lat > 90.0) lat = 90.0;
            if (lat < -90.0) lat = -90.0;
        }
        outs[pole]->lat = lat; outs[pole]->lon = lon; outs[pole]->alt = 0.0;
    }
    return 0;
}

/* L5: Earth dipole moment from n=1 Gauss coefficients */
double compute_dipole_moment(const IGRFModel *model) {
    double g10 = model->coeffs[schmidt_index(1, 0) - 1].g_nm;
    double g11 = model->coeffs[schmidt_index(1, 1) - 1].g_nm;
    double h11 = model->coeffs[schmidt_index(1, 1) - 1].h_nm;
    double a3 = GEOMAG_REF_RADIUS * GEOMAG_REF_RADIUS * GEOMAG_REF_RADIUS;
    double B0_T = sqrt(g10*g10 + g11*g11 + h11*h11) * 1e-9;
    return (4.0 * M_PI * a3 / (4.0 * M_PI * 1e-7)) * B0_T;
}

/* L5: McIlwain L-shell parameter (dipole approximation) */
double compute_l_shell(double g10, double g11, double h11,
                        const GeodeticCoord *loc) {
    double r, theta, phi;
    geodetic_to_geocentric_spherical(loc, &r, &theta, &phi);
    double dipole_colat = atan2(sqrt(g11*g11 + h11*h11), g10);
    double dipole_lon = atan2(h11, g11);
    double cos_alpha = cos(theta) * cos(dipole_colat)
                     + sin(theta) * sin(dipole_colat) * cos(phi - dipole_lon);
    double sin_alpha = sqrt(1.0 - cos_alpha * cos_alpha);
    if (sin_alpha < 1e-10) sin_alpha = 1e-10;
    return r / (GEOMAG_REF_RADIUS * sin_alpha * sin_alpha);
}

/* L5: Total magnetic field energy outside core (Lowes-Mauersberger spectrum) */
double compute_field_energy(const IGRFModel *model) {
    double a3 = GEOMAG_REF_RADIUS * GEOMAG_REF_RADIUS * GEOMAG_REF_RADIUS;
    double factor = a3 / 1e-7;
    double energy = 0.0;
    for (int n = 1; n <= model->nmax; n++) {
        double row_sum = 0.0;
        for (int m = 0; m <= n; m++) {
            int idx = schmidt_index(n, m);
            double g = model->coeffs[idx].g_nm;
            double h = model->coeffs[idx].h_nm;
            row_sum += (g * g + h * h) * 1e-18;
        }
        energy += factor * (n + 1) * row_sum;
    }
    return energy;
}

/* L6: South Atlantic Anomaly proximity indicator [0,1] */
double compute_saa_indicator(const IGRFModel *model, const GeodeticCoord *loc) {
    if (!model || !loc) return 0.0;
    MagVector B_full;
    MagneticElements elem_full;
    if (igrf_compute_field(model, loc, &B_full) != 0) return 0.0;
    compute_magnetic_elements(&B_full, &elem_full);

    double g10 = model->coeffs[schmidt_index(1, 0) - 1].g_nm;
    double g11 = model->coeffs[schmidt_index(1, 1) - 1].g_nm;
    double h11 = model->coeffs[schmidt_index(1, 1) - 1].h_nm;
    MagVector B_dipole;
    compute_dipole_field(g10, g11, h11, loc, &B_dipole);
    double F_dipole = mag_magnitude(&B_dipole);
    if (F_dipole < 1.0) return 0.0;

    double ratio = elem_full.total_intensity / F_dipole;
    double ind = 1.0 - ratio;
    if (ind < 0.0) ind = 0.0;
    if (ind > 1.0) ind = 1.0;
    return ind;
}

/* ========================================================================
 * IGRF-13 coefficients for epoch 2020.0 [nT] with SV [nT/yr]
 *
 * Source: IAGA Division V-MOD Working Group, IGRF-13, 2020.
 * Format: n, m, g_nm, h_nm, dg_nm, dh_nm
 * ======================================================================== */

#define IGRF13_NMAX 13

static const double igrf13_coeff_data[] = {
    1,0, -29404.5, 0.0, 5.7, 0.0,
    1,1, -1450.7, 4652.9, 7.4, -25.9,
    2,0, -2500.0, 0.0, -11.0, 0.0,
    2,1, 2982.0, -2991.6, -7.0, -30.2,
    2,2, 1676.8, -734.8, -1.9, -23.9,
    3,0, 1363.9, 0.0, 2.3, 0.0,
    3,1, -2381.0, -82.2, -5.9, 5.7,
    3,2, 1236.2, 241.8, 3.4, -1.0,
    3,3, 525.7, -542.9, -12.2, 1.1,
    4,0, 903.1, 0.0, -1.4, 0.0,
    4,1, 809.4, 282.0, -1.6, 0.6,
    4,2, 86.2, -158.4, -6.0, 6.9,
    4,3, -309.4, 199.8, 5.4, 3.7,
    4,4, 47.9, -350.1, -5.5, -5.6,
    5,0, -234.4, 0.0, -0.3, 0.0,
    5,1, 363.1, 47.7, 0.6, 0.0,
    5,2, 187.8, 208.4, -0.7, 2.5,
    5,3, -140.7, -121.3, 0.1, -0.9,
    5,4, -151.2, 32.2, 1.2, 3.0,
    5,5, 13.7, 99.1, 1.0, 0.5,
    6,0, 65.9, 0.0, -0.6, 0.0,
    6,1, 65.6, -19.1, -0.4, 0.2,
    6,2, 73.0, 25.0, 0.5, -1.8,
    6,3, -121.4, 52.7, 2.1, -1.4,
    6,4, -36.2, -64.4, -1.0, 1.0,
    6,5, 13.5, 9.0, -0.3, 0.1,
    6,6, -64.7, 68.1, 0.8, 0.7,
    7,0, 80.6, 0.0, -0.1, 0.0,
    7,1, -76.7, -51.4, -0.3, 0.5,
    7,2, -8.3, -16.8, 0.0, -0.2,
    7,3, 56.5, 2.3, 0.6, -0.7,
    7,4, 11.8, 23.5, 0.5, -0.2,
    7,5, 3.4, -2.0, 0.5, -0.2,
    7,6, -0.7, -27.1, -0.1, -0.1,
    7,7, 5.8, -8.3, -0.1, -0.1,
    8,0, 24.0, 0.0, -0.1, 0.0,
    8,1, 8.6, 10.2, 0.1, -0.3,
    8,2, -16.9, -18.1, -0.5, 0.2,
    8,3, -3.2, 13.3, 0.5, 0.2,
    8,4, -20.6, -14.6, -0.2, 0.3,
    8,5, 13.3, 16.2, 0.4, -0.1,
    8,6, 11.7, 5.7, 0.2, -0.2,
    8,7, -15.9, -9.1, 0.0, 0.0,
    8,8, -2.0, 2.1, 0.2, 0.0,
    9,0, 5.4, 0.0, -0.1, 0.0,
    9,1, 8.8, -21.6, -0.1, -0.2,
    9,2, 3.1, 10.8, 0.1, -0.1,
    9,3, -3.1, 11.7, 0.2, -0.2,
    9,4, 0.6, -6.8, 0.0, 0.1,
    9,5, -13.3, -6.9, -0.2, 0.0,
    9,6, -0.1, 7.8, 0.0, 0.1,
    9,7, 8.7, 1.0, 0.1, -0.1,
    9,8, -9.1, -3.9, 0.0, 0.0,
    9,9, -10.5, 8.5, 0.0, -0.1,
    10,0, -1.9, 0.0, 0.0, 0.0,
    10,1, -6.5, 3.3, 0.0, 0.0,
    10,2, 0.2, -0.4, 0.0, 0.0,
    10,3, 0.6, 4.6, 0.0, 0.0,
    10,4, -0.6, 4.4, 0.0, 0.0,
    10,5, 1.7, -7.9, 0.0, 0.0,
    10,6, -0.7, -0.6, 0.0, 0.0,
    10,7, 2.1, -4.1, 0.0, 0.0,
    10,8, 2.3, -2.8, 0.0, 0.0,
    10,9, -1.8, -1.1, 0.0, 0.0,
    10,10, -3.6, -8.7, 0.0, 0.0,
    11,0, 3.1, 0.0, 0.0, 0.0,
    11,1, -1.5, -0.1, 0.0, 0.0,
    11,2, -2.3, 2.1, 0.0, 0.0,
    11,3, 2.1, -0.6, 0.0, 0.0,
    11,4, -0.9, -1.1, 0.0, 0.0,
    11,5, 0.6, 0.7, 0.0, 0.0,
    11,6, -1.2, 0.6, 0.0, 0.0,
    11,7, 0.6, -0.2, 0.0, 0.0,
    11,8, -0.1, 0.2, 0.0, 0.0,
    11,9, 0.4, 0.3, 0.0, 0.0,
    11,10, 0.0, -0.5, 0.0, 0.0,
    11,11, -1.0, -0.7, 0.0, 0.0,
    12,0, -2.0, 0.0, 0.0, 0.0,
    12,1, -0.2, 0.1, 0.0, 0.0,
    12,2, 0.3, 1.1, 0.0, 0.0,
    12,3, -0.4, 0.3, 0.0, 0.0,
    12,4, -0.2, -0.3, 0.0, 0.0,
    12,5, 0.2, -0.1, 0.0, 0.0,
    12,6, 0.6, 0.3, 0.0, 0.0,
    12,7, -0.4, -0.3, 0.0, 0.0,
    12,8, -0.2, -0.1, 0.0, 0.0,
    12,9, 0.0, 0.2, 0.0, 0.0,
    12,10, 0.0, 0.0, 0.0, 0.0,
    12,11, 0.0, 0.1, 0.0, 0.0,
    12,12, 0.0, -0.1, 0.0, 0.0,
    13,0, 0.1, 0.0, 0.0, 0.0,
    13,1, -0.3, 0.0, 0.0, 0.0,
    13,2, -0.1, 0.5, 0.0, 0.0,
    13,3, 0.1, -0.5, 0.0, 0.0,
    13,4, -0.1, 0.1, 0.0, 0.0,
    13,5, 0.0, -0.1, 0.0, 0.0,
    13,6, 0.2, -0.3, 0.0, 0.0,
    13,7, 0.0, 0.0, 0.0, 0.0,
    13,8, 0.0, 0.2, 0.0, 0.0,
    13,9, 0.2, -0.1, 0.0, 0.0,
    13,10, -0.1, 0.1, 0.0, 0.0,
    13,11, 0.0, 0.0, 0.0, 0.0,
    13,12, 0.0, -0.1, 0.0, 0.0,
    13,13, 0.7, -0.2, 0.0, 0.0
};

static const int igrf13_ncoeffs = sizeof(igrf13_coeff_data) / (6 * sizeof(double));

static int load_igrf13_coeffs(IGRFModel *model) {
    int nmax = IGRF13_NMAX;
    int ncoeffs = nmax * (nmax + 3) / 2;
    assert(igrf13_ncoeffs == ncoeffs);

    model->coeffs = (GaussCoeff *)malloc(ncoeffs * sizeof(GaussCoeff));
    if (!model->coeffs) return -1;

    for (int i = 0; i < ncoeffs; i++) {
        int n = (int)igrf13_coeff_data[i * 6 + 0];
        int m = (int)igrf13_coeff_data[i * 6 + 1];
        int idx = schmidt_index(n, m) - 1;
        model->coeffs[idx].n = n;
        model->coeffs[idx].m = m;
        model->coeffs[idx].g_nm = igrf13_coeff_data[i * 6 + 2];
        model->coeffs[idx].h_nm = igrf13_coeff_data[i * 6 + 3];
        model->coeffs[idx].dg_nm = igrf13_coeff_data[i * 6 + 4];
        model->coeffs[idx].dh_nm = igrf13_coeff_data[i * 6 + 5];
    }

    model->igrf_version = 13;
    model->epoch = 2020.0;
    model->nmax = nmax;
    model->ncoeffs = ncoeffs;
    snprintf(model->model_name, sizeof(model->model_name), "IGRF-13");
    return 0;
}

int igrf_init_model(int version, IGRFModel *model) {
    if (!model) return -1;
    memset(model, 0, sizeof(IGRFModel));
    if (version == 13) return load_igrf13_coeffs(model);
    return -1;
}

void igrf_free_model(IGRFModel *model) {
    if (model && model->coeffs) {
        free(model->coeffs);
        model->coeffs = NULL;
        model->ncoeffs = 0;
    }
}
