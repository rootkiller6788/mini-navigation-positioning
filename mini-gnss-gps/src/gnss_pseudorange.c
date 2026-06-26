/* =========================================================================
 * gnss_pseudorange.c — Pseudorange models, atmospheric corrections, UERE
 *
 * Covers L1 (pseudorange eqn, error budget), L2 (Klobuchar ionosphere,
 * Saastamoinen & Hopfield troposphere), L5 (CMC analysis, error budgets).
 *
 * References:
 * - Klobuchar, J.A. (1987). IEEE Trans. Aerospace, AES-23(3), 325-331.
 * - Saastamoinen, J. (1972). Bulletin Géodésique, 105(1), 279-298.
 * - Hopfield, H.S. (1969). JGR, 74(18), 4487-4499.
 * - Misra, P. & Enge, P. (2011). Global Positioning System, 2nd ed. Ch.5-6.
 * ========================================================================= */

#include "gnss_pseudorange.h"
#include "gnss_common.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------
 * L5: Klobuchar ionospheric model
 *
 * Algorithm (Klobuchar 1987):
 *
 * 1. Compute Earth-centered angle ψ (geocentric angle between user IPP and sat)
 *    ψ = 0.0137 / (E + 0.11) - 0.022  [semi-circles]
 *
 * 2. Sub-ionospheric point geodetic coordinates:
 *    φ_I = φ_u + ψ·cos(A)
 *    λ_I = λ_u + ψ·sin(A)/cos(φ_I)
 *    Clamp: |φ_I| ≤ 0.416 (≈ 74.88°)
 *
 * 3. Geomagnetic latitude:
 *    φ_m = φ_I + 0.064·cos(λ_I - 1.617)  [semi-circles]
 *
 * 4. Local time at IPP:
 *    t = 4.32e4·λ_I + t_GPS  [seconds]
 *    t = t - 86400 if t ≥ 86400; t = t + 86400 if t < 0
 *
 * 5. Slant factor (obliquity):
 *    F = 1.0 + 16.0·(0.53 - E)³
 *
 * 6. Ionospheric delay amplitude and period:
 *    A_I = α₀ + α₁·φ_m + α₂·φ_m² + α₃·φ_m³  [semi-circles: seconds]
 *    P_I = β₀ + β₁·φ_m + β₂·φ_m² + β₃·φ_m³  [semi-circles: seconds]
 *    (A_I minimum = 0; P_I minimum = 72000 s)
 *
 * 7. Phase: X = 2π·(t - 50400) / P_I
 *
 * 8. Vertical delay: I_z = 5e-9 + A_I·cos(X) for |X| < π/2
 *                        I_z = 5e-9               for |X| ≥ π/2
 *
 * 9. Ionospheric delay = c · F · I_z  [m]
 *
 * Key insight: The broadcast α,β parameters are in seconds and
 * semi-circles. The cosine model peaks at 14:00 local time (50400 s).
 * ------------------------------------------------------------------------- */

static double semi_circle_to_rad(double sc) { return sc * M_PI; }

double gnss_iono_klobuchar(const gnss_klobuchar_params_t *params,
                            gnss_gpstime_t t,
                            gnss_lla_t lla,
                            gnss_ecef_t sat_pos,
                            double freq_hz) {
    /* Convert receiver LLA to semi-circles */
    double phi_u = lla.lat / M_PI;   /* semi-circles */
    double lambda_u = lla.lon / M_PI;
    double elev = gnss_sat_elevation(sat_pos,
                   gnss_lla_to_ecef(lla));
    double azim = gnss_sat_azimuth(sat_pos,
                   gnss_lla_to_ecef(lla));

    double E = elev / M_PI; /* elevation in semi-circles */

    /* Step 1: Earth-centered angle */
    double psi = 0.0137 / (E + 0.11) - 0.022;

    /* Step 2: Sub-ionospheric point */
    double phi_I = phi_u + psi * cos(azim);
    if (phi_I >  0.416) phi_I =  0.416;
    if (phi_I < -0.416) phi_I = -0.416;

    double lambda_I = lambda_u + psi * sin(azim)
                      / cos(semi_circle_to_rad(phi_I));

    /* Step 3: Geomagnetic latitude */
    double phi_m = phi_I + 0.064 * cos(semi_circle_to_rad(lambda_I - 1.617));

    /* Step 4: Local time at IPP */
    double local_t = 43200.0 * lambda_I + t.sow;
    if (local_t >= 86400.0) local_t -= 86400.0;
    if (local_t < 0.0)      local_t += 86400.0;

    /* Step 5: Slant factor */
    double F = 1.0 + 16.0 * pow(0.53 - E, 3.0);

    /* Step 6: Amplitude and period (polynomial in φ_m) */
    double A_I = params->alpha[0] + phi_m * (params->alpha[1]
               + phi_m * (params->alpha[2]
               + phi_m *  params->alpha[3]));
    if (A_I < 0.0) A_I = 0.0; /* non-negative physical constraint */

    double P_I = params->beta[0] + phi_m * (params->beta[1]
               + phi_m * (params->beta[2]
               + phi_m *  params->beta[3]));
    if (P_I < 72000.0) P_I = 72000.0; /* minimum 20 hours */

    /* Step 7: Phase angle */
    double X = 2.0 * M_PI * (local_t - 50400.0) / P_I;

    /* Step 8: Vertical ionospheric delay */
    double I_z_L1; /* L1 delay in seconds */
    if (fabs(X) < M_PI / 2.0) {
        I_z_L1 = 5.0e-9 + A_I * cos(X);
    } else {
        I_z_L1 = 5.0e-9; /* nighttime constant: 5 ns */
    }

    /* L1 delay in meters */
    double I_L1 = GNSS_C_LIGHT * F * I_z_L1;

    /* Scale to user frequency (ionospheric delay ∝ 1/f²) */
    double f1_sq = GNSS_L1_FREQ * GNSS_L1_FREQ;
    double f_sq  = freq_hz * freq_hz;
    double I = I_L1 * (f1_sq / f_sq);

    return I;
}

/* -------------------------------------------------------------------------
 * L5: Saastamoinen tropospheric model
 *
 * Zenith hydrostatic + wet delay, mapped by a continued-fraction
 * mapping function. The standard form used in GPS:
 *
 *   ZHD = 0.0022768 · P / (1 - 0.00266·cos(2φ) - 0.28e-6·h)
 *   ZWD = 0.002277 · (1255/T + 0.05) · e
 *
 * where P = pressure [hPa], T = temperature [K], e = vapor pressure [hPa],
 * φ = latitude, h = height [m] above ellipsoid.
 *
 * Mapping function m(E):
 *   m(E) = 1.001 / √(0.002001 + sin²(E))
 *
 * Total tropospheric delay: T = ZHD·m_d(E) + ZWD·m_w(E)
 * For simplicity, using combined mapping function.
 * ------------------------------------------------------------------------- */

double gnss_tropo_saastamoinen(const gnss_tropo_params_t *params,
                                double elevation_rad, double alt_m) {
    double P = params->pressure_hPa;
    double T = params->temperature_K;
    double e = params->humidity_hPa;

    /* Zenith hydrostatic delay (Saastamoinen 1972, Davis 1985) */
    double zhd = 0.0022768 * P / (1.0 - 0.00266 * cos(2.0 * 0.0) - 2.8e-7 * alt_m);

    /* Zenith wet delay */
    double zwd = 0.002277 * (1255.0 / T + 0.05) * e;

    /* Mapping function (simplified Niell-like) */
    double sin_el = sin(elevation_rad);
    double m_h = 1.001 / sqrt(0.002001 + sin_el * sin_el);
    double m_w = 1.001 / sqrt(0.002001 + sin_el * sin_el); /* same for simplicity */

    return zhd * m_h + zwd * m_w;
}

/* -------------------------------------------------------------------------
 * L5: Hopfield tropospheric model
 *
 * Two-quartic refractivity profile:
 *
 *   N(s) = N_s · (1 - s/H)⁴
 *
 * Dry component: H_d = 40136 + 148.72·(T - 273.16) [m]
 * Wet component: H_w = 11000 [m]
 *
 * Surface refractivity:
 *   N_d0 = 77.6·P/T             [N-units]
 *   N_w0 = 3.73e5·e/T² + 12.0·e/T  [N-units]
 *
 * Path delay in zenith:
 *   δ_d_zenith = 10⁻⁶ · N_d0 · H_d / 5
 *   δ_w_zenith = 10⁻⁶ · N_w0 · H_w / 5
 *
 * Mapping to elevation E:
 *   m(E) = 1 / sin(√(E² + 6.25°²))
 * ------------------------------------------------------------------------- */

double gnss_tropo_hopfield_dry(const gnss_hopfield_params_t *params,
                                double elevation_rad) {
    double T = params->T0;
    double P = params->P0;
    double H_d = 40136.0 + 148.72 * (T - 273.16);
    double N_d0 = 77.6 * P / T;
    double zenith_dry = 1.0e-6 * N_d0 * H_d / 5.0;

    /* Ifantadis mapping function (1986) */
    double el_deg = elevation_rad * 180.0 / M_PI;
    double sin_el = sin(sqrt(el_deg*el_deg + 6.25) * M_PI / 180.0);
    if (sin_el < 0.001) sin_el = 0.001;
    double map = 1.0 / sin_el;

    return zenith_dry * map;
}

double gnss_tropo_hopfield_wet(const gnss_hopfield_params_t *params,
                                double elevation_rad) {
    double T = params->T0;
    double RH = params->RH;
    /* Water vapor pressure from relative humidity (Tetens formula) */
    double e_sat = 6.11 * exp(17.27 * (T - 273.16) / (T - 35.86)); /* hPa */
    double e = RH * e_sat / 100.0;

    double H_w = 11000.0;
    double N_w0 = 3.73e5 * e / (T*T) + 12.0 * e / T;
    double zenith_wet = 1.0e-6 * N_w0 * H_w / 5.0;

    double el_deg = elevation_rad * 180.0 / M_PI;
    double sin_el = sin(sqrt(el_deg*el_deg + 6.25) * M_PI / 180.0);
    if (sin_el < 0.001) sin_el = 0.001;
    double map = 1.0 / sin_el;

    return zenith_wet * map;
}

double gnss_tropo_hopfield_total(const gnss_hopfield_params_t *params,
                                  double elevation_rad) {
    return gnss_tropo_hopfield_dry(params, elevation_rad)
         + gnss_tropo_hopfield_wet(params, elevation_rad);
}

/* -------------------------------------------------------------------------
 * L2: Simple tropospheric model (no meteo data)
 *
 * Uses empirical formula from Collins & Langley (UNB, 1997).
 * T = 2.3 / (sin(E) + 0.15) for temperate latitudes.
 * ------------------------------------------------------------------------- */

double gnss_tropo_simple(double elevation_rad) {
    double sin_el = sin(elevation_rad);
    return 2.3 / (sin_el + 0.15);
}

/* -------------------------------------------------------------------------
 * L5: Combined pseudorange corrections
 * ------------------------------------------------------------------------- */

gnss_range_corrections_t gnss_correct_pseudorange(
    const gnss_pseudorange_t *raw,
    const gnss_ephemeris_t *eph,
    const gnss_klobuchar_params_t *iono,
    const gnss_tropo_params_t *tropo,
    gnss_ecef_t rx_pos) {
    gnss_range_corrections_t corr;
    memset(&corr, 0, sizeof(corr));

    /* Satellite clock correction [m] */
    double dt = raw->time_rx.sow - eph->toe.sow;
    if (dt >  302400.0) dt -= 604800.0;
    if (dt < -302400.0) dt += 604800.0;
    double sat_clk = gnss_sat_clock_correction(eph, dt, 0.0);
    corr.sat_clock_corr = sat_clk * GNSS_C_LIGHT;

    /* Ionospheric correction */
    gnss_lla_t rx_lla = gnss_ecef_to_lla(rx_pos);
    if (iono) {
        corr.iono_corr = gnss_iono_klobuchar(iono, raw->time_rx, rx_lla,
                                              gnss_lla_to_ecef(rx_lla),
                                              GNSS_L1_FREQ);
    }

    /* Tropospheric correction */
    if (tropo) {
        corr.tropo_corr = gnss_tropo_saastamoinen(tropo, raw->sat_elevation, rx_lla.alt);
    } else {
        corr.tropo_corr = gnss_tropo_simple(raw->sat_elevation);
    }

    /* Sagnac effect */
    corr.sagnac_corr = gnss_sagnac_correction(gnss_lla_to_ecef(rx_lla), rx_pos);

    /* Total (note: corrections are subtracted from pseudorange) */
    corr.total_corr = corr.sat_clock_corr + corr.iono_corr
                    + corr.tropo_corr + corr.sagnac_corr;

    return corr;
}

/* -------------------------------------------------------------------------
 * L1: UERE Error Budget
 *
 * Typical GPS Standard Positioning Service (SPS) UERE budget:
 *   Satellite clock & ephemeris: 0.8 m (RMS)
 *   Ionospheric residual: 2.0-7.0 m (single freq, depending on activity)
 *   Tropospheric residual: 0.2 m (zenith), 0.5-2.0 m (low elevation)
 *   Multipath: 0.5-1.5 m (code), 0.5-1.0 cm (carrier)
 *   Receiver noise: 0.3-1.5 m (C/A code), 0.2-1.0 mm (carrier)
 *
 * Space Segment URE: ~0.5 m RMS (GPS III), ~1.5 m RMS (GPS IIR/IIR-M)
 * Total UERE (single freq, no SA): ~4-7 m RMS
 * ------------------------------------------------------------------------- */

gnss_error_budget_t gnss_uere_budget(double elevation_rad, int single_freq) {
    gnss_error_budget_t b;
    b.uere_sat_clock = 0.8; /* GPS broadcast clock & ephemeris (modern) */
    b.uere_iono_residual = single_freq ? (2.0 / sin(elevation_rad + 0.1)) : 0.1;
    b.uere_tropo_residual = 0.2 / sin(elevation_rad + 0.05);
    b.uere_multipath = 0.5;
    b.uere_receiver_noise = 0.3;

    double var_total = b.uere_sat_clock * b.uere_sat_clock
                     + b.uere_iono_residual * b.uere_iono_residual
                     + b.uere_tropo_residual * b.uere_tropo_residual
                     + b.uere_multipath * b.uere_multipath
                     + b.uere_receiver_noise * b.uere_receiver_noise;
    b.uere_total = sqrt(var_total);
    return b;
}

/* -------------------------------------------------------------------------
 * L2: Multipath error model
 *
 * Single-reflector multipath-induced pseudorange error:
 *
 *   δτ = α·D·cos(Δφ) / (1 + α·cos(Δφ))
 *
 * where α = relative amplitude of reflected signal,
 * D = correlator spacing [chips],
 * Δφ = carrier phase difference between direct and reflected.
 *
 * Maximum error at Δφ = 0: δτ_max = α·D / (1+α)
 * ------------------------------------------------------------------------- */

double gnss_multipath_error(const gnss_multipath_t *mp,
                             double chip_rate, double correlator_spacing) {
    double alpha = mp->relative_amplitude;
    double delta_phi = mp->phase_shift_rad;
    double wavelength = GNSS_C_LIGHT / chip_rate;
    double D_m = correlator_spacing * wavelength;

    double numerator = alpha * D_m * cos(delta_phi);
    double denominator = 1.0 + alpha * cos(delta_phi);

    if (fabs(denominator) < 1e-15) return D_m; /* limit case */
    return numerator / denominator;
}

/* -------------------------------------------------------------------------
 * L5: Code-Minus-Carrier (CMC) analysis
 *
 * CMC = P_code - λ·Φ_carrier
 *     = 2·I + M_code + ε_code - m_φ - ε_φ - λ·N
 *
 * Since I varies slowly (minutes to hours) and M_code is high-frequency,
 * CMC can be used to:
 *   - Characterize code multipath environment
 *   - Detect carrier-phase cycle slips (CMC jumps by λ·ΔN)
 *   - Estimate ionospheric TEC changes
 *
 * The ambiguity term λ·N acts as an unknown constant offset, so CMC is
 * typically de-trended or differenced.
 * ------------------------------------------------------------------------- */

int gnss_cmc_compute(const double *code_range, const double *carrier_phase,
                      size_t n_epochs, double wavelength, gnss_cmc_series_t *cmc) {
    if (!code_range || !carrier_phase || !cmc || n_epochs == 0) return -1;

    cmc->n_epochs = n_epochs;
    cmc->code_minus_carrier = (double*)malloc(n_epochs * sizeof(double));
    if (!cmc->code_minus_carrier) return -2;

    double sum = 0.0, sum_sq = 0.0;
    double min_val = 1e30, max_val = -1e30;
    size_t i;

    for (i = 0; i < n_epochs; i++) {
        double cmc_val = code_range[i] - wavelength * carrier_phase[i];
        cmc->code_minus_carrier[i] = cmc_val;
        sum += cmc_val;
        sum_sq += cmc_val * cmc_val;
        if (cmc_val < min_val) min_val = cmc_val;
        if (cmc_val > max_val) max_val = cmc_val;
    }

    cmc->mean = sum / (double)n_epochs;
    double variance = sum_sq / (double)n_epochs - cmc->mean * cmc->mean;
    cmc->std = sqrt(variance > 0.0 ? variance : 0.0);
    cmc->max_abs = fmax(fabs(min_val), fabs(max_val));

    return 0;
}

void gnss_cmc_free(gnss_cmc_series_t *cmc) {
    if (cmc && cmc->code_minus_carrier) {
        free(cmc->code_minus_carrier);
        cmc->code_minus_carrier = NULL;
    }
}
