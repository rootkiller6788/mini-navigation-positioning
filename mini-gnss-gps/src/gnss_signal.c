/* =========================================================================
 * gnss_signal.c — C/A code generation, Doppler, satellite position, Kepler
 *
 * Covers L1 (signal structure, frequencies, PRN codes),
 * L2 (Gold codes, Doppler shift), L4 (Kepler's laws applied to GPS),
 * L5 (satellite position from ephemeris, Newton-Raphson for Kepler's eqn).
 *
 * References:
 * - IS-GPS-200 §3.3.2 (C/A code generation)
 * - IS-GPS-200 §20.3.3.3.3.1 (satellite position computation)
 * - Tsui, J. (2005). Fundamentals of GPS Receivers, 2nd ed. (Ch.4)
 * - Gold, R. (1967). "Optimal binary sequences for SS comms." IEEE Trans IT.
 * ========================================================================= */

#include "gnss_signal.h"
#include "gnss_common.h"
#include <string.h>
#include <stdlib.h>

/* -------------------------------------------------------------------------
 * L1: Band descriptors
 * ------------------------------------------------------------------------- */

gnss_band_desc_t gnss_band_descriptor(gnss_band_t band) {
    gnss_band_desc_t d;
    d.band = band;
    switch (band) {
    case GNSS_BAND_L1:
        d.freq_hz = GNSS_L1_FREQ;
        d.chip_rate = GNSS_CA_CHIP_RATE;
        break;
    case GNSS_BAND_L2:
        d.freq_hz = GNSS_L2_FREQ;
        d.chip_rate = 511500.0; /* L2C chip rate */
        break;
    case GNSS_BAND_L5:
        d.freq_hz = GNSS_L5_FREQ;
        d.chip_rate = 10.23e6;
        break;
    case GNSS_BAND_E1:
        d.freq_hz = 1575.42e6;
        d.chip_rate = 1.023e6;
        break;
    case GNSS_BAND_E5a:
        d.freq_hz = 1176.45e6;
        d.chip_rate = 10.23e6;
        break;
    case GNSS_BAND_B1I:
        d.freq_hz = 1561.098e6;
        d.chip_rate = 2.046e6;
        break;
    case GNSS_BAND_B2I:
        d.freq_hz = 1207.14e6;
        d.chip_rate = 2.046e6;
        break;
    default:
        d.freq_hz = 0.0;
        d.chip_rate = 0.0;
        break;
    }
    d.wavelength_m = (d.freq_hz > 0.0) ? (GNSS_C_LIGHT / d.freq_hz) : 0.0;
    return d;
}

/* -------------------------------------------------------------------------
 * L2: C/A Code Generation (Gold codes)
 *
 * GPS C/A code is a Gold code of length 1023, generated from two
 * 10-bit maximal-length LFSRs:
 *
 * G1 polynominal: 1 + X³ + X¹⁰  → G1 register taps at bits 3 and 10
 * G2 polynominal: 1 + X² + X³ + X⁶ + X⁸ + X⁹ + X¹⁰
 *
 * The G2 output is delayed (phase-selected) by a specific pair of G2-taps
 * per satellite PRN (IS-GPS-200 Table 3-Ia).
 *
 * Initial state: all-ones (1111111111₂ = 1023).
 * Output is the modulo-2 sum of G1[10] and the G2 tapped bits.
 * ------------------------------------------------------------------------- */

/* Phase-selector taps per PRN (Table 3-Ia subset: PRN 1-32) */
static const int g2_phase_taps[37][2] = {
    {0,0}, /* PRN 0 — not used */
    {2,6},   {3,7},   {4,8},   {5,9},   /* PRN 1-4 */
    {1,9},   {2,10},  {1,8},   {2,9},   /* PRN 5-8 */
    {3,10},  {2,3},   {3,4},   {5,6},   /* PRN 9-12 */
    {6,7},   {7,8},   {8,9},   {9,10},  /* PRN 13-16 */
    {1,4},   {2,5},   {3,6},   {4,7},   /* PRN 17-20 */
    {5,8},   {6,9},   {1,3},   {4,6},   /* PRN 21-24 */
    {5,7},   {6,8},   {7,9},   {8,10},  /* PRN 25-28 */
    {1,6},   {2,7},   {3,8},   {4,9}    /* PRN 29-32 */
};
/* Extend for PRN 33-37 */
static const int g2_phase_taps_ext[5][2] = {
    {5,10}, {1,4}, {2,5}, {4,5}, {6,9}
};

static int get_g2_taps(int32_t prn, int taps[2]) {
    if (prn >= 1 && prn <= 32) {
        taps[0] = g2_phase_taps[prn][0];
        taps[1] = g2_phase_taps[prn][1];
        return 1;
    } else if (prn >= 33 && prn <= 37) {
        taps[0] = g2_phase_taps_ext[prn - 33][0];
        taps[1] = g2_phase_taps_ext[prn - 33][1];
        return 1;
    }
    return 0; /* Invalid PRN */
}

/* G1 feedback bit: tap at position 3 XOR position 10 */
static uint16_t g1_feedback_bit(uint16_t reg) {
    uint16_t b3  = (reg >> 2) & 1;   /* bit index 2 (0-based) = tap 3 */
    uint16_t b10 = (reg >> 9) & 1;   /* bit index 9 = tap 10 */
    return b3 ^ b10;
}

/* G2 feedback bit (polynomial: 1 + X² + X³ + X⁶ + X⁸ + X⁹ + X¹⁰) */
static uint16_t g2_feedback_bit(uint16_t reg) {
    uint16_t b2  = (reg >> 1) & 1;   /* bit index 1 = tap 2 */
    uint16_t b3  = (reg >> 2) & 1;   /* tap 3 */
    uint16_t b6  = (reg >> 5) & 1;   /* tap 6 */
    uint16_t b8  = (reg >> 7) & 1;   /* tap 8 */
    uint16_t b9  = (reg >> 8) & 1;   /* tap 9 */
    uint16_t b10 = (reg >> 9) & 1;   /* tap 10 */
    return b2 ^ b3 ^ b6 ^ b8 ^ b9 ^ b10;
}

/* Shift LFSR: shift right 1, insert feedback at MSB (bit 9) */
static uint16_t g1_shift(uint16_t reg) {
    uint16_t fb = g1_feedback_bit(reg);
    return ((reg >> 1) & 0x1FF) | (fb << 9);
}

static uint16_t g2_shift(uint16_t reg) {
    uint16_t fb = g2_feedback_bit(reg);
    return ((reg >> 1) & 0x1FF) | (fb << 9);
}

int gnss_ca_code_generate(int32_t prn, gnss_ca_code_t *code) {
    int taps[2];
    if (!get_g2_taps(prn, taps)) return -1;

    code->prn = prn;

    /* G1 and G2 initialized to all-ones (10-bit: 0x3FF = 1023) */
    uint16_t g1 = 0x3FF;
    uint16_t g2 = 0x3FF;

    int chip;
    for (chip = 0; chip < 1023; chip++) {
        /* G2 phase-selected output: XOR of two tapped bits */
        uint16_t g2_bit0 = (g2 >> (taps[0] - 1)) & 1;
        uint16_t g2_bit1 = (g2 >> (taps[1] - 1)) & 1;
        uint16_t g2_out = g2_bit0 ^ g2_bit1;

        /* G1 output is bit 10 (MSB) */
        uint16_t g1_out = (g1 >> 9) & 1;

        /* C/A chip = G1 XOR G2_out, mapped to ±1 */
        code->code[chip] = ((g1_out ^ g2_out) == 0) ? 1 : -1;

        /* Shift both registers */
        g1 = g1_shift(g1);
        g2 = g2_shift(g2);
    }

    return 0;
}

int gnss_ca_code_generate_shifted(int32_t prn, int32_t shift_chips,
                                   gnss_ca_code_t *code) {
    if (gnss_ca_code_generate(prn, code) != 0) return -1;

    /* Cyclic shift: move to desired code phase */
    shift_chips = ((shift_chips % 1023) + 1023) % 1023;
    if (shift_chips == 0) return 0;

    int temp[1023];
    memcpy(temp, code->code, sizeof(temp));
    int i;
    for (i = 0; i < 1023; i++) {
        code->code[i] = temp[(i + shift_chips) % 1023];
    }
    return 0;
}

/* -------------------------------------------------------------------------
 * L2: C/A code correlation properties
 * ------------------------------------------------------------------------- */

double gnss_ca_autocorrelation(const gnss_ca_code_t *code, int32_t lag) {
    lag = ((lag % 1023) + 1023) % 1023;
    int sum = 0, i;
    for (i = 0; i < 1023; i++) {
        sum += code->code[i] * code->code[(i + lag) % 1023];
    }
    return (double)sum / 1023.0;
}

double gnss_ca_crosscorrelation(const gnss_ca_code_t *a,
                                 const gnss_ca_code_t *b) {
    if (a->prn == b->prn) return gnss_ca_autocorrelation(a, 0);
    int sum = 0, i;
    for (i = 0; i < 1023; i++) {
        sum += a->code[i] * b->code[i];
    }
    return (double)sum / 1023.0;
}

/* -------------------------------------------------------------------------
 * L2: Doppler computation
 *
 * Δf = - (f_carrier / c) · v_rel
 * where v_rel = e·(v_sat - v_rx) is the radial relative velocity
 * (dot product of LOS unit vector with relative velocity vector)
 * ------------------------------------------------------------------------- */

gnss_doppler_t gnss_doppler_compute(gnss_ecef_t sat_pos, gnss_ecef_t rx_pos,
                                     const double sat_vel[3],
                                     const double rx_vel[3],
                                     double carrier_freq) {
    gnss_doppler_t dop;
    double dx = sat_pos.x - rx_pos.x;
    double dy = sat_pos.y - rx_pos.y;
    double dz = sat_pos.z - rx_pos.z;
    double range = sqrt(dx*dx + dy*dy + dz*dz);

    if (range < 1.0) {
        dop.range_rate = 0.0;
        dop.doppler_shift = 0.0;
        dop.relative_velocity[0] = 0.0;
        dop.relative_velocity[1] = 0.0;
        dop.relative_velocity[2] = 0.0;
        return dop;
    }

    /* LOS unit vector from receiver to satellite */
    double e_x = dx / range, e_y = dy / range, e_z = dz / range;

    /* Relative velocity vector */
    dop.relative_velocity[0] = sat_vel[0] - rx_vel[0];
    dop.relative_velocity[1] = sat_vel[1] - rx_vel[1];
    dop.relative_velocity[2] = sat_vel[2] - rx_vel[2];

    /* Radial component (range rate): positive = receding */
    dop.range_rate = e_x * dop.relative_velocity[0]
                   + e_y * dop.relative_velocity[1]
                   + e_z * dop.relative_velocity[2];

    /* Doppler shift: negative for approaching satellite */
    dop.doppler_shift = -(carrier_freq / GNSS_C_LIGHT) * dop.range_rate;

    return dop;
}

double gnss_doppler_to_range_rate(double doppler_hz, double carrier_freq) {
    return -(doppler_hz / carrier_freq) * GNSS_C_LIGHT;
}

/* -------------------------------------------------------------------------
 * L4: Kepler's equation solver (Newton-Raphson)
 *
 * M = E - e·sin(E)
 *
 * Newton iteration: E_{k+1} = E_k - (E_k - e·sin(E_k) - M) / (1 - e·cos(E_k))
 *
 * For GPS (e ≈ 0.01), convergence in 3-4 iterations to 1e-12 radians.
 * ------------------------------------------------------------------------- */

double gnss_kepler_solve(double M_rad, double eccentricity,
                          int max_iter, double tol) {
    /* Initial guess: E₀ = M for small eccentricity;
     * for moderate e, improve with E₀ = M + e·sin(M) */
    double E = M_rad + eccentricity * sin(M_rad);

    /* Special case: e = 0 → E = M directly */
    if (eccentricity < 1e-15) return M_rad;

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double sin_E = sin(E);
        double cos_E = cos(E);
        double f = E - eccentricity * sin_E - M_rad;
        double f_prime = 1.0 - eccentricity * cos_E;

        if (fabs(f_prime) < 1e-30) break;

        double delta = -f / f_prime;
        E += delta;

        if (fabs(delta) < tol) return E;
    }
    return E;
}

/* -------------------------------------------------------------------------
 * L4: Satellite clock correction (relativistic)
 *
 * Δt_sv = af0 + af1·Δt + af2·Δt² + Δt_rel
 *
 * Relativistic correction due to orbital eccentricity:
 *   Δt_rel = -2·√(GM·a)·e·sin(E) / c²
 *
 * For GPS: √(GM·a)·e·sin(E) ≈ 2290·e·sin(E) [ns] for typical semi-major axis.
 * The factor -2/c² converts to seconds: -2/c² ≈ -2.226×10⁻¹¹ s/(m²/s²)
 *
 * Actually: Δt_rel = -2·(r·v)/(c²) where r,v are ECEF position/velocity.
 * Equivalently: Δt_rel = -2·√(GM·a)·e·sin(E)/c² [seconds]
 * ------------------------------------------------------------------------- */

double gnss_sat_clock_correction(const gnss_ephemeris_t *eph,
                                  double delta_t, double E_rad) {
    double dt_clk = eph->af0 + eph->af1 * delta_t + eph->af2 * delta_t * delta_t;

    /* Relativistic correction */
    double a = eph->sqrt_a * eph->sqrt_a;
    double rel_corr = -2.0 * sqrt(GNSS_GM_EARTH * a)
                      * eph->e * sin(E_rad)
                      / (GNSS_C_LIGHT * GNSS_C_LIGHT);

    dt_clk += rel_corr - eph->tgd; /* Subtract TGD for single-frequency users */
    return dt_clk;
}

/* -------------------------------------------------------------------------
 * L5: Satellite ECEF position from broadcast ephemeris
 *
 * Follows IS-GPS-200 §20.3.3.3.3.1 exactly.
 *
 * Step-by-step:
 *  1. Δt = t - toe  (with week rollover handling)
 *  2. n = n₀ + Δn   where n₀ = √(GM / a³)
 *  3. M = M₀ + n·Δt
 *  4. E from Kepler's equation: M = E - e·sin(E)
 *  5. ν = 2·atan(√((1+e)/(1-e))·tan(E/2))
 *  6. φ = ν + ω
 *  7. Perturbations: δu = Cus·sin(2φ) + Cuc·cos(2φ)
 *                     δr = Crs·sin(2φ) + Crc·cos(2φ)
 *                     δi = Cis·sin(2φ) + Cic·cos(2φ)
 *  8. u = φ + δu, r = a·(1 - e·cos(E)) + δr, i = i₀ + i_dot·Δt + δi
 *  9. In-plane: x' = r·cos(u), y' = r·sin(u)
 * 10. Ω = Ω₀ + (Ω_dot - ω_e)·Δt - ω_e·toe
 * 11. ECEF: x = x'·cos(Ω) - y'·cos(i)·sin(Ω)
 *            y = x'·sin(Ω) + y'·cos(i)·cos(Ω)
 *            z = y'·sin(i)
 * ------------------------------------------------------------------------- */

int gnss_satpos_from_ephemeris(const gnss_ephemeris_t *eph,
                                gnss_gpstime_t t_rx,
                                gnss_ecef_t *sat_pos,
                                double sat_vel[3],
                                double *clock_bias) {
    if (!eph || !sat_pos) return -1;

    /* 1. Time from ephemeris (handle week boundary) */
    double dt = t_rx.sow - eph->toe.sow;
    if (dt >  302400.0) dt -= 604800.0;
    if (dt < -302400.0) dt += 604800.0;

    /* 2. Mean motion */
    double a = eph->sqrt_a * eph->sqrt_a;
    double n0 = sqrt(GNSS_GM_EARTH / (a * a * a));
    double n = n0 + eph->delta_n;

    /* 3. Mean anomaly */
    double M = eph->M0 + n * dt;

    /* 4. Solve Kepler's equation */
    double E = gnss_kepler_solve(M, eph->e, 20, 1e-12);

    /* 5. True anomaly */
    double sin_E = sin(E), cos_E = cos(E);
    double sin_nu = sqrt(1.0 - eph->e * eph->e) * sin_E / (1.0 - eph->e * cos_E);
    double cos_nu = (cos_E - eph->e) / (1.0 - eph->e * cos_E);
    double nu = atan2(sin_nu, cos_nu);

    /* 6. Argument of latitude */
    double phi = nu + eph->w;

    /* 7. Second harmonic perturbations */
    double sin_2phi = sin(2.0 * phi);
    double cos_2phi = cos(2.0 * phi);
    double du = eph->cus * sin_2phi + eph->cuc * cos_2phi;
    double dr = eph->crs * sin_2phi + eph->crc * cos_2phi;
    double di = eph->cis * sin_2phi + eph->cic * cos_2phi;

    /* 8. Corrected orbital elements */
    double u = phi + du;
    double r = a * (1.0 - eph->e * cos_E) + dr;
    double i = eph->i0 + eph->i_dot * dt + di;

    /* 9. In-plane coordinates */
    double cos_u = cos(u), sin_u = sin(u);
    double x_prime = r * cos_u;
    double y_prime = r * sin_u;

    /* 10. Corrected longitude of ascending node */
    double omega = eph->omega0
                 + (eph->omega_dot - GNSS_OMEGA_E) * dt
                 - GNSS_OMEGA_E * eph->toe.sow;

    double cos_omega = cos(omega), sin_omega = sin(omega);
    double cos_i = cos(i), sin_i = sin(i);

    /* 11. ECEF coordinates */
    sat_pos->x = x_prime * cos_omega - y_prime * cos_i * sin_omega;
    sat_pos->y = x_prime * sin_omega + y_prime * cos_i * cos_omega;
    sat_pos->z = y_prime * sin_i;

    /* Velocities (approximate: derivative of position) */
    if (sat_vel) {
        double dE_dt = n / (1.0 - eph->e * cos_E); /* dE/dt */
        double dnu_dt = sqrt(1.0 - eph->e*eph->e) * dE_dt
                       / (1.0 - eph->e * cos_E);
        double du_dt = dnu_dt;
        double dr_dt = a * eph->e * sin_E * dE_dt;
        double di_dt = eph->i_dot;
        double domega_dt = eph->omega_dot - GNSS_OMEGA_E;

        double dxp_dt = dr_dt * cos_u - r * sin_u * du_dt;
        double dyp_dt = dr_dt * sin_u + r * cos_u * du_dt;

        double sin_i_di = sin_i * di_dt;
        double cos_i_di = cos_i * di_dt;

        sat_vel[0] = dxp_dt * cos_omega - x_prime * sin_omega * domega_dt
                   - (dyp_dt * cos_i - y_prime * sin_i_di) * sin_omega
                   - y_prime * cos_i * cos_omega * domega_dt;
        sat_vel[1] = dxp_dt * sin_omega + x_prime * cos_omega * domega_dt
                   + (dyp_dt * cos_i - y_prime * sin_i_di) * cos_omega
                   - y_prime * cos_i * sin_omega * domega_dt;
        sat_vel[2] = dyp_dt * sin_i + y_prime * cos_i_di;
    }

    /* Clock correction */
    if (clock_bias) {
        *clock_bias = gnss_sat_clock_correction(eph, dt, E);
    }

    return 0;
}
