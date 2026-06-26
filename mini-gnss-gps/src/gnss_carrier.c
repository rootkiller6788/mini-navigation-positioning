/* =========================================================================
 * gnss_carrier.c — Carrier phase, cycle slip detection, smoothing, DGPS
 *
 * Covers L1 (carrier phase measurement model, integer ambiguity),
 * L2 (cycle slip concepts, double-differencing),
 * L5 (Hatch carrier smoothing, GF/MW linear combinations, cycle slip
 *     detection algorithms), L6 (DGPS corrections).
 *
 * References:
 * - Hatch, R. (1982). Proc. 3rd Int'l Geodetic Symp., 373-386.
 * - Hofmann-Wellenhof, B. et al. (2007). GNSS — GPS, GLONASS, Galileo, 2e.
 * - Teunissen, P.J.G. & Montenbruck, O. (2017). Springer Handbook of GNSS.
 * - Blewitt, G. (1990). "An automatic editing algorithm for GPS data."
 *   Geophysical Research Letters, 17(3), 199-202.
 * ========================================================================= */

#include "gnss_carrier.h"
#include "gnss_common.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* -------------------------------------------------------------------------
 * L5: Hatch carrier smoothing filter
 *
 * Recursive formulation (Hatch 1982):
 *
 *   P_sm[0] = P[0]                              (initialize)
 *   P_sm[n] = (1/N)·P[n] + (1 - 1/N)·(P_sm[n-1] + Φ[n] - Φ[n-1])
 *
 * This is a low-pass filter on code pseudorange, where the high-frequency
 * phase variation provides the dynamics and the code provides the absolute
 * reference. Effectively:
 *
 *   P_sm[n] = Φ[n] + (1-Z⁻¹·(1-1/N))⁻¹ · (1/N)·(P[n]-Φ[n])
 *
 * In the frequency domain, this is a first-order low-pass filter with
 * time constant τ = N·T_epoch. Noise reduction factor = √N.
 *
 * Ionospheric divergence (code-carrier divergence rate):
 *   d/dt (P - Φ) = 2 · dI/dt ≈ 0.02 m/s (daytime max)
 *
 * This limits the practical smoothing window to about N = 100 epochs
 * (10 m error after 1000 s). Dual-frequency iono-free combination
 * removes this constraint.
 * ------------------------------------------------------------------------- */

void gnss_hatch_init(gnss_hatch_filter_t *f, int window_n, int reset_on_slip) {
    f->window_width = window_n;
    f->n_samples = 0;
    f->prev_p_sm = 0.0;
    f->prev_phi = 0.0;
    f->initialized = 0;
    f->reset_on_slip = reset_on_slip;
}

double gnss_hatch_smooth(gnss_hatch_filter_t *f,
                          double pseudorange, double phase_range,
                          int cycle_slip) {
    if (cycle_slip && f->reset_on_slip) {
        gnss_hatch_reset(f);
    }

    double p_sm;

    if (!f->initialized) {
        /* First epoch: use raw pseudorange */
        p_sm = pseudorange;
        f->initialized = 1;
    } else {
        /* Hatch filter update: fixed window weight = 1/N */
        double N = (double)f->window_width;
        double weight = 1.0 / N;

        p_sm = weight * pseudorange
             + (1.0 - weight) * (f->prev_p_sm + phase_range - f->prev_phi);
    }

    /* Update history */
    f->prev_p_sm = p_sm;
    f->prev_phi = phase_range;
    f->n_samples++;

    return p_sm;
}

void gnss_hatch_reset(gnss_hatch_filter_t *f) {
    f->initialized = 0;
    f->n_samples = 0;
}

/* -------------------------------------------------------------------------
 * L2: Geometry-free cycle slip detection
 *
 * GF = λ₁·Φ₁ - λ₂·Φ₂
 *
 * GF = (1 - f₁²/f₂²)·I + (λ₁·N₁ - λ₂·N₂)
 *    = const · I + (constant ambiguity combination)
 *
 * Since I changes slowly (~0.01 cycles/s under quiet conditions,
 * up to ~0.1 cycles/s during ionospheric storms),
 * the GF combination changes by ~λ·ΔN when a cycle slip occurs
 * on either frequency.
 *
 * Detection: |ΔGF| > threshold (typ. 0.05 cycles = 10 cm on L1) → cycle slip.
 *
 * For single-frequency, use phase-rate comparison instead.
 * ------------------------------------------------------------------------- */

double gnss_slip_detect_gf(double phi1_cycles, double phi2_cycles,
                            double prev_phi1, double prev_phi2,
                            double lambda1, double lambda2) {
    /* Geometry-free phase combination in meters */
    double gf_prev = lambda1 * prev_phi1 - lambda2 * prev_phi2;
    double gf_curr = lambda1 * phi1_cycles - lambda2 * phi2_cycles;

    /* ΔGF = change in meters (or cycles w.r.t. L1) */
    double delta_gf = gf_curr - gf_prev;

    /* Normalize to L1 equivalents: ΔN_L1 · λ₁ */
    /* For single-frequency case (no phi2), just check Δφ vs. expected */
    return fabs(delta_gf);
}

/* -------------------------------------------------------------------------
 * L2: Cycle slip detection via Doppler (phase-rate consistency)
 *
 * Expected phase change: ΔΦ_expected = -Doppler · Δt  [cycles]
 * Actual phase change:   ΔΦ_actual   = Φ[n] - Φ[n-1]
 *
 * Slip metric: |ΔΦ_actual - ΔΦ_expected|  [cycles]
 *
 * A cycle slip of 1 cycle gives a mismatch of ~1 cycle = ~19 cm (L1).
 * Detection threshold: 2-3 cycles (to avoid false alarms from noise).
 * ------------------------------------------------------------------------- */

double gnss_slip_detect_doppler(double delta_phi_cycles,
                                 double doppler_hz, double dt) {
    double expected_change = -doppler_hz * dt;
    double delta = fabs(delta_phi_cycles - expected_change);
    return delta;
}

/* -------------------------------------------------------------------------
 * L5: Cycle slip detector initialization and processing
 * ------------------------------------------------------------------------- */

void gnss_slip_detector_init(gnss_slip_detector_t *detector) {
    if (!detector) return;
    detector->n_slips_detected = 0;
    detector->slip_epochs = NULL;
    detector->current_slip = 0;
    detector->metric = 0.0;
    detector->threshold = 0.0;
}

int gnss_slip_detect(gnss_slip_detector_t *detector,
                      double phi1, double phi2,
                      double prev_phi1, double prev_phi2,
                      double doppler_hz, double dt,
                      const gnss_slip_detector_config_t *cfg) {
    if (!detector || !cfg) return 0;

    /* Dual-frequency geometry-free detection */
    double gf_delta = gnss_slip_detect_gf(phi1, phi2, prev_phi1, prev_phi2,
                                           GNSS_L1_WAVELENGTH, GNSS_L2_WAVELENGTH);

    /* Phase-rate (Doppler) consistency */
    double delta_phi1 = phi1 - prev_phi1;
    double dop_delta = gnss_slip_detect_doppler(delta_phi1, doppler_hz, dt);

    /* Combined metric: max of both detectors */
    double metric = fmax(gf_delta / GNSS_L1_WAVELENGTH, dop_delta);

    detector->metric = metric;
    detector->threshold = cfg->jump_threshold;
    detector->current_slip = (metric > cfg->jump_threshold) ? 1 : 0;

    return detector->current_slip;
}

/* -------------------------------------------------------------------------
 * L5: Ionosphere-free linear combination
 *
 * f1²·P1 - f2²·P2
 * P_IF = ─────────────  =  ρ + c·(dt_r - dt_s) + T  (no I)
 *          f1² - f2²
 *
 * The first-order ionospheric term (∝ 1/f²) cancels exactly.
 * Remaining: 2nd and 3rd order iono terms (< 0.1% of total, < 1 cm).
 *
 * Noise amplification:
 *   σ(P_IF) = √((f1⁴ + f2⁴) / (f1²-f2²)²) · σ(P)
 *   For GPS L1/L2: σ(P_IF) ≈ 2.98·σ(P_code) — code noise tripled!
 * ------------------------------------------------------------------------- */

gnss_ionofree_combo_t gnss_ionofree_combination(double f1, double f2) {
    gnss_ionofree_combo_t ic;
    ic.f1 = f1;
    ic.f2 = f2;
    double f1sq = f1 * f1;
    double f2sq = f2 * f2;
    double den = f1sq - f2sq;

    if (fabs(den) < 1e-6) {
        /* Equal frequencies → degenerate (should not happen) */
        ic.coeff_p1 = 1.0;
        ic.coeff_p2 = 0.0;
        ic.noise_factor = 1.0;
    } else {
        ic.coeff_p1 = f1sq / den;
        ic.coeff_p2 = -f2sq / den;
        ic.noise_factor = sqrt((f1sq*f1sq + f2sq*f2sq) / (den*den));
    }

    return ic;
}

double gnss_ionofree_code(double P1, double P2, const gnss_ionofree_combo_t *ic) {
    return ic->coeff_p1 * P1 + ic->coeff_p2 * P2;
}

double gnss_ionofree_phase(double Phi1, double Phi2,
                            const gnss_ionofree_combo_t *ic) {
    return ic->coeff_p1 * Phi1 + ic->coeff_p2 * Phi2;
}

double gnss_ionofree_noise_factor(const gnss_ionofree_combo_t *ic) {
    return ic->noise_factor;
}

/* -------------------------------------------------------------------------
 * L5: Geometry-free linear combination
 *
 * Φ_GF = Φ₁ - Φ₂  (in meters)
 *       = I₁ - I₂ + λ₁·N₁ - λ₂·N₂
 *       = (1 - f₁²/f₂²)·I₁ + (λ₁·N₁ - λ₂·N₂)
 *
 * Contains only ionospheric delay and ambiguity (no geometry, no clock,
 * no troposphere, no orbit errors). Ideal for:
 *   - Cycle slip detection
 *   - Ionospheric TEC monitoring
 *   - Ambiguity estimation when I is known
 * ------------------------------------------------------------------------- */

gnss_geomfree_combo_t gnss_geomfree_combination(double f1, double f2) {
    gnss_geomfree_combo_t gc;
    gc.f1 = f1;
    gc.f2 = f2;
    gc.coeff_phi1 = 1.0;
    gc.coeff_phi2 = -1.0;
    return gc;
}

double gnss_geomfree_phase(double Phi1_m, double Phi2_m,
                            const gnss_geomfree_combo_t *gc) {
    return gc->coeff_phi1 * Phi1_m + gc->coeff_phi2 * Phi2_m;
}

/* -------------------------------------------------------------------------
 * L5: Melbourne-Wübbena (MW) linear combination
 *
 * MW = (L1_L - L2_L) - (P1/P2)·(f1-f2)/(f1+f2)
 *    = λ_WL·N_WL  (wide-lane ambiguity, geometry-free and iono-free)
 *
 * Wide-lane wavelength: λ_WL = c / (f1 - f2)
 * For GPS: λ_WL ≈ 86.2 cm  (L1-L2 wide-lane)
 *           λ_WL ≈ 5.9 m   (L2-L5 wide-lane)
 *
 * The MW combination is ideal for wide-lane ambiguity resolution because:
 *   - Geometry-free (satellite and receiver positions cancel)
 *   - Ionosphere-free (first-order term cancels)
 *   - Long wavelength makes integer rounding easy
 * ------------------------------------------------------------------------- */

gnss_mw_combo_t gnss_mw_combination(double f1, double f2) {
    gnss_mw_combo_t mc;
    mc.f1 = f1;
    mc.f2 = f2;
    mc.wide_lane_lambda = GNSS_C_LIGHT / (f1 - f2);

    /* P_narrow = (f1·P1 + f2·P2) / (f1 + f2) */
    /* MW = (f1·Φ1 - f2·Φ2)/(f1 - f2) - (f1·P1 + f2·P2)/(f1 + f2) */
    mc.coeff_p = 1.0;   /* simplified: MW = Φ_wl - P_nl,
                            Φ_wl = (f1·Φ1-f2·Φ2)/(f1-f2),
                            P_nl = (f1·P1+f2·P2)/(f1+f2) */

    return mc;
}

double gnss_mw_wide_lane(double P_m, double Phi1_m, double Phi2_m,
                          const gnss_mw_combo_t *mc) {
    /* Φ_wide-lane = (f1·Φ1 - f2·Φ2) / (f1-f2)  [meters] */
    double phi_wl = (mc->f1 * Phi1_m - mc->f2 * Phi2_m) / (mc->f1 - mc->f2);

    /* P_narrow-lane = (f1·P1 + f2·P2) / (f1+f2)  [meters]
     * (using one P measurement as approximation) */
    double p_nl = P_m; /* simplified: assume P1 ≈ P2 after corrections */

    /* MW = Φ_wl - P_nl  → divided by λ_wl gives wide-lane ambiguity */
    double mw_value = (phi_wl - p_nl) / mc->wide_lane_lambda;

    return mw_value; /* in wide-lane cycles */
}

/* -------------------------------------------------------------------------
 * L6: DGPS correction application
 *
 * Corrects a user pseudorange using a reference station correction:
 *
 *   P_corr = P_user - PRC - RRC · (t - t₀)
 *
 * Where PRC and RRC are broadcast from the reference station.
 * The UDRE quantifies the quality of the correction.
 * ------------------------------------------------------------------------- */

double gnss_dgps_apply(double pseudorange, const gnss_dgps_correction_t *corr,
                        double time_of_applicability) {
    double dt = time_of_applicability;
    return pseudorange - corr->prc - corr->rrc * dt;
}

int gnss_dgps_position_correction(gnss_ecef_t ref_known,
                                   gnss_ecef_t ref_computed,
                                   double *dx, double *dy, double *dz) {
    if (!dx || !dy || !dz) return -1;
    *dx = ref_known.x - ref_computed.x;
    *dy = ref_known.y - ref_computed.y;
    *dz = ref_known.z - ref_computed.z;
    return 0;
}
