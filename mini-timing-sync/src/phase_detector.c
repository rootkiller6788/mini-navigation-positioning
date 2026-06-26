/**
 * @file phase_detector.c
 * @brief Phase and frequency detector implementations, DPLL components
 *
 * Implements multiple phase detector types, PLL loop dynamics,
 * DPLL update algorithm, and lock detection.
 */

#include "phase_detector.h"
#include <math.h>
#include <string.h>

/* ===================================================================
 * L4: PLL Loop Dynamics
 * =================================================================== */

void pll_compute_dynamics(const LoopFilterCoeffs *coeffs, double K_gain,
                          PllDynamics *dyn)
{
    if (!coeffs || !dyn) return;

    double T = coeffs->T_update_s;
    double alpha = coeffs->alpha;
    double beta = coeffs->beta;

    /* Type-2 PLL (PI loop filter) dynamics:
     *
     * Open-loop transfer function:
     *   G(s) = K * (alpha + beta/s) * (1/s) = K*beta/s^2 + K*alpha/s
     *
     * Characteristic equation: s^2 + K*alpha*s + K*beta = 0
     *
     * Natural frequency: omega_n = sqrt(K * beta)
     * Damping factor:    zeta    = K*alpha / (2 * omega_n)
     *                            = alpha/2 * sqrt(K/beta)
     *
     * Noise bandwidth (for type-2):
     *   BL = omega_n/2 * (zeta + 1/(4*zeta))
     */

    if (K_gain <= 0.0 || beta <= 0.0) {
        dyn->natural_freq_hz = 0.0;
        dyn->damping_factor = 0.0;
        dyn->loop_bandwidth_hz = 0.0;
        dyn->phase_margin_deg = 0.0;
        return;
    }

    double omega_n = sqrt(K_gain * beta / T);
    double zeta = (alpha / (2.0 * T)) * sqrt(T / (K_gain * beta));

    dyn->natural_freq_hz = omega_n / (2.0 * M_PI);
    dyn->damping_factor = zeta;

    /* Noise bandwidth formula for type-2 PLL */
    if (zeta > 0.0) {
        dyn->loop_bandwidth_hz = (omega_n / (2.0 * M_PI))
                                 * (zeta + 1.0 / (4.0 * zeta));
    } else {
        dyn->loop_bandwidth_hz = 0.0;
    }

    /* Phase margin from damping factor:
     * PM = arctan(2*zeta / sqrt(sqrt(4*zeta^4 + 1) - 2*zeta^2))
     * Approximated for typical zeta values:
     */
    if (zeta <= 1.0) {
        dyn->phase_margin_deg = 100.0 * zeta;
    } else {
        dyn->phase_margin_deg = 65.0 + 20.0 * log10(zeta);
    }
}

double pll_lock_range(double pd_gain, double vco_gain,
                      PhaseDetectorType pd_type)
{
    /* Lock range (hold-in range) calculation */
    switch (pd_type) {
    case PD_TYPE_LINEAR:
        /* Linear (multiplier) PD: lock_range = K * pi/2 */
        return pd_gain * vco_gain * M_PI / 2.0;

    case PD_TYPE_XOR:
        /* XOR PD: lock_range = K * pi/2 (same as linear for square waves) */
        return pd_gain * vco_gain * M_PI / 2.0;

    case PD_TYPE_PFD:
        /* PFD: Theoretical lock range is infinite (type-2 PLL).
         * Return 0 to indicate "unlimited" */
        return 0.0;

    case PD_TYPE_HOGGE:
        /* Hogge PD for CDR: lock_range limited by data rate / 2 */
        return vco_gain * M_PI / 2.0;

    case PD_TYPE_ALEXANDER:
        /* Bang-bang PD: lock range = VCO tuning range */
        return vco_gain * M_PI;

    case PD_TYPE_TANLOCK:
        /* Tanlock: similar to linear but extended range */
        return pd_gain * vco_gain * 1.3 * M_PI;

    default:
        return 0.0;
    }
}

/* ===================================================================
 * L5: Phase Detector Implementations
 * =================================================================== */

double phase_error_from_samples(double sample1, double sample2, double K_pd)
{
    /* Product detector for sinusoidal inputs:
     *
     * s1 = A1 * sin(wt + phi1)
     * s2 = A2 * sin(wt + phi2)
     *
     * Product: s1*s2 = (A1*A2/2)*[cos(phi1-phi2) - cos(2wt+phi1+phi2)]
     *
     * After lowpass filtering (averaging), the 2wt term averages to zero:
     *   V_out ~= (A1*A2/2) * cos(phi1 - phi2)
     *
     * For small errors: V_out ~= K_pd * (phi1 - phi2)
     * where K_pd = -A1*A2/2 (around phi1-phi2 = pi/2 quadrature point)
     *
     * Normalized output: arcsin(V_out / (A1*A2/2)) = phase error
     *
     * Simplified: use product directly as phase error estimate.
     */
    double product = sample1 * sample2;

    /* Normalize by K_pd and compute phase */
    if (fabs(K_pd) < 1e-30) return 0.0;

    double normalized = product / K_pd;

    /* Clamp to [-1, 1] for arcsin domain */
    if (normalized > 1.0) normalized = 1.0;
    if (normalized < -1.0) normalized = -1.0;

    return asin(normalized);
}

double xor_phase_detector(int sig1_level, int sig2_level, double v_high)
{
    /* XOR phase detector:
     *
     * For two 50% duty cycle square waves with phase difference phi:
     *   XOR output duty cycle = |phi| / pi
     *   Average voltage = V_high * |phi| / pi
     *
     * Range: phi in [-pi/2, pi/2] (ambiguous beyond this)
     *
     * Output sign: positive if sig2 lags sig1.
     */
    int xor_out = (sig1_level != sig2_level) ? 1 : 0;

    /* In a real implementation, this value would be lowpass filtered.
     * Here we return the instantaneous value, which requires averaging.
     *
     * Scaled output: -V_high to +V_high for -pi/2 to +pi/2.
     *
     * To get bidirectional output, we need edge information.
     * Simplified: return +v_high when XOR high, -v_high when low.
     * The average after filtering gives phase error.
     */
    return xor_out ? v_high : -v_high;
}

double pfd_output(double ref_edge_time, double div_edge_time, double T_ref)
{
    /* Phase-Frequency Detector (PFD):
     *
     * Output states:
     *   UP = 1 when ref_edge leads div_edge
     *   DN = 1 when div_edge leads ref_edge
     *   Both clear after both edges detected.
     *
     * Net output over one cycle:
     *   Q = (Delta_t_lead) / T_ref   for ref leading
     *   Q = -(Delta_t_lag) / T_ref   for div leading
     *
     * Range: [-2*pi, 2*pi] theoretically (unlike [-pi/2, pi/2] for XOR)
     * This gives PFD its frequency discrimination capability.
     */
    if (T_ref <= 0.0) return 0.0;

    double delta_t = ref_edge_time - div_edge_time;

    /* Normalize to reference period */
    return delta_t / T_ref;
}

int hogge_phase_detect(double data_sample, double prev_data,
                       double edge_sample)
{
    /* Hogge phase detector for CDR:
     *
     * Uses two samplers:
     *   A: samples data at the center of the eye (decision point)
     *   B: samples data at the transition (edge)
     *
     * Output = A * (prev_B - B) or equivalently (A XOR prev_A) * (A XOR B)
     *
     * Simplified version:
     *   If data transition occurs (data != prev_data):
     *     error = sign(data) * edge_sample
     *
     * Positive error -> clock is early (edge samples before transition)
     * Negative error -> clock is late (edge samples after transition)
     */
    double transition = data_sample - prev_data;

    /* No transition -> no phase information */
    if (fabs(transition) < 1e-12) return 0;

    /* Phase error: transition sign * edge value */
    double error = (transition > 0 ? 1.0 : -1.0) * edge_sample;

    if (error > 0.5) return 1;
    if (error < -0.5) return -1;
    return 0;
}

/* ===================================================================
 * L6: DPLL Update and Lock Detection
 * =================================================================== */

void dpll_init(DpllState *dpll, double nominal_freq_hz,
               double sampling_freq_hz, int nco_bit_width)
{
    if (!dpll) return;
    (void)nco_bit_width; /* NCO bit width for future quantization model */

    double phase_step = nominal_freq_hz / sampling_freq_hz;

    dpll->phase_accumulator = 0.0;
    dpll->frequency_word = phase_step;
    dpll->integral_path = 0.0;
    dpll->double_integral_path = 0.0;
    dpll->lock_metric = 0.0;
    dpll->lock_count = 0;
}

void dpll_update(double phase_error_rad, const LoopFilterCoeffs *coeffs,
                 DpllState *dpll, double *freq_out_hz)
{
    if (!coeffs || !dpll) return;

    /* Loop filter: PI controller
     *
     * Filter output = alpha * phase_error + beta * integral(phase_error)
     *
     * For type-3 PLL, add gamma * double_integral(phase_error).
     */

    /* Update integral path */
    dpll->integral_path += coeffs->beta * phase_error_rad;
    dpll->double_integral_path += coeffs->gamma * phase_error_rad;

    /* Loop filter output = proportional + integral + double-integral */
    double filter_output = coeffs->alpha * phase_error_rad
                         + dpll->integral_path
                         + dpll->double_integral_path;

    /* Update NCO frequency word */
    dpll->frequency_word += filter_output;

    /* Prevent negative frequency */
    if (dpll->frequency_word < 0.0) {
        dpll->frequency_word = 0.0;
    }

    /* Update phase accumulator */
    dpll->phase_accumulator += dpll->frequency_word;

    /* Keep phase accumulator in [0, 1) range */
    dpll->phase_accumulator = dpll->phase_accumulator -
                              floor(dpll->phase_accumulator);

    /* Output frequency estimate */
    if (freq_out_hz) {
        *freq_out_hz = dpll->frequency_word;
    }

    /* Lock metric: exponential moving average of |phase_error| */
    double alpha_lock = 0.01;
    dpll->lock_metric = (1.0 - alpha_lock) * dpll->lock_metric
                       + alpha_lock * fabs(phase_error_rad);
}

int dpll_lock_detect(double phase_error_rad, double lock_threshold_rad2,
                     int lock_count_needed, DpllState *dpll)
{
    if (!dpll) return 0;

    /* Check if current phase error variance-like metric is below threshold */
    double var_est = phase_error_rad * phase_error_rad;

    if (var_est < lock_threshold_rad2) {
        dpll->lock_count++;
        if (dpll->lock_count >= lock_count_needed) {
            return 1; /* Locked */
        }
    } else {
        dpll->lock_count = 0;
    }

    return 0; /* Not locked */
}

int detect_cycle_slip(double current_phase_error, double prev_phase_error,
                      double threshold_rad)
{
    /* Cycle slip detection:
     *
     * A cycle slip occurs when the phase error jumps by ~2*pi.
     * Detection: |current - prev| > threshold (typically ~pi).
     *
     * This is common during:
     * - High dynamic stress on the tracking loop
     * - Signal fades (GPS)
     * - PLL unlock events
     */
    double delta = current_phase_error - prev_phase_error;

    /* Wrap delta to [-pi, pi] */
    while (delta > M_PI)  delta -= 2.0 * M_PI;
    while (delta < -M_PI) delta += 2.0 * M_PI;

    if (fabs(delta) > threshold_rad) {
        return 1;
    }
    return 0;
}
