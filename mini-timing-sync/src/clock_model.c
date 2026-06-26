/**
 * @file clock_model.c
 * @brief Clock modeling implementation: polynomial model, regression, drift
 */

#include "clock_model.h"
#include <math.h>
#include <string.h>

/* ===================================================================
 * L4: Polynomial Clock Model
 * =================================================================== */

double clock_model_time_error(const ClockParameters *params, double t)
{
    if (!params) return 0.0;
    /* x(t) = x0 + y0*t + (D*t^2)/2 + (A*t^3)/6
     *
     * Where:
     *   x0 = params->offset_ns         [ns]
     *   y0 = params->freq_offset_ppb * 1e-9 * 1e9 = params->freq_offset_ppb
     *        (ppb means 1e-9 fractional, so y0 dimensionless * 1e9 = ppb value)
     *   D  = params->drift_ppb_per_day / 86400 [ppb/s]
     *        Converted to ns/s: D_ns_per_s2 = D * 1e-9 / 86400 * 1e9 = D / 86400
     *   A  = params->aging_ppb_per_day2 / (86400^2) [ppb/s^2]
     *
     * Time error in ns:
     *   x(t) = offset_ns + freq_offset_ppb * t + 0.5 * drift_ns_per_s2 * t^2
     *          + (1/6) * aging_ns_per_s3 * t^3
     */
    double y0_ns_per_s = params->freq_offset_ppb; /* ppb -> ns/s */
    double D_per_day = params->drift_ppb_per_day;
    double A_per_day2 = params->aging_ppb_per_day2;

    /* D_per_s2 in ppb/s^2 = 1e-9 fractional per s^2.
     * To get ns: multiply by 1e9 => D_per_s2 * 1e9 ns/s^2 = D_per_day / 86400^2 * 1e9
     * Actually simpler: D in ppb/day means D*1e-9 fractional change per day.
     * For time error: integral of y is x.
     *   dx/dt = y (fractional frequency)
     *   y(t) = y0 + D*t + 0.5*A*t^2  (t in days)
     *   x(t) = x0 + integral(y dt)
     */

    double t_days = t / 86400.0;
    double x_ns = params->offset_ns
                + y0_ns_per_s * t
                + 0.5 * D_per_day * t_days * t_days * 1.0e9
                + (1.0 / 6.0) * A_per_day2 * t_days * t_days * t_days * 1.0e9;

    return x_ns;
}

double clock_offset_evolution(double initial_offset_ns,
                              double freq_offset_ppb,
                              double dt_s)
{
    /* offset(t+dt) = offset(t) + freq_offset * dt
     *
     * freq_offset_ppb is in parts-per-billion.
     * 1 ppb = 1e-9 fractional frequency.
     * Over 1 second, accumulated time error = 1e-9 * 1 second = 1 ns.
     * So freq_offset_ppb directly gives ns/s of time drift.
     */
    return initial_offset_ns + freq_offset_ppb * dt_s;
}

/* ===================================================================
 * L5: Least-Squares Linear Regression
 * =================================================================== */

int clock_linear_regression(const double *t_ref, const double *t_device,
                            int N, double *slope, double *intercept,
                            double *r_squared)
{
    if (!t_ref || !t_device || !slope || !intercept || !r_squared) return -1;
    if (N < 2) return -1;

    /* Compute means */
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < N; i++) {
        sum_x += t_ref[i];
        sum_y += t_device[i];
    }
    double mean_x = sum_x / (double)N;
    double mean_y = sum_y / (double)N;

    /* Compute sums for linear regression */
    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0;
    for (int i = 0; i < N; i++) {
        double dx = t_ref[i] - mean_x;
        double dy = t_device[i] - mean_y;
        Sxx += dx * dx;
        Sxy += dx * dy;
        Syy += dy * dy;
    }

    /* Avoid division by zero */
    if (fabs(Sxx) < 1e-30) return -1;

    /* Slope and intercept */
    *slope = Sxy / Sxx;
    *intercept = mean_y - (*slope) * mean_x;

    /* R-squared */
    if (fabs(Syy) < 1e-30) {
        *r_squared = 1.0; /* Perfect fit if all y equal */
    } else {
        double SS_res = 0.0;
        for (int i = 0; i < N; i++) {
            double y_pred = (*slope) * t_ref[i] + (*intercept);
            double residual = t_device[i] - y_pred;
            SS_res += residual * residual;
        }
        *r_squared = 1.0 - SS_res / Syy;
    }

    return 0;
}

/* ===================================================================
 * L5: Clock State Transition
 * =================================================================== */

void clock_state_transition(double state[3], double dt)
{
    if (!state) return;

    double x  = state[0]; /* offset [ns] */
    double y  = state[1]; /* frequency offset [ppb] */
    double d  = state[2]; /* drift [ppb/day] */

    /* Convert drift from ppb/day to ppb/s for dt in seconds */
    double d_per_s = d / 86400.0;

    /* State transition:
     *   x(t+dt) = x(t) + y*dt + 0.5*d*dt^2
     *   y(t+dt) = y(t) + d*dt
     *   d(t+dt) = d(t)
     */
    state[0] = x + y * dt + 0.5 * d_per_s * dt * dt;
    state[1] = y + d_per_s * dt;
    state[2] = d; /* Drift unchanged (constant drift model) */
}

/* ===================================================================
 * L5: Clock Comparison
 * =================================================================== */

ClockComparison clock_compare(const ClockParameters *a,
                              const ClockParameters *b,
                              double comparison_interval_s)
{
    ClockComparison result;
    memset(&result, 0, sizeof(result));

    if (!a || !b) return result;

    result.offset_ns = a->offset_ns - b->offset_ns;
    result.freq_offset_ppb = a->freq_offset_ppb - b->freq_offset_ppb;
    result.num_samples = 1;

    /* RMS residual approximation */
    double freq_diff_ns = result.freq_offset_ppb * comparison_interval_s;
    result.rms_residual_ns = sqrt(
        result.offset_ns * result.offset_ns + freq_diff_ns * freq_diff_ns
    );

    /* 95% confidence (k=2 coverage factor) */
    double combined_uncertainty = sqrt(
        a->standard_uncertainty_ns * a->standard_uncertainty_ns +
        b->standard_uncertainty_ns * b->standard_uncertainty_ns
    );
    result.confidence_95_ns = 2.0 * combined_uncertainty;

    return result;
}

/* ===================================================================
 * L5: Temperature-Compensated Frequency
 * =================================================================== */

double clock_temp_compensated_freq(double base_freq_ppb,
                                   double temp_coeff_ppb_per_c,
                                   double quad_coeff_ppb_per_c2,
                                   double temperature_c,
                                   double ref_temperature_c)
{
    /* Model: y(T) = y0 + alpha*(T-Tref) + beta*(T-Tref)^2
     *
     * This is the standard AT-cut crystal model.
     * For SC-cut crystals, the quadratic term dominates.
     * For TCXOs, the coefficients are pre-characterized and
     * this function is evaluated in real-time with temperature sensor input.
     */
    double delta_T = temperature_c - ref_temperature_c;
    return base_freq_ppb
           + temp_coeff_ppb_per_c * delta_T
           + quad_coeff_ppb_per_c2 * delta_T * delta_T;
}

/* ===================================================================
 * L6: Holdover Error Prediction
 * =================================================================== */

double clock_predict_holdover_error(const ClockParameters *params,
                                    double holdover_duration_s)
{
    if (!params) return 0.0;

    /* Model: U(t) = sqrt(U_x^2 + (U_y*t)^2 + (0.5*U_d*t^2)^2)
     *
     * Where:
     *   U_x = offset uncertainty (standard_uncertainty_ns) [ns]
     *   U_y = frequency offset uncertainty * t [ns]
     *   U_d = drift uncertainty contribution [ns]
     *
     * The drift contribution is from the drift term in the clock model
     * integrated twice (0.5 * drift * t^2).
     */

    double t_days = holdover_duration_s / 86400.0;

    /* Offset uncertainty (initial) */
    double U_x = params->standard_uncertainty_ns;

    /* Frequency uncertainty contribution: freq_offset * t */
    double U_y = fabs(params->freq_offset_ppb) * holdover_duration_s;

    /* Drift contribution: 0.5 * drift * t^2 */
    double U_d = 0.5 * fabs(params->drift_ppb_per_day) * t_days * t_days * 1.0e9;

    /* Aging contribution: (1/6) * aging * t^3 */
    double U_a = (1.0 / 6.0) * fabs(params->aging_ppb_per_day2)
                 * t_days * t_days * t_days * 1.0e9;

    /* RSS combination */
    double total = sqrt(U_x * U_x + U_y * U_y + U_d * U_d + U_a * U_a);

    return total;
}
