/**
 * @file clock_model.h
 * @brief Clock modeling: offset, skew, drift, aging models
 *
 * L1 Definitions: OscillatorSpec, ClockParameters, ClockComparison
 * L2 Concepts: Clock state evolution, polynomial clock model
 * L3 Math: Polynomial fitting, linear regression on timestamp series
 * L4 Laws: Clock offset equation, frequency stability bound
 *
 * Reference: IEEE 1588-2019 Annex B
 * Textbook: Riley "Handbook of Frequency Stability Analysis" (NIST SP 1065)
 */

#ifndef CLOCK_MODEL_H
#define CLOCK_MODEL_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double nominal_freq_hz;
    double freq_offset_ppb;
    double drift_ppb_per_day;
    double temperature_coeff_ppb_per_c;
} OscillatorSpec;

typedef struct {
    double offset_ns;
    double freq_offset_ppb;
    double drift_ppb_per_day;
    double aging_ppb_per_day2;
    double standard_uncertainty_ns;
} ClockParameters;

typedef struct {
    double offset_ns;
    double freq_offset_ppb;
    double rms_residual_ns;
    int    num_samples;
    double confidence_95_ns;
} ClockComparison;

/* L4: Polynomial clock model.
 * x(t) = x0 + y0*t + (D*t^2)/2 + (A*t^3)/6 + epsilon(t) */
double clock_model_time_error(const ClockParameters *params, double t);

/* L4: offset(t) = offset(t0) + y * (t - t0) * 1e9 */
double clock_offset_evolution(double initial_offset_ns,
                              double freq_offset_ppb, double dt_s);

/* L5: Least-squares linear regression on timestamp pairs */
int clock_linear_regression(const double *t_ref, const double *t_device,
                            int N, double *slope, double *intercept,
                            double *r_squared);

/* L5: Three-state clock model state transition */
void clock_state_transition(double state[3], double dt);

/* L5: Compare two clock state estimates */
ClockComparison clock_compare(const ClockParameters *a,
                              const ClockParameters *b,
                              double comparison_interval_s);

/* L5: Temperature-compensated frequency offset */
double clock_temp_compensated_freq(double base_freq_ppb,
                                   double temp_coeff_ppb_per_c,
                                   double quad_coeff_ppb_per_c2,
                                   double temperature_c,
                                   double ref_temperature_c);

/* L6: Predict clock offset after holdover duration */
double clock_predict_holdover_error(const ClockParameters *params,
                                    double holdover_duration_s);

#ifdef __cplusplus
}
#endif
#endif
