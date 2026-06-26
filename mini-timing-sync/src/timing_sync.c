/**
 * @file timing_sync.c
 * @brief Core timing synchronization implementation
 *
 * Implements timestamp arithmetic, two-way time transfer equations,
 * PI clock servo, Kalman clock tracking, holdover management,
 * and GPS 1PPS timing.
 */

#include "timing_sync.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * L1: Timestamp Utility Functions
 * =================================================================== */

double timing_timestamp_diff_ns(const Timestamp *ts1, const Timestamp *ts2)
{
    if (!ts1 || !ts2) return 0.0;
    double sec_diff = (double)(ts1->seconds - ts2->seconds);
    double ns_diff  = (double)(ts1->nanoseconds - ts2->nanoseconds);
    return sec_diff * 1.0e9 + ns_diff;
}

void timing_timestamp_add_ns(Timestamp *ts, int64_t ns_to_add)
{
    if (!ts) return;
    if (ns_to_add == 0) return;

    /* Handle negative addition */
    int64_t total_ns = (int64_t)ts->nanoseconds + ns_to_add;
    int64_t sec_adj = total_ns / 1000000000LL;
    int64_t ns_rem  = total_ns % 1000000000LL;

    /* Ensure ns_rem is non-negative */
    if (ns_rem < 0) {
        sec_adj -= 1;
        ns_rem  += 1000000000LL;
    }

    ts->seconds += sec_adj;
    ts->nanoseconds = (int32_t)ns_rem;

    /* Clamp to valid range */
    if (ts->seconds < 0) {
        ts->seconds = 0;
        ts->nanoseconds = 0;
    }
    if (ts->nanoseconds >= 1000000000) {
        ts->seconds += ts->nanoseconds / 1000000000;
        ts->nanoseconds = ts->nanoseconds % 1000000000;
    }
}

int timing_timestamp_cmp(const Timestamp *a, const Timestamp *b)
{
    if (!a || !b) return 0;
    if (a->seconds < b->seconds) return -1;
    if (a->seconds > b->seconds) return 1;
    if (a->nanoseconds < b->nanoseconds) return -1;
    if (a->nanoseconds > b->nanoseconds) return 1;
    return 0;
}

int timing_timestamp_valid(const Timestamp *ts)
{
    if (!ts) return 0;
    return (ts->seconds >= 0) &&
           (ts->nanoseconds >= 0) &&
           (ts->nanoseconds < 1000000000);
}

double timing_timestamp_to_double(const Timestamp *ts)
{
    if (!ts) return 0.0;
    return (double)ts->seconds * 1.0e9 + (double)ts->nanoseconds;
}

Timestamp timing_double_to_timestamp(double ns)
{
    Timestamp ts = {0, 0};
    if (ns < 0.0) return ts;

    /* Split into seconds and sub-second parts */
    double int_part;
    double frac_part = modf(ns / 1.0e9, &int_part);

    ts.seconds = (int64_t)int_part;
    ts.nanoseconds = (int32_t)(frac_part * 1.0e9 + 0.5);

    /* Handle rounding overflow */
    if (ts.nanoseconds >= 1000000000) {
        ts.seconds += 1;
        ts.nanoseconds -= 1000000000;
    }

    return ts;
}

/* ===================================================================
 * L4: Two-way Time Transfer Equations
 * =================================================================== */

int timing_compute_offset_delay(const PtpTimestamps *ts,
                                double *offset, double *delay)
{
    if (!ts || !offset || !delay) return -1;

    /* Validate timestamps */
    if (!timing_timestamp_valid(&ts->t1) ||
        !timing_timestamp_valid(&ts->t2) ||
        !timing_timestamp_valid(&ts->t3) ||
        !timing_timestamp_valid(&ts->t4)) {
        return -1;
    }

    /* Must have t1 < t2 (sync arrives after departure) */
    if (timing_timestamp_cmp(&ts->t1, &ts->t2) >= 0) return -1;
    /* Must have t3 < t4 (delay_resp arrives after delay_req) */
    if (timing_timestamp_cmp(&ts->t3, &ts->t4) >= 0) return -1;

    double t2_minus_t1 = timing_timestamp_diff_ns(&ts->t2, &ts->t1);
    double t4_minus_t3 = timing_timestamp_diff_ns(&ts->t4, &ts->t3);

    /* PTP two-way equations (IEEE 1588 Section 11.3):
     *   offset = ((t2 - t1) - (t4 - t3)) / 2
     *   delay  = ((t2 - t1) + (t4 - t3)) / 2
     */
    *offset = (t2_minus_t1 - t4_minus_t3) / 2.0;
    *delay  = (t2_minus_t1 + t4_minus_t3) / 2.0;

    return 0;
}

int timing_ntp_offset_delay(const Timestamp *T1, const Timestamp *T2,
                            const Timestamp *T3, const Timestamp *T4,
                            double *offset, double *delay)
{
    if (!T1 || !T2 || !T3 || !T4 || !offset || !delay) return -1;

    /* Validate all timestamps */
    if (!timing_timestamp_valid(T1) || !timing_timestamp_valid(T2) ||
        !timing_timestamp_valid(T3) || !timing_timestamp_valid(T4)) {
        return -1;
    }

    double T2_minus_T1 = timing_timestamp_diff_ns(T2, T1);
    double T3_minus_T4 = timing_timestamp_diff_ns(T3, T4);
    double T4_minus_T1 = timing_timestamp_diff_ns(T4, T1);
    double T3_minus_T2 = timing_timestamp_diff_ns(T3, T2);

    /* NTP equations (RFC 5905 Section 8):
     *   offset = ((T2 - T1) + (T3 - T4)) / 2
     *   delay  = (T4 - T1) - (T3 - T2)
     */
    *offset = (T2_minus_T1 + T3_minus_T4) / 2.0;
    *delay  = T4_minus_T1 - T3_minus_T2;

    /* Delay should be non-negative */
    if (*delay < 0.0) *delay = 0.0;

    return 0;
}

/* ===================================================================
 * L2: PI Clock Servo (Clock Discipline Algorithm)
 * =================================================================== */

void pi_servo_init(PiServoConfig *cfg, PiServoState *state)
{
    if (!cfg || !state) return;

    /* Default tuning: critically damped with 1000s time constant */
    cfg->Kp = 0.7;
    cfg->Ki = 0.3;
    cfg->integral_limit_ns = 1.0e6;   /* 1 ms max integral */
    cfg->max_correction_ns = 1.0e5;   /* 100 us max per step */

    state->integral_ns = 0.0;
    state->first_update = 1;
}

double pi_servo_update(PiServoConfig *cfg, PiServoState *state,
                       double measured_offset_ns)
{
    if (!cfg || !state) return 0.0;

    /* On first update, just measure without correcting */
    if (state->first_update) {
        state->first_update = 0;
        state->integral_ns = 0.0;
        return 0.0;
    }

    /* Update integral term */
    state->integral_ns += cfg->Ki * measured_offset_ns;

    /* Anti-windup: clamp integral */
    if (state->integral_ns > cfg->integral_limit_ns) {
        state->integral_ns = cfg->integral_limit_ns;
    } else if (state->integral_ns < -cfg->integral_limit_ns) {
        state->integral_ns = -cfg->integral_limit_ns;
    }

    /* PI control law: correction = Kp * offset + integral */
    double correction = cfg->Kp * measured_offset_ns + state->integral_ns;

    /* Clamp per-step correction */
    if (correction > cfg->max_correction_ns) {
        correction = cfg->max_correction_ns;
    } else if (correction < -cfg->max_correction_ns) {
        correction = -cfg->max_correction_ns;
    }

    return correction;
}

void pi_servo_reset(PiServoState *state)
{
    if (!state) return;
    state->integral_ns = 0.0;
    state->first_update = 1;
}

/* ===================================================================
 * L3: Kalman Filter for Clock Tracking
 * =================================================================== */

void clock_state_init(ClockState *state, double initial_offset_ns,
                      double initial_freq_ppb, double initial_drift_ppb_per_day,
                      double offset_uncertainty_ns,
                      double freq_uncertainty_ppb)
{
    if (!state) return;

    state->offset_ns = initial_offset_ns;
    state->freq_offset_ppb = initial_freq_ppb;
    state->drift_ppb_per_day = initial_drift_ppb_per_day;
    state->last_update_time = 0;

    /* Initialize covariance matrix as diagonal */
    memset(state->P, 0, sizeof(state->P));
    state->P[0][0] = offset_uncertainty_ns * offset_uncertainty_ns;
    state->P[1][1] = freq_uncertainty_ppb * freq_uncertainty_ppb;
    state->P[2][2] = 0.01; /* Small initial drift uncertainty */

    /* Small cross-correlation between offset and frequency */
    state->P[0][1] = 0.1 * offset_uncertainty_ns * freq_uncertainty_ppb;
    state->P[1][0] = state->P[0][1];
}

void clock_kalman_predict(ClockState *state, double dt,
                          const double Q[9])
{
    if (!state || !Q || dt <= 0.0) return;

    double F[3][3] = {
        {1.0, dt,    0.5 * dt * dt},
        {0.0, 1.0,   dt           },
        {0.0, 0.0,   1.0          }
    };

    /* State prediction: x_pred = F * x */
    double x_pred[3];
    x_pred[0] = F[0][0] * state->offset_ns +
                F[0][1] * state->freq_offset_ppb +
                F[0][2] * state->drift_ppb_per_day;
    x_pred[1] = F[1][0] * state->offset_ns +
                F[1][1] * state->freq_offset_ppb +
                F[1][2] * state->drift_ppb_per_day;
    x_pred[2] = F[2][0] * state->offset_ns +
                F[2][1] * state->freq_offset_ppb +
                F[2][2] * state->drift_ppb_per_day;

    /* Covariance prediction: P_pred = F * P * F^T + Q */
    double FPFt[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                double pf_sum = 0.0;
                for (int l = 0; l < 3; l++) {
                    pf_sum += state->P[k][l] * F[j][l];
                }
                sum += F[i][k] * pf_sum;
            }
            FPFt[i][j] = sum;
        }
    }

    /* Add process noise Q (row-major, 3x3) */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            state->P[i][j] = FPFt[i][j] + Q[i * 3 + j];
        }
    }

    /* Update state */
    state->offset_ns = x_pred[0];
    state->freq_offset_ppb = x_pred[1];
    state->drift_ppb_per_day = x_pred[2];
}

void clock_kalman_update(ClockState *state, double measured_offset_ns,
                         double measurement_noise_var_ns2)
{
    if (!state || measurement_noise_var_ns2 <= 0.0) return;

    /* Measurement model: z = [1, 0, 0] * x + v
     * H = [1, 0, 0], so H*P = [P00, P10, P20] */
    double R = measurement_noise_var_ns2;

    /* Innovation: y = z - H * x */
    double y = measured_offset_ns - state->offset_ns;

    /* Innovation covariance: S = H * P * H^T + R */
    double PHt[3];
    for (int i = 0; i < 3; i++) {
        PHt[i] = state->P[i][0]; /* H * P = [P00, P10, P20]^T since H=[1,0,0] */
    }
    double S = state->P[0][0] + R;

    /* Kalman gain: K = P * H^T / S */
    double K[3];
    for (int i = 0; i < 3; i++) {
        K[i] = PHt[i] / S;
    }

    /* State update: x = x + K * y */
    state->offset_ns += K[0] * y;
    state->freq_offset_ppb += K[1] * y;
    state->drift_ppb_per_day += K[2] * y;

    /* Covariance update: P = (I - K*H) * P */
    double newP[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            newP[i][j] = state->P[i][j] - K[i] * state->P[0][j];
        }
    }

    /* Enforce symmetry */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            state->P[i][j] = (newP[i][j] + newP[j][i]) / 2.0;
        }
    }
}

double clock_predict_offset(const ClockState *state, double dt_s)
{
    if (!state) return 0.0;
    /* x(t + dt) = offset + freq * dt + 0.5 * drift * dt^2 */
    return state->offset_ns +
           state->freq_offset_ppb * dt_s +
           0.5 * state->drift_ppb_per_day * (dt_s / 86400.0) * dt_s;
}

double clock_get_freq_ppb(const ClockState *state)
{
    if (!state) return 0.0;
    return state->freq_offset_ppb;
}

/* ===================================================================
 * L6: Holdover Management
 * =================================================================== */

int holdover_should_enter(double offset_ns, const HoldoverConfig *cfg)
{
    if (!cfg) return 0;
    return (fabs(offset_ns) > cfg->entry_threshold_ns) ? 1 : 0;
}

int holdover_can_exit(double offset_ns, const HoldoverConfig *cfg,
                      double holdover_duration_s)
{
    if (!cfg) return 0;

    /* Check if max duration exceeded */
    if (holdover_duration_s > cfg->max_holdover_duration_s) {
        return 0; /* Stay in holdover until manually recovered */
    }

    /* Check if offset is below exit threshold */
    if (fabs(offset_ns) < cfg->exit_threshold_ns) {
        return 1;
    }

    return 0;
}

double holdover_estimate_uncertainty(double holdover_duration_s,
                                     double drift_ppb_per_day,
                                     double aging_ppb_per_day2,
                                     double initial_error_ns)
{
    /* Convert drift to ns/s: drift_ppb_per_day * 1e-9 per day
     *   = drift_ppb_per_day * 1e-9 / 86400 per second
     * Aging: aging_ppb_per_day2 * 1e-9 / (86400^2) per second^2
     */
    double drift_ns_per_s = drift_ppb_per_day * 1.0e-9 / 86400.0;
    double aging_ns_per_s2 = aging_ppb_per_day2 * 1.0e-9 /
                             (86400.0 * 86400.0);

    /* Deterministic model:
     *   error(t) = initial + drift * t + 0.5 * aging * t^2
     */
    double t = holdover_duration_s;
    double deterministic_error =
        initial_error_ns +
        drift_ns_per_s * t * 1.0e9 +
        0.5 * aging_ns_per_s2 * t * t * 1.0e9;

    /* Add stochastic component: random-walk FM adds uncertainty ~ sqrt(t) */
    double stochastic_error = 1.0 * sqrt(t); /* 1 ns/sqrt(s) typical */

    return fabs(deterministic_error) + stochastic_error;
}

/* ===================================================================
 * L7: GPS Timing ? 1PPS Discipline
 * =================================================================== */

int gps_1pps_compute_offset(const Timestamp *local_timestamp,
                            const GpsTimingConfig *cfg,
                            double *offset)
{
    if (!local_timestamp || !cfg || !offset) return -1;
    if (!timing_timestamp_valid(local_timestamp)) return -1;

    /* 1PPS edge should align with integer second boundary.
     * The nanosecond field gives the offset from the second boundary.
     */
    *offset = (double)local_timestamp->nanoseconds;

    /* Compensate for cable and receiver delays */
    *offset -= cfg->antenna_cable_delay_ns;
    *offset -= cfg->receiver_processing_delay_ns;

    /* Wrap offset to [-0.5s, +0.5s] range */
    if (*offset > 5.0e8) {
        *offset -= 1.0e9;
    } else if (*offset < -5.0e8) {
        *offset += 1.0e9;
    }

    return 0;
}

double gps_sawtooth_correction(double raw_offset_ns,
                               double sawtooth_amplitude_ns,
                               double phase_radians)
{
    /* Sawtooth waveform: saw(t) = amplitude * (phase / (2*pi) - 0.5)
     * Correction subtracts the sawtooth from the raw offset.
     *
     * For clock quantization of T_q:
     *   sawtooth_amplitude = T_q (peak-to-peak)
     *   phase = 2*pi * (t_local mod T_q) / T_q
     */
    double phase_normalized = phase_radians / (2.0 * M_PI);
    /* Normalize to [0, 1) */
    phase_normalized = phase_normalized - floor(phase_normalized);
    /* Sawtooth: linear ramp from -amplitude/2 to +amplitude/2 */
    double sawtooth_value = sawtooth_amplitude_ns * (phase_normalized - 0.5);

    return raw_offset_ns - sawtooth_value;
}
